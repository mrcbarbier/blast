
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>



#include "specs.h"

//////////////////////////// DECOUPAGE EN BOITES //////////////



long int overlap=0 ;
double epsilon = 0.00000001;
int nbox[dim]=NBOX;

int fbox(double * r, double *r0, double Rmax)
{	
	double resp = pow(nbox[dim-1], 1./dim)+epsilon, temp, specoord[dim]; // coordonees dans une base speciale ex polaire
	int result=0, i,j,k ;
	double x[dim] , norm=0.;


	if(boxtyp == 0) {
		for(k=0;k<dim;k++)
		{
			x[k]=r[k];
			if(k>0)
				result +=  floor( (x[k]+0.5) *(nbox[k]/nbox[k-1])) *nbox[k-1];
			else
				result +=  floor( (x[k]+0.5) *nbox[k]) ;
			
		}
		
	}

	if(boxtyp==1)
	{
		for(i=0;i<dim;i++)
		{
			x[i] = (r[i]-r0[i])/Rmax ;
			norm += x[i]*x[i];

		}
		norm = sqrt(norm);

		if(norm>1.)
			return nbox[dim-1]+1;


		specoord[0]=0. ;
		for(k=0;k<dim ; k++)
			specoord[0] += x[k]*x[k] ;
		specoord[0]=sqrt(specoord[0]); // (sqrt(dim)/2.);
		specoord[0]=pow(specoord[0],dim) ;
		if(dim>2){
			for(i=1;i<dim-1;i++)
			{
				temp = 0. ;
				for(k=dim-1;k>=i;k--)
					temp+=x[i]*x[i];
				specoord[i] = atan2( sqrt(temp) , x[i-1] );
			}
		}
		specoord[dim-1]= atan2(x[dim-1],x[dim-2]);
		for(k=0;k<dim;k++)
		{
			if(k==dim-1)
				specoord[k] = specoord[k]/(2*PI)+0.5;
			if(k>0 && k<dim-1)
				specoord[k] /= PI ;
			
			if(k>0)
			{temp =nbox[k]/nbox[k-1] +epsilon;
				result +=  floor( specoord[k] *temp ) *nbox[k-1];}
			else
				result +=  floor(  specoord[k] *nbox[k]) ;
		}

	}

	//printf("%i %f %f\n",result,norm, Rmax);
	return result ;
}



double xbox(int b1, int k, double *r0, double Rmax)
{
	double resp = pow(nbox[dim-1], 1./dim)+epsilon;
	double scale, result, specoord[dim]  ,temp; // coordonees dans une base speciale ex polaire
	int i ;
	if(boxtyp == 0)
	{
		 
		if(k>0)
		{ 	temp =nbox[k]/nbox[k-1] +epsilon;
			return  ( (   b1/nbox[k-1]  )    % ( (int) temp )    +0.5) /(temp ) -0.5;
		}
		else
			return  (   b1    %  nbox[0]     +0.5) /nbox[0] -0.5;
	}
	if(boxtyp == 1) {
		for(i=0;i<dim;i++)
		{
		 
			if(i>0)
			{ 	temp =nbox[i]/nbox[i-1] +epsilon;
				specoord[i] =   ( (   b1/nbox[i-1]  )    % ( (int) temp )    +0.5) /temp  ;
			}
			else
				specoord[i] =   (   b1    %  nbox[0]     +0.5) /nbox[0] ;


			if(i==0)
				specoord[i] = pow(specoord[i],1./dim); // ( sqrt(dim)/2.) ;
			if(i==dim-1)
				specoord[i] = (specoord[i]-0.5)*2*PI;
			if(i>0 && i<dim-1)
				specoord[i] *= PI ;
			//printf("\n");
			//printf("%i : %f %f %i\n",i,specoord[i], b1/( (int) pow(resp,i)  )    % ( (int) resp )  +0.0 , b1 );
		}
		
		result = specoord[0] ;
		if(k>0 && dim>2){
			for(i=1;i<=k;i++)
				result *= sin(specoord[i]) ;
		}
		//if(k==dim-1)
	//		result *= sin(specoord[k]);
		if(k!=dim-1)
			result *= cos(specoord[k+1]);
		if(k==1 && dim ==2)
			result *= sin(specoord[k]);
	}
	//printf("%i : %f\n",k,result);
	return result *Rmax+ r0[k];
}

////////////////////////////////// CHAMPS BOITES ET RADIAUX //////////////////////////////////////////////////////////:

void r0calc(double ** r, double ** v, double * r0) // calcul de la position du centre de masse des particules en mouvement
{
	int i=0,j=0,N =0;
	for (i=0;i<dim;i++)
		r0[i]=0. ;
	for(i=0 ; i<nu ;i++)
	{	
		if(v[i][dim]>0){
			for(j=0 ; j<dim ; j++)	
			{
			
				r0[j]+= r[i][j];
				N +=1;
			}
		}
	}
	for (i=0;i<dim;i++)
		r0[i]/=N ;

}

double er(int j, double * r, double * r0)
{
	return (r[j]-r0[j])/r[dim];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//---------------------------------------- Autres fonctions --------------------------------------------------------------//////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


double FastArcTan(double x)
{
    return PI*x/4. - x*(fabs(x) - 1)*(0.2447 + 0.0663*fabs(x));
}


double FastArcTan2(double y,double x)
{
    
    if (x>0) return FastArcTan(y/x);
    if (x<0 && y>=0) return FastArcTan(y/x) +PI;
    if (x<0 && y<0) return FastArcTan(y/x) - PI;
    if (x==0 && y>0 ) return PI/2;
    else return -PI/2;
    
}

void convrad(double * r, double *r0)

{
	
	double specoord[dim]; // coordonees dans une base speciale ex polaire
	int i,k ;
	double x[dim],temp ;


	specoord[0]=0. ;
	for(i=0;i<dim ; i++)
	{
		x[i] = (r[i]-r0[i]) ;
		specoord[0] += x[i]*x[i] ;
	}
	specoord[0]=sqrt(specoord[0]);

	if(dim>2){
		for(i=1;i<dim-1;i++)
		{
			temp = 0. ;
			for(k=dim-1;k>=i;k--)
				temp+=x[i]*x[i];
			specoord[i] = FastArcTan2( sqrt(temp) , x[i-1] );
		}
	}
	specoord[dim-1]= FastArcTan2(x[dim-1],x[dim-2]);
	for(k=0;k<dim;k++)
		r[dim+k]=specoord[k];	

}


double rand01()
{
	// Function : random number between 0 and 1

	double scale=RAND_MAX+1.;
	double base=rand()/scale;
	double fine=rand()/scale;
	return base+fine/scale;
}

double timeval(double * r1, double * r2, double * v1, double * v2, double sigma)
{
	// Function : evaluates the time before a collision between two particles 1 and 2
	// given radius sigma, positions and velocities

	double SIGMA = 2*sigma ; // here we need the SUM of both radii !

	double r12[dim], v12[dim], b12=0., vv12=0.,rr12=0., disc, tij ;
	int i, j ;
	FILE * fout ;

	for(i=0;i<dim;i++)
	{
		if(walltype==0 || (walltype==3 && i!=0) ){
			r12[i]=r1[i]-r2[i]-floor((r1[i]-r2[i])+0.5);
			if(fabs(r12[i])>0.5)
				printf("Error : r12 too big \n");}
		else
			r12[i]=r1[i]-r2[i];

		v12[i]=v1[i]-v2[i];
		b12 += r12[i]*v12[i];
		vv12 += v12[i]*v12[i];
		rr12 += r12[i]*r12[i] ;
	}

	if (b12 <0)
	{
		disc = b12*b12 - vv12*(rr12-SIGMA*SIGMA);
		if (disc>0)
		{
			tij = (- b12 - sqrt(disc) )/vv12;
			if(tij<epsilon&&tij>0)
			{
				//printf("erreur tij=0\n");
				overlap+=1;
				//printf("overlap : %f\n",tij)
				return -1;
			}
			return tij ;
		}
		/*else
		{
			if(SIGMA - sqrt(rr12) > epsilon * SIGMA)
			{
				overlap += 1 ;
			}
		}*/


	}
	return -1 ;
}

/*
int main()
{
    double test[2*dim],blank[dim];
    int i;
    double sum=0;
    printf("Testing fonctions.cpp\n");
    for (i=0;i<1000;i++)
    {
        test[0]=rand01() ;
        test[1]=rand01() ;
        test[2]=0;
        test[3]=0;
        sum=test[0]+test[1];
        convrad(test,blank);
        if (sum!= test[0]+test[1])
            printf("Error in convrad\n");
        
    }
    return 1;
}

*/

