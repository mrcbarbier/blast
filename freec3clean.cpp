#include <iostream>
#include <sstream>
using std::cerr ;
using std::cout ;
#include <fstream>
using std::ifstream;
#include <cstdlib>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#include "specs.h"
#include "fonctions.h"

// MOLECULAR DYNAMICS 

// Global variables (some defined in fonctions.cpp)

double M2 = 1. ;

int * blastcol[nsys]; //is the particle part of the blast?

double timbig = pow(10,10); //  infinite collision time

/* Possible optimizations:
- Do not recompute global quantities every "evolpas" steps (search for Label01 below)
-

*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//----------------------------------------  Important variables summary ----------------------------------------------------////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*

Everything is indexed by system, then by particle

LOGIC:
Geometry
    r[nsys][nu][2dim]    positions (sys; particle; cartesian then polar coordinates [x,y,...,r,theta,...] )
    v[nsys][nu][2dim]    velocities (same)
    r0[nsys][dim], v0[nsys][dim]   center of mass properties (cartesian only)
    
Collision
    coltim[nsys][nu]     collision time
    partner[nsys][nu]    collision partner
    tfirst[nsys]         time before the next collision in the system
    nextp1,2[nsys]       identities of the next collision pair


MEASURED QUANTITIES:
Collisional

    dissipation[nsys][nu] energy dissipated per particle
    angmov[nsys][nu] angular motion
    collisions[nsys][nu] number of collisions per particle
    blastcol[nsys][nu] number of collisions *with blast particles* per particle

Global
    
    E[nsys] energy
    E0[nsys] initial energy (in principle always 1)
    momx[nsys] moment of order x of the velocity distribution along dimension dobs
    a2[nsys] Sonine moment a2

*/






////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//---------------------------------------- Initialisation ------------------------------------------------------------------////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double init(double **r[nsys], double **v[nsys], double sigma, double m[nu],
    double *coltim[nsys], int *partner[nsys], double tfirst[nsys] ,
    int *blastcol[nsys], int *collisions[nsys], double *dissipation[nsys],double *angmov[nsys],
    double *r0[nsys], double *v0[nsys], double E0[nsys],long int nextp1[nsys],long int nextp2[nsys] )  
{
    int * ptinit[nsys],i,j,k,p1,p2,sys ;
    double test = 0.,tij=0 ,psi,theta,ptot[dim],Rmix;
    
    double blank[dim];
    for(i=0;i<dim;i++)
        blank[i]=0. ;
    
    for(i=0;i<nsys;i++){
          if(partinit>=1)
            ptinit[i] = (int *) malloc ( (int)partinit * sizeof(int ) );  
    }

    //INITIALISATION DES MASSES, POSITIONS ET VITESSES
    if(stopr==0) printf("Rayon des particules : %.16g\n",sigma);
    string line;
    ifstream myfile ("initconfig.dat");

    //Initialize positions from file
    for (sys=0; sys<nsys ; sys++)
    {
        if(partinit>=1){
            for(i=0;i<(int) partinit;i++)
            {
                ptinit[sys][i]=i;
            }
        }
        Rmix=1. ;E0[sys]=0 ;
        for (i=0; i<nu ; i++)
        {
            v[sys][i] = (double *) malloc ( 2*dim * sizeof(double) );
            r[sys][i] = (double *) malloc ( 2*dim * sizeof(double) );
            
            blastcol[sys][i]=0;
            collisions[sys][i]=0 ;    
            dissipation[sys][i]=0 ;    
            angmov[sys][i]=0 ;    
            
            getline (myfile,line);
            istringstream iss(line);
            j=0 ;
            do
            {
                string sub;
                iss >> sub;
              
                r[sys][i][j]=atof(sub.c_str()) ;
                j+= 1;
                if(j==dim) break;
            } while (iss);
            convrad(r[sys][i],r0[sys]);
        }
        
    }
    //Look for particles initially in the blast (if given as a number)
    for (sys=0; sys<nsys ; sys++)
    {    for (i=0; i<nu ; i++)
        {
            if (partinit>=1){
                for(k=0;k<partinit;k++)
                {
                    if(r[sys][i][dim]<r[sys][ptinit[sys][k]][dim]){
                        ptinit[sys][k]=i;
                        break;
                    }
                }
            }
        }
    }

    //Initialize masses
    for (i=0; i<nu ; i++)
    {
            if (masslaw==0)
                m[i] = 1. ;
            else
            {    
                m[i]= pow(r[0][i][dim],masslaw) ; 
                //initial position is the same in all systems anyway
                test+=m[i];
            }
    }
    if (masslaw!=0)
    { 
        for (i=0; i<nu ; i++)
            m[i]/=test;
    }
    
    
    //Initialize velocities
    
    for (sys=0; sys<nsys ; sys++)
    {

        for (i=0; i<nu ; i++)
        {
            if(partinit<1.) //partinit given as radius (or width for planar shock)
            {
                if((planar==0 && r[sys][i][dim]<partinit) ||
                    (planar==1 && r[sys][i][0]< -0.5+ partinit))
                {
                    // If within the initial blast region
                    blastcol[sys][i]=1;
                    if(planar==0)
                    { //Radial shock
                        for(j=0;j<dim;j++)
                        {
                            if(radinit==0)
                            {
                                psi = rand01() ;
                                theta = rand01()*2*PI;
                                v[sys][i][j] = sqrt(-2*log(psi))*cos(theta);
                            }
                            else
                            {
                                v[sys][i][j]= r[sys][i][j]/r[sys][i][dim]*rand01();
                            }
                        }
                    }
                    else
                    { //Planar shock (velocity in dimension 0 only)
                        for(j=1;j<dim;j++)
                            v[sys][i][j]=0.;
                        v[sys][i][0]= 1.;
                    }
                }
                else
                {
                    // Outside the initial blast region
                    if( blastMach<0)
                    {   // Infinite Mach - zero temperature outside
                        for(j=0;j<dim;j++)
                           v[sys][i][j]=0.;
                    }
                    else
                    {
                        for(j=0;j<dim;j++)
                        {
                            psi = rand01() ;
                            theta = rand01()*2*PI;
                            v[sys][i][j] = sqrt(-2*log(psi))*cos(theta) / blastMach;
                        }
                    }
                }
            }
            else
            {
                //partinit given as a number of particles
                v[sys][i][j]=0.;
                for(k=0;k<partinit;k++)
                {
                    if(i==ptinit[sys][k]){
                        psi = rand01() ;
                        theta = rand01()*2*PI;
                        v[sys][i][j] = sqrt(-2*log(psi))*cos(theta);
                        blastcol[sys][i]=1;
                    }
                }
            }
            
            // Total momentum
            for(j=0;j<dim;j++)
               ptot[j] += v[sys][i][j]*m[i];
            
            convrad(v[sys][i],blank);
            if (blastcol[sys][i]==0 && r[sys][i][dim]<Rmix) Rmix=r[sys][i][dim];
            // Total energy
            E0[sys]+=m[i]*v[sys][i][dim]*v[sys][i][dim]/2;
        }
        if(partinit>1)
        {
            for(j=0;j<dim;j++)
                v[sys][ptinit[sys][0]][j]-=ptot[j]/m[ptinit[sys][0]] ; 
                //set momentum of center of mass to 0 if more than one particle
            convrad(v[sys][ptinit[sys][0]],blank);
        }
        
        
    }
    

    for(sys=0;sys<nsys;sys++)
    {
        for(i=0;i<nu;i++)
        {
            for(j=0;j<dim;j++)
            {
                v[sys][i][j]= v[sys][i][j]/ sqrt(E0[sys]);
                v0[sys][j]+= v[sys][i][j];
            }
        }

        // TEMPS DE PREMIERE COLLISION
        for (i= 0; i<nu ; i++)
        {
            coltim[sys][i] = timbig ;
            partner[sys][i] = nu+1;
        }


        for (i= 0; i<nu; i++)
        {
            for(j=i+1 ; j < nu ; j ++)
            {
                //if(v[sys][i][dim] > epsilon || v[sys][j][dim] > epsilon ){
                if( blastcol[sys][i]>0 ||  blastcol[sys][j]>0 ){
                    tij = timeval(r[sys][i],r[sys][j],v[sys][i],v[sys][j], sigma);
                    
                    if(tij>0.)
                    {
                        if (tij < coltim[sys][i])
                        {
                            coltim[sys][i] = tij ;
                            partner[sys][i] = j ;
                        }
                        if (tij < coltim[sys][j])
                        {
                            coltim[sys][j] = tij ;
                            partner[sys][j] = i ;
                        }
                        if(tij < tfirst[sys])
                        {
                            tfirst[sys] = tij;

                            p1 = i;
                            p2 = j;
                            nextp1[sys]=p1;
                            nextp2[sys]=p2;
                        }
                    }
                
                
                    double rdist=0 ;
                 }

            }

        }
        free(ptinit[sys]);
    }

    return Rmix;
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//---------------------------------------- Corps du programme --------------------------------------------------------------////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





int main()
{

    // VARIABLE DEFINITIONS =============================================================================
    
	double *coltim[nsys], *tcounter[nsys] ;
    int *partner[nsys], Nexp[nsys] ;
    double m[nu] ;
    long int i, j, k,l, p1, p2, pp1, pp2, coll=0 ;
    int sys ;
    double tij, tfirst[nsys]; 
    
    double ** v[nsys], ** r[nsys], *r0[nsys], *v0[nsys], *dissipation[nsys], *angmov[nsys];
    int *collisions[nsys],* ptinit[nsys] ;


    double V = pow(size, dim) ;
    double sigma ;
    if(dim ==1)
        sigma = rdens*V/nu ;
    if( dim == 2)
        sigma = sqrt(rdens*V/nu/PI);
    if (dim == 3)
        sigma = pow(rdens *V/nu/(4.*PI/3.) , 1./3.);
        
        
    double sig12[dim], v12[dim],signorm,norm,norm1,norm2,normoy, pscal;


    double t0[nsys], tphys =0. ;
    double E0[nsys], E[nsys], a2[nsys], mom[nsys], Rexp[nsys], mu1[nsys],mu2[nsys], mu3[nsys], mu4[nsys];
    double moyt, moytphys, moyRexp,moyNexp,moyE, moyErad, moymom ; // moyenne sur les systemes

    double  dE,dErad, dissip=0.,Ecm=0. ; //dissipation
    double vrelative=0, tmoy = 0, lasttime=0.,align=0., vec1=0., vec2=0., vec3=0., vec4=0.;
    double Erad[nsys], vrad[nsys] ; // composantes radiales
    long int nextp1[nsys],nextp2[nsys];

    double R,R1,R2,Rmix, dif, test, psi, theta , test2,test3, tc,tlast=0;
    int  nactiv;

    double gamma0 = (1.0-alpha*alpha)/(2.0*dim) ;
    if (alpha==1)
        gamma0=1.0 ;


    // Surface of unit sphere in dimension d
    double omega[10];
    for(i=0;i<10;i++)
        omega[i]=0.;
    omega[1]= 1. ;
    omega[2] = 2*PI ;
    omega[3] = 4*PI ;
    omega[4] = 2*PI*PI ;
    omega[5] = 8*PI*PI/3.0 ;
    omega[6] = PI*PI*PI ;

    // time scale = mean free path / sqrt( energy per particle)
    // mean free path = sigma * d-dimensional volume / rdens / khi (enskog factor)
    // hence time scale = ( sqrt(nu / E0[sys]) * sigma*omega[dim]/dim/rdens/khi)
    // For homogeneous elastic hard spheres, this gives back the collision frequency in Trizac (2 sqrt(pi)) for any nu and rdens
    

    double phi = rdens*PI/4. ;
    double khi =  (1.-7./16. * phi)/(1.-phi)/(1.-phi);    

    double mfp = sigma * omega[dim]/rdens/khi ;
    printf("Mean free path : %g\n",mfp);

    time_t start,end,rawtime;
    struct tm * ptm;

    srand( (unsigned)time( NULL ) );

    int nb= nbstart ;
    char affix[10] ;
    char root[30] ;
    char number[100];

    time (&start);

    FILE *fout;
    FILE *fin;
    fout=fopen("results/S.dat", "w");fclose(fout);        
    fout=fopen("results/evT.dat", "w");fclose(fout);
    fout=fopen("results/evdens.dat", "w");fclose(fout);

    fout=fopen("results/overlap.dat", "w");
    fprintf(fout,"\n");
    fclose(fout);
    if((fout=fopen("results/t.dat", "w"))==NULL) 
                        printf("Cannot open file.\n");
    fclose(fout);
    

    
    for(i=0;i<nsys;i++){

        tfirst[i] = timbig ;
        blastcol[i] = (int *) malloc ( nu * sizeof(int) );
        dissipation[i] = (double *) malloc ( nu * sizeof(double ) );
        angmov[i] = (double *) malloc ( nu * sizeof(double ) );
        collisions[i] = (int *) malloc ( nu * sizeof(int ) );
        r0[i] = (double *) malloc ( dim * sizeof(double) );
        v0[i] = (double *) malloc ( dim * sizeof(double) );
        v[i] = (double **) malloc ( nu * sizeof(double *) );
        r[i] = (double **) malloc ( nu * sizeof(double *) );

        tcounter[i] = (double *) malloc ( nu * sizeof(double) );
        coltim[i] = (double *) malloc ( nu * sizeof(double) );
        partner[i] = (int *) malloc ( nu * sizeof(int) );
        for(k=0;k<dim;k++){
            r0[i][k]=0.;
            v0[i][k]=0.;}
        t0[i]=0.;
        E0[i]=0.;
        E[i]=0.;
        mu1[i]=0.;
        mu2[i]=0.;
        mu3[i]=0. ;
        mu4[i]=0.;
        vrad[i]=0. ;
        Erad[i] = 0. ;
    }
    double blank[dim];
    for(i=0;i<dim;i++)
        blank[i]=0. ;


    // ===================== APPEL FONCTION D'INITIALISATION ====================================================
    
    Rmix=init(r, v,  sigma, m,coltim,partner, tfirst,blastcol,collisions,dissipation,angmov,r0,v0,E0,nextp1,nextp2);
    time (&end);
    dif = difftime (end,start);
    if (stopr==0) printf("Boucle (début : %g min) \n", dif/60);


    // ==================  INITIALISATION DES QUANTITES MESUREES DANS LA BOUCLE ==============================
    for(sys=0;sys<nsys;sys++)
    {
        for(i=0;i<nu;i++)
        {
            for(j=0;j<dim;j++)
            {
                E[sys]+=m[i]*v[sys][i][j]*v[sys][i][j] /2.;
            }
            mu1[sys] += (v[sys][i][dobs] ) /nu ;
            mu2[sys] += (v[sys][i][dobs]*v[sys][i][dobs])/nu ;
            mu3[sys] += (v[sys][i][dobs]*v[sys][i][dobs]*v[sys][i][dobs])/nu;
            mu4[sys] += (v[sys][i][dobs]*v[sys][i][dobs]*v[sys][i][dobs]*v[sys][i][dobs])/nu ;

            
        }
        a2[sys] = dim*dim/4.*mu4[sys]/mu2[sys]/mu2[sys] ;
        a2[sys] = a2[sys]*( 4./dim/(dim+2.) )-1. ;
        mom[sys] = 8*mu3[sys] - 12*mu1[sys]; 

        E0[sys]=E[sys];
        t0[sys] /= ( sqrt(nu / E0[sys]) * sigma*omega[dim]/dim/rdens/khi) ; //pour le rechargement, qui se fait a partir du temps physique

        //printf("tboltz : %g, gamma0 : %g, khi : %g\n",gamma0*colboltz, gamma0, khi);
    }
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////
    //---------------------------------------- BOUCLE --------------------------------------------////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////


    tphys = 0. ;
    while (coll < ncoll && tphys < tmax)
    {
       
        for(sys=0;sys<nsys;sys++)
        {

            i = 0 ;
            j = 0 ;
            pp1 = 0;
            pp2 = 0 ;
            p1 = nextp1[sys] ;
            p2 = nextp2[sys] ;
            pscal = 0.;
            for(i=0;i<dim;i++)
                sig12[i]=0. ;

            if(tfirst[sys]>=timbig-epsilon)
            {
                printf("Tij = TIMBIG : coll %li\n", coll);
                return 0;
            }

            signorm = 0. ; norm = 0. ; norm1=0. ; norm2=0.; normoy =0. ;

            for(i=0;i<dim;i++)
            {
                sig12[i]=r[sys][p2][i]-r[sys][p1][i] ;

                v12[i]=v[sys][p2][i]-v[sys][p1][i];
                signorm += sig12[i]*sig12[i] ;
                norm += v12[i]*v12[i] ;
                norm1 +=v[sys][p1][i]*v[sys][p1][i];
                norm2 +=v[sys][p2][i]*v[sys][p2][i];
                normoy += (v[sys][p1][i] + v[sys][p2][i])*(v[sys][p1][i] + v[sys][p2][i]) / 4. ;
            }
            j=0;


            //============================================== DEPLACEMENT ============================
            
            for(i=0;i<nu;i++)
            {                
                
                //printf(" %i %g %g %g %g\n",i,r[sys][i][0],r[sys][i][1],r0[sys][0]),r0[sys][1];
                convrad(r[sys][i],r0[sys]);
                if (coll == colwait) dissipation[sys][i]=0. ;
                coltim[sys][i] -= tfirst[sys] ;
                tcounter[sys][i]+= tfirst[sys] ;
                for(j=0;j<dim;j++)
                {
                    angmov[sys][i] += r[sys][i][dim]*v[sys][i][j]*(1-er(j,r[sys][i],r0[sys]));
                    r[sys][i][j] = r[sys][i][j] + v[sys][i][j]*tfirst[sys];
                    if (walltype == 3 && -0.5 + 0.05*t0[sys] > pow(nu,1./dim)*3*sigma) {coll = ncoll ;  printf("Piston au maximum\n");}
                    
                    if(r[sys][i][j] > 0.5 || r[sys][i][j]<-0.5 || (walltype ==3 && r[sys][i][j]<-0.5 + 0.05*t0[sys] ) ){
                        if (walltype==2) //sticky
                        {
                            for(k=0;k<dim;k++) v[sys][i][k]=0. ;
                        }
                        else
                        {
                            if (v[sys][i][j] != 0. && blastMach < 0 && planar != 1 ){
                                printf("Bord atteint !\n");//coll,i,v[sys][i][0],v[sys][i][1],v[sys][i][dim]);
                                coll = ncoll ; 
                            } 
                            else
                            {    if (walltype==0) //periodique
                                {
                                    r[sys][i][j]-=floor(r[sys][i][j]+0.5);
                                }
                                if (walltype == 1)//reflexion speculaire
                                {    v[sys][i][j]*=-1 ;
                                    if(r[sys][i][j] > 0.5 ) r[sys][i][j]=0.5-sigma ;
                                    else r[sys][i][j]=-0.5+sigma ;
                                }
                                if (walltype==3) //piston
                                {
                                    if (j!=0) r[sys][i][j]-=floor(r[sys][i][j]+0.5);
                                    else 
                                    {    if (v[sys][i][j] <0) v[sys][i][j]*=-1 ;
                                        if(r[sys][i][j] > 0.5 ) r[sys][i][j]=0.5-sigma ;
                                        else r[sys][i][j]=-0.5+ 0.05*t0[sys]+sigma ;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //================================================== COLLISION ========================

            pscal = 0. ; align = 0. ;
            for(i=0;i<dim;i++)
            {
                sig12[i] /= sqrt(signorm) ;
                pscal += v12[i] * sig12[i] ;}
            for(i=0;i<dim;i++)
                align += (1- (v[sys][p1][i]+v[sys][p2][i])/2./normoy )* pscal*sig12[i] ;
            
            R1 = r[sys][p1][dim];
            R2 =r[sys][p2][dim];

            if(R1<=0) R1 += epsilon ;
            if(R2 <= 0) R2 += epsilon;

            
            // Energy lost during the collision: total and radial (to see if there is anisotropy)
            dE=0.; 
            dErad = 0. ; // Unused for now

            for(i=0;i<dim;i++)
            {
                dE+= (m[p1]*v[sys][p1][i]*v[sys][p1][i]+m[p2]*v[sys][p2][i] *v[sys][p2][i]) /2. ;
                //dErad += m[p1]*( v[sys][p1][i] * (r[sys][p1][i]-r0[sys][i]) /R1 )*( v[sys][p1][i] *(r[sys][p1][i]-r0[sys][i]) /R1 ) /2.;
                //dErad += m[p2]*( v[sys][p2][i] * (r[sys][p2][i]-r0[sys][i]) /R2 )*( v[sys][p2][i] *(r[sys][p2][i]-r0[sys][i]) /R2 ) /2.;
        
                if ( (blastcol[sys][p1]>0 || blastcol[sys][p2]>0) && sqrt(norm)>cutoff ) {
                    // at least one blast particle involved, and the relative velocity is above threshold
                    v[sys][p1][i] += (1. + alpha)*m[p2]/(m[p1]+m[p2]) * pscal *sig12[i];
                    v[sys][p2][i] -= (1. + alpha)*m[p1]/(m[p2]+m[p1]) * pscal *sig12[i] ;
                    blastcol[sys][p1]+=1;
                    blastcol[sys][p2]+=1;
                }
                else{
                    v[sys][p1][i] += 2*m[p2]/(m[p1]+m[p2]) * pscal *sig12[i];
                    v[sys][p2][i] -= 2*m[p1]/(m[p2]+m[p1]) * pscal *sig12[i] ; }

                    //dErad -= m[p1]*( v[sys][p1][i] * (r[sys][p1][i]-r0[sys][i]) /R1 )*( v[sys][p1][i] *(r[sys][p1][i]-r0[sys][i]) /R1 ) /2.;
                    //dErad -= m[p2]*( v[sys][p2][i] * (r[sys][p2][i]-r0[sys][i]) /R2 )*( v[sys][p2][i] *(r[sys][p2][i]-r0[sys][i]) /R2 ) /2.;
                    dE-= (m[p1]*v[sys][p1][i]*v[sys][p1][i]+m[p2]*v[sys][p2][i] *v[sys][p2][i]) /2. ;
            }

            a2[sys] = dim*dim/4.*mu4[sys]/mu2[sys]/mu2[sys] ;
            a2[sys] = a2[sys]*( 4./dim/(dim+2.) )-1. ;
            mom[sys] = 8*mu3[sys]- 12*mu1[sys];
            //printf("%g %g %g %F\n",a2[sys], E0[sys]);
            t0[sys] += tfirst[sys] ;

            tphys=t0[sys]/( sqrt(nu / E0[sys]) * sigma*omega[dim]/dim/rdens/khi) ;

            Nexp[sys]=0.;
            convrad(v[sys][p1],blank);
            convrad(v[sys][p2],blank);

           // ======================================== TEMPS AVANT LA PROCHAINE COLLISION ============
           
            tfirst[sys]=timbig;
            for(i=0;i<nu;i++)
            {       
                if (blastcol[sys][i]>0)  Nexp[sys]+=1;
            
                if(i == p1 || i == p2 || partner[sys][i] == p1 || partner[sys][i] == p2 )
                {

                    coltim[sys][i]=timbig ;
                    for(j=i+1;j<nu;j++)
                    {

                        tij = timeval(r[sys][i],r[sys][j],v[sys][i],v[sys][j], sigma);
                        if(tij>0.)
                        {
                            if (tij < coltim[sys][i])
                            {
                                coltim[sys][i] = tij ;
                                partner[sys][i] = j ;
                            }
                            if (tij < coltim[sys][j])
                            {
                                coltim[sys][j] = tij ;
                                partner[sys][j] = i ;
                            }
                        }



                    }

                }
                if(coltim[sys][i] < tfirst[sys] &&coltim[sys][i]<timbig/2. && coltim[sys][i]>0.)
                {
                    tfirst[sys] = coltim[sys][i];
                    pp1 = i;
                    pp2 = partner[sys][i];
                }
            }



            j=p1;
            for(i=0;i<p1;i++)
            {
                tij = timeval(r[sys][i],r[sys][j],v[sys][i],v[sys][j], sigma);

                if(tij>0.)
                {
                    if (tij < coltim[sys][i])
                    {
                        coltim[sys][i] = tij ;
                        partner[sys][i] = j ;
                    }
                    if (tij < coltim[sys][j])
                    {
                        coltim[sys][j] = tij ;
                        partner[sys][j] = i ;
                    }
                    if(tij < tfirst[sys])
                    {
                        tfirst[sys] = tij;
                        pp1 = i;
                        pp2 = j;
                    }
                }
                else
                {
                    if(tij>-1){
                    overlap+=1;
                    fout=fopen("results/overlap.dat", "a");
                    fprintf(fout,"%li %li %g\n",i, j, tij);
                    fclose(fout);}
                }
            }
            j=p2;
            for(i=0;i<p2;i++)
            {
                tij = timeval(r[sys][i],r[sys][j],v[sys][i],v[sys][j], sigma);
                if(tij>0.)
                {
                    if (tij < coltim[sys][i])
                    {
                        coltim[sys][i] = tij ;
                        partner[sys][i] = j ;
                    }
                    if (tij < coltim[sys][j])
                    {
                        coltim[sys][j] = tij ;
                        partner[sys][j] = i ;
                    }
                    if(tij < tfirst[sys])
                    {
                        tfirst[sys] = tij;
                        pp1 = i;
                        pp2 = j;
                    }
                }
            }

            
            for(k=0;k<dim;k++){
                if (v0[sys][k]>epsilon)
                    r0[sys][k]+=tfirst[sys] * v0[sys][k]/Nexp[sys];
                }

            if(alpha==1)
                r0calc(r[sys],v[sys],r0[sys]);

            nextp1[sys] = pp1;
            nextp2[sys] = pp2 ;

            //printf("sys %i tfirst %g (t = %g) %li %li %g %g\n",sys,tfirst[sys], t0[sys],p1,p2, v[sys][p1][dim],v[sys][p2][dim]);
            mom[sys]=8*mu3[sys] - 12 * mu1[sys];
            //Rexp[sys] /= Nexp[sys] ;
            vrad[sys] /= Nexp[sys];


            // QUANTITES COLLISIONNELLES ///

            collisions[sys][p1]+=1;
            collisions[sys][p2]+=1;

            dissip += dE ;
            dissipation[sys][p1] += dE/2.;
            dissipation[sys][p2] += dE/2.;
            tmoy += tfirst[sys] ;
            vrelative += pscal ;    


        }
        coll += 1 ;

        // -----------------------------------------OUTPUT --------------------------------------------------------

            
        // -----------------------------------------Outputs rapides (perode=evolpas) --------------------------------------------
        if ( ncoll <30000000 || coll%evolpas==0)
        {
        
            // On recalcule les quantites par systeme pour eviter les problemes de precision 
            // Peut-etre a enlever 
            // (Label01)
            for(sys=0;sys<nsys;sys++)
            {
                E[sys] = 0. ;
                Erad[sys] = 0. ;
                mu1[sys] = 0. ;
                mu2[sys] = 0. ;
                mu3[sys] = 0. ;
                mu4[sys] = 0. ;
                Nexp[sys] = 0 ;
                Rexp[sys]=0. ;
                vrad[sys]=0. ; 
                //if(coll%20==0)
                //    r0calc(r,v,r0);

                for(i=0;i<nu;i++)
                {


                    if (blastcol[sys][i]==0 && r[sys][i][dim]<Rmix) Rmix=r[sys][i][dim];
                    if( blastcol[sys][i]>0 ){
                        R = r[sys][i][dim];
                        for(k=0;k<dim;k++)
                        {    
                            E[sys]+= m[i]*v[sys][i][k]*v[sys][i][k]/2. ;
                            Ecm+= m[i]*v[sys][i][k];
                            vrad[sys] += fabs( v[sys][i][k] * (r[sys][i][k]-r0[sys][k]) /R );
                        }
                        mu1[sys] += v[sys][i][dobs]  /nu ;
                        mu2[sys] += v[sys][i][dobs]*v[sys][i][dobs] /nu ;
                        mu3[sys] += v[sys][i][dobs]*v[sys][i][dobs]*v[sys][i][dobs]/nu;
                        mu4[sys]+= v[sys][i][dobs]*v[sys][i][dobs]*v[sys][i][dobs]*v[sys][i][dobs]/nu ;

                        Nexp[sys] +=1 ;
                        if(R>Rexp[sys])
                            Rexp[sys]=R ;
                    for(k=0;k<dim;k++)
                    Erad[sys]+= m[i]*( v[sys][i][k] * (r[sys][i][k]-r0[sys][k]) /R )*( v[sys][i][k] *(r[sys][i][k]-r0[sys][k]) /R ) /2. ;
                }
            }
            
            
            // QUANTITES MOYENNES
        
            moyt=0 ;
            moytphys =0;
            moyRexp=0;
            moyNexp =0;
            moyE =0;
            moyErad =0;
            moymom =0;
            for(sys=0;sys<nsys;sys++)
            {
                moyt += t0[sys] ;
                moytphys += t0[sys]/( sqrt(nu / E0[sys]) * sigma*omega[dim]/dim/rdens/khi) ;
                moyRexp += Rexp[sys] ;
                moyNexp += Nexp[sys] ;
                moyE += E[sys] ;
                moyErad += Erad[sys] ;
                moymom += mom[sys];
            }
            moyt /=nsys ;
            moytphys /=nsys;
            moyRexp/=nsys;
            moyNexp /=nsys;
            moyE /=nsys;
            moyErad /=nsys;
            moymom /=nsys;
        
            if((fout=fopen("results/t.dat", "a"))==NULL) 
                                printf("Cannot open file.\n");
            Ecm = sqrt(Ecm*Ecm)/moyNexp/nsys ;    
            fprintf(fout,"%g %g %g %g %g\n",moyt,moyRexp,moyNexp,moyE, Ecm ); //mom[i], a2[i]
            Ecm = 0. ;
            fclose(fout);
            }

        }
        
        // -----------------------------------------Outputs lents (periode=printpas) --------------------------------------------
        if( (tstep <epsilon && ((coll >= colwait  && (coll-colwait) % printpas == 0 ) || coll == ncoll)) ||(tstep>0 && (moytphys-tlast > tstep) ))
        {
            tlast = moytphys ;
            if(stopr==0)
                printf("Affichage (%li overlaps)\n", overlap);
            nb = nb +1 ;
            //if (stopr==0)printf("Erad / E = %g \n",moyErad/moyE);


            // Collisional quantities
            fout=fopen("results/evT.dat", "a");
            fprintf(fout,"%g %g %g %g ",moyt, dissip/(moytphys-lasttime), vrelative/(moytphys-lasttime), tmoy/(moytphys-lasttime));
            fclose(fout);
            
            lasttime=moytphys;
             dissip=0; nactiv=0 ; vrelative=0 ; tmoy = 0;


            // POSITIONS ET VITESSES

            if(nsys<nslim){
                sprintf(affix,"%i", nb);
                sprintf(root, "results/r") ;
                strcat(root,affix) ;
                strcat(root,".dat") ;
                if((fout=fopen(root, "w"))==NULL) 
                   printf("Cannot open file.\n");
            }

            for (sys = 0; sys<nsys ;sys++)
            {

                fprintf(fout,"# %g  %g    %g    ",t0[sys],rdens,alpha);
                for(k=0; k<dim ; k++)
                    fprintf(fout,"%g    ",r0[sys][k]);
                fprintf(fout,"\n## ");
                for(k=0; k<dim ; k++)
                    fprintf(fout,"r    ");
                for(k=0; k<dim ; k++)
                    fprintf(fout,"v    ");
                fprintf(fout,"m sigma  blastcol  dissipation  angmov \n");
                for (i = 0; i<nu ;i++)
                {
                    if(nsys<nslim){
                        for(k=0; k<dim ; k++)
                            fprintf(fout,"%.24g    ",r[sys][i][k]);
                        for(k=0; k<dim ; k++)
                            fprintf(fout,"%g    ",v[sys][i][k]);
                        fprintf(fout,"%g %g %i    %g    %g   ",m[i],sigma,blastcol[sys][i], dissipation[sys][i],angmov[sys][i]);
                        //collisions[sys][i]=0. ; //set to 0 to count only collisions during interval
                    
                        fprintf(fout,"\n");
                    
                    }
       
                }
            
            }
        
            if(nsys<nslim) fclose(fout);

        }
    }

    time (&end);
    dif = difftime (end,start);
    if (stopr==0)printf("Temps écoulé : %g min\n", dif/60);

    time ( &rawtime );
    ptm = gmtime ( &rawtime );
    if (stopr==0)printf ("Heure de fin :  %2d:%02d \n", (ptm->tm_hour+2)%24, ptm->tm_min);
    
    for(sys=0 ; sys<nsys ;sys++)
    {
        for(j=0;j<nu;j++)
        {
            free(r[sys][j]);
            free(v[sys][j]);
        }
        free(r[sys]) ;
        free(v[sys]) ;
        free(r0[sys]) ;
        free(v0[sys])  ;
        free(coltim[sys]);
        free(tcounter[sys]);
        free(partner[sys]);
        free(dissipation[sys]);
        free(angmov[sys]);
        free(collisions[sys]);
    }
    return 0 ;
}


