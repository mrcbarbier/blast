import os
from matplotlib import pyplot as plt, colors, axes 
import numpy as np
from numpy import array
import pylab
import brewer2mpl
import re
from math import floor,cos,pi,sin
from datatools import *


nrad=80 #number of boxes in radial direction
nang=80 #number of angular slices
dim=2 #not finished in dim>2

rescale='xy'

folder='./results/'
indices=[7] #Indices of files to plot (list or range, empty=all)

#======================================================

listdir= os.listdir(folder)
index={i:int(i.split('.')[0][1:] ) for i in listdir if re.findall("r.*?.dat",i) and i.split('.')[0][1:].isdigit() }

npart=0
files=[]
for i in sorted(index, key=lambda e:index[e] ):
    if i[0]=='r' and (not indices or True in ['r{}.'.format(x) in i for x in indices]) :
        files.append((i,open(folder+i,'r')))
    


def dic_to_mat(dic):
    amax=[0]
    eshape=[]
    for i in dic:
        try:
            while len(amax)<len(i):
                amax.append(0)
            for idx,j in enumerate(i):
                if amax[idx]<=j:
                    amax[idx]=j+1
        except:
            amax[0]=max(amax[0],i+1)
        if hasattr(dic[i],'__iter__'):
            try:
                eshape=list( dic[i].shape)
            except:
                eshape=[len(dic[i]) ]
    mat=np.zeros(amax+eshape )
    for i in dic:
        try:
            j=list(i)
            elem=mat
            while len(j)>1:
                elem=elem[j[0]]
                j.pop(0)
            elem[j[0]]=dic[i]
        except:
            mat[i]=dic[i]
    return mat

def convrad(r, r0=None):
    '''Convert vector r to radial coordinates, using r0 as origin'''
    r=r.T
    specoord=np.zeros(r.shape)
    dim=len(r)
    if r0 is None:
        r0=np.zeros(r.shape)
    x=np.zeros(r.shape)
    for i in range(dim):
        x[i] = (r[i]-r0[i]) 
        specoord[0] += x[i]**2
    specoord[0]=np.sqrt(specoord[0])
    for i in range(1,dim-1):
        temp = 0. ;
        for k in range(dim-1,i-1,-1): 
            temp+=x[k]**2
        specoord[i] = np.arctan2( sqrt(temp) , x[i-1] )
    specoord[dim-1]= np.arctan2(x[dim-1],x[dim-2])
    return specoord.T

def batch_convrad(vec,*args):
    newvec=vec.copy()
    for i in range(newvec.shape[0]):
        for j in range(newvec.shape[1]):
            newvec[i,j]=convrad(vec[i,j],*args)
    return newvec

def extract(f):
    '''Extract data from file'''
    rv=[]
    m=[]
    sigma=[]
    blastcol=[]
    global npart
    for line in f:
        if not '#' in line:
            ls=[float(x) for x in line.split()]
            ll= [ls.pop(0)  for l in range(2*dim) ]
            m.append(ls.pop(0))
            sigma.append(ls.pop(0))
            blastcol.append(ls.pop(0))
            rv.append(ll) 
        else:
            tags=line[1:]
    rv=array(rv).T
    npart=max(npart,rv.shape[1] )
    r,v=rv[:dim].T,rv[dim:].T
    f.seek(0)
    return tags,r,v,m,sigma,blastcol
    
def treat(f):
    '''Extract hydro fields from file'''
    tags,r,v,m,sigma,blastcol=extract(f)
    radr=convrad(r)
    radv=convrad(v)
    angbox={}
    slices={}
    Rslice={}
    in_box={}
    boxes={}
    boxcenter={}
    radcenter={}
    area={}
    for n in range(npart):
        angbox[n]=floor((radr[n][1]/(2*pi)+.5)*nang)
        slices.setdefault(angbox[n],[]).append(n)
    for th in range(nang):
        Rmin=min(radr[n][0] if not blastcol[n] else 
            1 for n in slices[th] )
        Rmax= max(radr[n][0] if blastcol[n] else Rmin for n in slices[th] )
        Rslice[th]= [Rmin,Rmax,min(0.5,Rmax*1.4)]
    for n in range(npart):
        th= int(floor(radr[n][1]/2/pi*nang))
        if radr[n][1]<0:
            th+=nang
        box=floor(radr[n][0]/Rslice[th][2]*nrad)
        if box>=nrad:
            continue
        box=tuple(int(z) for z in (box,angbox[n]))
        in_box[n]=box
        boxes.setdefault(box,[]).append(n)
    rho,rho0,volfrac,u,E,Th={},{},{},{},{},{}
    for i in range(nrad):
        for j in range(nang):
            box=(i,j)
            boxcenter[box]=(box[0]+.5)*np.array((cos(box[1]*2*pi/nang),sin(box[1]*2*pi/nang)))/nrad *Rslice[box[1]][2]
            area[box]=pi* ( 2*box[0]+1)*Rslice[box[1]][2]**2 /nang/nrad**2
    for box in boxes:
        #radcenter[box]=np.array( ( (box[0]+.5)/nrad*Rslice[box[1]], box[1]*2*pi/nang ) )
        rho[box]=sum(m[p] for p in boxes[box])
        rho0[box]=sum(m[p] for p in boxes[box] if not blastcol[p])
        volfrac[box]=sum(sigma[p]**2*pi for p in boxes[box])/area[box]
        u[box]=sum(m[p]*v[p] for p in boxes[box])/rho[box]
        E[box]=sum(m[p]*np.dot(v[p],v[p])/2 for p in boxes[box])/rho[box]
        test=((m[p]*v[p]-u[box])/np.sqrt(2*m[p] ) for p in boxes[box])
        Th[box]=sum(np.dot(t,t) for t in test) /rho[box]
        rho[box]/=area[box]
        rho0[box]/=area[box]
    boxcenter,Rslice,rho,rho0,volfrac,u,E,Th=map(dic_to_mat, (boxcenter,Rslice,rho,rho0,volfrac,u,E,Th))
    return boxes,boxcenter,Rslice,rho,rho0,volfrac,u,E,Th

if __name__=='__main__':
        
    for z,f in files:
        print 'File:', z
        boxes,box,R,rho,rho0,phi,u,E,Th=treat(f)
        Rmin,Rmax,Rslice=np.mean(R,0)
        Rmin=min(R[:,0])
        urad=batch_convrad(u)
        boxrad=batch_convrad(box)
        r=np.mean(boxrad,axis=1)[:,0]
        if 'x' in rescale:
            r/=Rmin
        n=np.mean(phi,axis=1)
        n0=np.mean(rho0,axis=1)
        ur=np.mean(urad,axis=1)[:,0]
        temp=np.mean(Th,axis=1)
        if 'y' in rescale:
            rmax= int(round((Rmin/Rslice)*nrad))-1
            n/=n[rmax]
            ur/=ur[rmax]
            temp/=temp[rmax]
            n0/=max(n0)
        plot(r,n,xlabel='r',ylabel=r'$\phi$',hold=1 )
        plot(r,n0,xlabel='r',ylabel=r'$\rho_0$',hold=1 )
        plot(r,ur,xlabel='r',ylabel=r'$u$',hold=1)
        plot(r,temp,xlabel='r',ylabel=r'$\Theta$')
        fout=open("profiles.dat",'w')
        for i in range(boxrad.shape[0]):
            fout.write('{} {} {} {}\n'.format(r[i],n[i],ur[i],temp[i]) )
        fout.close()
