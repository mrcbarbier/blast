import os
from matplotlib import pyplot as plt, colors, axes 
import numpy as np
from numpy import array
import pylab
import brewer2mpl
import re
from math import floor,cos,pi,sin
from datatools import *


nrad=100 #number of boxes in radial direction
nang=20 #number of angular slices
dim=2 #not finished in dim>2


folder='./results/'
indices=[] #Indices of files to plot (list or range, empty=all)

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
        Rslice[th]= min(0.5,(min(radr[n][0] if not blastcol[n] else 
            1 for n in slices[th] )+ max(radr[n][0] if blastcol[n] else 0 for n in slices[th] ))/2 *1.4)
    for n in range(npart):
        th= int(floor(radr[n][1]/2/pi*nang))
        if radr[n][1]<0:
            th+=nang
        box=floor(radr[n][0]/Rslice[th]*nrad)
        if box>=nrad:
            continue
        box=tuple(int(z) for z in (box,angbox[n]))
        in_box[n]=box
        boxes.setdefault(box,[]).append(n)
    rho,volfrac,u,E,Th={},{},{},{},{}
    for i in range(nrad):
        for j in range(nang):
            box=(i,j)
            boxcenter[box]=(box[0]+.5)*np.array((cos(box[1]*2*pi/nang),sin(box[1]*2*pi/nang)))/nrad *Rslice[box[1]]
            area[box]=pi* ( 2*box[0]+1)*Rslice[box[1]]**2 /nang/nrad**2
    for box in boxes:
        #radcenter[box]=np.array( ( (box[0]+.5)/nrad*Rslice[box[1]], box[1]*2*pi/nang ) )
        rho[box]=sum(m[p] for p in boxes[box])
        volfrac[box]=sum(sigma[p]**2*pi for p in boxes[box])/area[box]
        u[box]=sum(m[p]*v[p] for p in boxes[box])/rho[box]
        E[box]=sum(m[p]*np.dot(v[p],v[p])/2 for p in boxes[box])/rho[box]
        rho[box]/=area[box]
        test=((m[p]*v[p]-u[box])/np.sqrt(2*m[p] ) for p in boxes[box])
        Th[box]=sum(np.dot(t,t) for t in test) 
    boxcenter,Rslice,rho,volfrac,u,E,Th=map(dic_to_mat, (boxcenter,Rslice,rho,volfrac,u,E,Th))
    return boxcenter,Rslice,rho,volfrac,u,E,Th

    
for z,f in files:
    print z
    box,R,rho,phi,u,E,Th=treat(f)
    urad=batch_convrad(u)
    boxrad=batch_convrad(box)
    plot(np.mean(boxrad,axis=1)[:,0],np.mean(phi,axis=1),xlabel='r',ylabel=r'$\phi$' )
    plot(np.mean(boxrad,axis=1)[:,0],np.mean(urad,axis=1)[:,0] ,xlabel='r',ylabel=r'$u$')
    plot(np.mean(boxrad,axis=1)[:,0],np.mean(Th,axis=1),xlabel='r',ylabel=r'$\Theta$')

