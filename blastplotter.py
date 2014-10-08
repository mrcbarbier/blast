import os
from matplotlib import pyplot as plt, colors, axes 
import numpy as np
from numpy import array
import pylab
import brewer2mpl
import re
from math import floor,cos,pi,sin


nrad=20
nang=20
dim=2

npart=0
files=[]

folder='../newsimu/'

listdir= os.listdir(folder)
index={i:int(i.split('.')[0][1:] ) for i in listdir if re.findall("r.*?.dat",i) and i.split('.')[0][1:].isdigit() }


for i in sorted(index, key=lambda e:index[e] ):
    if i[0]=='r':
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
    
def extract(f):
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
    tags,r,v,m,sigma,blastcol=extract(f)
    radr=convrad(r)
    radv=convrad(v)
    angbox={}
    slices={}
    Rslice={}
    in_box={}
    boxes={}
    boxcenter={}
    for n in range(npart):
        angbox[n]=floor((radr[n][1]/(2*pi)+.5)*nang)
        slices.setdefault(angbox[n],[]).append(n)
    for th in range(nang):
        Rslice[th]= (min(radr[n][0] if not blastcol[n] else 1 for n in slices[th] )+ max(radr[n][0] if blastcol[n] else 0 for n in slices[th] ))/2
    for n in range(npart):
        box=floor(radr[n][0]/Rslice[th]*nrad)
        if box>=nrad:
            continue
        box=tuple(int(z) for z in (box,angbox[n]))
        in_box[n]=box
        boxes.setdefault(box,[]).append(n)
    rho,u,E,Th={},{},{},{}
    for box in boxes:
        
        boxcenter[box]=(box[0]+.5)*np.array((cos(box[1]*2*pi/nang),sin(box[1]*2*pi/nang)))/nrad
        rho[box]=sum(m[p] for p in boxes[box])
        u[box]=sum(m[p]*v[p] for p in boxes[box])/rho[box]
        E[box]=sum(m[p]*np.dot(v[p],v[p])/2 for p in boxes[box])/rho[box]
        test=((m[p]*v[p]-u[box])/np.sqrt(2*m[p] ) for p in boxes[box])
        Th[box]=sum(np.dot(t,t) for t in test) 
    boxcenter,Rslice,rho,u,E,Th=map(dic_to_mat, (boxcenter,Rslice,rho,u,E,Th))
    return boxcenter,Rslice,rho,u,E,Th

    
for z,f in files:
    print z
    box,R,rho,u,E,Th=treat(f)
    plt.plot(np.mean(Th,axis=1))
    plt.show()
