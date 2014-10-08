import os
from matplotlib import pyplot as plt, colors, axes 
import numpy as np
from numpy import array
import pylab
import brewer2mpl
import re

files=[]

folder='./results/'

listdir= os.listdir(folder)
index={i:int(i.split('.')[0][1:] ) for i in listdir if re.findall("r.*?.dat",i) and i.split('.')[0][1:].isdigit() }


for i in sorted(index, key=lambda e:index[e] ):
    if i[0]=='r':
        files.append((i,open(folder+i,'r')))
    

for z,f in files:
    print z
    rf=[]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(aspect=1)
    sigma=[]
    for line in f:
        if not '#' in line:
            ll= [float(l)  for l in line.split()[:4] ]
            sigma.append(float(line.split()[5]))
            rf.append(ll) 
    ar = array(rf)
    x,y,vx,vy=ar.T[:4]
    sigma_in_points =  ax.transData.transform((0,sigma[0]) )-ax.transData.transform((0,0) )
    plt.scatter(x,y,alpha=.1,s=3.14*sigma_in_points**2)
    x,y,vx,vy=ar[ [i for i in range(ar.shape[0]) if ar[i].all() ]  ] .T[:4]


    plt.quiver(x,y,vx,vy )
    plt.savefig(folder+z.split('.')[0]+'.png' )
    plt.close()
    
