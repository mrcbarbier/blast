import os
from matplotlib import pyplot as plt, colors, axes 
import numpy as np
from numpy import array
import pylab
import brewer2mpl
import re
from math import hypot
from subprocess import call


folder='./results/'
indices=[30]        #indices of the files to plot (empty=all)
save=0              #Save plot to file
imgformat='jpg'
movie=0             #Make a movie from the pictures

window=[-.5,.5],[-.5,.5]  #Default window
#window=[-0.06,0.0],[-0.26,-0.18] 
#window=[0.166,0.246],[0.015,0.035]

remove_overlaps=0  #Remove overlapping particles (very slow if window is big)


listdir= os.listdir(folder)
index={i:int(i.split('.')[0][1:] ) for i in listdir if re.findall("r.*?.dat",i) and i.split('.')[0][1:].isdigit() }


files=[]
for i in sorted(index, key=lambda e:index[e] ):
    if i[0]=='r' and (not indices or True in ['r{}.'.format(x) in i for x in indices]) :
        files.append((i,open(folder+i,'r')))

for z,f in files:
    print z
    rf=[]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set(aspect=1)

    plt.xlim(window[0])
    plt.ylim(window[1])

    sigma=[]
    blastcols=[]
    nb=0
    for line in f:
        if not '#' in line:
            LINE=line.split()
            ll= [float(l)  for l in LINE[:4] ]
            sigma.append(float(LINE[5]))
            rf.append(ll) 
            blastcols.append(int(LINE[6]))
            nb+=1
    ar = array(rf)
    x,y,vx,vy=ar.T[:4]
    sigma=np.array(sigma)
    if remove_overlaps: #Remove overlaps
        print("Removing overlaps")
        nonover=[]
        blastpart=[]
        for i in range(x.shape[0]):
            if (window[0][0]-sigma[i]< x[i] < window[0][1]+sigma[i] and window[1][0]-sigma[i]< y[i] < window[1][1]+sigma[i] ):
                nonover.append(i)
                if blastcols[i]>0:
                    blastpart.append(i)
        non=tuple(nonover)
        for idx,i in enumerate(non):
            for j in non[idx+1:]:
                if hypot(x[i]-x[j],y[i]-y[j])<(sigma[i]+sigma[j])*0.9:
                    try:
                        nonover.remove(i)
                    except:
                        pass
#        print("Removed {} overlaps.".format(len(overlaps)))
        x=x[nonover]
        y=y[nonover]
        vx=vx[nonover]
        vy=vy[nonover]
        sigma=sigma[nonover]
                
    sigma_in_points =  ax.transData.transform((sigma[0],sigma[0]) )-ax.transData.transform((0,0) )
    plt.scatter(x,y,alpha=.1,s=3.14*min(sigma_in_points)**2)
    
    x,y,vx,vy=x[vx!=0 or vy!=0],y[vx!=0 or vy!=0],vx[vx!=0 or vy!=0],vy[vx!=0 or vy!=0]
    plt.quiver(x,y,vx,vy ,units='x',width=0.00042, scale=.6, alpha=0.8)
    if save:
        plt.savefig(folder+z.split('.')[0]+'.{}'.format(imgformat) )
    else:
        plt.show()
    plt.close()
    
if movie:
    try:
        call("ffmpeg -r 18  -i ./results/r%d.{} ./results/movie.mp4".format(imgformat),shell=True)
    except:
        print("Please install FFmpeg (or copy FFmpeg.exe in this folder).")
