import os
from matplotlib import pyplot as plt, colors, axes 
import numpy as np
from numpy import array
import pylab
import brewer2mpl
import re
from math import hypot
from subprocess import call

from plot_profiles import treat

folder='./results/'
indices=[]
            #indices of the files to plot (empty=all)
save=1
             #Save plot to file
imgformat='png'
movie=0
             #Make a movie from the pictures


### AESTHETIC DETAILS
fields={
    'show':0,
    'smooth':1.4
    }

particles={
    'show':True,
    'cmap':plt.cm.afmhot,# plt.cm.Greys,#plt.cm.hot
    'arrows':[(0.00032,2.6),(0.0042,1.6)][1],#
    'alpha':.5,
    'restcolor':(.8,.3,.3),#(1.,0.5,0.5),(1.,1.,1.),#
    }

hold=0
background=True #Fill background
bgcmap=plt.cm.afmhot

mycmap = plt.cm.get_cmap('seismic')

### WINDOW FOR PLOTTING

dftwindow=np.array([[-.5,.5],[-.5,.5]])  #Default window (whole system)
window=dftwindow[:]
#window='rescale' #rescales the window to zoom out as blast expands

remove_overlaps=0  #Remove overlapping particles (very slow if window is big)


listdir= os.listdir(folder)
index={i:int(i.split('.')[0][1:] ) for i in listdir if re.findall("r.*?.dat",i) and i.split('.')[0][1:].isdigit() }



files=[]
for i in sorted(index, key=lambda e:index[e] ):
    if i[0]=='r' and (not indices or True in ['r{}.'.format(x) in i for x in indices]) :
        files.append((i,open(folder+i,'r')))

if hold:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set(aspect=1)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)

idx=0
for z,f in files:
    print z
    rf=[]
    
    boxes,box,R,rho,rho0,phi,u,E,Th=treat(f)
    
    if not hold:
        plt.clf()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set(aspect=1)
        plt.setp(ax.get_xticklabels(), visible=False)
        plt.setp(ax.get_yticklabels(), visible=False)
        
    plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)

    if background:
        window=np.array(window)
        size=100
        a,b=np.linspace(window[0,0],window[0,1],size),np.linspace(window[1,1],window[1,0],size).reshape(1,size).T
        bg=np.zeros((2,size,size))
        bg[0,:,:]+=a
        bg[1,:,:]+=b

        #import matplotlib.colors as pltcolors
        #cmblue = pltcolors.LinearSegmentedColormap('BLUE', cdict)
        
        decal=.1,0#.5#.5,.5
        
        mp=decal[0]+np.exp(-np.sum(np.abs(bg)**3,axis=0))
        mp[-1,-1]=np.max(mp)+decal[1]
        mp[0,0]=0
        mp=mp[-1,-1]-mp
        plt.imshow(mp ,interpolation='bicubic',extent=window.ravel(),cmap=bgcmap)

        
    plt.xlim(window[0])
    plt.ylim(window[1])
    
    if fields['smooth']:
        smooth=fields['smooth']
        import scipy.ndimage as ndimage
        for field in [rho,rho0,phi,E,Th]:
            field[:]=ndimage.gaussian_filter(field, sigma=(smooth, smooth), order=0)
        
    
    if fields['show']:
        field=Th
        
        def circularize(arr):
            #Add initial point at end of array so that pcolor makes a full circle
            return np.concatenate((arr,arr[:,:1]),1)
        box=-circularize(box)
        field=circularize(field)
        plt.pcolor(box[:,:,0],box[:,:,1],field,cmap=mycmap)
        plt.savefig(folder+'f'+z.split('.')[0]+'.{}'.format(imgformat) ,frameon=False,dpi=350,tight_layout=True)
        fig = plt.figure(frameon=False)
        ax = fig.add_subplot(111)
        ax.set(aspect=1)
        plt.xlim(window[0])
        plt.ylim(window[1])
        
    if particles['show']:
            

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
        if window=='rescale':
            xb=x[vx!=0 or vy!=0]
            yb=y[vx!=0 or vy!=0]
            window=[np.min(xb),np.max(xb)],[np.min(yb),np.max(yb)]
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
        #colors=np.sqrt(vx**2+vy**2)
        #colors/=np.max(colors)
        colors=[.5 for particle in range(len(x))]
        colors = np.zeros(len(x))
        for b in boxes:
            for part in boxes[b]:
                colors[part]=(vx[part]**2 + vy[part]**2)
                colors[part]=.3+E[b]/np.max(E)
        
        alp=particles['alpha']
        nonzero=(vx!=0) + (vy!=0)
        if hold:
            alp/=len(files)
        plt.scatter(x,y,s=3.14*min(sigma_in_points)**2,color=[particles['restcolor']  if ix>0 else (0,0,0) for ix in range(x.shape[0])],alpha=.2,cmap=plt.cm.bone)
    
        if hold:
            alp=alp*(idx+1)
        x,y,vx,vy=x[nonzero],y[nonzero],vx[nonzero],vy[nonzero]
        colors=colors[nonzero]
        colors[0]=0
        vx*=10
        vy*=10
        plt.quiver(x,y,vx,vy,colors,units='x',width=particles['arrows'][0], scale=particles['arrows'][1], alpha=alp,cmap=particles['cmap'])
        #img2 = np.frombuffer(fig.canvas.buffer_rgba(), np.uint8).reshape(h, w, -1).copy()
    if not hold:
        if save:
            plt.savefig(folder+z.split('.')[0]+'.{}'.format(imgformat) ,frameon=False,dpi=350,tight_layout=True)
        else:
            plt.show()
        plt.close()
    idx+=1

if hold:   
        if save:
            plt.savefig(folder+'held'+'.{}'.format(imgformat) ,frameon=False,dpi=350,tight_layout=True)
        else:
            plt.show()
if movie:
    try:
        call("ffmpeg -r 18  -i ./results/r%d.{} ./results/movie.mp4".format(imgformat),shell=True)
    except:
        print("Please install FFmpeg (or copy FFmpeg.exe in this folder).")
