from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os



##Setup plotting environment
plt.style.use('science')


fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1), (0,0))


#Load data
path = '/Users/tomkimpson/Data/ThesisData/MPD/'
MPDFile= path +'trajectory.txt'


def plot(f):
    data = np.loadtxt(f)

    t = data[:,0]
    x = data[:,4] /1e3
    y = data[:,5] /1e3
    z = data[:,6] /1e3

    
    xB = data[:,10] /1e3
    yB = data[:,11] /1e3
    zB = data[:,12] /1e3


    x0 = x[0]
    y0 = y[0]
    z0 = z[0]

    x = x-x0
    y = y-y0
    z = z-z0


    xB = xB-x0
    yB = yB-y0
    zB = zB-z0


    ax1.plot(x,y,c='C0')
    ax1.plot(xB,yB,c='C2', linestyle='--')


  #  for i in range(len(x)):
   #     if (i % 2) == 0:
    #        xx = [x[i],xB[i]]
     #       yy = [y[i],yB[i]]
      #      ax1.plot(xx,yy,c='C1')


plot(MPDFile)

fs = 20
#Label the axes
ax1.set_xlabel(r'$x [r_{\rm g}]$',fontsize=fs)
ax1.set_ylabel(r'$y [r_{\rm g}]$',fontsize=fs)


ax1.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
ax1.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers




savefile = '/Users/tomkimpson/Data/ThesisData/RVM.png'
plt.savefig(savefile,dpi=100,bbox='tight')
plt.show()

