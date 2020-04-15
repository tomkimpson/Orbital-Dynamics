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

width = 30
height = width/2


fig = plt.figure(figsize=(width,height))
ax1 = plt.subplot2grid((3,3), (0,0), rowspan=2)
ax2 = plt.subplot2grid((3,3), (2,0))
ax3 = plt.subplot2grid((3,3), (0,1), rowspan=3,colspan=2)


#Load data
path = '/Users/tomkimpson/Data/ThesisData/MPD/'
f = path + 'trajectory.txt'


#some parameters
c = 3e8
ThetaObs = np.pi/4

def get_data(f):
    data = np.loadtxt(f)




    t = data[:,0] * 1e3 #ms
    x = data[:,1]
    y = data[:,2]
    z = data[:,6] #meters
    phi= data[:,-1] / (2*np.pi)


    #plot the trajectory
    ax1.plot(x,y)
    ax2.plot(x,z/1e3)

    #and the black hole
    ax1.scatter(0,0,c='r')


    deltaT = z*np.cos(ThetaObs) / c
    deltaT = deltaT  *1e6 #microseconds
    ax3.plot(phi,deltaT,c='C2')

    return t,z


get_data(f)

fs = 20
#Label the axes
plt.setp(ax1.get_xticklabels(),visible=False)

ax1.set_ylabel(r'$y [r_{\rm g}]$',fontsize=fs)
ax2.set_xlabel(r'$x [r_{\rm g}]$',fontsize=fs)
ax2.set_ylabel(r'$z$ [km]',fontsize=fs)


ax3.set_xlabel(r'$\phi / 2 \pi$',fontsize=fs)
ax3.set_ylabel(r'$\Delta t [\mu s]$',fontsize=fs)

plt.subplots_adjust(hspace=0.01)


ax1.locator_params(axis='both', nbins=5)
ax1.tick_params(axis='both', which='major', labelsize=fs-4)

ax2.locator_params(axis='both', nbins=5)
ax2.tick_params(axis='both', which='major', labelsize=fs-4)

ax3.locator_params(axis='both', nbins=10)
ax3.tick_params(axis='both', which='major', labelsize=fs-4)


savefile = '/Users/tomkimpson/Data/ThesisData/spin_curvature_time_delay.png'
plt.savefig(savefile,dpi=100,bbox='tight')
plt.show()

