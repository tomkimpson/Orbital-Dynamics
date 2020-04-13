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



#fig, ax1,= plt.subplots(1, 1, figsize=(10,10))

fig = plt.figure(figsize=(10,15))
ax1 = plt.subplot2grid((3,1), (0,0), rowspan=2)
ax2 = plt.subplot2grid((3,1), (2,0))


#Load data
path = '/Users/tomkimpson/Data/ThesisData/MPD/'
MPDFile= path +'trajectory.txt'


def plot(f):
    data = np.loadtxt(f)

    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]

    xKM = data[:,4] / 1e3
    yKM = data[:,5] /1e3
    zKM = data[:,6] /1e3

    ax1.plot(x,y,c='C0')
    ax2.plot(x,zKM,c='C0')



def Format2D(ax):

    fs = 20

    #Label the axes
    ax.set_ylabel(r'$y [r_{\rm g}]$',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH
    ax.scatter(0,0,c='r')

    #axes limits
    sq = 25
    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)



plot(MPDFile)
Format2D(ax1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

ax2.set_xlabel(r'$x [r_{\rm g}]$',fontsize=fs)
ax2.set_ylabel(r'$z$ [km]',fontsize=fs)
ax2.locator_params(axis='both', nbins=5) #set number of xticks
ax2.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
plt.setp(ax1.get_xticklabels(),visible=False)

plt.subplots_adjust(hspace=0.01)


savefile = '/Users/tomkimpson/Data/ThesisData/orbitaldynamics2d.png'
plt.savefig(savefile,dpi=100,bbox='tight')
plt.show()

plt.show()
