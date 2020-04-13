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
ax1 = plt.subplot2grid((1,1), (0,0),projection='3d')


#Load data
path = '/Users/tomkimpson/Data/ThesisData/MPD/'
MPDFile= path +'trajectory.txt'


def plot(f):
    data = np.loadtxt(f)

    t = data[:,0]
    x = data[:,1]
    y = data[:,2]
    z = data[:,3]


    ax1.plot(x,y,z)
    ax1.scatter(x[0],y[0],z[0])



def Format3D(ax):

    fs = 20

    #Label the axes
    ax.set_xlabel(r'$x [r_g]$',fontsize=fs)
    ax.set_ylabel(r'$y [r_g]$',fontsize=fs)

    ax.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers



    #Draw the BH
    ax.scatter(0,0,c='r')

    #axes limits
    sq = 20
    ax.set_xlim(-sq,sq)
    ax.set_ylim(-sq,sq)
    ax.set_zlim(-sq,sq)

    ax.set_axis_off()

plot(MPDFile)
Format3D(ax1)
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
fs = 20

savefile = '/Users/tomkimpson/Data/ThesisData/orbitaldynamics3d.png'
plt.savefig(savefile,dpi=100,bbox='tight')
plt.show()





plt.show()
