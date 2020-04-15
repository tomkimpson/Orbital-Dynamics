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
MPDFile2= path +'trajectory2.txt'


def plot(f,ID):
    data = np.loadtxt(f)

    t = data[:,0] * 1e3 #ms
    pitch_coordinate = data[:,13]
    pitch_tetrad = data[:,14] 

    if ID == 1:
        ls = 'solid'
        c1 = 'C0'
        c2 = 'C1'

    else:
        ls = '--'
        c1 = 'C2'
        c2 = 'C3'




    ax1.plot(t,pitch_coordinate,linestyle='solid',c=c1)
    ax1.plot(t,pitch_tetrad,linestyle='--',c=c1)



plot(MPDFile,1)
plot(MPDFile2,2)

fs = 20
#Label the axes
ax1.set_xlabel(r'$t$ [ms]',fontsize=fs)
ax1.set_ylabel(r'$\tilde{\omega}$ [rad]',fontsize=fs)


ax1.locator_params(axis='both', nbins=5) #set number of xticks                                                                                                                                                                           
ax1.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers




savefile = '/Users/tomkimpson/Data/ThesisData/pitch_angle.png'
plt.savefig(savefile,dpi=100,bbox='tight')
plt.show()

