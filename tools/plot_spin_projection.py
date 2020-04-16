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

width = 20
height = 10


fig = plt.figure(figsize=(width,height))
ax1 = plt.subplot2grid((1,2), (0,0))
ax2 = plt.subplot2grid((1,2), (0,1),sharey=ax1)


#Load data
path = '/Users/tomkimpson/Data/ThesisData/MPD/spin_projection/'

files = glob.glob(path+'spin/*.txt')



#f = path + 'spin.txt'
#files = [f]


def get_data(f,ax):
    data = np.loadtxt(f)

    print (f)


    tau = data[:,0]
    thetaSL = data[:,1]
    phiSL = data[:,2]
    sx = data[:,3] 
    sy = data[:,4] 


    sx = sx / sx[0]
    sy = sy / sy[0]

    ax.plot(sx,sy,c='C2')



axes = plt.gcf().get_axes()
i = 0

efiles = ['e01', 'e09']

#for f in files:


for e in efiles:
    fname = path+'spin/'+e+'.txt'
    get_data(fname,axes[i])
    i +=1


#formatting

#Fontsize
fs = 20

#adjust wspace between plots
plt.subplots_adjust(wspace=0.01)


#remove some labels
plt.setp(ax2.get_yticklabels(),visible=False)

for ax in axes:
    ax.locator_params(axis='both', nbins=5)
    ax.tick_params(axis='both', which='major', labelsize=fs-4)
    ax.set_xlabel(r'$S^x$',fontsize=fs)



ax1.set_ylabel(r'$S^y$',fontsize=fs)



savefile = '/Users/tomkimpson/Data/ThesisData/spin_projection.png'
plt.savefig(savefile,dpi=100,bbox='tight')





plt.show()

