from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
import os
from scipy import interpolate


##Setup plotting environment
plt.style.use('science')


fig = plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((2,1), (0,0))
ax2 = plt.subplot2grid((2,1), (1,0),sharex=ax1)


#Load data
path = '/Users/tomkimpson/Data/ThesisData/MPD/beamwidth/'
files = glob.glob(path+'*.txt')


t0 = 0.7
t1 = 1.4
gamma = 20.0
global_x = np.linspace(t0,t1,1000)



def plot(f,c):
    data = np.loadtxt(f)

    t = data[:,0] * 1e3 #ms
    pitch_tetrad = data[:,14]  * 180 / np.pi
 #   ax1.plot(t,pitch_tetrad)

    func = interpolate.interp1d(t, pitch_tetrad)
    ynew = func(global_x)
    

    ax1.plot(global_x,ynew,c=c)

    #Get the minimum
    idx = np.argmin(ynew)
    ax1.scatter(global_x[idx], ynew[idx],c=c, marker = 'x')
    ax1.axvline(global_x[idx],c=c, linestyle='--')
    ax2.axvline(global_x[idx],c=c, linestyle='--')

    print (f)
    print ('Minimum time = ', global_x[idx])


    #Get the lower gamma crossing
    xx = global_x[:idx]
    yy = ynew[:idx]

    func = interpolate.interp1d(yy,xx) #given yy, what is xx?
    xLower = func(gamma)

    ax1.scatter(xLower,gamma,c=c)

    #Get the upper gamma crossing
    xx = global_x[idx:]
    yy = ynew[idx:]

    func = interpolate.interp1d(yy,xx) #given yy, what is xx?
    xUpper = func(gamma)

    ax1.scatter(xUpper,gamma,c=c)


    tree_x = [xLower, xLower]
    tree_y = [0,1]
    ax2.plot(tree_x,tree_y,c=c)
    ax2.scatter(tree_x[-1], tree_y[-1], c=c)


    
    tree_x = [xUpper, xUpper]
    ax2.plot(tree_x,tree_y,c=c)
    ax2.scatter(tree_x[-1], tree_y[-1], c=c)



i = 0
for f in files:
    col = 'C'+str(i)
    plot(f,col)
    i += 1
fs = 20
#Label the axes
ax1.set_ylabel(r'$\tilde{\omega}$ [deg]',fontsize=fs)

axes = [ax1,ax2]
for ax in axes:
    ax.locator_params(axis='both', nbins=5)
    ax.tick_params(axis='both', which='major', labelsize=fs-4)
    ax.set_xlabel(r'$t$ [ms]',fontsize=fs)



plt.subplots_adjust(hspace=0.01)
plt.setp(ax1.get_xticklabels(),visible=False)
plt.setp(ax2.get_yticklabels(),visible=False)


ax1.set_xlim(0.8,1.3)
ax1.set_ylim(0,30)
ax2.set_ylim(0,1.1)

ax1.axhline(gamma, linestyle='--', c='0.5')


ax1.set_xlabel(r'$t$ [ms]',fontsize=fs)

#labels = [item.get_text() for item in ax1.get_xticklabels()]
#labels[1] = 'Testing'
#ax1.set_xticklabels(labels)



savefile = '/Users/tomkimpson/Data/ThesisData/beamwidth_numerical.png'
plt.savefig(savefile,dpi=100,bbox='tight')

plt.show()

