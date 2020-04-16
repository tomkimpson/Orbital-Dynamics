from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D


#Setup plotting environment
plt.style.use('science')
fig, ax1 =plt.subplots(1, 1, sharex=True, figsize=(10,10))



path = '/Users/tomkimpson/Data/ThesisData/MPD/precession_time_delay/'



gamma = 10.0 * np.pi/180.0 
psi = np.pi/12


def cot(x):

    return np.cos(x) / np.sin(x)


def csc(x):
    return 1/np.sin(x)


def beamwidth(stheta):
    out = np.arccos((-np.cos(stheta)*cot(psi) + 1.4142135623730951*np.cos(gamma)*csc(psi) - cot(psi)*np.sin(stheta))/(np.cos(stheta) - np.sin(stheta)))
    return 2*out



def get_data(f):

    data = np.loadtxt(f)
    tau = data[:,6] #normalised w.r.t orbital period
 
    
    thetaSL = data[:,1]
    w = beamwidth(thetaSL)
    w = w / w[0]

    ax1.plot(tau,w)





eStrings = ['e06','e07','e08', 'e09']
#eStrings = ['e08','e09']
for e in eStrings:
    fname = path+e+'.txt'
    get_data(fname)



#Format
fs = 20



all_axes = plt.gcf().get_axes()
for ax in all_axes:
    ax.locator_params(axis='both', nbins=5) #set number of xticks
    ax.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers
    

ax1.set_xlabel(r'$ \tau / P $',fontsize=fs)
ax1.set_ylabel(r'$w/w_0$',fontsize=fs)

savefile = '/Users/tomkimpson/Data/ThesisData/beamwidth.png'
plt.savefig(savefile,dpi=100,bbox='tight')
plt.show()



