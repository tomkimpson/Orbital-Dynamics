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

left, bottom, width, height = [0.2, 0.65, 0.2, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])


path = '/Users/tomkimpson/Data/ThesisData/MPD/precession_time_delay/'


def get_data(f):

    data = np.loadtxt(f)
    tau = data[:,6] #normalised w.r.t orbital period
    Xc = data[:,5]


    #Calculate the time delay
    delta_Xc = Xc-Xc[0]

    Ps = 1e-3 #seconds
    delta_t = Ps *delta_Xc / (2*np.pi)
    
    delta_t = delta_t * 1e6 #microseconds

    ax1.plot(tau,delta_t)
    ax2.plot(tau,delta_t)





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
    

#configure inset plot
ax2.set_xlim(2.05,2.25)
ax2.set_ylim(105,150)
plt.setp(ax2.get_xticklabels(),visible=False)
plt.setp(ax2.get_yticklabels(),visible=False)
ax2.locator_params(axis='both', nbins=5) #set number of xticks
ax2.tick_params(axis='both', which='major', labelsize=fs-4) #set size of numbers


ax1.set_xlabel(r'$ \tau / P $',fontsize=fs)
ax1.set_ylabel(r'$\Delta t [\mu s]$',fontsize=fs)

savefile = '/Users/tomkimpson/Data/ThesisData/spin_time_delay.png'
plt.savefig(savefile,dpi=100,bbox='tight')
plt.show()



