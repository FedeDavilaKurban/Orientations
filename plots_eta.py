#%%
import sys
import numpy as np
from scipy import spatial
import scipy.stats
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import random
from config import writePath, units
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
#%%

minradV = 7.
maxradV = 0.

lowMcut = -1


plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15


for sec in [3]:

    fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=False)

    if sec==1: title='High Spin - Low Mass Galaxies'
    elif sec==3: title='High Spin - High Mass Galaxies'
    elif sec==4: title='Low Spin - Low Mass Galaxies'
    elif sec==6: title='Low Spin - High Mass Galaxies'
    elif sec==123: title='High Spin Galaxies'
    elif sec==456: title='Low Spin Galaxies'
    elif sec==0: title='All Galaxies'
    else: title=sec

    axs[0].set_title(title)

    for ax, vtype in zip(axs,['a','r','s']):

        etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,sec,vtype),\
                names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax'])
                
        y = etaTable['eta'].data
        yerr = etaTable['eta_std'].data

        yran_mean = etaTable['eta_random_mean'].data
        yran_err = etaTable['eta_random_std'].data

        x = (etaTable['rmin'].data+etaTable['rmax'].data)/2
        
        # Theoretical n1/n2 value for random spins
        ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
        ax.plot(x,yran_mean,c='k',alpha=.7)

        ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.3, color='k')

        ax.errorbar(x,y,yerr=yerr)

        #ax.legend()
        ax.set_ylabel(r'$\eta$')
        if sec==3:  
            ax.set_ylim([1.5,3.5])
        else: ax.set_ylim([2.,2.8])

        # plt.ylabel(r'$n_1/n_2$')
        # plt.xlabel('R/Rv')
        # plt.legend()

    # x_ticks_labels = []
    # for i in range(nbins):
    #     x_ticks_labels.append( '{:.1f}-{:.1f}'.format(r1[i],r2[i]) )
    axs[2].set_xticks(np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]))
    #axs[2].set_xticklabels(x_ticks_labels)
    axs[2].set_xlabel('R/Rv')

    axs[0].text(1.25,axs[0].get_ylim()[1]*.9,'All Voids')
    axs[1].text(1.25,axs[0].get_ylim()[1]*.9,'Rising Voids')
    axs[2].text(1.25,axs[0].get_ylim()[1]*.9,'Shell Voids')


    plt.savefig('../plots/eta_vs_r/Eta_allvtypes_sec{}.jpg'.format(sec))

# %%
#%%
"""
2nd version -- No "all voids"
"""

minradV = 7.
maxradV = 0.

lowMcut = -1


plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15


for sec in [456]:

    fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=False)

    if sec==1: title='High Spin - Low Mass Galaxies'
    elif sec==3: title='High Spin - High Mass Galaxies'
    elif sec==4: title='Low Spin - Low Mass Galaxies'
    elif sec==6: title='Low Spin - High Mass Galaxies'
    elif sec==123: title='High Spin Galaxies'
    elif sec==456: title='Low Spin Galaxies'
    elif sec==0: title='All Galaxies'
    else: title=sec

    axs[0].set_title(title)

    for ax, vtype in zip(axs,['r','s']):

        etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,sec,vtype),\
                names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax'])
                
        y = etaTable['eta'].data
        yerr = etaTable['eta_std'].data

        yran_mean = etaTable['eta_random_mean'].data
        yran_err = etaTable['eta_random_std'].data

        x = (etaTable['rmin'].data+etaTable['rmax'].data)/2
        
        # Theoretical n1/n2 value for random spins
        ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
        ax.plot(x,yran_mean,c='k',alpha=.7)

        ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.15, color='k')
        ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.15, color='k')
        ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.15, color='k')

        ax.errorbar(x,y,yerr=yerr)

        #ax.legend()
        ax.set_ylabel(r'$\eta$')
        if sec==3:  
            ax.set_ylim([1.5,3.5])
        else: ax.set_ylim([2.,2.8])

        # plt.ylabel(r'$n_1/n_2$')
        # plt.xlabel('R/Rv')
        # plt.legend()


    axs[1].set_xticks(np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]))
    axs[1].set_xlabel('R/Rv')

    #axs[0].text(1.25,axs[0].get_ylim()[1]*.9,'All Voids')
    axs[0].text(1.25,axs[0].get_ylim()[1]*.9,'Rising Voids')
    axs[1].text(1.25,axs[0].get_ylim()[1]*.9,'Shell Voids')


    plt.savefig('../plots/eta_vs_r/Eta_sec{}.jpg'.format(sec))
# %%
