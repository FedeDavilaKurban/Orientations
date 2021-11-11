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
#%%

minradV = 7.
maxradV = 0.

lowMcut = -1


plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15


for sec in [14,36]:

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


for sec in [14,36]:

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
"""
3rd version -- (eta-eta0)/sigma(eta) vtype r y s
"""
minradV = 7.
maxradV = 0.

lowMcut = -1


plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15

fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=False)

for sec in [14,36]:

    if sec==14: label = "Low Mass"
    if sec==36: label = "High Mass"

    if sec==1: title='High Spin - Low Mass Galaxies'
    elif sec==3: title='High Spin - High Mass Galaxies'
    elif sec==4: title='Low Spin - Low Mass Galaxies'
    elif sec==6: title='Low Spin - High Mass Galaxies'
    elif sec==123: title='High Spin Galaxies'
    elif sec==456: title='Low Spin Galaxies'
    elif sec==0: title='All Galaxies'
    else: title=sec

    #axs[0].set_title(title)

    for ax, vtype in zip(axs,['r','s']):

        etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,sec,vtype),\
                names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax'])
                
        yran_mean = etaTable['eta_random_mean'].data
        yran_err = etaTable['eta_random_std'].data

        y = (etaTable['eta'].data-yran_mean)/yran_err
        yerr = etaTable['eta_std'].data/yran_err

        x = (etaTable['rmin'].data+etaTable['rmax'].data)/2
        
        # Theoretical n1/n2 value for random spins
        ax.hlines(0,x[0],x[-1],linestyles=':')
        #ax.plot(x,yran_mean,c='k',alpha=.7)

        #ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.15, color='k')
        #ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.15, color='k')
        #ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.15, color='k')

        ax.errorbar(x,y,yerr=yerr,label=label,fmt='o-',capsize=3,ms=5)

        #ax.legend()
        ax.set_ylabel(r'$\eta$')
        #if sec==3:  
        #    ax.set_ylim([1.5,3.5])
        #else: ax.set_ylim([2.,2.8])

        # plt.ylabel(r'$n_1/n_2$')
        # plt.xlabel('R/Rv')
        # plt.legend()


axs[1].set_xticks(np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]))
axs[1].set_xlabel('R/Rv')

axs[0].text(1.25,axs[0].get_ylim()[1]*.9,'Rising Voids')
axs[1].text(1.25,axs[0].get_ylim()[1]*.9,'Shell Voids')

axs[0].legend(loc='upper left')

plt.savefig('../plots/eta_vs_r/Eta_mass.jpg'.format(sec))

# %%
"""
3rd version -- (eta-eta0)/sigma(eta) vtype all
"""
import seaborn as sns
colors = sns.color_palette()
blue, oran = colors[0], colors[1]

minradV = 7.
maxradV = 0.
vtype = 'a'

lowMcut = -1


plt.rcParams['figure.figsize'] = (9, 6)
plt.rcParams['font.size'] = 15

fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=False)

#plt.fill_between(x, -1, 1, alpha=.025, color='k')
#plt.fill_between(x, -3, 3, alpha=.03, color='k')

plt.text(1.2,4.45,'All Voids')

for sec in [1,2,3,4,5,6]:

    if sec==14: label = "Low Mass"
    if sec==36: label = "High Mass"
    if sec==123: label = 'High Spin'
    if sec==456: label = 'Low Spin'

    if sec==1: 
        label = 'H-L Gxs'
        c = blue
        fmt = 'x:'
    if sec==2: 
        label = 'H-I Gxs'
        c = blue
        fmt = '^--'
    if sec==3: 
        label = 'H-H Gxs'
        c = blue
        fmt = 'o-'
    if sec==4: 
        label = 'L-L Gxs'
        c = oran
        fmt = 'x:'
    if sec==5: 
        label = 'L-I Gxs'
        c = oran
        fmt = '^--'
    if sec==6: 
        label = 'L-H Gxs'
        c = oran
        fmt = 'o-'


    etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype),\
            names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax'])
            
    yran_mean = etaTable['eta_random_mean'].data
    yran_err = etaTable['eta_random_std'].data

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    x = (etaTable['rmin'].data+etaTable['rmax'].data)/2
    
    plt.hlines(0,x[0],x[-1],linestyles=':')

    plt.errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=5,color=c)

    plt.ylabel(r'$\eta$')


plt.xticks(np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]))
plt.xlabel('R/Rv')

plt.legend(loc='upper left',ncol=2)

#plt.savefig('../plots/eta_vs_r/Eta_allsections_allvoids.jpg'.format(sec))
#%%
"""
Plot No Subsamples
"""
import seaborn as sns
#colors = sns.color_palette()
#blue, oran = colors[0], colors[1]

minradV = 7.
maxradV = 0.
sec=0

plt.rcParams['figure.figsize'] = (9, 6)
plt.rcParams['font.size'] = 18

r1=[.8,.9,1.,1.1,1.2,1.3,1.4]
r2=[.9,1.,1.1,1.2,1.3,1.4,1.5]
x = (np.array(r1)+np.array(r2))/2

plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')
plt.hlines(0,.8,1.5,linestyles=':')

for vtype, label, fmt in zip(['a','r','s'],\
                            ['All Voids', 'R-Voids', 'S-Voids'],\
                            ['o-','x--','^:']):

    etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype),\
            names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax','N'])

    yran_mean = etaTable['eta_random_mean'].data
    yran_err = etaTable['eta_random_std'].data

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    

    plt.errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color='k')

plt.ylabel(r'$\zeta$')
plt.xlabel(r'$\mathrm{r/R_v}$')
plt.xlim([.8,1.5])
plt.legend()

plt.savefig('../plots/eta_nosubsample.pdf')
# %%
"""
Plot High Spin, High Mass, High/Low Velocity
"""
import seaborn as sns
#colors = sns.color_palette()
#blue, oran = colors[0], colors[1]

minradV = 7.
maxradV = 0.
sec=3

plt.rcParams['figure.figsize'] = (9, 6)
plt.rcParams['font.size'] = 18

r1=[.8,.9,1.,1.1,1.2,1.3,1.4]
r2=[.9,1.,1.1,1.2,1.3,1.4,1.5]
x = (np.array(r1)+np.array(r2))/2

plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')
plt.hlines(0,.8,1.5,linewidth=.6,color='k',alpha=.7)

for vtype, label, fmt in zip(['a','r','s'],\
                            ['All Voids', 'R-Type', 'S-Type'],\
                            ['o-','x--','^:']):

    filename = '../data/eta/eta_vfilterhi_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])

    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    

    plt.errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color='k')

plt.ylabel(r'$\zeta$')
plt.xlabel(r'$\mathrm{r/R_v}$')
plt.xlim([.8,1.5])
plt.ylim([-4,8])

plt.legend(loc='upper left')
# %%
"""
Plot High Spin, High Mass, High/Low Velocity
"""
import seaborn as sns
#colors = sns.color_palette()
#blue, oran = colors[0], colors[1]

minradV = 7.
maxradV = 0.
sec=3

plt.rcParams['figure.figsize'] = (9, 6)
plt.rcParams['font.size'] = 18

r1=[.8,.9,1.,1.1,1.2,1.3,1.4]
r2=[.9,1.,1.1,1.2,1.3,1.4,1.5]
x = (np.array(r1)+np.array(r2))/2

plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')
plt.hlines(0,.8,1.5,linewidth=.6,color='k',alpha=.7)

for vtype, label, fmt in zip(['a','r','s'],\
                            ['All Voids', 'R-Type', 'S-Type'],\
                            ['o-','x--','^:']):

    filename = '../data/eta/eta_vfilterhi_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])

    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    

    plt.errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color='k')

plt.legend(loc='upper left')

for vtype, label, fmt in zip(['a','r','s'],\
                            ['All Voids', 'R-Type', 'S-Type'],\
                            ['o-','x--','^:']):

    filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])

    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    

    plt.errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color='C00')


plt.ylabel(r'$\zeta$')
plt.xlabel(r'$\mathrm{r/R_v}$')
plt.xlim([.8,1.5])
plt.ylim([-4,8])

#plt.legend(loc='upper left')
plt.savefig('../plots/hSpin-hMass-hlVel.jpg')

 # %%
"""
Plot High Spin, Low Velocity, High/Low Mass
"""
import seaborn as sns
#colors = sns.color_palette()
#blue, oran = colors[0], colors[1]

minradV = 7.
maxradV = 0.
plt.xscale
plt.rcParams['figure.figsize'] = (9, 6)
plt.rcParams['font.size'] = 18

r1=[.8,.9,1.,1.1,1.2,1.3,1.4]
r2=[.9,1.,1.1,1.2,1.3,1.4,1.5]
x = (np.array(r1)+np.array(r2))/2

plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
plt.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')
plt.hlines(0,.8,1.5,linewidth=.6,color='k',alpha=.7)

# sec=1
# for vtype, label, fmt in zip(['a','r','s'],\
#                             ['All Voids', 'R-Type', 'S-Type'],\
#                             ['o-','x--','^:']):

#     filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
#             .format(minradV,maxradV,sec,vtype)
#     etaTable = ascii.read(filename,\
#             names=['eta','eta_std','rmin','rmax','N'])

#     yran_mean = 1./(np.sqrt(2)-1)
#     yran_err = np.sqrt(28.1421/etaTable['N'])

#     y = (etaTable['eta'].data-yran_mean)/yran_err
#     yerr = etaTable['eta_std'].data/yran_err

    

#     plt.errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color='k')

# plt.legend(loc='upper left')

sec=3
for vtype, label, fmt in zip(['a','r','s'],\
                            ['All Voids', 'R-Type', 'S-Type'],\
                            ['o-','x--','^:']):

    filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])

    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    

    plt.errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color='C00')

plt.legend(loc='upper left')

plt.ylabel(r'$\zeta$')
plt.xlabel(r'$\mathrm{r/R_v}$')
plt.xlim([.8,1.5])
plt.ylim([-4,8])

plt.savefig('../plots/hSpin-lVel-hlMass.jpg')

# %%
"""
Trying Something
"""
import seaborn as sns
#colors = sns.color_palette()
#blue, oran = colors[0], colors[1]


colors = sns.color_palette()
plt.rcParams['figure.figsize'] = (15, 10)
plt.rcParams['font.size'] = 18

fig, axs = plt.subplots(2, 2, constrained_layout=True, sharex=True, sharey=True)

minradV = 7.
maxradV = 0.
plt.xscale
plt.rcParams['figure.figsize'] = (9, 6)
plt.rcParams['font.size'] = 18

r1=[.8,.9,1.,1.1,1.2,1.3,1.4]
r2=[.9,1.,1.1,1.2,1.3,1.4,1.5]
x = (np.array(r1)+np.array(r2))/2

for ax in axs:
    for ax_ in ax:
        ax_.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
        ax_.fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')
        ax_.hlines(0,.8,1.5,linewidth=.6,color='k',alpha=.7)
        ax_.set_ylim([-4,8])

sec=1
# for vtype, label, fmt in zip(['a','r','s'],\
#                             ['All Voids', 'R-Type', 'S-Type'],\
#                             ['o-','x--','^:']):

#     filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
#             .format(minradV,maxradV,sec,vtype)
#     etaTable = ascii.read(filename,\
#             names=['eta','eta_std','rmin','rmax','N'])

#     yran_mean = 1./(np.sqrt(2)-1)
#     yran_err = np.sqrt(28.1421/etaTable['N'])

#     y = (etaTable['eta'].data-yran_mean)/yran_err
#     yerr = etaTable['eta_std'].data/yran_err

    

#     axs[0,0].errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color='k')

for sec, c, label in zip([1,3],['C00','C03'],['Low Mass', 'High Mass']):
    vtype='a'
    filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])
    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])
    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err
    axs[0,0].errorbar(x,y,yerr=yerr,label=label,fmt='o-',capsize=3,ms=7,color=c)
    axs[0,0].legend(framealpha=.6)

for sec, c in zip([1,3],['C00','C03']):
    for vtype, label, fmt in zip(['r','s'], ['R-type', 'S-type'], ['x--','^:']):
        filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,sec,vtype)
        etaTable = ascii.read(filename,\
                names=['eta','eta_std','rmin','rmax','N'])
        yran_mean = 1./(np.sqrt(2)-1)
        yran_err = np.sqrt(28.1421/etaTable['N'])
        y = (etaTable['eta'].data-yran_mean)/yran_err
        yerr = etaTable['eta_std'].data/yran_err
        axs[0,1].errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color=c)
        #axs[0,1].legend(framealpha=.6)

sec=3
for vfilter, c, label in zip(['hi','lo'],['C00','C03'],['Hi V', 'Low V']):
    vtype='a'
    filename = '../data/eta/eta_vfilter{}_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(vfilter,minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])
    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])
    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err
    axs[1,0].errorbar(x,y,yerr=yerr,label=label,fmt='o-',capsize=3,ms=7,color=c)
    axs[1,0].legend(framealpha=.6)

for vfilter, c in zip(['hi','lo'],['C00','C03']):
    for vtype, label, fmt in zip(['r','s'], ['R-type', 'S-type'], ['x--','^:']):
        filename = '../data/eta/eta_vfilter{}_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                .format(vfilter,minradV,maxradV,sec,vtype)
        etaTable = ascii.read(filename,\
                names=['eta','eta_std','rmin','rmax','N'])
        yran_mean = 1./(np.sqrt(2)-1)
        yran_err = np.sqrt(28.1421/etaTable['N'])
        y = (etaTable['eta'].data-yran_mean)/yran_err
        yerr = etaTable['eta_std'].data/yran_err
        axs[1,1].errorbar(x,y,yerr=yerr,label=label,fmt=fmt,capsize=3,ms=7,color=c)
        #axs[1,1].legend(framealpha=.6)


#plt.legend(loc='upper left')

# %%
