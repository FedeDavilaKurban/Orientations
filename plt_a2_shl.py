#%%
import sys
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import random
from config import writePath, units
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
#%%

fig , ax = plt.subplots(nrows = 3, ncols = 1, sharex=False, sharey=True, figsize=(8,10))

a2_mean = []
a2_std = []
a2_ran_mean = []
a2_ran_std = []

exp_name = 'shl'

"""
PLOT 1
"""
exp_ids = ["{0:03}".format(i) for i in range(1,6)]
my_xticks = []
pvalues = []
for exp_id in exp_ids:
    exp = exp_name+'_'+exp_id
    filename = writePath+'Proyectos/Orientations/data/'+exp+'_a2.dat'
    names = ['a2_mean','a2_std','a2_ran_mean','a2_ran_std','pvalue']
    a2Table = ascii.read(filename,names=names)

    a2_mean.append( a2Table['a2_mean'].data[0] )
    a2_std.append( a2Table['a2_std'].data[0] )
    a2_ran_mean.append( a2Table['a2_ran_mean'].data[0] )
    a2_ran_std.append( a2Table['a2_ran_std'].data[0] )

    pvalues.append(a2Table['pvalue'].data[0])

    exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(exp)
    my_xticks.append(r'${}R_v$'.format(str(rmin)))

a2_mean = np.array(a2_mean)
a2_std = np.array(a2_std)
a2_ran_mean = np.array(a2_ran_mean)
a2_ran_std = np.array(a2_ran_std)

x=[int(i) for i in exp_ids]

#plt.figure(figsize=(15, 5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax[0].fill_between(x,a2_ran_mean-a2_ran_std,a2_ran_mean+a2_ran_std,alpha=0.4,color=cycle[0],label=r'$\sigma_{\langle a2 \rangle_{Ran}}$')
ax[0].fill_between(x,a2_ran_mean-2*a2_ran_std,a2_ran_mean+2*a2_ran_std,alpha=0.4,color=cycle[0])
ax[0].fill_between(x,a2_ran_mean-3*a2_ran_std,a2_ran_mean+3*a2_ran_std,alpha=0.4,color=cycle[0])
ax[0].plot(x,a2_ran_mean,color=cycle[0],label=r'$\langle a2 \rangle_{Ran}$',marker='o',markerfacecolor='none')
ax[0].fill_between(x,a2_mean-a2_std,a2_mean+a2_std,alpha=0.6,color=cycle[1],label=r'$\sigma_{\langle a2 \rangle}$')
ax[0].plot(x,a2_mean,color=cycle[1],label=r'$\langle a2 \rangle$',marker='o',markerfacecolor='none')

for i in range(len(x)):
    ax[0].text(x[i]-.1,-0.005,'p='+str(pvalues[i]))

ax[0].set_xticks(x)
ax[0].set_xticklabels(my_xticks)

"""
PLOT 2
"""
exp_ids = ["{0:03}".format(i) for i in range(6,11)]
my_xticks = []
a2_mean = []
a2_std = []
a2_ran_mean = []
a2_ran_std = []
pvalues = []
for exp_id in exp_ids:
    exp = exp_name+'_'+exp_id
    filename = writePath+'Proyectos/Orientations/data/'+exp+'_a2.dat'
    names = ['a2_mean','a2_std','a2_ran_mean','a2_ran_std','pvalue']
    a2Table = ascii.read(filename,names=names)

    a2_mean.append( a2Table['a2_mean'].data[0] )
    a2_std.append( a2Table['a2_std'].data[0] )
    a2_ran_mean.append( a2Table['a2_ran_mean'].data[0] )
    a2_ran_std.append( a2Table['a2_ran_std'].data[0] )

    pvalues.append(a2Table['pvalue'].data[0])

    exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(exp)
    my_xticks.append(r'${}R_v$'.format(str(rmax)))

a2_mean = np.array(a2_mean)
a2_std = np.array(a2_std)
a2_ran_mean = np.array(a2_ran_mean)
a2_ran_std = np.array(a2_ran_std)

x=[int(i) for i in exp_ids]

#plt.figure(figsize=(15, 5))
cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax[1].fill_between(x,a2_ran_mean-a2_ran_std,a2_ran_mean+a2_ran_std,alpha=0.4,color=cycle[0])
ax[1].fill_between(x,a2_ran_mean-2*a2_ran_std,a2_ran_mean+2*a2_ran_std,alpha=0.4,color=cycle[0])
ax[1].fill_between(x,a2_ran_mean-3*a2_ran_std,a2_ran_mean+3*a2_ran_std,alpha=0.4,color=cycle[0])
ax[1].plot(x,a2_ran_mean,color=cycle[0],marker='o',markerfacecolor='none')
ax[1].fill_between(x,a2_mean-a2_std,a2_mean+a2_std,alpha=0.6,color=cycle[1])
ax[1].plot(x,a2_mean,color=cycle[1],marker='o',markerfacecolor='none')
ax[1].set_xticks(x)
ax[1].set_xticklabels(my_xticks)

for i in range(len(x)):
    ax[1].text(x[i]-.1,-0.005,'p='+str(pvalues[i]))

ax[1].set_xticks(x)
ax[1].set_xticklabels(my_xticks)

"""
PLOT 3
"""
exp_ids = ["{0:03}".format(i) for i in range(11,16)]
my_xticks = []
a2_mean = []
a2_std = []
a2_ran_mean = []
a2_ran_std = []
pvalues = []
for exp_id in exp_ids:
    exp = exp_name+'_'+exp_id
    filename = writePath+'Proyectos/Orientations/data/'+exp+'_a2.dat'
    names = ['a2_mean','a2_std','a2_ran_mean','a2_ran_std','pvalue']
    a2Table = ascii.read(filename,names=names)

    a2_mean.append( a2Table['a2_mean'].data[0] )
    a2_std.append( a2Table['a2_std'].data[0] )
    a2_ran_mean.append( a2Table['a2_ran_mean'].data[0] )
    a2_ran_std.append( a2Table['a2_ran_std'].data[0] )

    pvalues.append(a2Table['pvalue'].data[0])

    exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(exp)
    my_xticks.append(r'${}R_v$'.format(str(rmin)))

a2_mean = np.array(a2_mean)
a2_std = np.array(a2_std)
a2_ran_mean = np.array(a2_ran_mean)
a2_ran_std = np.array(a2_ran_std)

x=np.array([int(i) for i in exp_ids])
xp=x+(.5)
#plt.figure(figsize=(15, 5))

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax[2].fill_between(xp,a2_ran_mean-a2_ran_std,a2_ran_mean+a2_ran_std,alpha=0.4,color=cycle[0])
ax[2].fill_between(xp,a2_ran_mean-2*a2_ran_std,a2_ran_mean+2*a2_ran_std,alpha=0.4,color=cycle[0])
ax[2].fill_between(xp,a2_ran_mean-3*a2_ran_std,a2_ran_mean+3*a2_ran_std,alpha=0.4,color=cycle[0])
ax[2].plot(xp,a2_ran_mean,color=cycle[0],marker='o',markerfacecolor='none')
ax[2].fill_between(xp,a2_mean-a2_std,a2_mean+a2_std,alpha=0.6,color=cycle[1])
ax[2].plot(xp,a2_mean,color=cycle[1],marker='o',markerfacecolor='none')
for i in range(len(xp)):
    ax[2].text(xp[i]-.1,0.025,'p='+str(pvalues[i]))

ax[2].set_xticks(x)
ax[2].set_xticklabels(my_xticks)

ax[0].set_title('Increasing inner radius', fontsize=14)
ax[1].set_title('Increasing outer radius', fontsize=14)
ax[2].set_title('Increasing both radii', fontsize=14)

ax[0].text(1., .04, r'$R_{outer}=1.1R_v$', fontsize=14)
ax[1].text(6., .04, r'$R_{inner}=0.9R_v$', fontsize=14)
#ax[2].text(11., .04, r'$R_{outer}=R_{inner}+0.2R_v$', fontsize=14)

#plt.text(6, .005, 'Increasing outer radius', fontsize=14)
#plt.text(10.5, .005, 'Increasing both radii', fontsize=14)

ax[0].set_xlabel(r'$\mathrm{R_{inner}}$', fontsize=14)
ax[1].set_xlabel(r'$\mathrm{R_{outer}}$', fontsize=14)
ax[2].set_xlabel(r'$\mathrm{R_{inner}}$', fontsize=14)
for i in [0,1,2]:
    ax[i].set_ylabel('a2', fontsize=14)
    #ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

ax[0].legend(fontsize=14,ncol=2)

plt.tight_layout()
plt.savefig('../plots/a2_shl1.png')


# %%
