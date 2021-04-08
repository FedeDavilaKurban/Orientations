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

#%%
fig , ax = plt.subplots(nrows = 2, ncols = 1, sharex=False, sharey=True, figsize=(8,6))

a2_mean = []
a2_std = []
a2_ran_mean = []
a2_ran_std = []

exp_name = 's5'

"""
PLOT 1
"""
exp_ids = ["{0:03}".format(i) for i in range(1,4)]
my_xticks = ['Very isolated'+r'$\,(\sigma_5\geq P_{90\%})$','Very crowded'+r'$\,(\sigma_5\leq P_{10\%})$','No subsample']
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
    #my_xticks.append(r'${}$'.format(str(s5)))

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
ax[0].plot(x,a2_ran_mean,color=cycle[0],label=r'$\langle a2 \rangle_{Ran}$')
ax[0].fill_between(x,a2_mean-a2_std,a2_mean+a2_std,alpha=0.6,color=cycle[1],label=r'$\sigma_{\langle a2 \rangle}$')
ax[0].plot(x,a2_mean,color=cycle[1],label=r'$\langle a2 \rangle$')

for i in range(len(x)):
    ax[0].text(x[i],-0.0075,str(pvalues[i]))

ax[0].set_xticks(x)
ax[0].set_xticklabels([])

"""
PLOT 2
"""
exp_ids = ["{0:03}".format(i) for i in range(4,7)]
#my_xticks = []
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
    #my_xticks.append(r'${}$'.format(str(s5)))

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
ax[1].plot(x,a2_ran_mean,color=cycle[0])
ax[1].fill_between(x,a2_mean-a2_std,a2_mean+a2_std,alpha=0.6,color=cycle[1])
ax[1].plot(x,a2_mean,color=cycle[1])

ax[0].set_title(r'$\sigma_5$'+' Subsamples', fontsize=14)
ax[0].text(1,.0125,r'$R_v\geq7\mathrm{Mpc}$', fontsize=12)
ax[1].text(4,.0125,r'$6\mathrm{Mpc}\geq R_v\geq7\mathrm{Mpc}$', fontsize=12)

for i in range(len(x)):
    ax[1].text(x[i]-.05,-0.0075,'p='+str(pvalues[i]))
ax[1].set_xticks(x)
ax[1].set_xticklabels(my_xticks)

for i in [0,1]:
    ax[i].set_ylabel('a2', fontsize=14)
    #ax[i].xaxis.set_major_locator(MaxNLocator(integer=True))

ax[0].legend(fontsize=14,ncol=2)

plt.tight_layout()
plt.savefig('../plots/a2_s5.png')

# %%
