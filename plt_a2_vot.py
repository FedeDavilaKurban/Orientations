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
a2_mean = []
a2_std = []
a2_ran_mean = []
a2_ran_std = []

plt.figure(figsize=(7,4)) 

exp_name = 'vot'
exp_ids = ["{0:03}".format(i) for i in range(1,4)]

my_xticks = ['All Voids','Rising Voids','Shell Voids']
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
    #my_xticks.append(r'${}$'.format(str(vtype)))

a2_mean = np.array(a2_mean)
a2_std = np.array(a2_std)
a2_ran_mean = np.array(a2_ran_mean)
a2_ran_std = np.array(a2_ran_std)

x=[int(i) for i in exp_ids]

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.fill_between(x,a2_ran_mean-a2_ran_std,a2_ran_mean+a2_ran_std,alpha=0.4,color=cycle[0],label=r'$\sigma_{\langle a2 \rangle_{Ran}}$')
plt.fill_between(x,a2_ran_mean-2*a2_ran_std,a2_ran_mean+2*a2_ran_std,alpha=0.4,color=cycle[0])
plt.fill_between(x,a2_ran_mean-3*a2_ran_std,a2_ran_mean+3*a2_ran_std,alpha=0.4,color=cycle[0])

plt.plot(x,a2_ran_mean,color=cycle[0],label=r'$\langle a2 \rangle_{Ran}$')
plt.fill_between(x,a2_mean-a2_std,a2_mean+a2_std,alpha=0.6,color=cycle[1],label=r'$\sigma_{\langle a2 \rangle}$')
plt.plot(x,a2_mean,color=cycle[1],label=r'$\langle a2 \rangle$')

plt.title('Void Types', fontsize=14)
plt.text(1,.006,r'$R_v\geq7\mathrm{Mpc}$', fontsize=12)

for i in range(len(x)):
    plt.text(x[i]-.05,-0.01,'p='+str(pvalues[i]))

plt.xticks(x, my_xticks)

plt.legend(fontsize=14,ncol=2)

plt.savefig('../plots/a2_vot.png')
# %%

