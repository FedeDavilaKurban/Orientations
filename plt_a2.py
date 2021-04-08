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

exp_name = 'voR'
exp_ids = ["{0:03}".format(i) for i in range(2,8)]

for exp_id in exp_ids:
    exp = exp_name+'_'+exp_id
    filename = writePath+'Proyectos/Orientations/data/'+exp+'a2.dat'
    names = ['a2_mean','a2_std','a2_ran_mean','a2_ran_std']
    a2Table = ascii.read(filename,names=names)

    a2_mean.append( a2Table['a2_mean'].data[0] )
    a2_std.append( a2Table['a2_std'].data[0] )
    a2_ran_mean.append( a2Table['a2_ran_mean'].data[0] )
    a2_ran_std.append( a2Table['a2_ran_std'].data[0] )

a2_mean = np.array(a2_mean)
a2_std = np.array(a2_std)
a2_ran_mean = np.array(a2_ran_mean)
a2_ran_std = np.array(a2_ran_std)

x=[int(i) for i in exp_ids]

cycle = plt.rcParams['axes.prop_cycle'].by_key()['color']
plt.fill_between(x,a2_ran_mean-a2_ran_std,a2_ran_mean+a2_ran_std,alpha=0.6,color=cycle[0])
plt.plot(x,a2_ran_mean,color=cycle[0])
plt.fill_between(x,a2_mean-a2_std,a2_mean+a2_std,alpha=0.6,color=cycle[1])
plt.plot(x,a2_mean,color=cycle[1])
# %%
