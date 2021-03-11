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

# %%

minradV = 7
maxradV = 8
rmin = 0.9 
rmax = 1.2
sec = 0
#fxa =  
vtype = 'a'

voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)
print('Num of voids:',len(voids))

#I DO THIS FOR DISK SPACE REASONS
if len(voids)>500:
    print('Too many voids!')
    voids = voids[random.choices(list(range(len(voids))),k=300)]

    print('Num of voids: ',len(voids))

gxs = readTNG()

if units=='Mpc': 
    lbox=205.
    for col in ['x','y','z']:
        gxs[col]/=1000.
if units=='kpc': 
    lbox=205000.

#Applying some filters
gxs.remove_row(np.where(gxs['y']==lbox)[0][0])
gxs.remove_row(np.where(gxs['x']<0.)[0][0])

gxs = gxs[:1000]
data = np.column_stack((gxs['x'],gxs['y'],gxs['z']))
tree = spatial.cKDTree(data=data, boxsize=lbox)

nearest_dist, nearest_ind = tree.query(data, k=2)  # k=2 nearest neighbors where k1 = identity
dists = nearest_dist[:, 1]    # drop id; assumes sorted -> see args!
inds = nearest_ind[:, 1]
# %%
p10=np.percentile(dists,q=10)
p90=np.percentile(dists,q=90)
plt.hist(dists)
plt.axvline(p10,ymin=0,ymax=250,ls='--',color='r')
plt.axvline(p90,ymin=0,ymax=250,ls='--',color='r')

print(len(dists[dists>p90]))
print(len(dists[dists<p10]))
# %%
