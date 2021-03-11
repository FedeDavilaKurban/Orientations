#%%
import sys
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import cosCalc, readVoids, JvsM, readExp
import random
from config import writePath, illustrisPath, units
import matplotlib.pyplot as plt
#%%
#nvs = range(len(readVoids(7.)))

id_int = 0
fig , ax = plt.subplots(nrows = 2, ncols = 7, sharex=True, sharey=True, figsize=(15,3))

for i in range(0,2):
    for j in range(0,7):

        id_int+=1
        id_str = str('{:03d}'.format(id_int))
        if id_str=='014': break
        exp, minradV, maxradV, rmin, rmax, sec, fxa, vtype = readExp('sfr_'+id_str)

        data = ascii.read(writePath+'Proyectos/Orientations/data/{}/sfr_cos.dat'.format(exp),names=['cos','sfr'])
        hist=ax[i][j].hist2d(data['cos'],np.log10(data['sfr']),bins=[10,10])
        ax[i][j].set_ylim(-2.5,1.)
plt.colorbar(hist[3])

# for ax in axes.flat:
#     #ax.set_axis_off()

#     id_int+=1
#     id_str = str('{:03d}'.format(id_int))
#     if id_str=='014': break
#     exp, minradV, maxradV, rmin, rmax, sec, fxa, vtype = readExp('sfr_'+id_str)

#     data = ascii.read(writePath+'Proyectos/Orientations/data/{}/sfr_cos.dat'.format(exp),names=['cos','sfr'])

#     im = ax.imshow(np.column_stack((data['cos'],np.log10(data['sfr']))), cmap='viridis', vmin=0, vmax=100)
# cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.95)

#%%
data = ascii.read(writePath+'Proyectos/Orientations/data/sfr_001/sfr_cos.dat'.format(exp),names=['cos','sfr'])
cos = data['cos']
sfr = data['sfr']

fig, ax = plt.subplots() 
plt.hist2d(cos,np.log10(sfr),bins=[10,10])
plt.xlabel(r'$\mathrm{cos}\theta$')
plt.ylabel('SFR')
plt.colorbar() 
#plt.savefig('../plots/{}.png'.format(exp))
# %%
