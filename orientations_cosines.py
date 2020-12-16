#%%

import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *

# %%
"""
Read Voids, Galaxies, and create cKDTree
"""
units = 'kpc'

voids = readVoids(7.)
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

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

# %%
"""
Orentations (cosines) per void
"""
nvs = 0
nvs = range(len(voids))
rmaxs = [1.4]
rmins = [0.7]

#Choose one of the following as parameter in orientations()
secs = [3]
sections = [4,5,45,56,456]
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

ran_iter = 100
s5 = 0 #Sigma5; considera sólo galaxias con sigma5 <= s5. Poner s5=0 para desactivar

for nv in nvs:
    for rmin in rmins:
        for rmax in rmaxs:
            for sec in allsections:
                print('nv, rmin, rmax, sec:',nv,rmin,rmax,sec)

                cos = orientations(gxs,tree,units,voids,nv,rmin,rmax,sec,s5)
                ascii.write(Table(np.reshape(cos,(len(cos),1))),'../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'],overwrite=True)

