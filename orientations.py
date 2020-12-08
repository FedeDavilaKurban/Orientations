#%%

import numpy as np
from scipy import spatial
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
Orentations
"""

#nv = 0
rmaxs = [1.2]
rmins = [0.7]

sec = [6]
sections = [4,5,45,56,456]
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

ran_iter = 100
s5 = 0 #Sigma5; considera sólo galaxias con sigma5 <= s5. Poner s5=0 para desactivar

for s5 in []:
    orientations(rmins,rmaxs,sections,s5,ran_iter,write=True,ranfits=True)

