#%%
"""
Gets the cosines of the orientations of the galaxies in void shells 
"""
import sys
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import random
from config import writePath, units

# %%
"""
Read Voids, Galaxies, and create cKDTree
"""
print("""
##################
Cosines
##################
""")

exp, minradV, rmin, rmax, sec, fxa = readExp(sys.argv[1])

voids = readVoids(minrad=minradV)
print('Num of voids:',len(voids))

#I DO THIS FOR DISK SPACE REASONS
if len(voids)>500:
    print('Too many voids!')
    voids = voids[random.choices(list(range(len(voids))),k=200)]

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

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

"""
Orentations (cosines) per void
"""
#nvs = 0
nvs = range(len(voids))

ran_iter = 100

for nv in nvs:
    cos = orientations(gxs,tree,units,voids,nv,rmin,rmax,sec,fxa)
    filename = writePath+'Proyectos/Orientations/data/'+exp+'/cos.dat'
    ascii.write(Table(np.reshape(cos,(len(cos),1))),filename,names=['cos'],overwrite=True)
print('Written up to', filename)

print('N_total of gxs: ',len(cos)/2)