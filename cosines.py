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

exp, minradV, maxradV, rmin, rmax, sec, fxa, vtype = readExp(sys.argv[1])
print('Codename of experiment:', exp)
print('minradVoid = {}Mpc'.format(minradV))
print('maxradVoid = {}Mpc'.format(maxradV))
print('rmin = {}Rvoid'.format(rmin))
print('rmax = {}Rvoid'.format(rmax))
print('sec =',sec)
print('fxa =',fxa)
print('vtype =',vtype)

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

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

"""
Orentations (cosines) per void
"""
#nvs = 0
nvs = range(len(voids))

n_total = 0 
for nv in nvs:
    cos = orientations(gxs,tree,units,voids,nv,rmin,rmax,sec,fxa)

    if exp=='shl_019':
        cos = random.choices(cos,k=int(len(cos)/10))


    filename = writePath+'Proyectos/Orientations/data/'+exp+'/cos_void{}_{}.dat'.format(nv,vtype)
    ascii.write(Table(np.reshape(cos,(len(cos),1))),filename,names=['cos'],overwrite=True)
    
    print('Num of gxs in Void{}:'.format(nv),len(cos))

    n_total += len(cos)
print('Written up to', filename)

print('N_total of gxs: ',n_total)

print('N_avg of gxs: ',float(n_total)/len(voids))
