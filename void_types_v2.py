"""
Determines void types and writes 'minradV{}_vtype.dat'
"""
#%%
from astropy.io import ascii
from astropy.table import Table
from config import writePath
import matplotlib.pyplot as plt
from orientationsTools import readVoids
import numpy as np

minradV = 7.
voids = readVoids(minradV)

list_voidtype = []
for nv in range(len(voids)):
    r3 = voids['maxdeltaint_2-3r'][nv]

    if r3>0:
        list_voidtype.append('s')
    if r3<0:
        list_voidtype.append('r')

data = np.column_stack([range(len(voids)),list_voidtype])
filename = writePath+'Proyectos/Orientations/data/vtype/minradV{}_vtype.dat'.format(minradV,nv)
names = ['id','type']
ascii.write(data,filename,names=names)
    




# %%
#plt.plot(profile['r'],profile['rho'])

# %%
