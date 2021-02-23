#%%
from astropy.io import ascii
from astropy.table import Table
from config import writePath
import matplotlib.pyplot as plt
from orientationsTools import readVoids
import numpy as np

minradV=7.0
nv=0


filename = writePath+'Proyectos/Orientations/data/vtype/minradV{}_vtype.dat'.format(minradV)
names = ['id','type']
vtypes = ascii.read(filename,names=names)

voids = readVoids(minrad=minradV)
svoids_id = np.where(vtypes['type']=='s')

voids[svoids_id]
# %%
