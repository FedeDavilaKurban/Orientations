#%%
"""
Calculates ECDF, the residues, and their fits, using the cosines from orientations_cosines.py
"""
import sys
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import matplotlib.pyplot as plt
import random
from config import writePath, units

print("""
####################################
ECDF and Residues for single Voids
####################################
""")

exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(sys.argv[1])
print('Codename of experiment:', exp)
print('minradVoid = {}Mpc'.format(minradV))
print('maxradVoid = {}Mpc'.format(maxradV))
print('rmin = {}Rvoid'.format(rmin))
print('rmax = {}Rvoid'.format(rmax))
print('sec =',sec)
print('s5 =',s5)
print('vtype =',vtype)

voids = readVoids(minradV,maxradV,vtype=vtype)

#I DO THIS FOR DISK SPACE REASONS
if len(voids)>500:
    print('Too many voids!')
    voids = voids[random.choices(list(range(len(voids))),k=300)]

    print('Num of voids: ',len(voids))

nvs = range(len(voids))

for nv in nvs:
    cosTable = ascii.read(writePath+'Proyectos/Orientations/data/'+exp+'/cos_void{}_{}.dat'.format(nv,vtype),names=['cos'])
    cosArray = cosTable['cos'].data 

    #N = len(cosArray) 

    #ECDF
    cos,ecdf,y = ecdf_residues(cosArray)

    filename = writePath+'Proyectos/Orientations/data/'+exp+'/ecdf_void{}_{}'.format(nv,vtype)
    #print(filename)
    ascii.write(Table(np.column_stack([cos,ecdf(cos),y])),filename,names=['cos','ecdf','y'],overwrite=True)


print('Written up to', filename)
