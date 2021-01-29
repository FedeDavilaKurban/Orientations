#%%
"""
Calculates ECDF, the residues, and their fits, using the cosines from orientations_cosines.py
"""

import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import matplotlib.pyplot as plt


voids = readVoids(7.)

nvs = range(len(voids))
rmax = 1.4
rmin = 0.7
#sec = 3

#Choose one of the following as parameter in orientations()
section = [3]
sections = [4,5,45,56,456]
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

secs = [1,2,3]

ran_iter = 100
s5 = 0 #Sigma5; considera sólo galaxias con sigma5 <= s5. Poner s5=0 para desactivar

for sec in allsections:
    print(sec)
    for nv in nvs:
        cosTable = ascii.read('../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'])
        cosArray = cosTable['cos'].data 

        #cos_flattened = np.concatenate(cos).ravel()

        N = len(cosArray) #N/2 is the number of galaxies of interest

        #ECDF
        cos,ecdf,y = ecdf_residues(cosArray)

        #plt.plot(cos,y,alpha=.01,color='k')

        filename = '../data/ecdf_sec{}_rmin{}_rmax{}_void{}'.format(sec,rmin,rmax,nv)
        print(filename)
        ascii.write(Table(np.column_stack([cos,ecdf(cos),y])),filename,names=['cos','ecdf','y'],overwrite=True)


# %%
