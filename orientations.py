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
Orentations per void
"""
nvs = 0
rmaxs = [1.2]
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
            for sec in secs:
                print('nv, rmin, rmax, sec:',nv,rmin,rmax,sec)

                cos = orientations(gxs,tree,units,voids,nv,rmin,rmax,sec,s5)
                ascii.write(Table(np.reshape(cos,(len(cos),1))),'../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'],overwrite=True)

# %%
"""
ECDF and residue fits
"""
#Read cosines of galaxies in shells of voids of interest
cos = []
for nv in nvs:
    for rmin in rmins:
        for rmax in rmaxs:
            for sec in secs:
                cosTable = ascii.read('../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'])
                cos.append( cosTable['cos'].data )

cos_flattened = np.concatenate(cos).ravel()
print(cos_flattened)
N = len(cos_flattened) #N/2 is the number of galaxies of interest

#ECDF, fits
cos,y,ecdf,yfit,d_yfit,a2 = fits(cos_flattened)
print(a2)
#a2string="{:.3f}".format(a2)
ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_sec{}_rmin{}_rmax{}_s5:{}'.format(sec,rmin,rmax,s5),names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)

#Random Sample
xmean, ymean, ystd = ranOrientations(ran_iter,N)
ascii.write(Table(np.column_stack([xmean,ymean,ystd])),'../data/randomfits_sec{}_rmin{}_rmax{}_s5:{}'.format(sec,rmin,rmax,s5),names=['xmean','ymean','ystd'],overwrite=True)
