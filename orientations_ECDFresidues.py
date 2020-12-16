#%%
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *

voids = readVoids(7.)

nvs = range(len(voids))
rmaxs = [1.4]
rmins = [0.7]

#Choose one of the following as parameter in orientations()
secs = [3]
sections = [4,5,45,56,456]
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

ran_iter = 100
s5 = 0 #Sigma5; considera sólo galaxias con sigma5 <= s5. Poner s5=0 para desactivar

"""
ECDF and residue fits
"""
#Read cosines of galaxies in shells of voids of interest

for jk in nvs:
    print('jk: {}/{}'.format(jk,len(nvs)-1))
    cos = []
    for nv in nvs:
        if jk==nv: continue
        for rmin in rmins:
            for rmax in rmaxs:
                for sec in secs:
                    cosTable = ascii.read('../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'])
                    cos.append( cosTable['cos'].data )

    cos_flattened = np.concatenate(cos).ravel()
    #print(cos_flattened)
    N = len(cos_flattened) #N/2 is the number of galaxies of interest

    #ECDF, fits
    cos,y,ecdf,yfit,d_yfit,a2 = fits(cos_flattened)
    #print(a2)
    #a2string="{:.3f}".format(a2)
    ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_sec{}_rmin{}_rmax{}_jk{}'.format(sec,rmin,rmax,jk),names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)

    #Random Sample
    xmean, ymean, ystd = ranOrientations(ran_iter,N)
    ascii.write(Table(np.column_stack([xmean,ymean,ystd])),'../data/randomfits_sec{}_rmin{}_rmax{}_jk{}'.format(sec,rmin,rmax,jk),names=['xmean','ymean','ystd'],overwrite=True)

# %%
