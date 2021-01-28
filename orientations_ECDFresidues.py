#%%
"""
Calculates ECDF, the residues, and their fits, using the cosines from orientations_cosines.py
"""

import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *

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
    for jk in nvs:
        #print('jk: {}/{}'.format(jk,len(nvs)-1))
        cos = []
        for nv in nvs:
            if jk==nv: continue
            cosTable = ascii.read('../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'])
            cos.append( cosTable['cos'].data )

        cos_flattened = np.concatenate(cos).ravel()

        N = len(cos_flattened) #N/2 is the number of galaxies of interest

        #ECDF
        cos,ecdf,y = ecdf_residues(cos_flattened)

        filename = '../data/ecdf_sec{}_rmin{}_rmax{}_jk{}'.format(sec,rmin,rmax,jk)
        #print(filename)
        ascii.write(Table(np.column_stack([cos,ecdf(cos),y])),filename,names=['cos','ecdf','y'],overwrite=True)

    dataList, x_ran, y_var, y_mean, y_sd = jk_mean_sd(1000,sec,rmin,rmax) #El argumento es el numero de valores random del eje x para evaluar mean y SD de las curvas JK

    #Fits the mean values of the JK curves corresponding to the x_ran
    yfit,d_yfit,a2 = fits(x_ran,y_mean)
    #OJO!!! calculo el a2 pero no lo estoy escribiendo

    ascii.write(Table(np.column_stack([yfit,d_yfit])),'../data/fits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['yfit','d_yfit'],overwrite=True)
    ascii.write(Table(np.column_stack([x_ran,y_var,y_mean,y_sd])),'../data/jackknifeValues_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['x_ran','y_var','y_mean','y_sd'],overwrite=True)



    #Random Sample
    n_forsum = []
    for nv in nvs:
        #print(nv)
        cosTable = ascii.read('../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'])
        n_forsum.append( len(cosTable['cos'].data) )

    N_forRandom = np.sum(n_forsum) #Equals to the number of gxs in all the voids considered 
    xmean, ymean, sigma = ranOrientations(ran_iter,N_forRandom)
    ascii.write(Table(np.column_stack([xmean,ymean,sigma])),'../data/randomfits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['xmean','ymean','ystd'],overwrite=True)


# %%
