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
import random
from config import writePath, units

print("""
####################################
ECDF and Residues for stacked Voids
####################################
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


voids = readVoids(minradV,maxradV,vtype=vtype)

#I DO THIS FOR DISK SPACE REASONS
if len(voids)>500:
    print('Too many voids!')
    voids = voids[random.choices(list(range(len(voids))),k=300)]

    print('Num of voids:',len(voids))

nvs = range(len(voids))

ran_iter = 500

print('Jackknife...')
for jk in nvs:
    cos = []
    for nv in nvs:
        if jk==nv: continue
        cosTable = ascii.read(writePath+'Proyectos/Orientations/data/'+exp+'/cos_void{}_{}.dat'.format(nv,vtype),names=['cos'])
        cos.append( cosTable['cos'].data )

    cos_flattened = np.concatenate(cos).ravel()

    #N = len(cos_flattened) #number of galaxies of interest

    #ECDF
    cos,ecdf,y = ecdf_residues(cos_flattened)

    filename = writePath+'Proyectos/Orientations/data/'+exp+'/ecdf_jk{}_{}.dat'.format(jk,vtype)
    
    ascii.write(Table(np.column_stack([cos,ecdf(cos),y])),filename,names=['cos','ecdf','y'],overwrite=True)

print('Written files up to:',filename)

dataList, x_ran, y_var, y_mean, y_sd = jk_mean_sd(1000,sec,rmin,rmax,exp,len(voids),vtype) #El argumento es el numero de valores random del eje x para evaluar mean y SD de las curvas JK

#Fits the mean values of the JK curves corresponding to the x_ran
yfit,d_yfit,a2 = fits(x_ran,y_mean)

print("""
------------
a2 = {}
-----------
""".format(a2))
#OJO!!! calculo el a2 pero no lo estoy escribiendo

fitsfilename = writePath+'Proyectos/Orientations/data/'+exp+'/fits_{}.dat'.format(vtype)
print('Writing',fitsfilename)
ascii.write(Table(np.column_stack([yfit,d_yfit])),fitsfilename,names=['yfit','d_yfit'],overwrite=True)

jkfilename = writePath+'Proyectos/Orientations/data/'+exp+'/jackknifeValues_{}.dat'.format(vtype)
print('Writing',jkfilename)
ascii.write(Table(np.column_stack([x_ran,y_var,y_mean,y_sd])),jkfilename,names=['x_ran','y_var','y_mean','y_sd'],overwrite=True)


#Random Sample
print('Random sample...')
#Esto lo hago solo para saber cuantas galaxias de interes hay (en el stacking de voids)
#Seguro hay otra forma mejor de sacar ese numero porque ya lo tengo de antes solo que no lo escribo en ningun lado
n_forsum = []
for nv in nvs:
    #print(nv)
    cosTable = ascii.read(writePath+'Proyectos/Orientations/data/'+exp+'/cos_void{}_{}.dat'.format(nv,vtype),names=['cos'])
    n_forsum.append( len(cosTable['cos'].data) )

N_forRandom = np.sum(n_forsum) #Equals to the number of gxs in all the voids considered 
xmean, ymean, sigma = ranOrientations(ran_iter,N_forRandom)

randomfitsfilename = writePath+'Proyectos/Orientations/data/'+exp+'/randomfits_{}.dat'.format(vtype)
print('Writing',randomfitsfilename)
ascii.write(Table(np.column_stack([xmean,ymean,sigma])),randomfitsfilename,names=['xmean','ymean','ystd'],overwrite=True)

print('END')
