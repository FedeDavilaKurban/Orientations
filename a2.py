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
import matplotlib.pyplot as plt

#%%
print("""
####################################
ECDF and Residues for stacked Voids
####################################
""")

exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp('voR_001')
print('Codename of experiment:', exp)
print('minradVoid = {}Mpc'.format(minradV))
print('maxradVoid = {}Mpc'.format(maxradV))
print('rmin = {}Rvoid'.format(rmin))
print('rmax = {}Rvoid'.format(rmax))
print('sec =',sec)
print('s5 =',s5)
print('vtype =',vtype)


voids = readVoids(minradV,maxradV,vtype=vtype)
nvs = range(len(voids))

ran_iter = 500

#%%
print('Jackknife...')
a2List=[]
for jk in nvs:

    filename = writePath+'Proyectos/Orientations/data/'+exp+'/ecdf_jk{}_{}.dat'.format(jk,vtype)
    ascii.read(filename,names=['cos','ecdf','y'])

    dataList, x_ran, y_var, y_mean, y_sd = jk_mean_sd(1000,sec,rmin,rmax,exp,len(voids),vtype) #El argumento es el numero de valores random del eje x para evaluar mean y SD de las curvas JK

    #Fits the mean values of the JK curves corresponding to the x_ran
    yfit,d_yfit,a2 = fits(x_ran,y_mean)

    a2List.append(a2)

a2_var = np.var(a2List,ddof=1) 
a2_mean = np.mean(a2List) 
a2_sd = np.sqrt(a2_var)


# %%
def ranOrientations_(n_iter,N):
    """
    Generate random orientations, calculate mean and standard deviations
    Returns: xmean (the random "cosines") , ymean (mean value of the fits), ystd (stand dev of fits)

    iter: integer number of iterations
    """
    import numpy as np

    a2_ran = []
    yran_fit = []
    xran_fit = []
    for _ in range(n_iter):
        np.random.seed(_)
        rancos = np.random.uniform(0.,1.,int(N/2))

        cos_,ecdf_,y_ = ecdf_residues(rancos)
        yfit_,d_yfit_,a2_ = fits(cos_,y_)
        
        a2_ran.append(a2_)
        yran_fit.append(yfit_)
        xran_fit.append(cos_)

    xmean = np.mean(xran_fit,axis=0)
    ymean = np.mean(yran_fit,axis=0)
    ystd = np.std(yran_fit,axis=0)

    a2_ran_mean = np.mean(a2_ran)
    a2_ran_std = np.std(a2_ran)

    del cos_, ecdf_, y_, yfit_, d_yfit_, a2_, yran_fit, xran_fit, a2_ran

    return xmean,ymean,ystd,a2_ran_mean,a2_ran_std

#%%

#Random Sample
print('Random sample...')
#Esto lo hago solo para saber cuantas galaxias de interes hay (en el stacking de voids)
#Seguro hay otra forma mejor de sacar ese numero porque ya lo tengo de antes solo que no lo escribo en ningun lado
# n_forsum = []
# for nv in nvs:
#     #print(nv)
#     cosTable = ascii.read(writePath+'Proyectos/Orientations/data/'+exp+'/cos_void{}_{}.dat'.format(nv,vtype),names=['cos'])
#     n_forsum.append( len(cosTable['cos'].data) )

N_forRandom = 5000 #Equals to the number of gxs in all the voids considered 'np.sum(n_forsum)'
xmean, ymean, sigma, a2_ran_mean, a2_ran_std = ranOrientations_(ran_iter,N_forRandom)
#%%
plt.plot(xmean,ymean)
plt.plot(xmean,sigma)
plt.hlines(a2_ran_mean, -1., 1.)
plt.hlines(a2_ran_mean+a2_ran_std, -1., 1.,linestyles='--')


#%%
plt.plot(xmean,ymean)
plt.fill_between(xmean,-sigma,sigma,color='k',alpha=.4,label=r'$1\sigma$')
plt.fill_between(xmean,-2*sigma,2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
plt.fill_between(xmean,-3*sigma,3*sigma,color='k',alpha=.2,label=r'$3\sigma$')
plt.hlines(a2_ran_mean, -1., 1.)

plt.hlines(a2_ran_std, -1., 1.,linestyles='--')
plt.hlines(-a2_ran_std, -1., 1.,linestyles='--')

plt.hlines(2*a2_ran_std, -1., 1.,linestyles='--')
plt.hlines(-2*a2_ran_std, -1., 1.,linestyles='--')

plt.hlines(3*a2_ran_std, -1., 1.,linestyles='--')
plt.hlines(-3*a2_ran_std, -1., 1.,linestyles='--')
# %%
