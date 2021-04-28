#%%
import sys
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import random
from config import writePath, units
import matplotlib.pyplot as plt

# %%
def ranOrientations_(n_iter,N, a2):
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

    pvalue =  float(len([x for x in a2_ran if x < a2])) / float(len(a2_ran)) 

    a2_ran_mean = np.mean(a2_ran)
    a2_ran_std = np.std(a2_ran)

    del cos_, ecdf_, y_, yfit_, d_yfit_, a2_, yran_fit, xran_fit, a2_ran

    return xmean,ymean,ystd,a2_ran_mean,a2_ran_std, pvalue
#%%
print("""
####################################
ECDF and Residues for stacked Voids
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
    voids = voids[random.choices(list(range(len(voids))),k=500)]

    print('Num of voids:',len(voids))
nvs = range(len(voids))

ran_iter = 500

#%%
print('Jackknife...')
dataList=[]
a2List=[]
for jk in nvs:
    #print(jk,nvs[-1])
    filename = writePath+'Proyectos/Orientations/data/'+exp+'/ecdf_jk{}_{}.dat'.format(jk,vtype)
    names = ['cos','ecdf','y']  
    dataList.append(ascii.read(filename, names=names))
    #plt.plot(dataList[-1]['cos'],dataList[-1]['y'],alpha=.1,color='k')
    yfit,d_yfit,a2 = fits(dataList[-1]['cos'],dataList[-1]['y'])
    a2List.append(a2)

a2_var = np.var(a2List,ddof=1) 
a2_mean = np.mean(a2List) 
a2_std = np.sqrt(a2_var) 

del dataList, a2List

#%%

#Random Sample
print('Random sample...')
#Esto lo hago solo para saber cuantas galaxias de interes hay (en el stacking de voids)
#Seguro hay otra forma mejor de sacar ese numero porque ya lo tengo de antes solo que no lo escribo en ningun lado
n_forsum = []
for nv in nvs:
    #print(nv,nvs[-1])
    cosTable = ascii.read(writePath+'Proyectos/Orientations/data/'+exp+'/cos_void{}_{}.dat'.format(nv,vtype),names=['cos'])
    n_forsum.append( len(cosTable['cos'].data) )

N_forRandom = np.sum(n_forsum) #Equals to the number of gxs in all the voids considered 'np.sum(n_forsum)'
xmean, ymean, sigma, a2_ran_mean, a2_ran_std, pvalue = ranOrientations_(ran_iter,N_forRandom,a2_mean)

#%%
# """
# ######################################################################################
# PLOTS FOR RANDOM SAMPLE
# ######################################################################################
# """
# plt.plot(xmean,ymean)
# plt.plot(xmean,sigma)
# plt.hlines(a2_ran_mean, -1., 1.)
# plt.hlines(a2_ran_mean+a2_ran_std, -1., 1.,linestyles='--')

# #%%
# plt.plot(xmean,ymean)
# plt.fill_between(xmean,-sigma,sigma,color='k',alpha=.4,label=r'$1\sigma$')
# plt.fill_between(xmean,-2*sigma,2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
# plt.fill_between(xmean,-3*sigma,3*sigma,color='k',alpha=.2,label=r'$3\sigma$')
# plt.hlines(a2_ran_mean, -1., 1.)

# plt.hlines(a2_ran_std, -1., 1.,linestyles='--')
# plt.hlines(-a2_ran_std, -1., 1.,linestyles='--')

# plt.hlines(2*a2_ran_std, -1., 1.,linestyles='--')
# plt.hlines(-2*a2_ran_std, -1., 1.,linestyles='--')

# plt.hlines(3*a2_ran_std, -1., 1.,linestyles='--')
# plt.hlines(-3*a2_ran_std, -1., 1.,linestyles='--')
# """
# ######################################################################################
# """

filename = writePath+'Proyectos/Orientations/data/'+exp+'_a2.dat'
dataTable = Table(np.column_stack((a2_mean,a2_std,a2_ran_mean,a2_ran_std,pvalue)))
names = ['a2_mean','a2_std','a2_ran_mean','a2_ran_std','p_value']
ascii.write(dataTable,filename,names=names)

print("Created file ", filename)
print("END")



# %%
