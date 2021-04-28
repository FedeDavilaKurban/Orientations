#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from orientationsTools import readVoids, readExp
import random
from config import writePath

#%%
fig , ax = plt.subplots(6, sharex=True, sharey=True, figsize=(8,10))

id_int = 0

for i in range(0,6):
    #id_int+=1
    id_str = str('{:03d}'.format(i+1))

    exp, minradV, maxradV, rmin, rmax, sec, fxa, vtype = readExp('voR_'+id_str)

    #Control sample for Sigma Error
    random = ascii.read(writePath+'Proyectos/Orientations/data/voR_{}/randomfits_a.dat'.format(id_str),names=['xmean','ymean','ystd'])
    xmean, ymean, sigma = random['xmean'],random['ymean'],random['ystd']
    ax[i].fill_between(xmean,ymean-sigma,ymean+sigma,color='k',alpha=.4,label=r'$1\sigma$')
    ax[i].fill_between(xmean,ymean-2*sigma,ymean+2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
    ax[i].fill_between(xmean,ymean-3*sigma,ymean+3*sigma,color='k',alpha=.2,label=r'$3\sigma$')

    # #Individual voids
    # nvs = range(len(readVoids(minradV)))
    # for nv in nvs:
    #     filename = '/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/voR_{}/ecdf_void{}'.format(id_str,nv)
    #     void_curve = ascii.read(filename,names=['cos','ecdf','y'])
    #     plt.plot(void_curve['cos'],void_curve['y'],alpha=.025,color='k')

    #JK mean and var
    data = ascii.read(writePath+'Proyectos/Orientations/data/voR_{}/jackknifeValues_a.dat'.format(id_str))
    x_ran,y_var,y_mean,y_sd = data['x_ran'],data['y_var'],data['y_mean'],data['y_sd']
    ax[i].fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
    ax[i].set_title(r'${}\leq R_v \leq {}$'.format(minradV,maxradV))



#plt.suptitle('MinradV={}Mpc, sec={}, fxa={}'.format(minradV,sec,fxa))
plt.tight_layout()
plt.savefig('../plots/voR3.png')
plt.show(block=False)

# %%
