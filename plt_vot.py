#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from orientationsTools import readVoids, readExp
import random
from config import writePath

#%%
fig , ax = plt.subplots(3, sharex=True, sharey=True, figsize=(8,10))

nvs = range(len(readVoids(7.)))

id_int = 0

for i in range(0,3):
    id_int+=1
    id_str = str('{:03d}'.format(id_int))

    exp, minradV, rmin, rmax, sec, fxa, vtype = readExp('vot_'+id_str)

    #Control sample for Sigma Error
    random = ascii.read(writePath+'Proyectos/Orientations/data/vot_{}/randomfits_{}.dat'.format(id_str,vtype),names=['xmean','ymean','ystd'])
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
    data = ascii.read(writePath+'Proyectos/Orientations/data/vot_{}/jackknifeValues_{}.dat'.format(id_str,vtype))
    x_ran,y_var,y_mean,y_sd = data['x_ran'],data['y_var'],data['y_mean'],data['y_sd']
    ax[i].fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
    #ax[i].set_title('Rv = {}'.format(minradV))



#plt.suptitle('MinradV={}Mpc, sec={}, fxa={}'.format(minradV,sec,fxa))
plt.tight_layout()
plt.savefig('../plots/vot_2.png')
plt.show(block=False)

# %%
