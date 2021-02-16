#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from orientationsTools import readVoids, readExp
import random
from config import writePath




#%%
fig , ax = plt.subplots(nrows = 5, ncols = 3, sharex=True, sharey=True, figsize=(8,12))

id_int = 0

for i in range(0,5):
    for j in range(0,3):
        id_int+=1
        id_str = str('{:03d}'.format(id_int))

        exp, minradV, rmin, rmax, sec, fxa = readExp('sec_'+id_str)

        #JK mean and var
        data = ascii.read(writePath+'Proyectos/Orientations/data/sec_{}/jackknifeValues.dat'.format(id_str))
        x_ran,y_var,y_mean,y_sd = data['x_ran'],data['y_var'],data['y_mean'],data['y_sd']
        ax[i][j].fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
        ax[i][j].set_title('Section {}'.format(sec))

        #Control sample for Sigma Error
        random = ascii.read(writePath+'Proyectos/Orientations/data/sec_{}/randomfits.dat'.format(id_str),names=['xmean','ymean','ystd'])
        xmean, ymean, sigma = random['xmean'],random['ymean'],random['ystd']
        ax[i][j].fill_between(xmean,ymean-sigma,ymean+sigma,color='k',alpha=.4,label=r'$1\sigma$')
        ax[i][j].fill_between(xmean,ymean-2*sigma,ymean+2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
        ax[i][j].fill_between(xmean,ymean-3*sigma,ymean+3*sigma,color='k',alpha=.2,label=r'$3\sigma$')

plt.suptitle('MinradV={}Mpc, rmin={}Rv, rmax={}Rv, fxa={}'.format(minradV,rmin,rmax,fxa))
plt.savefig('../plots/sec.png')
plt.show(block=False)

# %%
