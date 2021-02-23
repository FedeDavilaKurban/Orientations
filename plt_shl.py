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

        exp, minradV, rmin, rmax, sec, fxa, vtype = readExp('shl_'+id_str)

        #Control sample for Sigma Error
        random = ascii.read(writePath+'Proyectos/Orientations/data/shl_{}/randomfits_a.dat'.format(id_str),names=['xmean','ymean','ystd'])
        xmean, ymean, sigma = random['xmean'],random['ymean'],random['ystd']
        ax[i][j].fill_between(xmean,ymean-sigma,ymean+sigma,color='k',alpha=.4,label=r'$1\sigma$')
        ax[i][j].fill_between(xmean,ymean-2*sigma,ymean+2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
        ax[i][j].fill_between(xmean,ymean-3*sigma,ymean+3*sigma,color='k',alpha=.2,label=r'$3\sigma$')

        #Individual voids
        # nvs = range(len(readVoids(minradV)))
        # for nv in nvs:
        #     filename = '/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/shl_{}/ecdf_void{}'.format(id_str,nv)
        #     void_curve = ascii.read(filename,names=['cos','ecdf','y'])
        #     ax[i][j].plot(void_curve['cos'],void_curve['y'],alpha=.025,color='k')

        #JK mean and var
        data = ascii.read(writePath+'Proyectos/Orientations/data/shl_{}/jackknifeValues_a.dat'.format(id_str))
        x_ran,y_var,y_mean,y_sd = data['x_ran'],data['y_var'],data['y_mean'],data['y_sd']
        ax[i][j].fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
        ax[i][j].set_title('rmin={}Rv, rmax={}Rv, sec={}'.format(rmin,rmax,sec))


#plt.suptitle('MinradV={}Mpc, sec={}, fxa={}'.format(minradV,sec,fxa))
plt.ylim(-.015,.015)
#plt.tight_layout()
plt.savefig('../plots/shl_1.png')
plt.show(block=False)

# %%
# fig , ax = plt.subplots(nrows = 2, ncols = 2, sharex=True, sharey=True, figsize=(8,6))

# #id_int=15
# for i in range(0,2):
#     for j in range(0,2):
#         id_int+=1
#         id_str = str('{:03d}'.format(id_int))

#         exp, minradV, rmin, rmax, sec, fxa, vtype = readExp('shl_'+id_str)

#         #Control sample for Sigma Error
#         random = ascii.read(writePath+'Proyectos/Orientations/data/shl_{}/randomfits_a.dat'.format(id_str),names=['xmean','ymean','ystd'])
#         xmean, ymean, sigma = random['xmean'],random['ymean'],random['ystd']
#         ax[i][j].fill_between(xmean,ymean-sigma,ymean+sigma,color='k',alpha=.4,label=r'$1\sigma$')
#         ax[i][j].fill_between(xmean,ymean-2*sigma,ymean+2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
#         ax[i][j].fill_between(xmean,ymean-3*sigma,ymean+3*sigma,color='k',alpha=.2,label=r'$3\sigma$')

#         #Individual voids
#         # nvs = range(len(readVoids(minradV)))
#         # for nv in nvs:
#         #     filename = '/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/shl_{}/ecdf_void{}'.format(id_str,nv)
#         #     void_curve = ascii.read(filename,names=['cos','ecdf','y'])
#         #     ax[i][j].plot(void_curve['cos'],void_curve['y'],alpha=.025,color='k')

#         #JK mean and var
#         data = ascii.read(writePath+'Proyectos/Orientations/data/shl_{}/jackknifeValues_a.dat'.format(id_str))
#         x_ran,y_var,y_mean,y_sd = data['x_ran'],data['y_var'],data['y_mean'],data['y_sd']
#         ax[i][j].fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
#         ax[i][j].set_title('rmin={}Rv, rmax={}Rv, sec={}'.format(rmin,rmax,sec))

        



# #plt.suptitle('MinradV={}Mpc, sec={}, fxa={}'.format(minradV,sec,fxa))
# plt.ylim(-.015,.015)
# plt.tight_layout()
# plt.savefig('../plots/shl_2.png')
# plt.show(block=False)