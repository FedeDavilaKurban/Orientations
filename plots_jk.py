#%%
"""
Quiero plotear los residuos de una sola seccion (generalmente va a ser la 3)
con los errores jackknife

Eventualmente añadire los errores sigma 1 y 2.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from orientationsTools import readVoids
import random

#%%

minradV = 7.
#nvs = range(81)
rmin, rmax = .7, 1.4
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

for minradV in [5.,6.,7.,8.,9.,10.]:
    voids = readVoids(minradV)
    #I DO THIS FOR DISK SPACE REASONS
    if len(voids)>500:
        print('Too many voids!')
        voids = voids[random.choices(list(range(len(voids))),k=200)]
    nvs = range(len(voids))

    for sec in allsections:

        random = ascii.read('/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/minradVoid{}/randomfits_sec{}_rmin{}_rmax{}'.format(int(minradV),sec,rmin,rmax),names=['xmean','ymean','ystd'])
        xmean, ymean, sigma = random['xmean'],random['ymean'],random['ystd']

        data = ascii.read('/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/minradVoid{}/jackknifeValues_sec{}_rmin{}_rmax{}'.format(int(minradV),sec,rmin,rmax),names=['x_ran','y_var','y_mean','y_sd'])
        x_ran,y_var,y_mean,y_sd = data['x_ran'],data['y_var'],data['y_mean'],data['y_sd']

        fits = ascii.read('/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/minradVoid{}/fits_sec{}_rmin{}_rmax{}'.format(int(minradV),sec,rmin,rmax),names=['yfit','d_yfit'])
        yfit, d_yfit = fits['yfit'],fits['d_yfit']

        for nv in nvs:
            filename = '/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/minradVoid{}/ecdf_sec{}_rmin{}_rmax{}_void{}'.format(int(minradV),sec,rmin,rmax,nv)
            void_curve = ascii.read(filename,names=['cos','ecdf','y'])
            plt.plot(void_curve['cos'],void_curve['y'],alpha=.025,color='k')

        #dataList=[]
        #for jk in range(81):
        #    dataList.append( ascii.read('/media/fede/TOSHIBA EXT/Proyectos/Orientations/data/minradVoid/ecdf_sec{}_rmin{}_rmax{}_jk{}'.format(int(minradV),rmin,rmax,jk),names=['cos','ecdf','y']) )
        #    plt.plot(dataList[-1]['cos'],dataList[-1]['y'],alpha=.01,color='k')

        plt.fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
        #plt.plot(x_ran,yfit)

        plt.fill_between(xmean,ymean-sigma,ymean+sigma,color='k',alpha=.4,label=r'$1\sigma$')
        plt.fill_between(xmean,ymean-2*sigma,ymean+2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
        plt.fill_between(xmean,ymean-3*sigma,ymean+3*sigma,color='k',alpha=.2,label=r'$3\sigma$')
        plt.legend(loc=1)
        plotname = '../plots/minradVoid{}/minradVoid{}_sec{}_rmin{}_rmax{}.pdf'.format(int(minradV),int(minradV),sec,rmin,rmax)
        print(plotname)
        plt.savefig(plotname,dpi=800)
        #plt.show()
        plt.close()
