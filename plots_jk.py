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
#%%
rmin, rmax = .7, 1.4
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

for sec in allsections:

    random = ascii.read('../data/randomfits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['xmean','ymean','ystd'])
    xmean, ymean, sigma = random['xmean'],random['ymean'],random['ystd']

    data = ascii.read('../data/jackknifeValues_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['x_ran','y_var','y_mean','y_sd'])
    x_ran,y_var,y_mean,y_sd = data['x_ran'],data['y_var'],data['y_mean'],data['y_sd']

    fits = ascii.read('../data/fits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['yfit','d_yfit'])
    yfit, d_yfit = fits['yfit'],fits['d_yfit']

    #dataList=[]
    #for jk in range(81):
    #    dataList.append( ascii.read('../data/ecdf_sec{}_rmin{}_rmax{}_jk{}'.format(sec,rmin,rmax,jk),names=['cos','ecdf','y']) )
    #    plt.plot(dataList[-1]['cos'],dataList[-1]['y'],alpha=.01,color='k')

    plt.fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
    plt.plot(x_ran,yfit)

    plt.fill_between(xmean,ymean-sigma,ymean+sigma,color='k',alpha=.4,label=r'$1\sigma$')
    plt.fill_between(xmean,ymean-2*sigma,ymean+2*sigma,color='k',alpha=.3,label=r'$2\sigma$')
    plt.fill_between(xmean,ymean-3*sigma,ymean+3*sigma,color='k',alpha=.2,label=r'$3\sigma$')
    plt.legend()
    plotname = 'sec{}_rmin{}_rmax{}.png'.format(sec,rmin,rmax)
    print(plotname)
    plt.savefig(plotname)
    plt.show()
    plt.close()
