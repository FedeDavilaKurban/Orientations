
#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table

allsections = [1,2,3,123,12,23,45,56,4,5,6,456]
rmin, rmax = .7, 1.2
sec=6
for sec in [4,5,45,56,456]:
    fig, axes = plt.subplots(2, 3,figsize=(8,6),sharex=True,sharey=True)
    axs = []
    for i in range(2):
        for j in range(3):
            axs.append( axes[i,j] )
    for s5,ax in zip([50,100,200,300,400,0],axs):

        if s5!=0:
            data = ascii.read('../data/fits_sec{}_rmin{}_rmax{}_s5:{}'.format(sec,rmin,rmax,s5),names=['cos','ecdf','y','yfit','d_yfit'])
            random = ascii.read('../data/randomfits_sec{}_rmin{}_rmax{}_s5:{}'.format(sec,rmin,rmax,s5),names=['xmean','ymean','ystd'])

        else:
            data = ascii.read('../data/fits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax,s5),names=['cos','ecdf','y','yfit','d_yfit'])
            random = ascii.read('../data/randomfits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax,s5),names=['xmean','ymean','ystd'])


        cos = data['cos']
        #ecdf = data['ecdf']
        y = data['y']
        yfit = data['yfit']
        d_yfit = data['d_yfit']

        xmean = random['xmean']
        ymean = random['ymean']
        ystd = random['ystd']

        ax.set_title('Sigma5<={}'.format(s5))
        ax.step(cos,y)
        ax.plot(cos,yfit,'r--')
        ax.fill_between(xmean,ymean-ystd,ymean+ystd,color='k',alpha=.6,label=r'$1\sigma$')
        ax.fill_between(xmean,ymean-2*ystd,ymean+2*ystd,color='k',alpha=.4,label=r'$2\sigma$')
        #ax.legend()

    axes[0,0].legend(loc=2)
    plt.tight_layout()
    plt.savefig('../plots/cosines_sec{}_rmin{}_rmax{}_s5.png'.format(sec,rmin,rmax),dpi=300)
    plt.show()
# %%
