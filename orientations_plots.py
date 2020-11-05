
#%%
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table

allsections = [1,2,3,123,12,23,45,56,4,5,6,456]
rmin, rmax = .7, 1.2
for rmin in [0.5,.6,.8,.9]:
    print(rmin)
    for rmax in [1.2,]:
        print(rmax)
        fig, axes = plt.subplots(3, 4,figsize=(8,6),sharex=True,sharey=True)

        axs = []
        for i in range(3):
            for j in range(4):
                axs.append( axes[i,j] )

        for sec,ax in zip(allsections,axs):

            data = ascii.read('../data/fits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['cos','ecdf','y','yfit','d_yfit'])
            random = ascii.read('../data/randomfits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['xmean','ymean','ystd'])

            cos = data['cos']
            #ecdf = data['ecdf']
            y = data['y']
            yfit = data['yfit']
            d_yfit = data['d_yfit']

            xmean = random['xmean']
            ymean = random['ymean']
            ystd = random['ystd']

            ax.set_title('Sec={}'.format(sec))
            ax.step(cos,y)
            ax.plot(cos,yfit,'r--')
            ax.fill_between(xmean,ymean-ystd,ymean+ystd,color='k',alpha=.6,label=r'$1\sigma$')
            ax.fill_between(xmean,ymean-2*ystd,ymean+2*ystd,color='k',alpha=.4,label=r'$2\sigma$')
            #ax.legend()

        axes[0,0].legend(loc=2)
        plt.tight_layout()
        plt.savefig('../plots/cosines_allsections_rmin{}_rmax{}.png'.format(rmin,rmax),dpi=300)
        plt.show()
# %%
