#%%
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from astropy.io import ascii
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15
import seaborn as sns
#%%

minradV=7.0
maxradV=0.0
r1 = [0.8,0.9,1.0,1.1,1.2,1.3,1.4]
r2 = [0.9,1.0,1.1,1.2,1.3,1.4,1.5]

sec = 3

fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=True)

for vtype, ax in zip(['r','s'],axs):
    for rmin, rmax in zip(r1,r2):

        vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data
        sigma5 = ascii.read('../data/sigma5/sigma5_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                            .format(minradV,maxradV,rmin,rmax,sec,vtype),names=['sigma5'])['sigma5'].data

        ax.scatter(vrad,np.log10(sigma5),color='k')

        ax.set_ylabel(r'$log_{10}\Sigma_5$')

axs[1].set_xlabel('Vrad')

plt.savefig('../plots/vel_vs_sigma.png')
# %%
