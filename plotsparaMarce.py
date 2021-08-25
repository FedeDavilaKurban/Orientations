#%%
import sys
import numpy as np
from scipy import spatial
import scipy.stats
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import random
from config import writePath, units
import matplotlib.pyplot as plt
import seaborn as sns

def get_eta_bs(x,Nbs=1000):
    bs_eta = []
    for _ in range(Nbs):  
        bs_x = np.random.choice(x, size=len(x))
    
        n_perp = len(np.where(bs_x>0.)[0])
        n_prll = len(np.where(bs_x<0.)[0])
        bs_eta.append( n_perp / n_prll )

    eta = np.mean(bs_eta) 
    eta_std = np.std(bs_eta, ddof=1)

    return eta, eta_std, bs_eta
#%%
minradV, maxradV = 7.0, 0.0
rmin, rmax = 1.2, 1.3
sec = 123
vtype = 'a'

beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta'].data
N = len(beta)

# etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
#                     .format(minradV,maxradV,sec,vtype),\
#                     names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax','N'])
x = np.log10(beta)
eta, eta_std, bs_eta = get_eta_bs(x)
#%%
plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['font.size'] = 15

plt.hist(x,bins=30,density=True)
plt.vlines(0,0,.82,ls='--')
plt.xlabel(r'$\mathrm{log_{10}(\beta)}$')
plt.savefig('../plots/beta_paraMarce.pdf')
# %%
eta0 = 1./(np.sqrt(2)-1)

plt.hist(bs_eta,bins=30,density=True)
plt.vlines(eta0,0,18,ls='--')
plt.vlines(eta,0,18,ls='-')
plt.xlabel(r'$\eta=n(\beta>1)/n(\beta<1)$')
plt.savefig('../plots/eta_paraMarce.pdf')
# %%
etaran_std = np.sqrt(28.1421/N)
zeta = (eta-eta0)/etaran_std
zeta_err = eta_std/etaran_std

plt.errorbar(zeta,1,xerr=zeta_err,fmt='o-',ms=5,capsize=5)
plt.xlabel(r'$\zeta=(\eta-\eta_0)/\sigma_n$')
plt.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
plt.savefig('../plots/zeta_paraMarce.pdf')
# %%
fig, axs = plt.subplots(1, 1, constrained_layout=True, sharex=True, sharey=True)


etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
        .format(minradV,maxradV,sec,vtype),\
        names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax','N'])
etaTable.remove_row(0)
eta0 = 1./(np.sqrt(2)-1)
eta_ran_std = np.sqrt(28.1421/etaTable['N'])
y = (etaTable['eta']-eta0)/eta_ran_std
yerr = etaTable['eta_std']/eta_ran_std
x = (etaTable['rmin'].data+etaTable['rmax'].data)/2

axs.hlines(0,x[0],x[-1],linestyles=':')
axs.fill_between(x, -1, 1, alpha=.035, color='k')
axs.fill_between(x, -3, 3, alpha=.035, color='k')

axs.errorbar(x,y,yerr=yerr,label='High Spin',fmt='o-',capsize=3,ms=5)
axs.set_xticks(np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]))

axs.set_ylabel(r'$\zeta$')
axs.set_xlabel('R/Rv')

axs.legend(loc='upper left',framealpha=.4)
plt.savefig('../plots/highspin_paraMarce.pdf')
# %%
