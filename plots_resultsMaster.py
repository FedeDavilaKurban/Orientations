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

def get_random_vec(R, N):

    phi = 2*np.pi*np.random.random(N)
    costheta = 1-2*np.random.random(N)
    u = np.random.random(N)

    theta = np.arccos( costheta )
    r = R * np.cbrt( u )

    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    
    return x,y,z

def get_beta_random(N_beta):

    sx, sy, sz = get_random_vec(1.,N_beta)

    s_randvec = np.column_stack((sx,sy,sz))

    s_norms = np.linalg.norm(s_randvec, axis=1)

    sx /= s_norms
    sy /= s_norms
    sz /= s_norms

    sp = np.sqrt(sx**2+sy**2)

    b = sp/np.abs(sz)

    return b

def get_eta_random(N_eta,N_beta):
    
    eta_ran = np.zeros(N_eta)

    for i in range(N_eta):

        b = get_beta_random(N_beta)
        eta_ran[i] = len(b[b>1]) / len(b[b<1])

    return eta_ran
    
def get_eta_bs(x,Nbs=1000):
    bs_eta = []
    for _ in range(Nbs):  
        bs_x = np.random.choice(x, size=len(x))
    
        n_perp = len(np.where(bs_x>0.)[0])
        n_prll = len(np.where(bs_x<0.)[0])
        bs_eta.append( n_perp / n_prll )

    eta = np.mean(bs_eta) 
    eta_std = np.std(bs_eta, ddof=1)

    return eta, eta_std



#%%
minradV = 7.
maxradV = 0.

lowMcut = -1

colors = sns.color_palette()
plt.rcParams['figure.figsize'] = (20, 15)
plt.rcParams['font.size'] = 18

fig, axs = plt.subplots(4, 3, constrained_layout=True, sharex=True, sharey=True)
axs[0,0].set_ylim([-3,7])

r1 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.9,1.0,1.1,1.2,1.3,1.4,1.5])
x = (r1+r2)/2

for ax_ in axs:
    for ax in ax_:
            ax.hlines(0,r1[0],r2[-1],linestyles=':')
            ax.fill_between([0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
            ax.fill_between([0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')

for ax_,secs in zip(axs[:2],([25,36,14],[123,456])):

    for ax, vtype in zip(ax_,['a','r','s']):

        for sec in secs:

            if sec==14: label = "Low Mass"
            if sec==36: label = "High Mass"
            if sec==25: label = "Intermediate Mass"
            if sec==456: label = "Low Spin"
            if sec==123: label = "High Spin"

            etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,sec,vtype),\
                    names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax','N'])
            etaTable.remove_row(0)

            yran_mean = etaTable['eta_random_mean'].data
            yran_err = etaTable['eta_random_std'].data

            y = (etaTable['eta'].data-yran_mean)/yran_err
            yerr = etaTable['eta_std'].data/yran_err

            x = (etaTable['rmin'].data+etaTable['rmax'].data)/2

            if sec==25:
                ax.errorbar(x,y,yerr=yerr,label=label,fmt='--',capsize=3,ms=5,\
                    color='k')
            else: ax.errorbar(x,y,yerr=yerr,label=label,fmt='o-',capsize=3,ms=5)

            ax.set_xticks(np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]))

    ax_[0].legend(loc='upper left',ncol=2,framealpha=.4)
    ax_[0].set_ylabel(r'$\zeta$')

###########################################################

ax_ = axs[2]
print('Sigma5')
for sfilter in ['hi','lo']:

    for sec in [0]:

        for vtype,ax in zip(['a','r','s'],ax_):

            etaT = ascii.read('../data/eta/eta_sfilter{}_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                                    .format(sfilter,minradV,maxradV,sec,vtype))
        
            eta0 = 1./(np.sqrt(2)-1)
            eta_ran_std = np.sqrt(28.1421/etaT['N'].data)
            y = (etaT['eta'].data-eta0)/eta_ran_std
            yerr = etaT['eta_std'].data/eta_ran_std

            if sfilter=='hi': label='High '+r'$\Sigma_5$'
            if sfilter=='lo': label='Low '+r'$\Sigma_5$'

            x = (etaTable['rmin'].data+etaTable['rmax'].data)/2
            print(x)
            ax.errorbar(x,y,yerr=yerr,capsize=3,fmt='o-',ms=5,label=label)

ax_[0].legend(loc='upper left',ncol=2,framealpha=.4)
ax_[0].set_ylabel(r'$\zeta$')
###########################################################


###########################################################
ax_ = axs[3]
print('Vrad')
for vfilter in ['hi','lo']:

    for sec in [0]:

        for vtype,ax in zip(['a','r','s'],ax_):
            
            etaT = ascii.read('../data/eta/eta_vfilter{}_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
                                    .format(vfilter,minradV,maxradV,sec,vtype))
        
            eta0 = 1./(np.sqrt(2)-1)
            eta_ran_std = np.sqrt(28.1421/etaT['N'].data)
            y = (etaT['eta'].data-eta0)/eta_ran_std
            yerr = etaT['eta_std'].data/eta_ran_std

            x = (etaTable['rmin'].data+etaTable['rmax'].data)/2
            
            if vfilter=='hi': label='High '+r'$\mathrm{V_{rad}}$'
            if vfilter=='lo': label='Low '+r'$\mathrm{V_{rad}}$'

            ax.errorbar(x,y,yerr=yerr,capsize=3,fmt='o-',ms=5,label=label)

ax_[0].legend(loc='upper left',ncol=2,framealpha=.4)
ax_[0].set_ylabel(r'$\zeta$')

###########################################################
for ax in axs[3]:
    ax.set_xlabel(r'$\mathrm{r/R_v}$')

axs[0,0].set_ylim([-4,7])
axs[0,0].set_xlim([.8,1.5])
axs[0,0].set_title('All Voids')
axs[0,1].set_title('R-Voids')
axs[0,2].set_title('S-Voids')
#plt.tight_layout()
#plt.savefig('../plots/allvoidsallgalaxies_results.jpg')
plt.savefig('../plots/allvoidsallgalaxies_results.pdf')


#%%
plt.rcParams['figure.figsize'] = (8, 12)
plt.rcParams['font.size'] = 18

minradV, maxradV = 7.0, 0.0
#vtype = 'a'
r1 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.9,1.0,1.1,1.2,1.3,1.4,1.5])
x = (r1+r2)/2

fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

for vtype,ax in zip(['a','r','s'],axs):
    
    ax.hlines(0,r1[0],r2[-1],linestyles=':')
    ax.fill_between([0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
    ax.fill_between([0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')


    for sec in [3,1]:

        eta=[]
        eta_std=[]
        eta_random_mean=[]
        eta_random_std=[]
        n_gal=[]

        for rmin,rmax in zip(r1,r2):

            beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

            vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

            sigma5 = ascii.read('../data/sigma5/sigma5_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                            .format(minradV,maxradV,rmin,rmax,sec,vtype),names=['sigma5'])['sigma5'].data

            vp50 = np.percentile(vrad,50)
            sp50 = np.percentile(sigma5,50)

            if sec==3:
                beta = beta[vrad<vp50]
                #sigma5 = sigma5[vrad<vp50]

                #beta = beta[sigma5<sp50]
                label = 'Low {}, H-H Galaxies'.format(r'$\mathrm{V_{rad}}$')

            if sec==1:
                beta = beta[vrad>vp50]
                #sigma5 = sigma5[vrad>vp50]

                #beta = beta[sigma5>sp50]
                label = 'High {}, L-L Galaxies'.format(r'$\mathrm{V_{rad}}$')

            logb = np.log10(beta.data)

            eta_, eta_std_ = get_eta_bs(logb)

            eta.append( eta_ )
            eta_std.append( eta_std_ )

            N = len(beta)
            n_gal.append(N)

            eta_random_mean.append( 1./(np.sqrt(2)-1) )
            eta_random_std.append( np.sqrt(28.1421/N) )

        eta=np.array(eta)
        eta_std=np.array(eta_std)
        eta_random_mean=np.array(eta_random_mean)
        eta_random_std=np.array(eta_random_std)

        y = (eta-eta_random_mean)/eta_random_std
        yerr = eta_std/eta_random_std

        ax.errorbar(x,y,yerr=yerr,capsize=3,fmt='o-',ms=5,label=label)

    ax.set_xlabel(r'$\mathrm{r/R_v}$')
    ax.legend(loc='upper left',framealpha=.4,ncol=1)
    ax.set_ylabel(r'$\zeta$')

axs[0].set_ylim([-4,8])
axs[0].set_xlim([.8,1.5])
axs[0].set_title('All Voids')
axs[1].set_title('R Voids')
axs[2].set_title('S Voids')

#plt.savefig('../plots/bestsignal.png')

# %%
