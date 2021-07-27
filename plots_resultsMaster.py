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


plt.rcParams['figure.figsize'] = (15, 13)
plt.rcParams['font.size'] = 15

fig, axs = plt.subplots(4, 3, constrained_layout=True, sharex=True, sharey=True)

for ax_,secs in zip(axs[:2],([36,14],[123,456])):

    for ax, vtype in zip(ax_,['a','r','s']):

        for sec in secs:

            if sec==14: label = "Low Mass"
            if sec==36: label = "High Mass"
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
            
            ax.hlines(0,x[0],x[-1],linestyles=':')
            ax.fill_between(x, -1, 1, alpha=.035, color='k')
            ax.fill_between(x, -3, 3, alpha=.035, color='k')

            ax.errorbar(x,y,yerr=yerr,label=label,fmt='o-',capsize=3,ms=5)

            ax.set_xticks(np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]))

    ax_[0].legend(loc='upper left')
    ax_[0].set_ylabel(r'$\zeta$')

###########################################################
r1 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.9,1.0,1.1,1.2,1.3,1.4,1.5])

ax_ = axs[2]
print('Sigma5')
for sfilter in ['hi','lo']:
    print(sfilter)

    for sec in [0]:
        print(sec)
        eta = []
        eta_std = []

        eta_random_mean = []
        eta_random_std = []
        
        for vtype in ['a','r','s']:
            print(vtype)
            
            n_gal = []

            for rmin,rmax in zip(r1,r2):
                
                beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

                sigma5 = ascii.read('../data/sigma5/sigma5_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype),names=['sigma5'])['sigma5'].data

                p50 = np.percentile(sigma5,50)
                
                if sfilter=='hi': beta = beta[sigma5>p50]
                elif sfilter=='lo': beta = beta[sigma5<p50]

                x = np.log10(beta.data)

                eta_, eta_std_ = get_eta_bs(x)
    
                eta.append( eta_ )
                eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

                # Obtain mean and var of control samples
                N = len(beta)
                n_gal.append(N)

                eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
                eta_random_mean.append( np.mean(eta_random) )
                #eta_random_var.append( np.var(eta_random,ddof=1) )
                eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )

            print('N =',n_gal)

        nbins = len(r1)
        x = (r1+r2)/2

        ####################
        eta_random_mean = np.array(eta_random_mean)
        eta_random_std = np.array(eta_random_std)

        yran_mean_a = eta_random_mean[:nbins]
        yran_mean_r = eta_random_mean[nbins:2*nbins]
        yran_mean_s = eta_random_mean[-nbins:]

        yran_std_a = eta_random_std[:nbins]
        yran_std_r = eta_random_std[nbins:2*nbins]
        yran_std_s = eta_random_std[-nbins:]
        ####################

        ya = (eta[:nbins]-yran_mean_a)/yran_std_a
        yr = (eta[nbins:2*nbins]-yran_mean_r)/yran_std_r
        ys = (eta[-nbins:]-yran_mean_s)/yran_std_s

        ya_err = eta_std[:nbins]/yran_std_a
        yr_err = eta_std[nbins:2*nbins]/yran_std_r
        ys_err = eta_std[-nbins:]/yran_std_s

        ####################


        for ax, y, yerr in zip(ax_,(ya,yr,ys),(ya_err,yr_err,ys_err)):

            ax.hlines(0,x[0],x[-1],linestyles=':')

            ax.fill_between(x, -1, 1, alpha=.035, color='k')
            ax.fill_between(x, -3, 3, alpha=.035, color='k')


            if sfilter=='hi': label='High '+r'$\Sigma_5$'
            if sfilter=='lo': label='Low '+r'$\Sigma_5$'

            ax.errorbar(x,y,yerr=yerr,capsize=3,fmt='o-',ms=5,label=label)

ax_[0].legend(loc='upper left')
ax_[0].set_ylabel(r'$\zeta$')
###########################################################


###########################################################
ax_ = axs[3]
print('Vrad')
for vfilter in ['hi','lo']:
    print(vfilter)

    for sec in [0]:
        print(sec)
        eta = []
        eta_std = []

        eta_random_mean = []
        eta_random_std = []
        
        for vtype in ['a','r','s']:
            print(vtype)
            
            n_gal = []

            for rmin,rmax in zip(r1,r2):
                
                beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

                vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

                p50 = np.percentile(vrad,50)
                
                if vfilter=='hi': beta = beta[vrad>p50]
                elif vfilter=='lo': beta = beta[vrad<p50]

                x = np.log10(beta.data)

                eta_, eta_std_ = get_eta_bs(x)
    
                eta.append( eta_ )
                eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

                # Obtain mean and var of control samples
                N = len(beta)
                n_gal.append(N)

                eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
                eta_random_mean.append( np.mean(eta_random) )
                #eta_random_var.append( np.var(eta_random,ddof=1) )
                eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )

            print('N =',n_gal)

        nbins = len(r1)
        x = (r1+r2)/2

        ####################
        eta_random_mean = np.array(eta_random_mean)
        eta_random_std = np.array(eta_random_std)

        yran_mean_a = eta_random_mean[:nbins]
        yran_mean_r = eta_random_mean[nbins:2*nbins]
        yran_mean_s = eta_random_mean[-nbins:]

        yran_std_a = eta_random_std[:nbins]
        yran_std_r = eta_random_std[nbins:2*nbins]
        yran_std_s = eta_random_std[-nbins:]
        ####################

        ya = (eta[:nbins]-yran_mean_a)/yran_std_a
        yr = (eta[nbins:2*nbins]-yran_mean_r)/yran_std_r
        ys = (eta[-nbins:]-yran_mean_s)/yran_std_s

        ya_err = eta_std[:nbins]/yran_std_a
        yr_err = eta_std[nbins:2*nbins]/yran_std_r
        ys_err = eta_std[-nbins:]/yran_std_s

        ####################


        for ax, y, yerr in zip(ax_,(ya,yr,ys),(ya_err,yr_err,ys_err)):

            ax.hlines(0,x[0],x[-1],linestyles=':')

            ax.fill_between(x, -1, 1, alpha=.035, color='k')
            ax.fill_between(x, -3, 3, alpha=.035, color='k')

            if vfilter=='hi': label='High '+r'$\mathrm{V_{rad}}$'
            if vfilter=='lo': label='Low '+r'$\mathrm{V_{rad}}$'

            ax.errorbar(x,y,yerr=yerr,capsize=3,fmt='o-',ms=5,label=label)

ax_[0].legend(loc='upper left')
ax_[0].set_ylabel(r'$\zeta$')

###########################################################
for ax in axs[3]:
    ax.set_xlabel('R/Rv')

axs[0,0].set_title('All Voids')
axs[0,1].set_title('R-Voids')
axs[0,2].set_title('S-Voids')
plt.tight_layout()
plt.savefig('../plots/allvoidsallgalaxies_results.jpg')
#%%

