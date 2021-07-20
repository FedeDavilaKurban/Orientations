#%%
"""
Calculate Beta (J_perp/J_prll) for different bins in Radius

With this "_writing" version I want to write the betas first, 
and do calculations and plots later

Plots:
- Eta vs R for a, r, and s voids
- B vs Mass
- B vs Vrad&Vtra
- B vs Gas Mass
- Eta vs Low/High Mass
- Eta vs Low/High Vrad&Vtra
- Eta vs Low/High Gas Mass
"""

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
from mpl_toolkits.mplot3d import axes3d

# %%

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


########################################################
########################################################
#%%

# rinner_i = np.float64(0.7)
# rinner_f = np.float64(1.5)
# rstep = np.float64(0.1)
# r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
# r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)[-1]

r1 = np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
#%%%%%%%%
#############################################################################
#############################################################################
########################PLOTS################################################
#############################################################################
#############################################################################
#############################################################################
# %%
"""
B vs Mass
"""

minradV = 7.
maxradV = 0.
sec = 3
rmin, rmax = 1.1,1.2

plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    mass = ascii.read('../data/mass/-1/mass_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['mass'].data

    x = np.log10(mass)
    y = np.log10(beta)
    m,b,rvalue,pvalue,std = scipy.stats.linregress(x,y) 

    ax.scatter(x, y, s=.5)

    ax.plot(x,m*x+b,ls=':',color='k',label='r_lr={:.2f}'.format(rvalue))
    ax.legend()

    ax.set_title('section {} - vtype {}'.format(sec,vtype))
    ax.set_ylabel(r'$log_{10}(\beta)$')

axs[2].set_xlabel(r'$log_{10}(M[10^{10}M_{\odot}])$')
plt.savefig('../plots/BetaEta_vs_MassVel/BvsMass.jpg')

#%%
"""
B vs Gmass
"""
import matplotlib.colors as colors

minradV = 7.
maxradV = 0.
sec = 3
rmin, rmax = 1.1,1.2

plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    gmass = ascii.read('../data/mass/-1/mass_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['gmass'].data

    # smass = ascii.read('../data/mass/smass_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
    #                 .format(minradV,maxradV,rmin,rmax,sec,vtype))['smass'].data

    x = np.log10(beta[np.where(gmass!=0.)[0]])
    y = np.log10(gmass[np.where(gmass!=0.)[0]])

    ax.scatter(x,y,s=.5)

    m,b,rvalue,pvalue,std = scipy.stats.linregress(x,y) 
    ax.plot(x,m*x+b,ls=':',color='k',label='r_lr={:.2f}'.format(rvalue))
    ax.set_title('section {} - vtype {}'.format(sec,vtype))
    ax.set_ylabel(r'$log_{10}(Gass Mass)$')
    ax.legend()

axs[2].set_xlabel(r'$log_{10}(\beta)$')

plt.savefig('../plots/BetaEta_vs_MassVel/BvsGasMass.jpg')


#%%
"""
Eta vs Low/High Mass
"""
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


plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    mass = ascii.read('../data/mass/-1/mass_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['mass'].data

    logmass = np.log10(mass)
    p50 = np.percentile(logmass,50)

    lmass1 = logmass[logmass<p50]
    lmass2 = logmass[logmass>p50]

    beta1 = beta[logmass<p50]
    beta2 = beta[logmass>p50]

    logbeta1 = np.log10(beta1)
    logbeta2 = np.log10(beta2)

    eta1, eta1_std = get_eta_bs(logbeta1)
    eta2, eta2_std = get_eta_bs(logbeta2)

    ax.errorbar(1,eta1,eta1_std,fmt='o',capsize=5)
    ax.errorbar(2,eta2,eta2_std,fmt='o',capsize=5)
    ax.set_title('section {} - vtype {}'.format(sec,vtype))
    ax.set_ylabel(r'$\eta$')

axs[2].set_xlim([0.,3.])
axs[2].set_xlabel('Mass')
axs[2].set_xticks([1,2])
axs[2].set_xticklabels(['Low Mass','High Mass'])

plt.savefig('../plots/BetaEta_vs_MassVel/EtavsLoHiMass.jpg')


# %%

"""
Eta vs Low/High GasMass
"""
minradV=7.
maxradV=0.
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

rmin, rmax = 1.1,1.2

plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    gmass = ascii.read('../data/mass/-1/mass_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['gmass'].data

    #loggmass = np.log10(gmass)
    loggmass = gmass
    p50 = np.percentile(loggmass,50)

    lgmass1 = loggmass[loggmass<p50]
    lgmass2 = loggmass[loggmass>p50]

    beta1 = beta[loggmass<p50]
    beta2 = beta[loggmass>p50]

    logbeta1 = np.log10(beta1)
    logbeta2 = np.log10(beta2)

    eta1, eta1_std = get_eta_bs(logbeta1)
    eta2, eta2_std = get_eta_bs(logbeta2)

    ax.errorbar(1,eta1,eta1_std,fmt='o',capsize=5)
    ax.errorbar(2,eta2,eta2_std,fmt='o',capsize=5)
    ax.set_title('section {} - vtype {}'.format(sec,vtype))
    ax.set_ylabel(r'$\eta$')

axs[2].set_xlim([0.,3.])
axs[2].set_xlabel('GasMass')
axs[2].set_xticks([1,2])
axs[2].set_xticklabels(['Low Gas Mass','High Gas Mass'])

plt.savefig('../plots/BetaEta_vs_MassVel/EtavsLoHiGasMass.jpg')

#%%
"""
Vrad vs Vtra para gxs en la sec 3, voids R, 0.9-1.2Rv

Hay correlacion entre la vrad y vtra en esta muestra de galaxias?
"""
plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['font.size'] = 15

sec = 3
vtype = 'a'
rmin, rmax = 1.1,1.2

vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
    .format(minradV,maxradV,rmin,rmax,sec,vtype))

plt.scatter(vel['vtra'],vel['vrad'],s=0.5)

# m,b,rvalue,pvalue,std = scipy.stats.linregress(vel['vtra'],vel['vrad']) 
# plt.plot(vel['vtra'],m*vel['vtra']+b,ls='-',color='k',label='All voids, r_lr={:.2f}'.format(rvalue))

# vtype = 'r'
# vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
#     .format(minradV,maxradV,rmin,rmax,sec,vtype))
# m,b,rvalue,pvalue,std = scipy.stats.linregress(vel['vtra'],vel['vrad']) 
# plt.plot(vel['vtra'],m*vel['vtra']+b,ls=':',color='k',label='R voids, r_lr={:.2f}'.format(rvalue))

# vtype = 's'
# vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
#     .format(minradV,maxradV,rmin,rmax,sec,vtype))
# m,b,rvalue,pvalue,std = scipy.stats.linregress(vel['vtra'],vel['vrad']) 
# plt.plot(vel['vtra'],m*vel['vtra']+b,ls='-.',color='k',label='S voids, r_lr={:.2f}'.format(rvalue))

plt.legend()
plt.xlabel('Vtra')
plt.ylabel('Vrad')
plt.xscale('log')
plt.yscale('log')
plt.savefig('../plots/BetaEta_vs_MassVel/Vtra_vs_Vrad.jpg')
