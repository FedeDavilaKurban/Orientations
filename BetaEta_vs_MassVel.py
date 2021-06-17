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
    


########################################################
########################################################
#%%

rinner_i = np.float64(0.6)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)


#%%%%%%%%
#############################################################################
#############################################################################
########################PLOTS################################################
#############################################################################
#############################################################################
#############################################################################
"""
Eta vs R for a, r, and s voids
"""

minradV = 7.
maxradV = 0.

lowMcut = -1

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

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

for sec in [3,4,5,6,123,456]:
    print(sec)
    eta = []
    eta_std = []

    eta_random_mean = []
    eta_random_var = []
    eta_random_std = []
    
    n_gal = []
    for vtype in ['a','r','s']:
        print(vtype)
        
        for rmin,rmax in zip(r1,r2):
            
            beta = ascii.read('../data/beta/{}/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                .format(str(lowMcut),minradV,maxradV,rmin,rmax,sec,vtype))['beta']

            N = len(beta)
            #print(N)

            x = np.log10(beta.data)

            # Obtain mean and var of eta with Bootstrap
            # bs_eta = []
            # for _ in range(1000):  #so B=1000
            #     bs_x = np.random.choice(x, size=len(x))
            
            #     n_perp = len(np.where(bs_x>0.)[0])
            #     n_prll = len(np.where(bs_x<0.)[0])
            #     bs_eta.append( n_perp / n_prll )
            eta_, eta_std_ = get_eta_bs(x)

            # eta.append( np.mean(bs_eta) )
            # eta_std.append( np.std(bs_eta, ddof=1))#/np.sqrt(len(bs_eta)) )
 
            eta.append( eta_ )
            eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

            # Obtain mean and var of control samples

        
            eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
            eta_random_mean.append( np.mean(eta_random) )
            eta_random_var.append( np.var(eta_random,ddof=1) )
            eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )



    nbins = len(r1)
    x = r1


    ya = eta[:nbins]
    yr = eta[nbins:2*nbins]
    ys = eta[-nbins:]

    ya_err = eta_std[:nbins]
    yr_err = eta_std[nbins:2*nbins]
    ys_err = eta_std[-nbins:]

    ####################
    eta_random_mean = np.array(eta_random_mean)
    eta_random_var = np.array(eta_random_var)
    eta_random_std = np.array(eta_random_std)

    yran_mean_a = eta_random_mean[:nbins]
    yran_mean_r = eta_random_mean[nbins:2*nbins]
    yran_mean_s = eta_random_mean[-nbins:]

    yran_var_a = eta_random_var[:nbins]
    yran_var_r = eta_random_var[nbins:2*nbins]
    yran_var_s = eta_random_var[-nbins:]

    yran_std_a = eta_random_std[:nbins]
    yran_std_r = eta_random_std[nbins:2*nbins]
    yran_std_s = eta_random_std[-nbins:]
    ####################

    fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=False)

    axs[0].set_title('section {}'.format(sec))


    for ax, y, yerr, yran_mean, yran_err, label, in zip(axs,\
                                    (ya,yr,ys),\
                                    (ya_err,yr_err,ys_err),\
                                    (yran_mean_a,yran_mean_r,yran_mean_s),\
                                    (yran_std_a,yran_std_r,yran_std_s),('All','Rising','Shell')):

        # Theoretical n1/n2 value for random spins
        ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
        ax.plot(x,yran_mean,c='k',alpha=.7)

        ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.3, color='k')

        ax.errorbar(x,y,yerr=yerr,label=label)

        ax.legend()

        # plt.ylabel(r'$n_1/n_2$')
        # plt.xlabel('R/Rv')
        # plt.legend()

    x_ticks_labels = []
    for i in range(nbins):
        x_ticks_labels.append( '{:.1f}-{:.1f}'.format(r1[i],r2[i]) )
    axs[2].set_xticks(x)
    axs[2].set_xticklabels(x_ticks_labels)

    plt.savefig('../plots/BetaEta_vs_MassVel/Eta_allvtypes_sec{}.jpg'.format(sec))
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
# %%
"""
B vs Vrad & Vtra
"""

minradV = 7.
maxradV = 0.
sec = 3
rmin, rmax = 1.1,1.2

plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 2, constrained_layout=True, sharex=False, sharey=True)

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))

    for v, ax_i, in zip([vel['vrad'],vel['vtra']],ax):

        m,b,rvalue,pvalue,std = scipy.stats.linregress(v,np.log10(beta)) 

        ax_i.scatter(v,np.log10(beta),s=.5)
        ax_i.set_title('section {} - vtype {}'.format(sec,vtype))

        ax_i.plot(v,m*v+b,ls=':',color='k',label='r_lr={:.2f}'.format(rvalue))
        ax_i.legend()

    ax[0].set_ylabel(r'$log_{10}(\beta)$')

axs[2,0].set_xlabel(r'$V_{rad}$')
axs[2,1].set_xlabel(r'$V_{tra}$')

plt.savefig('../plots/BetaEta_vs_MassVel/BvsVel.jpg')

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
Eta vs Low/High Vrad and Vtra
"""
minradV=7.
maxradV=0.
sec=3
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
fig, axs = plt.subplots(3, 2, constrained_layout=True, sharex=True, sharey=True)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))

    for v, color, fmt, ax_i, in zip([vel['vrad'],vel['vtra']],colors[:3],['o','x'],ax):

        p30 = np.percentile(v,30)
        p50 = np.percentile(v,50)
        p70 = np.percentile(v,70)

        lomask = v<p50
        himask = v>p50

        v1 = v[lomask]
        v2 = v[himask]
        #print(np.mean(v1),np.mean(v2))

        beta1 = beta[lomask]
        beta2 = beta[himask]

        logbeta1 = np.log10(beta1)
        logbeta2 = np.log10(beta2)

        eta1, eta1_std = get_eta_bs(logbeta1)
        eta2, eta2_std = get_eta_bs(logbeta2)
        print(eta1,eta2)
        ax_i.errorbar(1,eta1,eta1_std,fmt=fmt, color=color, capsize=5)
        ax_i.errorbar(2,eta2,eta2_std,fmt=fmt, color=color, capsize=5)
        ax_i.hlines(2.41,.5,2.5,ls=':',color='k')
        ax_i.set_title('section {} - vtype {}'.format(sec,vtype))

    ax[0].set_ylabel(r'$\eta$')

axs[2,1].set_xlim([0.,3.])
axs[2,0].set_xlabel('V_rad')
axs[2,1].set_xlabel('V_tra')
axs[2,1].set_xticks([1,2])
axs[2,1].set_xticklabels(['Low','High'])

plt.savefig('../plots/BetaEta_vs_MassVel/EtavsLoHiVel.jpg')

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

# %%
"""
Eta vs R for LowVrad
"""
rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

for sec in [3]:
    print(sec)
    eta = []
    eta_std = []

    eta_random_mean = []
    eta_random_var = []
    eta_random_std = []
    
    n_gal = []
    for vtype in ['a','r','s']:
        print(vtype)
        
        for rmin,rmax in zip(r1,r2):
            
            beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

            vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

            p50 = np.percentile(vrad,50)
            beta = beta[vrad<p50]

            x = np.log10(beta.data)

            eta_, eta_std_ = get_eta_bs(x)
 
            eta.append( eta_ )
            eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

            # Obtain mean and var of control samples
            N = len(beta)
        
            eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
            eta_random_mean.append( np.mean(eta_random) )
            eta_random_var.append( np.var(eta_random,ddof=1) )
            eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )



    nbins = len(r1)
    x = r1


    ya = eta[:nbins]
    yr = eta[nbins:2*nbins]
    ys = eta[-nbins:]

    ya_err = eta_std[:nbins]
    yr_err = eta_std[nbins:2*nbins]
    ys_err = eta_std[-nbins:]

    ####################
    eta_random_mean = np.array(eta_random_mean)
    eta_random_var = np.array(eta_random_var)
    eta_random_std = np.array(eta_random_std)

    yran_mean_a = eta_random_mean[:nbins]
    yran_mean_r = eta_random_mean[nbins:2*nbins]
    yran_mean_s = eta_random_mean[-nbins:]

    yran_var_a = eta_random_var[:nbins]
    yran_var_r = eta_random_var[nbins:2*nbins]
    yran_var_s = eta_random_var[-nbins:]

    yran_std_a = eta_random_std[:nbins]
    yran_std_r = eta_random_std[nbins:2*nbins]
    yran_std_s = eta_random_std[-nbins:]
    ####################

    fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=False)

    axs[0].set_title('section {}'.format(sec))


    for ax, y, yerr, yran_mean, yran_err, label, in zip(axs,\
                                    (ya,yr,ys),\
                                    (ya_err,yr_err,ys_err),\
                                    (yran_mean_a,yran_mean_r,yran_mean_s),\
                                    (yran_std_a,yran_std_r,yran_std_s),('All','Rising','Shell')):

        # Theoretical n1/n2 value for random spins
        ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
        ax.plot(x,yran_mean,c='k',alpha=.7)

        ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.3, color='k')

        ax.errorbar(x,y,yerr=yerr,label=label)
        print(y)
        ax.legend()

        # plt.ylabel(r'$n_1/n_2$')
        # plt.xlabel('R/Rv')
        # plt.legend()

    x_ticks_labels = []
    for i in range(nbins):
        x_ticks_labels.append( '{:.1f}-{:.1f}'.format(r1[i],r2[i]) )
    axs[2].set_xticks(x)
    axs[2].set_xticklabels(x_ticks_labels)

    #plt.savefig('../plots/BetaEta_vs_MassVel/Eta_LoVrad_allvtypes_sec{}_.jpg'.format(sec))
# %%
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
# %%

"""
Vrad vs R
"""
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15
import seaborn as sns

vtype = 's'
sec = 3
colors = sns.color_palette()

for vtype, c in zip(['a','r','s'],colors[:3]):
    vrad_mean=[]
    vtra_mean=[]
    r_mean=[]

    for rmin,rmax in zip(r1,r2):
        vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                        .format(minradV,maxradV,rmin,rmax,sec,vtype))
        vrad_mean.append( np.mean(vel['vrad'].data))
        vtra_mean.append(np.mean(vel['vtra'].data))
        r_mean.append( (rmin+rmax)/2)

    plt.plot(r_mean,vrad_mean,marker='o',color=c,label='Vrad, {} voids'.format(vtype))
    plt.plot(r_mean,vtra_mean,ls='--',marker='x',color=c,label='Vtra, {} voids'.format(vtype))



plt.legend(ncol=3)
plt.xlabel('Rv')
plt.ylabel('Vrad, Vtra')
plt.savefig('../plots/BetaEta_vs_MassVel/VradVtra_vs_R_allvtypes.jpg')
# %%
