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
# %%
"""

"""
plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15

minradV=7.
maxradV=0.

# rinner_i = np.float64(0.8)
# rinner_f = np.float64(1.5)
# rstep = np.float64(0.1)
# r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
# r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

r1 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.9,1.0,1.1,1.2,1.3,1.4,1.5])

fig, axs = plt.subplots(1, 1, constrained_layout=True, sharex=True, sharey=True)

for sec in [3]:
    print(sec)
    #eta = []
    #eta_std = []

    eta_random_mean = []
    eta_random_std = []
    n_gal = []
    
    for vtype in ['s']:
        print(vtype)
        
        for rmin,rmax in zip(r1,r2):
            
            beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

            p50 = np.percentile(sigma5,50)
            
            #x = np.log10(beta.data)

            # Obtain mean and var of control samples
            N = len(beta)
            print(rmin,rmax,N)
            n_gal.append(N)

            eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
            eta_random_mean.append( np.mean(eta_random) )
            #eta_random_var.append( np.var(eta_random,ddof=1) )
            eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )


    nbins = len(r1)
    x = (r1+r2)/2

    eta_theo_mean = 1/(np.sqrt(2)-1)

    ####################
    n_gal = np.array(n_gal)
    eta_random_mean = np.array(eta_random_mean)
    eta_random_std = np.array(eta_random_std)

    yran_mean_r = eta_random_mean[:nbins]
    yran_mean_s = eta_random_mean[-nbins:]

    yran_std_r = eta_random_std[:nbins]
    yran_std_s = eta_random_std[-nbins:]
    ####################

    #yr = (eta[:nbins]-yran_mean_r)/yran_std_r
    #ys = (eta[-nbins:]-yran_mean_s)/yran_std_s

    #yr_err = eta_std[:nbins]/yran_std_r
    #ys_err = eta_std[-nbins:]/yran_std_s

    ####################

    if sec==3: 
        title='H-H Galaxies'

        axs.set_title(title)

        axs.set_xlabel('R/Rv')



        for yran_mean, yran_err in zip((yran_mean_r,yran_mean_s),\
                                        (yran_std_r,yran_std_s)):

            # Theoretical n1/n2 value for random spins
            axs.hlines(eta_theo_mean,x[0],x[-1],linestyles=':')
            axs.plot(x,yran_mean,c='k',alpha=.7)

            axs.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.1, color='k')
            axs.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.1, color='k')
            axs.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.1, color='k')


            #ax.errorbar(x,y,yerr=yerr,label=label)
            #print(y)
            
            axs.set_ylabel(r'$\eta$')

            axs.set_ylim([1.5,3.5])

        axs.scatter(x,eta_theo_mean+np.sqrt(28.1421/n_gal))

        axs.set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs.set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        #axs[2].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

        #plt.savefig('../plots/sigma5/Eta_{}Sigma5_sec{}.jpg'.format(sec))
# %%
