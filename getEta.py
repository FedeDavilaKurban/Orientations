"""
Calculate and write Eta
"""
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

r1 = np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
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


for sec in [0,1,2,3,4,5,6,123,456]:
    print(sec)

    
    for vtype in ['a','r','s']:
        print(vtype)
        
        eta = []
        eta_std = []

        eta_random_mean = []
        eta_random_var = []
        eta_random_std = []

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

        ascii.write(np.column_stack([eta,eta_std,eta_random_mean,eta_random_std,r1,r2]),\
            '../data/eta/eta_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype),\
            names=['eta','eta_std','eta_random_mean','eta_random_std','rmin','rmax'],\
            overwrite=True)
# %%
