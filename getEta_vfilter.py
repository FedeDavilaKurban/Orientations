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

minradV, maxradV = 7.0, 0.0
r1 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.9,1.0,1.1,1.2,1.3,1.4,1.5])
r1 = [.9]
r2 = [1.4]

for vfilter in ['hi','lo']:
    print(vfilter)

    for sec in [0,123,456,1,3,4,6,14,36]:
        print(sec)
        
        for vtype in ['a','r','s']:
            print(vtype)

            eta = []
            eta_std = []           
            n_gal = []

            for rmin,rmax in zip(r1,r2):
                
                beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']
                
                vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

                p50 = np.percentile(vrad,50)
                
                if vfilter=='hi': 
                    #p67 = np.percentile(vrad,67)
                    beta = beta[vrad>p50]
                elif vfilter=='lo': 
                    #p33 = np.percentile(vrad,33)
                    beta = beta[vrad<p50]

                n_gal.append(len(beta))

                x = np.log10(beta.data)

                eta_, eta_std_ = get_eta_bs(x)
    
                eta.append( eta_ )
                eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

            ascii.write(np.column_stack([eta,eta_std,r1,r2,n_gal]),\
            '../data/eta/eta_vfilter{}_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(vfilter,minradV,maxradV,sec,vtype),\
            names=['eta','eta_std','rmin','rmax','N'],\
            overwrite=True)


# %%
