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
    #x should be log10(beta)
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


#%%%%%%%%
"""
Eta vs R for a, r, and s voids
"""
r1 = np.array([0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.9,1.0,1.1,1.2,1.3,1.4,1.5])

forTable=False

if forTable:
    r1 = [.9]
    r2 = [1.4]

minradV = 7.
maxradV = 0.

lowMcut = -1

for vfilter in ['lo']:

    for sec in [3]:
        print('sec:',sec)

        for vtype in ['a','r']:
            print('vtype:',vtype)

            voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)

            # eta = []
            # eta_std = []

            # eta_random_mean = []
            # eta_random_var = []
            # eta_random_std = []
            
            # n_gal=[]
            for rmin,rmax in zip(r1,r2):

                zeta_jk = []
                #n_gal_jk = []
                for jk in range(len(voids)):

                    etaTable = ascii.read('../data/eta/eta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}_jk{}.txt'\
                                        .format(minradV,maxradV,rmin,rmax,sec,vtype,jk))['eta','rmin','rmax','N']


                    eta_ran = 1./(np.sqrt(2)-1)
                    eta_ran_std = np.sqrt(28.1421/etaTable['N'])

                    zeta_jk.append( (etaTable['eta']-eta_ran)/eta_ran_std )

                    #N = len(beta)
                    #n_gal_jk.append(N)

                #print(zeta_jk)

                if forTable: 
                    filename = '../data/zeta/zeta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}_jk_forTable.txt'
                else: 
                    filename = '../data/zeta/zeta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}_jk.txt'
                #print(filename)
                ascii.write(np.column_stack([zeta_jk]),\
                    filename\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype,jk),\
                    names=['zeta_jk'],\
                    overwrite=True)
# %%
