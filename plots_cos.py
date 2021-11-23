#%%
#import sys
import numpy as np
#from scipy import spatial
#import scipy.stats
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
#import random
#from config import writePath, units
import matplotlib.pyplot as plt
#%%
########################################################################
"""
Histograma de cosenos para zeta=0,1,2,3
"""
########################################################################

minradV = 7.
maxradV = 0.
sec = 0

plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['font.size'] = 15

plt.hlines(1.,-0.05,1.05,ls='-',color='k',linewidth=3,alpha=.6)

for rmin,rmax,vtype,zeta,ls in ([.8,.9,'a',0,'-'],\
                                [1.,1.1,'s',1,'-.'],\
                                [1.4,1.5,'r',2,'--'],\
                                [1.4,1.5,'a',3,'-']):
    print(rmin,rmax,vtype)

    cosT = ascii.read('../data/cos/cos_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,rmin,rmax,sec,vtype),\
                names=['cos'])

    cos = cosT['cos'].data
    cosm = np.median(cos)
    print('median: ',cosm)

    label = r'$\zeta\simeq${}'.format(zeta)

    plt.hist(cos,histtype='step',ls=ls,\
                density=True,bins=6,label=label,linewidth=2)

    #plt.vlines(cosm,.85,1.15,ls=ls)

plt.xlabel(r'$cos(\lambda)$')
plt.ylim([.85,1.15])
plt.legend(ncol=2,loc='upper right')
plt.show()
#%%
########################################################################
"""
Histograma de arccos(lambda) para zeta=0,1,2,3
"""
########################################################################

minradV = 7.
maxradV = 0.
sec = 0

plt.rcParams['figure.figsize'] = (8, 6)
plt.rcParams['font.size'] = 15

#plt.hlines(1.,-0.05,1.05,ls='-',color='k',linewidth=3,alpha=.6)

for rmin,rmax,vtype,zeta,ls in ([.8,.9,'a',0,'-'],\
                                [1.,1.1,'s',1,'-.'],\
                                [1.4,1.5,'r',2,'--'],\
                                [1.4,1.5,'a',3,'-']):
    print(rmin,rmax,vtype)

    cosT = ascii.read('../data/cos/cos_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,rmin,rmax,sec,vtype),\
                names=['cos'])

    cos = cosT['cos'].data
    lmbda = np.arccos(cos)*180./np.pi
    lmbdaM = np.median(lmbda)
    print('mean and median: ',np.mean(lmbda),lmbdaM)

    label = r'$\zeta\simeq${}'.format(zeta)

    plt.hist(lmbda,histtype='step',ls=ls,\
                density=True,bins=50,label=label,linewidth=2)

    #plt.vlines(cosm,.85,1.15,ls=ls)

plt.xlabel(r'$\lambda$')
#plt.ylim([.85,1.15])
plt.legend(ncol=1,loc='upper left')
plt.show()
# %%
########################################################################
"""
Histograma de log(Beta) para zeta=0,1,2,3
"""
########################################################################

plt.vlines(0,0.,1,linewidth=3,alpha=.6,color='k')

for rmin,rmax,vtype,zeta,ls in ([.8,.9,'a',0,'-'],\
                                [1.,1.1,'s',1,'-.'],\
                                [1.4,1.5,'r',2,'--'],\
                                [1.4,1.5,'a',3,'-']):
    print(rmin,rmax,vtype)

    betaT = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,rmin,rmax,sec,vtype))

    beta = betaT['beta'].data
    np.where(beta==0)
    lbeta = np.log10(beta)
    lbetaM = np.median(np.log10(beta))
    print('median: ',lbetaM)

    label = r'$\zeta\simeq${}'.format(zeta)

    plt.hist(lbeta,histtype='step',ls=ls,\
                density=True,bins=15,label=label,linewidth=2)

#    plt.vlines(lbetaM,0.,0.7,ls=ls)

plt.xlabel(r'$log_{10}(\beta)$')
plt.ylim([0,.9])
plt.legend(ncol=1,loc='upper right')
plt.show()

# beta = ascii.read('../data/beta/-1/beta_minradV7.0_maxradV0.0_rmin1.1_rmax1.2_sec3_vtyper.txt')['beta'].data
# vrad = ascii.read('../data/vel/-1/vel_minradV7.0_maxradV0.0_rmin1.1_rmax1.2_sec3_vtyper.txt')['vrad'].data
# p50 = np.percentile(vrad,50)
# beta = beta[vrad<p50]
# print('median: ',np.mean(np.log10(beta)))
# plt.hist(np.log10(beta),histtype='step',ls=ls,\
#             density=True,bins=30,label=label,linewidth=2,color='k')
# plt.show()

# %%
########################################################################
"""
Histograma de arctg(Beta) para zeta=0,1,2,3
"""
########################################################################

#plt.vlines(0,0.,1,linewidth=3,alpha=.6,color='k')

for rmin,rmax,vtype,zeta,ls in ([.8,.9,'a',0,'-'],\
                                [1.,1.1,'s',1,'-.'],\
                                [1.4,1.5,'r',2,'--'],\
                                [1.4,1.5,'a',3,'-']):
    print(rmin,rmax,vtype)

    betaT = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,rmin,rmax,sec,vtype))

    beta = betaT['beta'].data
    lmbda = np.arctan(beta)*180/np.pi
    print('mean and median: ',np.mean(lmbda),np.median(lmbda))

    label = r'$\zeta\simeq${}'.format(zeta)

    plt.hist(lmbda,histtype='step',ls=ls,\
                density=True,bins=15,label=label,linewidth=2)

#    plt.vlines(lbetaM,0.,0.7,ls=ls)

plt.xlabel(r'$\lambda$')
#plt.ylim([0,.9])
plt.legend(ncol=1,loc='upper left')
plt.show()

# %%
