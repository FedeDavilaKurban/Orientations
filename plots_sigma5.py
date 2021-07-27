#%%
"""

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
# %%
"""
B vs Sigma5
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

    sigma5 = ascii.read('../data/sigma5/sigma5_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['sigma5'].data

    y = np.log10(beta)
    x = np.log10(sigma5)
    m,b,rvalue,pvalue,std = scipy.stats.linregress(x,y) 


    p50 = np.percentile(sigma5,50)
    x1 = x[sigma5<p50]
    x2 = x[sigma5>p50]
    y1 = y[sigma5<p50]
    y2 = y[sigma5>p50]

    ax.scatter(x1, y1, s=.5)
    ax.scatter(x2, y2, s=.5)

    ax.plot(x,m*x+b,ls=':',color='k',label='r_lr={:.2f}'.format(rvalue))
    #ax.legend()

    ax.set_title('section {} - vtype {}'.format(sec,vtype))
    ax.set_ylabel(r'$log_{10}(\beta)$')


axs[2].set_xlabel(r'$log_{10}(\Sigma_5)\,[kpc^{-2}]$')


plt.savefig('../plots/sigma5/BvsSigma5.jpg')

#%%
"""
Eta vs Low/High Sigma5
"""
minradV=7.
maxradV=0.

rmin, rmax = 1.1,1.2
sec = 3

plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    sigma5 = ascii.read('../data/sigma5/sigma5_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype),names=['sigma5'])['sigma5'].data

    p50 = np.percentile(sigma5,50)

    sigma5_1 = sigma5[sigma5<p50]
    sigma5_2 = sigma5[sigma5>p50]

    beta1 = beta[sigma5<p50]
    beta2 = beta[sigma5>p50]

    logbeta1 = np.log10(beta1)
    logbeta2 = np.log10(beta2)

    eta1, eta1_std = get_eta_bs(logbeta1)
    eta2, eta2_std = get_eta_bs(logbeta2)
    print(eta1, eta1_std, eta2, eta2_std)

    ax.hlines(1/(np.sqrt(2)-1),.5,2.5,linestyles=':')

    ax.errorbar(1,eta1,eta1_std,fmt='o',capsize=5)
    ax.errorbar(2,eta2,eta2_std,fmt='o',capsize=5)
    ax.set_title('section {} - vtype {}'.format(sec,vtype))
    ax.set_ylabel(r'$\eta$')

axs[2].set_xlim([0.,3.])
axs[2].set_xlabel(r'$\Sigma_5$')
axs[2].set_xticks([1,2])
axs[2].set_xticklabels(['Low '+r'$\Sigma_5$','High '+r'$\Sigma_5$'])

#plt.savefig('../plots/sigma5/EtavsSigma5.jpg')

#%%
"""
Eta vs R for Lo/Hi Sigma5
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

for sfilter in ['hi','lo']:
    print(sfilter)

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

                sigma5 = ascii.read('../data/sigma5/sigma5_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype),names=['sigma5'])['sigma5'].data

                p50 = np.percentile(sigma5,50)
                
                if sfilter=='hi': beta = beta[sigma5>p50]
                elif sfilter=='lo': beta = beta[sigma5<p50]

                x = np.log10(beta.data)

                eta_, eta_std_ = get_eta_bs(x)
    
                eta.append( eta_ )
                eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

                if vtype=='r' and rmin==1.1 and rmax==1.2: print(eta_,eta_std_) 

                # Obtain mean and var of control samples
                N = len(beta)
            
                eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
                eta_random_mean.append( np.mean(eta_random) )
                eta_random_var.append( np.var(eta_random,ddof=1) )
                eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )



        nbins = len(r1)
        x = (r1+r2)/2


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

        if sec==3: 
            if sfilter=='hi': title='H-H Galaxies - High {}'.format(r'$\Sigma_5$')
            elif sfilter=='lo': title='H-H Galaxies - Low {}'.format(r'$\Sigma_5$')

        axs[0].set_title(title)

        axs[2].set_xlabel('R/Rv')

        axs[0].text(1.25,3.15,'All Voids')
        axs[1].text(1.25,3.15,'Rising Voids')
        axs[2].text(1.25,3.15,'Shell Voids')

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
            #print(y)
            
            ax.set_ylabel(r'$\eta$')

            ax.set_ylim([1.5,3.5])

        axs[0].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[1].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[2].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

        #plt.savefig('../plots/sigma5/Eta_{}Sigma5_allvtypes_sec{}_.jpg'.format(sfilter,sec))

# %%
#%%
"""
Eta vs R for Lo/Hi Sigma5 -- 2nd Version (no "all voids")
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

for sfilter in ['hi','lo']:
    print(sfilter)

    for sec in [3]:
        print(sec)
        eta = []
        eta_std = []

        eta_random_mean = []
        eta_random_var = []
        eta_random_std = []
        
        n_gal = []
        for vtype in ['r','s']:
            print(vtype)
            
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

                if vtype=='r' and rmin==1.1 and rmax==1.2: print(eta_,eta_std_) 

                # Obtain mean and var of control samples
                N = len(beta)
            
                eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
                eta_random_mean.append( np.mean(eta_random) )
                eta_random_var.append( np.var(eta_random,ddof=1) )
                eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )



        nbins = len(r1)
        x = (r1+r2)/2


        #ya = eta[:nbins]
        yr = eta[:nbins]
        ys = eta[-nbins:]

        #ya_err = eta_std[:nbins]
        yr_err = eta_std[:nbins]
        ys_err = eta_std[-nbins:]

        ####################
        eta_random_mean = np.array(eta_random_mean)
        eta_random_var = np.array(eta_random_var)
        eta_random_std = np.array(eta_random_std)

        #yran_mean_a = eta_random_mean[:nbins]
        yran_mean_r = eta_random_mean[:nbins]
        yran_mean_s = eta_random_mean[-nbins:]

        #yran_var_a = eta_random_var[:nbins]
        yran_var_r = eta_random_var[:nbins]
        yran_var_s = eta_random_var[-nbins:]

        #yran_std_a = eta_random_std[:nbins]
        yran_std_r = eta_random_std[:nbins]
        yran_std_s = eta_random_std[-nbins:]
        ####################

        fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=True)

        if sec==3: 
            if sfilter=='hi': title='H-H Galaxies - High {}'.format(r'$\Sigma_5$')
            elif sfilter=='lo': title='H-H Galaxies - Low {}'.format(r'$\Sigma_5$')

        axs[0].set_title(title)

        axs[1].set_xlabel('R/Rv')

        #axs[0].text(1.25,3.15,'All Voids')
        axs[0].text(1.25,3.15,'Rising Voids')
        axs[1].text(1.25,3.15,'Shell Voids')

        for ax, y, yerr, yran_mean, yran_err, label, in zip(axs,\
                                        (yr,ys),\
                                        (yr_err,ys_err),\
                                        (yran_mean_r,yran_mean_s),\
                                        (yran_std_r,yran_std_s),('Rising','Shell')):

            # Theoretical n1/n2 value for random spins
            ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
            ax.plot(x,yran_mean,c='k',alpha=.7)

            ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.15, color='k')
            ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.15, color='k')
            ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.15, color='k')

            ax.errorbar(x,y,yerr=yerr,label=label)
            #print(y)
            
            ax.set_ylabel(r'$\eta$')

            ax.set_ylim([1.5,3.5])

        axs[0].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[1].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        #axs[2].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

        plt.savefig('../plots/sigma5/Eta_{}Sigma5_sec{}.jpg'.format(sfilter,sec))

# %%
"""
Eta vs R for Lo/Hi Sigma5 -- 3rd Version (eta/sigma(eta))
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

fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=True)

for sfilter in ['hi','lo']:
    print(sfilter)

    for sec in [3]:
        print(sec)
        eta = []
        eta_std = []

        eta_random_mean = []
        #eta_random_var = []
        eta_random_std = []
        
        for vtype in ['r','s']:
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

                #if vtype=='r' and rmin==1.1 and rmax==1.2: print(eta_,eta_std_) 

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

        yran_mean_r = eta_random_mean[:nbins]
        yran_mean_s = eta_random_mean[-nbins:]

        yran_std_r = eta_random_std[:nbins]
        yran_std_s = eta_random_std[-nbins:]
        ####################

        yr = (eta[:nbins]-yran_mean_r)/yran_std_r
        ys = (eta[-nbins:]-yran_mean_s)/yran_std_s

        yr_err = eta_std[:nbins]/yran_std_r
        ys_err = eta_std[-nbins:]/yran_std_s

        ####################

        if sec==3: 
            title='H-H Galaxies'

        axs[0].set_title(title)

        axs[1].set_xlabel('R/Rv')

        #axs[0].text(1.25,3.15,'All Voids')
        axs[0].text(1.25,4.45,'Rising Voids')
        axs[1].text(1.25,4.45,'Shell Voids')

        for ax, y, yerr in zip(axs,(yr,ys),(yr_err,ys_err)):

            # Theoretical n1/n2 value for random spins
            #ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
            ax.hlines(0,x[0],x[-1],linestyles=':')
            #ax.plot(x,yran_mean,c='k',alpha=.7)

            ax.fill_between(x, -1, 1, alpha=.025, color='k')
            ax.fill_between(x, -1, 3, alpha=.03, color='k')


            if sfilter=='hi': label='High '+r'$\Sigma_5$'
            if sfilter=='lo': label='Low '+r'$\Sigma_5$'

            ax.errorbar(x,y,yerr=yerr,capsize=3,fmt='o-',ms=5,label=label)

            ax.set_ylabel(r'$(\eta-\eta_0)/\sigma(\eta)$')


            #ax.set_ylim([1.5,3.5])

        axs[0].legend(loc='upper left',framealpha=.6)

        axs[0].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[1].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

    #plt.savefig('../plots/sigma5/Eta_Sigma5_sec{}_v3.jpg'.format(sec))
# %%
