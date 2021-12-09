"""
Quiero hacer un plot que sea zeta_cos vs zeta_eta
Para linkear nuestro método con algo tradicional como el coseno
"""
#%%
import sys
from config import illustrisPath
sys.path.append(illustrisPath+'param_tools')
# use param_tools: https://github.com/maxkapur/param_tools
from param_tools import r_surface, surface_area

import numpy as np
from matplotlib import pyplot as plt
from functools import partial
#%%
"""
Definitions
"""
def ellipsoid(t, u, a=1, b=1, c=1):
    return np.array([a*np.sin(u)*np.cos(t), b*np.sin(u)*np.sin(t), c*np.cos(u)])
domain_t = [0, 2*np.pi]
domain_u = [0, np.pi]

def get_cos(xs,ys,zs):

    cos = np.zeros(len(xs))
    for j in range(len(xs)):
        norm = np.sqrt(xs[j]**2 + ys[j]**2 + zs[j]**2)
        cos[j] = abs(zs[j]) / norm
        
    return cos
        
def get_beta(xs,ys,zs):

    beta = np.zeros(len(xs)) 
    for j in range(len(xs)):
        para = zs[j]
        perp = np.sqrt(xs[j]**2+ys[j]**2)
        beta[j] = perp / abs(para)

    return beta

def get_eta_bs(x,Nbs=1000):
    bs_eta = []
    for _ in range(Nbs):  
        bs_x = np.random.choice(x, size=len(x))

        n_perp = np.sum(bs_x>1.)
        n_prll = np.sum(bs_x<1.)
        bs_eta.append( n_perp / n_prll )
        
        # n_perp = np.sum(bs_x<0.)
        # n_prll = np.sum(bs_x>0.)
        # bs_eta.append( n_prll / n_perp )
    
    eta = np.mean(bs_eta) 
    eta_std = np.std(bs_eta, ddof=1)
    
    return eta, eta_std, bs_eta

#%%
"""
Code parameters
"""
Nran = 1000
Nbs = 100
nseed = 5

eta_ran = 1/(np.sqrt(2)-1)
#eta_ran = eta_ran**(-1)
eta_ran_std = np.sqrt(28.1421/Nran) #beta = perp/prll
#eta_ran_std = np.sqrt(.8284/Nran) #beta = prll/perp

cc = np.array([.6,.8,1])
aa = np.array([1,1,1])
bb = aa
   
#%%      
"""
Plot Parameters
"""   
fig = plt.figure(figsize=(10,12))
fig.subplots_adjust(hspace=0.4, wspace=0.2, bottom=.2)
ncol = 2
nrow = 4

ax0 = fig.add_subplot(nrow,ncol,(1,2),projection='polar')
ax1 = fig.add_subplot(nrow,ncol,3) #cos
ax2 = fig.add_subplot(nrow,ncol,4) #betas
ax3 = fig.add_subplot(nrow,ncol,6) #etas 
ax4 = fig.add_subplot(nrow,ncol,5) #eta vs cos
#ax5 = fig.add_subplot(nrow,ncol,(5,6)) #zetas

ax1.axhline(1,ls='--',lw=1,color='slategrey')
ax2.axvline(0,ls='--',lw=1,color='slategrey')
ax3.axvline(eta_ran,ls='--',lw=1,color='slategrey')
ax4.axvline(0.5,ls='--',lw=1,color='slategrey')
ax4.axhline(eta_ran,ls='--',lw=1,\
    color='slategrey')

bins=6
clrs = ['#361e15',\
        #'#402218', \
        '#865439', \
        #'#C68B59',\
        '#D7B19D']

#from pylab import *
#clrs = cm.get_cmap('seismic', len(cc))  # matplotlib color palette name, n colors

zeta_eta = []
zeta_cos = []

cos = []
beta = []
eta = []
eta_arrays = []
eta_stds = []

rseeds = np.arange(0,0+nseed)
for rseed in rseeds:

    for k, (a, b, c) in enumerate(zip(aa,bb,cc)):
        #print(k,a,b,c)

        f = partial(ellipsoid, a=a, b=b, c=c)

        np.random.seed(rseed)
        x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
        xs = x[0,:]
        ys = x[1,:]
        zs = x[2,:]
        
        cos.append( get_cos(xs,ys,zs) )
        beta.append( get_beta(xs,ys,zs) )
        eta_m, eta_std, eta_array = get_eta_bs(beta[-1],Nbs=Nbs) 
        eta.append(eta_m)
        eta_arrays.append(eta_array)
        eta_stds.append(eta_stds)

zeta_cos = np.zeros(nseed)
zeta_eta = np.zeros(nseed)
for i in range(nseed):
    zeta_cos[i] = (np.mean(cos[i])-0.5)/np.sqrt(1./12)
    zeta_eta[i] = (np.mean(eta[i])-eta_ran)/eta_ran_std


#######################################################
#EYE OF SAURON
def ellipsp(t, a, c):
    e2 = 1 - (c/a)**2
    r = np.sqrt(c**2 / (1-e2*np.cos(t)**2))
    return r
for k, (a, b, c) in enumerate(zip(aa,bb,cc)):
    t_ = np.linspace(0, 2*np.pi, 100)
    r_ = ellipsp(t_, c, a)
    if c<=a: e2=1-(c/a)**2
    if c>a: e2=1-(a/c)**2 
    ax0.plot(t_, r_, clrs[k],label=r'$e^2=$'+f'{e2:.1f}')
ax0.set_theta_zero_location("N")
ax0.legend(bbox_to_anchor=(-.2, 1.05))
#######################################################

for i in range(len(cc)):
    ax1.hist(cos[i], bins=bins, histtype='step', \
        color=clrs[i], density=True, lw=2)
    ax2.hist(np.log10(beta[i]), bins=bins*2, histtype='step', \
        color=clrs[i], density=True, lw=2)
    ax3.hist(eta_arrays[i], bins=bins*3, histtype='step', \
        color=clrs[i], density=True, lw=2)
    ax3.axvline(np.mean(eta_arrays[i]),ls=':',color=clrs[i])

for i in range(len(cc)):
    
    xx = np.mean(cos,axis=1)[i::len(cc)]
    yy = eta[i::len(cc)]
    
    ax4.scatter(xx, yy, color=clrs[i])
    ax4.errorbar(np.mean(xx), np.mean(yy), xerr=np.std(xx), yerr=np.std(yy),\
        color='k',capsize=3)
    
ax1.text(.7,.9,'Nran={}'.format(Nran),\
    fontsize=10, alpha=.7, color='k', transform = ax1.transAxes)
ax4.text(.3,.9,'{} random samples for each c/a'.format(nseed),\
    fontsize=10,alpha=.7,color='k',transform = ax4.transAxes)

ax1.set_xlabel(r'$\cos(\lambda)$')
ax2.set_xlabel(r'$log_{10}(\beta)$')
ax3.set_xlabel(r'$\eta$')
ax4.set_xlabel(r'$\langle cos(\lambda)\rangle$')
ax4.set_ylabel(r'$\eta$')

# %%
