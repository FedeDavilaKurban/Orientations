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
nseed = 50

eta_ran = 1/(np.sqrt(2)-1)
#eta_ran = eta_ran**(-1)
eta_ran_std = np.sqrt(28.1421/Nran) #beta = perp/prll
#eta_ran_std = np.sqrt(.8284/Nran) #beta = prll/perp

cc = np.array([.6,.8,1])
aa = np.array([1,1,1])
bb = aa
   
#%%      

zeta_eta = []
zeta_cos = []

cos = []
beta = []
eta = []
eta_arrays = []
eta_stds = []

rseeds = np.arange(4,4+nseed)
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

#%%
"""
Plot Parameters
"""   
fig = plt.figure(figsize=(13,17))
fig.subplots_adjust(hspace=0.6, wspace=0.5, bottom=.2)
plt.rcParams['font.size'] = 16
fs = 18

ncol = 2
nrow = 3

ax0 = fig.add_subplot(nrow,ncol,2,projection='polar')
ax1 = fig.add_subplot(nrow,ncol,1) #cos
ax2 = fig.add_subplot(nrow,ncol,3) #betas
ax3 = fig.add_subplot(nrow,ncol,4) #etas 
#ax4 = fig.add_subplot(nrow,ncol,5) #eta vs cos
ax5 = fig.add_subplot(nrow,ncol,(5,6)) #zetas
#ax5 = ax4.twinx().twiny()

ax1.axhline(1,ls='--',lw=1,color='slategrey')
ax2.axvline(0,ls='--',lw=1,color='slategrey')
ax3.axvline(eta_ran,ls='--',lw=1,color='slategrey')
#ax4.axvline(0.5,ls='--',lw=1,color='slategrey')
#ax4.axhline(eta_ran,ls='--',lw=1,\
#    color='slategrey')
ax5.axvline(0,ls=':',color='slategrey')
ax5.axhline(0,ls=':',color='slategrey')

bins=8
clrs = ['#361e15',\
        #'#402218', \
        '#865439', \
        #'#C68B59',\
        '#D7B19D']

#from pylab import *
#clrs = cm.get_cmap('seismic', len(cc))  # matplotlib color palette name, n colors

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
ax0.legend(bbox_to_anchor=(0.05, 1.))
#######################################################

for i in range(len(cc)):
    ax1.hist(cos[i], bins=bins, histtype='step', \
        color=clrs[i], density=True, lw=2)
    ax2.hist(np.log10(beta[i]), bins=bins*3, histtype='step', \
        color=clrs[i], density=True, lw=2)
    ax3.hist(eta_arrays[i], bins=bins*3, histtype='step', \
        color=clrs[i], density=True, lw=2)
    ax3.axvline(np.mean(eta_arrays[i]),ls=':',color=clrs[i])

    ytext=np.arange(len(cc))/10+.1
    s = r'$\langle cos(\lambda)\rangle=$'+f'{cos[i].mean():.2f}'+\
        r'$\pm$'+f'{cos[i].std():.2f}'
    ax1.text(.1,ytext[i],s,\
        color=clrs[i],transform = ax1.transAxes)

    ytext=np.arange(10)[-len(cc):]/10
    s = r'$\langle \eta\rangle=$'+f'{np.mean(eta_arrays[i]):.2f}'+\
        r'$\pm$'+f'{np.std(eta_arrays[i],ddof=1):.2f}'
    ax3.text(.5,ytext[i],s,\
        color=clrs[i],transform = ax3.transAxes)
    ax3.text(.5,np.min(ytext)-.1,r'$\eta_{ran}\simeq$'\
        +f'{eta_ran:.2f}'+r'$\pm$'+f'{eta_ran_std:.2f}',\
        color='slategrey',transform = ax3.transAxes)
    ax3.text(.05,.8,f'{Nbs} bootstrap \n resamplings',\
        color='k',transform = ax3.transAxes)

zeta_cos=[]
zeta_eta=[]
for i in range(len(cc)):    
    xx = np.mean(cos,axis=1)[i::len(cc)]
    yy = eta[i::len(cc)]
    
    #ax4.scatter(xx, yy, color=clrs[i], alpha=.5)
    #ax4.errorbar(np.mean(xx), np.mean(yy), \
    #    xerr=np.std(xx), yerr=np.std(yy), lw=2,\
    #    color='k',capsize=5)
    
    zc = (xx-0.5)/np.sqrt(1/12)
    zeta_cos.append(zc)
    ze = (yy-eta_ran)/eta_ran_std
    zeta_eta.append(ze)
    
    
xp = []
yp = []
for i in range(len(aa)):
    ax5.scatter(zeta_cos[i::len(aa)],zeta_eta[i::len(aa)],\
        color=clrs[i],alpha=.5)

    zc_m = np.mean(zeta_cos[i::len(aa)])
    zc_std = np.std(zeta_cos[i::len(aa)],ddof=1)
    ze_m = np.mean(zeta_eta[i::len(aa)])
    ze_std = np.std(zeta_eta[i::len(aa)],ddof=1)

    xp.append(zc_m)
    yp.append(ze_m)

    ax5.errorbar(zc_m,ze_m,xerr=zc_std,yerr=ze_std,lw=2,\
        color='black',capsize=5,label=f'c/a={cc[i]/aa[i]:.1f}')

    
ax1.text(.6,.9,'Nran={}'.format(Nran),\
    fontsize=fs, alpha=.7, color='k', transform = ax1.transAxes)
ax2.text(.65,.9,r'$\beta=S_\perp/S_\parallel$'.format(Nran),\
    fontsize=fs+4, alpha=.7, color='k', transform = ax2.transAxes)
ax5.text(.75,.8,'{} random samples \n for each {}'.format(nseed,r'$e^2$'),\
    fontsize=fs,alpha=.7,color='k',transform = ax5.transAxes)

ax1.set_xlabel(r'$\cos(\lambda)$')
ax2.set_xlabel(r'$log_{10}(\beta)$')
ax3.set_xlabel(r'$\eta$')
#ax4.set_xlabel(r'$\langle cos(\lambda)\rangle$')
#ax4.set_ylabel(r'$\eta$')
ax5.set_xlabel(r'$\zeta_{cos}=(\langle cos\rangle -cos_{ran})/\sigma_{cos}}$')
ax5.set_ylabel(r'$\zeta_{\eta}=(\eta-\eta_{ran})/\sigma_{\eta}}$')
#ax5.set_xlim([-.35,.14])
plt.tight_layout()
plt.savefig(f'../plots/zetaCos_v_zetaEta_v3/eta_Nran{Nran}_Nbs{Nbs}_nseed{nseed}.pdf')
# %%
