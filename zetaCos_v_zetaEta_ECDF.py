#%%
"""
Quiero hacer un plot que sea zeta_cos vs zeta_eta
Para linkear nuestro método con algo tradicional como el coseno
"""
import sys
from config import illustrisPath
sys.path.append(illustrisPath+'param_tools')
# use param_tools: https://github.com/maxkapur/param_tools
from param_tools import r_surface, surface_area

import numpy as np
from matplotlib import pyplot as plt
from functools import partial
#%%

def ellipsoid(t, u, a=1, b=1, c=1):
    return np.array([a*np.sin(u)*np.cos(t), b*np.sin(u)*np.sin(t), c*np.cos(u)])
domain_t = [0, 2*np.pi]
domain_u = [0, np.pi]

#%%
# for rseed in [2]:
#     fig = plt.figure(figsize=(10,15))
#     fig.subplots_adjust(hspace=0.3, wspace=0.1, bottom=.2)
#     ax1 = fig.add_subplot(3,2,(3,4)) #etas
#     ax2 = fig.add_subplot(3,2,1) #cos
#     #ax4 = fig.add_subplot(2,2,3) #beta
#     ax3 = fig.add_subplot(3,2,2,projection='polar')
#     ax5 = fig.add_subplot(3,2,(5,6))

#     bins = np.linspace(1.9, 4.8, 20)

#     clrs = ['#361e15','#402218', '#865439', '#C68B59', '#D7B19D']

#     cs = [.6,.8,1.,1.,1.]
#     bs = [1.,1.,1.,1-(1/.8-1),1-(1/.6-1)]

#     Nran = 500
#     Netas = 300*2

#     np.random.seed(rseed)
#     for k, (b, c) in enumerate(zip(bs,cs)):
#         print(b,c)

#         # if c>0.:
#         #     f = partial(ellipsoid, c=c)
#         # else:
#         #     b=-c
#         #     c=1
#         #     f = partial(ellipsoid, b=b)
#         f = partial(ellipsoid, a=b, c=c)

#         eta = np.zeros(Netas)
#         for i in range(Netas):
#             x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
#             xs = x[0,:]
#             ys = x[1,:]
#             zs = x[2,:]
            
#             beta = np.zeros(Nran)
#             for j in range(len(xs)):
#                 para = zs[j]
#                 perp = np.sqrt(xs[j]**2+ys[j]**2)
#                 beta[j] = perp / abs(para)
#             eta[i] = sum(beta>1)/sum(beta<1)

#         eta.sort()
#         ax1.step(eta, np.arange(0,Netas)/Netas, linewidth=2, 
#                 color=clrs[k], label=f'b={b:.1f},c={c:.1f}')
#         ax1.axvline(eta.mean(), color=clrs[k], linestyle=':')

#         print(eta.mean(), eta.mean()-1/(np.sqrt(2)-1) )

#         x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
#         xs = x[0,:]
#         ys = x[1,:]
#         zs = x[2,:]   

#         cos = np.zeros(Nran)
#         beta = np.zeros(Nran)
#         for j in range(len(xs)):
#             para = zs[j]
#             perp = np.sqrt(xs[j]**2 + ys[j]**2)
#             beta[j] = perp / abs(para)

#             norm = np.sqrt(xs[j]**2 + ys[j]**2 + zs[j]**2)
#             cos[j] = abs(para) / norm

#         cos.sort()


#         ###################
#         theta_m = np.mean(np.arccos(cos))*180./np.pi
#         print('coseno promedio:',cos.mean())
#         print('angulo promedio:',theta_m)

#         ####################
#         #ax4.hist(np.log10(beta), histtype='step', color=clrs[k], linewidth=2, density=True)
#         ax2.step(cos, np.arange(0.,Nran)/Nran,\
#                 color=clrs[k], linewidth=2)

#         #################
#         zeta_eta = (eta.mean()-(1./(np.sqrt(2)-1)))/eta.std()
#         zeta_cos = (cos.mean()-.5)/cos.std()
#         #coef = np.polyfit(zeta_cos,zeta_eta,1)
#         #poly1d_fn = np.poly1d(coef) 
#         ax5.scatter(zeta_cos,zeta_eta,color=clrs[k])
#         #ax5.plot(poly1d_fn(zeta_cos), '--k')

#     ax1.axvline(2.41, color='slategrey', linestyle='--',label='Random')
#     #ax2.step(np.arange(0,1), np.arange(0,1),color='slategrey', linestyle='--')
#     #ax4.axvline(0, color='slategrey', linestyle='--')

#     ax1.legend()

#     ax1.set_xlabel(r'$\eta$')
#     ax2.set_xlabel(r'$cos(\lambda)$')
#     ax5.set_xlabel(r'$\zeta_{cos}$')
#     ax5.set_ylabel(r'$\zeta_{\eta}$')

#     ax5.set_xlim([-.3,.3])

#     #ax4.set_xlabel(r'$log_{10}(\beta)$')
#     ################


#     def ellipsp(t, a, b):
#         e2 = 1 - (b/a)**2
#         r = np.sqrt(b**2 / (1-e2*np.cos(t)**2))
#         return r

#     t = np.linspace(0, 2*np.pi, 100)

#     for k, (b, c) in enumerate(zip(bs,cs)):

#         # if c>0:
#         #     r = ellipsp(t, c, 1)
#         # else:
#         #     r = ellipsp(t, 1, -c)
#         r = ellipsp(t, c, b)
#         ax3.plot(t, r, clrs[k])

#     ax3.set_theta_zero_location("N")
#     ax3.grid(True)

#     plt.show()
#     #plt.savefig('../plots/zetaCos_v_zetaEta_ECDF/zetaCos_v_zetaEta_ECDF_ranseed{}.jpg'.format(rseed))
 # %%

"""
Quiero hacer el plot zeta vs zeta con varios random seed
"""
#%%
def get_eta_bs(x,Nbs=1000):
    bs_eta = []
    #bs_cos_mean = []
    for _ in range(Nbs):  
        bs_x = np.random.choice(x, size=len(x))
        #bs_cos = np.random.choice(cos, size=len(cos))

        n_perp = np.sum(bs_x>0.)
        n_prll = np.sum(bs_x<0)
        bs_eta.append( n_perp / n_prll )
    
        #bs_cos_mean.append( np.mean(bs_cos) )

    eta = np.mean(bs_eta) 
    eta_std = np.std(bs_eta, ddof=1)
    
    return eta, eta_std, bs_eta

fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(hspace=0.3, wspace=0.2, bottom=.2)
ax3 = fig.add_subplot(3,2,1) #cos
ax2 = fig.add_subplot(3,2,2,projection='polar')
ax4 = fig.add_subplot(3,2,3) #etas 
ax1 = fig.add_subplot(3,2,4) 

ax5 = fig.add_subplot(3,2,(5,6)) #zetas

clrs = ['#361e15',\
        #'#402218', \
        '#865439', \
        #'#C68B59',\
        '#D7B19D']

cc = np.array([.6,1,1])
aa = np.array([1,1,.6])
bb = aa

#from pylab import *
#clrs = cm.get_cmap('seismic', len(cc))  # matplotlib color palette name, n colors

Nran = 500
Nbs = 1000
#Netas = 300*2

zeta_eta = []
zeta_cos = []

all_etas = []
all_cos_m = []
diff_etas=[]

nseed = 50
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
        
        beta = np.zeros(Nran)
        cos = np.zeros(Nran)

        for j in range(len(xs)):
            para = zs[j]
            perp = np.sqrt(xs[j]**2+ys[j]**2)
            beta[j] = perp / abs(para)
    
            #beta[j] = abs(para)/perp

            norm = np.sqrt(xs[j]**2 + ys[j]**2 + zs[j]**2)
            cos[j] = abs(para) / norm
        
        eta_bs, eta_std, eta_array = \
            get_eta_bs(np.log10(beta),Nbs=Nbs)

        #Esto lo hice porque queria ver un
        #histograma de la diferencia entre estas dos etas
        #a.k.a.: diff_etas
        #eta = sum(beta>1)/sum(beta<1)
        #diff_etas.append(eta-eta_bs)

        eta_ran = 1/(np.sqrt(2)-1)
        eta_ran_std = 28.1421/Nran

        all_etas.append(eta_bs)
        all_cos_m.append(cos.mean())
        #################
        ze = (eta_bs-eta_ran)/eta_ran_std
        zc = (cos.mean()-.5)/np.sqrt(1./12)

        zeta_eta.append( ze )
        zeta_cos.append( zc )
        # ax1.scatter(zc,ze,\
        #         alpha=.6,color=clrs[k],label=f'c/a={c/a:.1f}')
        # QUIERO PROBAR PLOTTEAR SOLO ETA VS COS
        ax1.axvline(0.5,ls='--',lw=1,color='slategrey')
        ax1.axhline(eta_ran,ls='--',lw=1,\
            color='slategrey')
        ax1.scatter(cos.mean(),eta_bs,\
                alpha=.5,color=clrs[k],label=f'c/a={c/a:.1f}')

#################################################################
        def ellipsp(t, a, c):
            e2 = 1 - (c/a)**2
            r = np.sqrt(c**2 / (1-e2*np.cos(t)**2))
            return r
        t_ = np.linspace(0, 2*np.pi, 100)
        r_ = ellipsp(t_, c, a)
        ax2.plot(t_, r_, clrs[k],label=f'c/a={c/a:.1f}')
#################################################################

        if rseed==rseeds[0]:
            bins = 6
            ax3.hist(cos,bins=bins, histtype='step',\
                 color=clrs[k], density=True, lw=2)
            ax4.hist(eta_array, bins=bins*3, histtype='step',\
                 color=clrs[k], density=True, lw=2)
            ax4.axvline(eta_ran,color='slategrey',ls='--')
            ax4.axvline(eta_bs,color=clrs[k],ls=':')

            ytext=np.arange(len(cc))/10+.1
            s = r'$\langle cos(\lambda)\rangle=$'+f'{cos.mean():.2f}'+\
                r'$\pm$'+f'{cos.std():.2f}'
            ax3.text(.1,ytext[k],s,\
                color=clrs[k],transform = ax3.transAxes)

            ytext=np.arange(10)[-len(cc):]/10
            s = r'$\langle \eta\rangle=$'+f'{np.mean(eta_array):.2f}'+\
                r'$\pm$'+f'{np.std(eta_array,ddof=1):.2f}'
            ax4.text(.65,ytext[k],s,\
                color=clrs[k],transform = ax4.transAxes)
            ax4.text(.65,np.min(ytext)-.1,r'$\eta_{ran}\simeq$'\
                +f'{eta_ran:.2f}'+r'$\pm$'+f'{eta_ran_std:.2f}',\
                color='slategrey',transform = ax4.transAxes)
            ax4.text(.67,.1,f'{Nbs} bootstrap \n resamplings',\
                color='k',transform = ax4.transAxes)

        ax3.text(.75,.9,'Nran={}'.format(Nran),\
            fontsize=10,alpha=.7,color='k', transform = ax3.transAxes)

        ax1.text(.3,.9,'{} random samples for each c/a'.format(nseed),\
            fontsize=10,alpha=.7,color='k',transform = ax1.transAxes)

##############################################################

    ax2.set_theta_zero_location("N")
    #ax2.grid(True)
 
    if rseed==rseeds[0]: ax2.legend(bbox_to_anchor=(-.1, 1.05))


    ax1.set_xlabel(r'$\langle cos(\lambda)\rangle$')
    ax1.set_ylabel(r'$\eta$')

    #ax1.set_xlim([-.4,.3])

    ax3.set_xlabel(r'$\cos(\lambda)$')
    ax4.set_xlabel(r'$\eta$')
    ################



# for i in range(len(aa)):
#     zc_m = np.mean(zeta_cos[i::len(aa)])
#     zc_std = np.std(zeta_cos[i::len(aa)],ddof=1)
#     ze_m = np.mean(zeta_eta[i::len(aa)])
#     ze_std = np.std(zeta_eta[i::len(aa)],ddof=1)

#     ax1.errorbar(zc_m,ze_m,xerr=zc_std,yerr=ze_std,color='k',capsize=3)

# QUIERO PROBAR PLOTTEAR SOLO ETA VS COS
for i in range(len(aa)):
    cos_m = np.mean(all_cos_m[i::len(aa)])
    cos_m_std = np.std(all_cos_m[i::len(aa)],ddof=1)
    eta_m = np.mean(all_etas[i::len(aa)])
    eta_m_std = np.std(all_etas[i::len(aa)],ddof=1)

    ax1.errorbar(cos_m,eta_m,xerr=cos_m_std,yerr=eta_m_std,\
        color='k',capsize=3)

#plt.show()

###################################

for i in range(len(aa)):
    ax5.scatter(zeta_cos[i::len(aa)],zeta_eta[i::len(aa)],\
        color=clrs[i],alpha=.5)


    zc_m = np.mean(zeta_cos[i::len(aa)])
    zc_std = np.std(zeta_cos[i::len(aa)],ddof=1)
    ze_m = np.mean(zeta_eta[i::len(aa)])
    ze_std = np.std(zeta_eta[i::len(aa)],ddof=1)

    ax5.errorbar(zc_m,ze_m,xerr=zc_std,yerr=ze_std,alpha=2,lw=3,\
        color=clrs[i],capsize=5,label=f'c/a={cc[i]/aa[i]:.1f}')

#ax5.set_yscale('log')
#ax5.set_xscale('log')
ax5.axvline(0,ls=':',color='slategrey')
ax5.axhline(0,ls=':',color='slategrey')
ax5.legend()
ax5.set_xlabel(r'$\zeta_{cos}=(\langle cos\rangle -cos_{ran})/\sigma_{cos}}$')
ax5.set_ylabel(r'$\zeta_{\eta}=(\eta-\eta_{ran})/\sigma_{\eta}}$')
plt.show()

#plt.savefig('../plots/zetaCos_v_zetaEta_ECDF/zetaCos_v_zetaEta_ECDF_v2.jpg'.format(rseed))

# %%
