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
for rseed in [2]:
    fig = plt.figure(figsize=(10,15))
    fig.subplots_adjust(hspace=0.3, wspace=0.1, bottom=.2)
    ax1 = fig.add_subplot(3,2,(3,4)) #etas
    ax2 = fig.add_subplot(3,2,1) #cos
    #ax4 = fig.add_subplot(2,2,3) #beta
    ax3 = fig.add_subplot(3,2,2,projection='polar')
    ax5 = fig.add_subplot(3,2,(5,6))

    bins = np.linspace(1.9, 4.8, 20)

    clrs = ['#361e15','#402218', '#865439', '#C68B59', '#D7B19D']

    cs = [.6,.8,1.,1.,1.]
    bs = [1.,1.,1.,1-(1/.8-1),1-(1/.6-1)]

    Nran = 500
    Netas = 300*2

    np.random.seed(rseed)
    for k, (b, c) in enumerate(zip(bs,cs)):
        print(b,c)

        # if c>0.:
        #     f = partial(ellipsoid, c=c)
        # else:
        #     b=-c
        #     c=1
        #     f = partial(ellipsoid, b=b)
        f = partial(ellipsoid, a=b, c=c)

        eta = np.zeros(Netas)
        for i in range(Netas):
            x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
            xs = x[0,:]
            ys = x[1,:]
            zs = x[2,:]
            
            beta = np.zeros(Nran)
            for j in range(len(xs)):
                para = zs[j]
                perp = np.sqrt(xs[j]**2+ys[j]**2)
                beta[j] = perp / abs(para)
            eta[i] = sum(beta>1)/sum(beta<1)

        eta.sort()
        ax1.step(eta, np.arange(0,Netas)/Netas, linewidth=2, 
                color=clrs[k], label=f'b={b:.1f},c={c:.1f}')
        ax1.axvline(eta.mean(), color=clrs[k], linestyle=':')

        print(eta.mean(), eta.mean()-1/(np.sqrt(2)-1) )

        x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
        xs = x[0,:]
        ys = x[1,:]
        zs = x[2,:]   

        cos = np.zeros(Nran)
        beta = np.zeros(Nran)
        for j in range(len(xs)):
            para = zs[j]
            perp = np.sqrt(xs[j]**2 + ys[j]**2)
            beta[j] = perp / abs(para)

            norm = np.sqrt(xs[j]**2 + ys[j]**2 + zs[j]**2)
            cos[j] = abs(para) / norm

        cos.sort()


        ###################
        theta_m = np.mean(np.arccos(cos))*180./np.pi
        print('coseno promedio:',cos.mean())
        print('angulo promedio:',theta_m)

        ####################
        #ax4.hist(np.log10(beta), histtype='step', color=clrs[k], linewidth=2, density=True)
        ax2.step(cos, np.arange(0.,Nran)/Nran,\
                color=clrs[k], linewidth=2)

        #################
        zeta_eta = (eta.mean()-(1./(np.sqrt(2)-1)))/eta.std()
        zeta_cos = (cos.mean()-.5)/cos.std()
        #coef = np.polyfit(zeta_cos,zeta_eta,1)
        #poly1d_fn = np.poly1d(coef) 
        ax5.scatter(zeta_cos,zeta_eta,color=clrs[k])
        #ax5.plot(poly1d_fn(zeta_cos), '--k')

    ax1.axvline(2.41, color='slategrey', linestyle='--',label='Random')
    #ax2.step(np.arange(0,1), np.arange(0,1),color='slategrey', linestyle='--')
    #ax4.axvline(0, color='slategrey', linestyle='--')

    ax1.legend()

    ax1.set_xlabel(r'$\eta$')
    ax2.set_xlabel(r'$cos(\lambda)$')
    ax5.set_xlabel(r'$\zeta_{cos}$')
    ax5.set_ylabel(r'$\zeta_{\eta}$')

    ax5.set_xlim([-.3,.3])

    #ax4.set_xlabel(r'$log_{10}(\beta)$')
    ################


    def ellipsp(t, a, b):
        e2 = 1 - (b/a)**2
        r = np.sqrt(b**2 / (1-e2*np.cos(t)**2))
        return r

    t = np.linspace(0, 2*np.pi, 100)

    for k, (b, c) in enumerate(zip(bs,cs)):

        # if c>0:
        #     r = ellipsp(t, c, 1)
        # else:
        #     r = ellipsp(t, 1, -c)
        r = ellipsp(t, c, b)
        ax3.plot(t, r, clrs[k])

    ax3.set_theta_zero_location("N")
    ax3.grid(True)

    plt.show()
    #plt.savefig('../plots/zetaCos_v_zetaEta_ECDF/zetaCos_v_zetaEta_ECDF_ranseed{}.jpg'.format(rseed))
 # %%

"""
Quiero hacer el plot zeta vs zeta con varios random seed
"""
#%%
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
#from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize=(6,10))
fig.subplots_adjust(hspace=0.3, wspace=0.1, bottom=.2)
ax1 = fig.add_subplot(2,1,2) 
ax2 = fig.add_subplot(2,1,1,projection='polar')

bins = np.linspace(1.9, 4.8, 20)

clrs = ['#361e15','#402218', '#865439', '#C68B59', '#D7B19D']

cc = np.array([.5,.8,1,1,1])**1
aa = np.array([1,1,1,.8,.4])**1
bb = aa

#from pylab import *
#clrs = cm.get_cmap('seismic', len(cc))  # matplotlib color palette name, n colors

Nran = 500
#Netas = 300*2

zeta_eta = []
zeta_cos = []

nseed = 30

for rseed in np.arange(0,nseed):

    np.random.seed(rseed)
    for k, (a, b, c) in enumerate(zip(aa,bb,cc)):
        #print(b,c)

        # if c>0.:
        #     f = partial(ellipsoid, c=c)
        # else:
        #     b=-c
        #     c=1
        #     f = partial(ellipsoid, b=b)
        f = partial(ellipsoid, a=a, b=b, c=c)

        # eta = np.zeros(Netas)
        # for i in range(Netas):
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

            norm = np.sqrt(xs[j]**2 + ys[j]**2 + zs[j]**2)
            cos[j] = abs(para) / norm
        
        eta_m, eta_std = get_eta_bs(np.log10(beta),Nbs=10)

        #################
        zeta_eta.append( (eta_m-(1./(np.sqrt(2)-1)))/eta_std )
        zeta_cos.append( (cos.mean()-.5)/cos.std() )
        ax1.scatter(zeta_cos[-1],zeta_eta[-1],\
                alpha=.6,color=clrs[k],label=f'c/a={c/a:.1f}')

#################################################################
        # if rseed==0:
        #     ax2.scatter(ys,zs,color=clrs[k])

        def ellipsp(t, a, c):
            e2 = 1 - (c/a)**2
            r = np.sqrt(c**2 / (1-e2*np.cos(t)**2))
            return r
        t_ = np.linspace(0, 2*np.pi, 100)
        r_ = ellipsp(t_, a, c)
        ax2.plot(t_, r_, clrs[k])

        # # Set of all spherical angles:
        # u = np.linspace(0, 2 * np.pi, 100)
        # v = np.linspace(0, np.pi, 100)

        # # Cartesian coordinates that correspond to the spherical angles:
        # # (this is the equation of an ellipsoid):
        # ex = a * np.outer(np.cos(u), np.sin(v))
        # ey = b * np.outer(np.sin(u), np.sin(v))
        # ez = c * np.outer(np.ones_like(u), np.cos(v))

        # # Plot:
        # ax2.plot_surface(ex, ey, ez,  rstride=10, cstride=4,\
        #      color=clrs[k],alpha=.4)

        # # Adjustment of the axes, so that they all have the same span:
        # max_radius = max(1/np.sqrt((a,b,c)))
        # for axis in 'xyz':
        #     getattr(ax2, 'set_{}lim'.format(axis))((-max_radius, max_radius))

    ax2.set_theta_zero_location("N")
    #ax2.grid(True)

    #ax1.axvline(2.41, color='slategrey', linestyle='--',label='Random')
 
    if rseed==0: ax1.legend()

    ax1.set_xlabel(r'$\zeta_{cos}$')
    ax1.set_ylabel(r'$\zeta_{\eta}$')

    #ax1.set_xlim([-.4,.3])

    #ax4.set_xlabel(r'$log_{10}(\beta)$')
    ################



for i in range(len(aa)):
    zc_m = np.mean(zeta_cos[i::len(aa)])
    zc_std = np.std(zeta_cos[i::len(aa)],ddof=1)
    ze_m = np.mean(zeta_eta[i::len(aa)])
    ze_std = np.std(zeta_eta[i::len(aa)],ddof=1)

    ax1.errorbar(zc_m,ze_m,xerr=zc_std,yerr=ze_std,color='k',capsize=3)

plt.show()
# %%
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=plt.figaspect(1))  # Square figure
ax = fig.add_subplot(111, projection='3d')


a, b, c = (1,1,1.2)

# Set of all spherical angles:
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)

# Cartesian coordinates that correspond to the spherical angles:
# (this is the equation of an ellipsoid):
x = a * np.outer(np.cos(u), np.sin(v))
y = b * np.outer(np.sin(u), np.sin(v))
z = c * np.outer(np.ones_like(u), np.cos(v))

# Plot:
ax.plot_surface(x, y, z,  rstride=10, cstride=4, color='b',alpha=.4)

# Adjustment of the axes, so that they all have the same span:
max_radius = max(1/np.sqrt((a,b,c)))
for axis in 'xyz':
    getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))

plt.show()

# %%
