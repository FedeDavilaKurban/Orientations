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

def ellipsoid(t, u, a=1, b=1, c=1.1):
    return np.array([a*np.sin(u)*np.cos(t), b*np.sin(u)*np.sin(t), c*np.cos(u)])
domain_t = [0, 2*np.pi]
domain_u = [0, np.pi]

#%%

fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(hspace=0.3, wspace=0.1, bottom=.2)
ax1 = fig.add_subplot(2,2,(3,4)) #etas
ax2 = fig.add_subplot(2,2,1) #cos
#ax4 = fig.add_subplot(2,2,3) #beta
ax3 = fig.add_subplot(2,2,2,projection='polar')

bins = np.linspace(1.9, 4.8, 20)

clrs = ['#402218', '#865439', '#C68B59', '#D7B19D']

Nran = 500
Netas = 300*10

np.random.seed(0)
for k, c in enumerate(np.arange(0.7, 1, 0.1)):
    print(c)

    f = partial(ellipsoid, c=c)

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

    ax1.hist(eta, bins=bins, histtype='step', linewidth=2, 
            color=clrs[k], density=True, label=f'{c:.1f}')
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

    # cos = np.concatenate((-cos, cos))
    #lmbda = np.arctan(beta)
    #lmbda = np.concatenate((lmbda, lmbda+np.pi))

    #ax4.hist(np.log10(beta), histtype='step', color=clrs[k], linewidth=2, density=True)
    ax2.hist(cos, histtype='step', color=clrs[k], linewidth=2, density=True)

ax1.axvline(2.41, color='slategrey', linestyle='--',label='Random')
ax2.axhline(1, color='slategrey', linestyle='--')
#ax4.axvline(0, color='slategrey', linestyle='--')

ax1.legend()

ax1.set_xlabel(r'$\eta$')
ax2.set_xlabel(r'$cos(\lambda)$')
#ax4.set_xlabel(r'$log_{10}(\beta)$')
################


def ellipsp(t, a, b):
    e2 = 1 - (b/a)**2
    r = np.sqrt(b**2 / (1-e2*np.cos(t)**2))
    return r

t = np.linspace(0, 2*np.pi, 100)

for k, c in enumerate([0.7, 0.8, 0.9, 1.]):

    r = ellipsp(t, c, 1)
    ax3.plot(t, r, clrs[k])

ax3.set_theta_zero_location("N")
ax3.grid(True)

plt.show()
# %%
