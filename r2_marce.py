
#%%
import numpy as np
from matplotlib import pyplot as plt


def f1(N):
    u = np.random.uniform(0, 1, N)
    b = np.tan(np.arccos(u))
    return b


def f2(N):
    X = np.random.normal(loc=0, scale=1, size=(N, 3))
    n =  np.linalg.norm(X, axis=1)
    s_par = X[:,2]
    s_perp = np.linalg.norm(X[:,0:2], axis=1) 
    b = s_perp / s_par
    return b


def f3(N):
    X = np.random.uniform(0, 1, (3*N, 3))
    n =  np.linalg.norm(X, axis=1)
    filtro = n<1
    s_par = X[filtro][:,2]
    s_perp = np.linalg.norm(X[filtro][:,0:2], axis=1) 
    b = s_perp / s_par
    return b

def fb(b):
    return b*pow(1+b**2, -3/2)
fb=np.vectorize(fb)


def Fb(b):
    return 1-pow(1+b**2, -3/2)
Fb=np.vectorize(Fb)

# ===================================================================
 
N = 10000
b = f1(N)

# factor de corrección:
bmin = 0.1
bmax = 10
fac = 1 - Fb(0.1) - (1-Fb(10))
fac = 1/fac

brange_lin = np.linspace(0.1, 10, 500)
fbteo_lin = fb(brange_lin)
brange_log = np.linspace(0.1, 10, 500)
fbteo_log = fb(brange_log)

fig, ax = plt.subplots(2, 2, figsize=(8, 6))


blin = np.linspace(bmin, bmax, 30)
blog = np.logspace(np.log10(bmin), np.log10(bmax), 30)

ax[0,0].hist(b, bins=blin, density=True, color='thistle',
             label='histogram of the\nrandom sample')
ax[0,0].plot(brange_lin, fbteo_lin*fac, color='navy', 
             label='probability density')
ax[0,0].set_xlabel(r'$\beta$', fontsize=14)
ax[0,0].set_ylabel(r'1/N dN/d$\beta$,  f$_B(\beta)$', fontsize=14)
ax[0,0].legend()


ax[1,0].hist(b, bins=blog, density=True, color='thistle')
ax[1,0].plot(brange_lin, fbteo_lin*fac, color='navy')
ax[1,0].set_xlabel(r'$\beta$', fontsize=14)
ax[1,0].set_ylabel(r'1/N dN/d$\beta$,  f$_B(\beta)$', fontsize=14)
                            
ax[0,1].hist(b, bins=blin, density=True, color='thistle')
ax[0,1].plot(brange_log, fbteo_log*fac, color='navy')
ax[0,1].set_xscale('log')
ax[0,1].set_xlabel(r'$\beta$', fontsize=14)
ax[0,1].set_ylabel(r'1/N dN/dlog($\beta$),  f$_B(\beta)$', fontsize=14)

ax[1,1].hist(b, bins=blog, density=True, color='thistle')
ax[1,1].plot(brange_log, fbteo_log*fac, color='navy')
ax[1,1].set_xscale('log')
ax[1,1].set_xlabel(r'$\beta$', fontsize=14)
ax[1,1].set_ylabel(r'1/N dN/dlog($\beta$),  f$_B(\beta)$', fontsize=14)

fig.tight_layout()
plt.show()
#fig.savefig('b2.pdf')











# Histogramas
# ax.hist(np.log10(f1(N)), bins=nbins, ec='cadetblue', fc='none', lw=1.5,
#         histtype='step', density=True, label='Método 1')
# ax.hist(np.log10(f2(N)), bins=nbins, ec='coral', fc='none', lw=1.5,
#         histtype='step', density=True, label='Método 2')
# ax.hist(np.log10(f3(N)), bins=nbins, ec='orchid', fc='none', lw=1.5, 
#         histtype='step', density=True, label='Método 3')
# 
# 
# # teorica
# b = np.logspace(-3.5, 2, 1000)
# fb_teo1 = fb1(b)
# fb_teo2 = fb2(b)
# 
# # randon variable simulation
# bins = np.logspace(-3, 1, 50)
# b_rnd = np.log10(f1(N))
# Hy, _ = np.histogram(b_rnd, bins=bins, density=True)
# deltas = bins[1:]-bins[:-1]
# #A = np.dot(H[0], deltas)
# #Hx = (H[1][1:]+H[1][:-1])/2
# #Hy = H[0]/A
# ax.step(bins[:-1], Hy, where='pre', color='lightseagreen', 
#         label='histogram of simulated beta')
# 
# ax.plot(np.log10(b), fb_teo1, color='k', label='distribution of beta1')
# ax.plot(np.log10(b), fb_teo2, color='green', label='distribution of beta2')
# 
# 
# ax.set_xlim(-2, 3)
# ax.set_xlabel(r'$log_{10}(\beta)$')
# ax.set_ylabel(r'$1/N~~dN/dlog_{10}(\beta),\,f(\beta)$')
# #ylabel(r'$\beta$ probability distribution')
# ax.legend(loc='upper left')    
# ax.grid()
# 
# fig.savefig('b.pdf')
 

# %%
