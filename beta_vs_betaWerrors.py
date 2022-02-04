#%%
import numpy as np
from functools import partial
from param_tools import r_surface, surface_area
from matplotlib import pyplot as plt

def get_beta(xs,ys,zs):

    beta = np.zeros(len(xs)) 
    for j in range(len(xs)):
        para = zs[j]
        perp = np.sqrt(xs[j]**2+ys[j]**2)
        beta[j] = perp / abs(para)

    return beta

def get_beta_wErr(xs,ys,zs,stdErr):
    """
    Calculate beta with gaussian errors 
    in the components "para" and "perp"

    Args:
        xs (array): x-component of vector
        ys (array): y-component of vector
        zs (array): z-component of vector
        stdErr (scalar): Stand. dev. of gaussian distribution for error

    Returns:
        beta [array]: perp / abs(para)
    """


    beta = np.zeros(len(xs)) 
    for j in range(len(xs)):
        para = zs[j] + np.random.normal(0, stdErr)
        perp = np.sqrt(xs[j]**2+ys[j]**2) + np.random.normal(0, stdErr)
        beta[j] = abs(perp) / abs(para)

    return beta

def ellipsoid(t, u, a=1, b=1, c=1):
    return np.array([a*np.sin(u)*np.cos(t), b*np.sin(u)*np.sin(t), c*np.cos(u)])
domain_t = [0, 2*np.pi]
domain_u = [0, np.pi]


#%%
Nran = 1000
stdErr = .1
nseed = 10

cc = np.array([.6,.8,1.])
aa = np.array([1,1,1])
bb = aa

#%%



rseeds = np.arange(4,4+nseed)
for rseed in range(nseed):
    
    beta = []
    beta_wErr = []
    
    for k, c in enumerate(cc):

        f = partial(ellipsoid, c=c)

        np.random.seed(rseed)
        x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
        xs = x[0,:]
        ys = x[1,:]
        zs = x[2,:]

        beta.append( get_beta(xs,ys,zs) )
        beta_wErr.append( get_beta_wErr(xs,ys,zs,stdErr) )
        # plt.hist(np.log10(beta[-1]),histtype='step')
        # plt.hist(np.log10(beta_wErr[-1]),histtype='step')
        # plt.title(f'c={c}')
        # plt.show()
        
        np.save(f'../data/beta_vs_betaWerrors/beta_Nran{Nran}_stdErr{stdErr}_rseed{rseed}',\
            beta)
        np.save(f'../data/beta_vs_betaWerrors/beta_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}',\
            beta_wErr)


# %%
fig = plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = 16
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax = fig.add_subplot(111)

rseeds = np.arange(4,4+nseed)
for rseed in range(nseed)[:1]:
    
    beta = np.load(f'../data/beta_vs_betaWerrors/beta_Nran{Nran}_stdErr{stdErr}_rseed{rseed}.npy')
    beta_wErr = np.load(f'../data/beta_vs_betaWerrors/beta_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}.npy')

    for i in range(len(cc)):
        if rseed==0: label=f'c={cc[i]}'
        else: label=None
        
        ax.hist(np.log10(beta[i]),histtype='step',label=label,color=clrs[i])
        ax.hist(np.log10(beta_wErr[i]),histtype='step',label=label,color=clrs[i],lw=2)

ax.legend()
plt.show()
# %%

fig = plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = 16
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax = fig.add_subplot(111)

rseeds = np.arange(4,4+nseed)
for rseed in range(nseed):
    
    beta = np.load(f'../data/beta_vs_betaWerrors/beta_Nran{Nran}_stdErr{stdErr}_rseed{rseed}.npy')
    beta_wErr = np.load(f'../data/beta_vs_betaWerrors/beta_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}.npy')

    for i in range(len(cc)):
        if rseed==0: label=f'c={cc[i]}'
        else: label=None
        
        x = beta[i]
        y = beta_wErr[i]
        percErr = abs(y-x)/x
      
        logx = np.log10(beta[i])
        logy = np.log10(beta_wErr[i])
        
        ax.scatter(logx,percErr,color=clrs[i],label=label)
        #ax.errorbar(logx.mean(),percErr.mean(),xerr=logx.std(),yerr=percErr.std(),\
        #    label=label,color=clrs[i])

ax.legend()
plt.show()
# %%
    