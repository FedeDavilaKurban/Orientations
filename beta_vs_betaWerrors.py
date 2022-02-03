"""
Quiero ver el bias entre los betas normales y los betas
calculados metiendo ruido gaussiano en las componentes de los vectores
generados aleatoriamente.
Esto lo hago para simular posibles errores observacionales.

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
    err = np.random.normal(0, stdErr)
    for j in range(len(xs)):
        para = zs[j] + err
        perp = np.sqrt(xs[j]**2+ys[j]**2) + err
        beta[j] = abs(perp) / abs(para)

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
Nrans = [1000]
Nbs = 100
nseed = 50
stdErr = .01

cc = np.array([.6,.8,1.])
aa = np.array([1,1,1])
bb = aa

eta_ran = 1/(np.sqrt(2)-1)

#print(f'Nran={Nran}, Nbs={Nbs}, nseed={nseed}, cosErr={cosErr}')
#%%      
fig = plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = 16
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax = fig.add_subplot(111)

for Nran in Nrans:
    
    print(f'N={Nran}')

    beta = []
    beta_wErr = []

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
    
            beta.append( get_beta(xs,ys,zs) )
            beta_wErr.append( get_beta_wErr(xs,ys,zs,stdErr) )

            x = np.log10(beta[-1])
            y = np.log10(beta_wErr[-1])
            ratio = np.array(y/x)
            
            ax.hist(np.log10(np.array(beta)))
            #ax.scatter(x,ratio)
            #ax.errorbar(x.mean(),ratio.mean(),xerr=x.std(),yerr=ratio.std(),\
            #    color=clrs[k],capsize=5)#,label=label)
            
            
    



# %%

plt.hist(np.log10(beta))

# %%
