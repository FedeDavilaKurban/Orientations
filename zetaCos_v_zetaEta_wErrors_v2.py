"""
Quiero encontrar un bias 0 entre los etas normales y los etas
calculados metiendo ruido gaussiano en las componentes de los vectores
generados aleatoriamente.
Esto lo hago para simular posibles errores observacionales.

Calculo los betas y etas con las componentes originales y las "borroneadas".
El cociente entre estos etas es el bias que introducen los errores observacionales.

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
    for j in range(len(xs)):
        para = zs[j] + np.random.normal(0, stdErr)
        perp = np.sqrt(xs[j]**2+ys[j]**2) + np.random.normal(0, stdErr)
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

# def get_eta_bs_wErrors(beta,cosErr,Nbs=1000):
#     bs_eta = []
#     for _ in range(Nbs):  
#         bs_beta = np.random.choice(beta, size=len(beta))

#         theta = np.arctan(bs_beta) 
#         x = np.cos(theta)
#         dx = abs(np.random.normal(0,cosErr,len(theta))) #The std of the gaussian
#                                                         #is 10% of the cosine mean
#         dbeta = x**(-2)*np.array(dx)

#         #print(np.sum(bs_beta-dbeta<0))
#         n_perp = np.sum(bs_beta-dbeta>1.)
#         n_prll = np.sum(bs_beta+dbeta<1.)
#         bs_eta.append( n_perp / n_prll )
        
#         # n_perp = np.sum(bs_x<0.)
#         # n_prll = np.sum(bs_x>0.)
#         # bs_eta.append( n_prll / n_perp )
    
#     eta = np.mean(bs_eta) 
#     eta_std = np.std(bs_eta, ddof=1)
    
#     return eta, eta_std, bs_eta

#%%
"""
Code parameters
"""
Nrans = [500]
Nbs = 100
nseed = 50
stdErr = .1

cc = np.array([.6,.8,1.])
aa = np.array([1,1,1])
bb = aa

eta_ran = 1/(np.sqrt(2)-1)

#print(f'Nran={Nran}, Nbs={Nbs}, nseed={nseed}, cosErr={cosErr}')
#%%      

for Nran in Nrans:
    
    print(f'N={Nran}')

    zeta_eta = []
    zeta_cos = []

    cos = []
    beta = []
    eta = []
    eta_arrays = []
    eta_stds = []

    cos_wErr = []
    beta_wErr = []

    #These arrays are the ones that take into account
    #the gaussian noise introduced in the components
    eta_wErr = []
    eta_arrays_wErr = []
    eta_stds_wErr = []

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
            
            #cos = get_cos(xs,ys,zs) 
            beta.append( get_beta(xs,ys,zs) )
            #beta = np.tan(np.arccos(cos)) 

            eta_m, eta_std, eta_array = get_eta_bs(beta[-1],Nbs=Nbs) 
            eta.append(eta_m)
            #eta_arrays.append(eta_array)
            eta_stds.append(eta_stds)

            #beta_wErr = np.tan(np.arccos(cos_wErr)) 
            beta_wErr.append( get_beta_wErr(xs,ys,zs,stdErr) )

            eta_m_wErr, eta_std_wErr, eta_array_wErr = get_eta_bs(beta_wErr[-1],Nbs=Nbs) 
            eta_wErr.append(eta_m_wErr)
            #eta_arrays_wErr.append(eta_array_wErr)
            eta_stds_wErr.append(eta_stds_wErr)

    np.save(f'../data/zetaCos_v_zetaEta_wErrors_v2/eta_Nran{Nran}_stdErr{stdErr}',\
        eta)
    np.save(f'../data/zetaCos_v_zetaEta_wErrors_v2/eta_wErr_Nran{Nran}_stdErr{stdErr}',eta_wErr)

    np.save(f'../data/zetaCos_v_zetaEta_wErrors_v2/beta_Nran{Nran}_stdErr{stdErr}',\
        beta)
    np.save(f'../data/zetaCos_v_zetaEta_wErrors_v2/beta_wErr_Nran{Nran}_stdErr{stdErr}',beta_wErr)

#%%

fig = plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = 16
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax = fig.add_subplot(111)

Nrans = [500,1000,5000]

for j, Nran in enumerate(Nrans):
    
    eta = np.load(f'../data/zetaCos_v_zetaEta_wErrors_v2/eta_Nran{Nran}_stdErr{stdErr}.npy')
    eta_wErr = np.load(f'../data/zetaCos_v_zetaEta_wErrors_v2/eta_wErr_Nran{Nran}_stdErr{stdErr}.npy')

    for i in range(len(aa)):
        
        x = np.array(eta[i::len(aa)])
        y = np.array(eta_wErr[i::len(aa)])

        ratio = y/x

        ax.axhline(1,color='k',ls='--')
        ax.axvline(eta_ran,color='k',ls='--')
        
        
        if i==0: label=f'N={Nran}'
        else: label=None
        ax.errorbar(x.mean(),ratio.mean(),xerr=x.std(),yerr=ratio.std(),\
            color=clrs[j],capsize=5,label=label)
        
    ax.legend()
    ax.set_xlabel(r'$\eta$')
    ax.set_ylabel(r'$\eta_\mathrm{with\,errors} / \eta$')
    
    
#ax.text(.025,.95,r'$\langle \eta / \eta_\mathrm{w/error} \rangle =$'+f'{np.mean(ratio):.4f}',\
#                    transform = ax.transAxes)
#plt.savefig(f'../plots/zetaCos_v_zetaEta_wErrors/ratioVSeta_stdErr{stdErr}.jpg')
        

# %%
fig = plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = 16
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax = fig.add_subplot(111)

Nrans = [500]

for j, Nran in enumerate(Nrans):
    
    beta = np.load(f'../data/zetaCos_v_zetaEta_wErrors_v2/beta_Nran{Nran}_stdErr{stdErr}.npy')
    beta_wErr = np.load(f'../data/zetaCos_v_zetaEta_wErrors_v2/beta_wErr_Nran{Nran}_stdErr{stdErr}.npy')

    #beta = np.log10(beta)
    #beta_wErr = np.log10(beta_wErr)

    for i in range(len(aa)):
        
        x = np.array(beta[i::len(aa)])
        y = np.array(beta_wErr[i::len(aa)])
        ratio = y/x
        percErr = abs(y-x)
        
        #ax.axhline(1,color='k',ls='--')
        #ax.axvline(eta_ran,color='k',ls='--')
        
        
        if i==0: label=f'N={Nran}'
        else: label=None
        
        #Probando algunas cosas:
        #ax.hist(np.log10(x))
        #ax.scatter(x,percErr,color='k')
        
        ax.errorbar(x.mean(),percErr.mean(),xerr=x.std(),yerr=percErr.std(),\
            color=clrs[j],capsize=5,label=label)
      
    #ax.set_xscale('log')  
    #ax.set_yscale('log')  

    ax.legend()
    ax.set_xlabel(r'$log_{10}(\beta)$')
    #ax.set_ylabel(r'$log_{10}(\beta_\mathrm{with\,errors}) / log_{10}(\beta)$')
    
    
#ax.text(.025,.95,r'$\langle \eta / \eta_\mathrm{w/error} \rangle =$'+f'{np.mean(ratio):.4f}',\
#                    transform = ax.transAxes)
#plt.savefig(f'../plots/zetaCos_v_zetaEta_wErrors/ratioVSbeta_stdErr{stdErr}.jpg')
        

# %%
