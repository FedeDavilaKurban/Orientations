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

def get_eta_bs(x,Nbs=1000):
    """
    Calculate eta from beta

    Args:
        x (array): log10(beta)
        Nbs (int, optional): Number of bootstrap resamplings. 
                                            Defaults to 1000.

    Returns:
       eta (float): mean of boostrap etas
       eta_std (float): std of boostrap etas
       bs_eta (array): all boostrap etas

    """
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

def ellipsoid(t, u, a=1, b=1, c=1):
    return np.array([a*np.sin(u)*np.cos(t), b*np.sin(u)*np.sin(t), c*np.cos(u)])
domain_t = [0, 2*np.pi]
domain_u = [0, np.pi]


#%%
Nran = 500
stdErr = .1
nseed = 10
Nbs = 100

cc = np.array([.6,.8,1.])
aa = np.array([1,1,1])
bb = aa

#%%

for rseed in range(nseed):
    
    for k, c in enumerate(cc):

        beta = []
        beta_wErr = []
        f = partial(ellipsoid, c=c)

        np.random.seed(rseed)
        x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
        xs = x[0,:]
        ys = x[1,:]
        zs = x[2,:]

        beta.append( get_beta(xs,ys,zs) )
        beta_wErr.append( get_beta_wErr(xs,ys,zs,stdErr) )
        
        eta_m, eta_std, eta_array = get_eta_bs(beta[-1],Nbs=Nbs) 
        eta_m_wErr, eta_std_wErr, eta_array_wErr = get_eta_bs(beta_wErr[-1],Nbs=Nbs) 

        
        np.save(f'../data/errorsOfBetaAndEta/beta_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}',\
            beta)
        np.save(f'../data/errorsOfBetaAndEta/beta_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}',\
            beta_wErr)
        np.save(f'../data/errorsOfBetaAndEta/eta_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}',\
            eta_m)
        np.save(f'../data/errorsOfBetaAndEta/eta_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}',\
            eta_m_wErr)


# %%

fig = plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = 16
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax = fig.add_subplot(111)

c=.6
for rseed in range(nseed):
    
    beta = np.load(f'../data/errorsOfBetaAndEta/beta_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}.npy')
    beta_wErr = np.load(f'../data/errorsOfBetaAndEta/beta_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}.npy')

    print(len(beta[0]),len(beta_wErr[0]))

    x = beta
    y = beta_wErr
    percErr = abs(y-x)/x
    
    logx = np.log10(beta)
    logy = np.log10(beta_wErr)
    
    ax.scatter(logx,percErr,color=clrs[0])
# %%

fig = plt.figure(figsize=(8,6))
plt.rcParams['font.size'] = 16
clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
ax = fig.add_subplot(111)

Nran = 5000
for i, c in enumerate(cc):
    for rseed in range(nseed):
    
        eta = np.load(f'../data/errorsOfBetaAndEta/eta_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}.npy')
        eta_wErr = np.load(f'../data/errorsOfBetaAndEta/eta_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}.npy')

        #print(len(beta[0]),len(beta_wErr[0]))

        x = eta
        y = eta_wErr
        percErr = abs(y-x)/x
        
        # logx = np.log10(beta)
        # logy = np.log10(beta_wErr)
        
        if rseed==0: label=f'c={c}'
        else: label=None
        ax.scatter(x,percErr,color=clrs[i],label=label)

ax.set_title(f'N={Nran}')  
ax.legend()
ax.set_xlabel(r'$\eta$')
ax.set_ylabel('Percentual Error')
plt.savefig(f'../plots/errorsOfBetaAndEta/N{Nran}.jpg')

# %%
