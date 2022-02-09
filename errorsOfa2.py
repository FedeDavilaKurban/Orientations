#%%
import numpy as np
from functools import partial

from numpy.core.fromnumeric import size
from param_tools import r_surface, surface_area
from matplotlib import pyplot as plt

def fits(cos1,y,verbose=False):
    "Perform fits on the residues of the ecdf"
    
    from scipy.optimize import curve_fit

    def func(x,a,b,c,d,e):
        return a \
            + b*np.sin( np.pi*x ) \
            + c*np.sin( 2.*np.pi*x ) \
            + d*np.sin( 3.*np.pi*x ) \
            + e*np.sin( 4.*np.pi*x )
    def dfunc(x,b,c,d,e):
        return np.pi*b*np.cos( np.pi*x ) \
            + 2*np.pi*c*np.cos( 2.*np.pi*x ) \
            + 3*np.pi*d*np.cos( 3.*np.pi*x ) \
            + 4*np.pi*e*np.cos( 4.*np.pi*x )

    x = np.array(cos1)

    coeffs, cov = curve_fit(func, x, y)
    yfit = func(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
    d_yfit = dfunc(x,coeffs[1],coeffs[2],coeffs[3],coeffs[4])

    a1 = coeffs[1]
    a2 = coeffs[2]
    a3 = coeffs[3]
    a4 = coeffs[4]
    if verbose==True: print('a1=',a1,'; a2=',a2,'; a3=',a3,'; a4=',a4)

    del x, cov
    return yfit,d_yfit,a1

def ecdf_residues(cos_list,verbose=False):
    """
    Calculates ECDF and residues. Fits the residues.
    Returns: ecdf, fit, derivative of fit, and coefficient a2
    """
    #import numpy as np
    from statsmodels.distributions.empirical_distribution import ECDF

    #cos_neg=-1*np.array(cos_list)
    #cos1 = np.sort(np.concatenate([cos_list,cos_neg]))
    cos1 = np.sort(cos_list)
    ecdf = ECDF(cos1) # Empirical cumulated distribution function
    y = ecdf(cos1)-(cos1)#+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 

    return cos1,ecdf,y

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

def get_cos_wErr(xs,ys,zs,stdErr):
    
    cos = np.zeros(len(xs))
    for j in range(len(xs)):
        
        para = zs[j] + np.random.normal(0, stdErr)
        perp = np.sqrt(xs[j]**2+ys[j]**2) + np.random.normal(0, stdErr)    
        norm = np.sqrt(para**2+perp**2)
        
        cos[j] = abs(para) / norm
        
    return cos
# %%
"""
Code parameters
"""
Nrans = [5000]
Nbs = 100
nseed = 100
stdErr = .1

cc = np.array([.6,.8,1])
aa = np.array([1,1,1])
bb = aa
#%%

for Nran in Nrans:
    print(Nran)
    for rseed in range(nseed):

        for k, c in enumerate(cc):
            #print(k,a,b,c)

            f = partial(ellipsoid, c=c)

            np.random.seed(rseed)
            x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
            xs = x[0,:]
            ys = x[1,:]
            zs = x[2,:]

            cos = get_cos(xs,ys,zs) 
            newcos, ecdf, residues = ecdf_residues(cos)
            
            a2_bs=[]
            for i in range(Nbs):
                c_bs = np.random.choice(cos,size=len(cos))
                nc_, e_, y_ = ecdf_residues(c_bs)
                yf_, d_yf, a2_ = fits(nc_,y_)
                a2_bs.append(a2_)
            
            a2 = np.mean(a2_bs)
        
            np.save(f'../data/errorsOfa2/a2_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}',\
                a2)
            
            #####################
            #With Erors
            #####################

            cos = get_cos_wErr(xs,ys,zs,stdErr) 
            newcos, ecdf, residues = ecdf_residues(cos)
            
            a2_bs=[]
            for i in range(Nbs):
                c_bs = np.random.choice(cos,size=len(cos))
                nc_, e_, y_ = ecdf_residues(c_bs)
                yf_, d_yf, a2_ = fits(nc_,y_)
                a2_bs.append(a2_)
            
            a2_wErr = np.mean(a2_bs)
            #a2_std_wErr = np.std(a2_bs,ddof=1)        
            
            np.save(f'../data/errorsOfa2/a2_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}',\
                a2_wErr)
        
#%%
Nrans=[1000,5000]
nseed = 100

for Nran in Nrans:   
    
    fig = plt.figure(figsize=(8,6))
    plt.rcParams['font.size'] = 16
    clrs = plt.rcParams['axes.prop_cycle'].by_key()['color']
    ax = fig.add_subplot(111)

    if Nran==1000: db=.03
    elif Nran==5000: db=.01
    ax.plot(np.linspace(-.05,.15,100),np.linspace(-.05,.15,100),color='k',ls='--')
    ax.fill_between(np.linspace(-.05,.15,100),np.linspace(-.05,.15,100)-db,\
        np.linspace(-.05,.15,100)+db,color='k',alpha=.5)
    
    for i, c in enumerate(cc):
        for rseed in range(nseed):
        
            a2 = np.load(f'../data/errorsOfa2/a2_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}.npy')
            a2_wErr = np.load(f'../data/errorsOfa2/a2_wErr_Nran{Nran}_stdErr{stdErr}_rseed{rseed}_c{c}.npy')

            x = a2
            y = a2_wErr
            percErr = abs(y-x)/abs(x)
        
            if rseed==0: label=f'c={c}'
            else: label=None
            ax.scatter(x,y,color=clrs[i],label=label)

    ax.set_title(f'N={Nran}')  
    ax.legend()
    ax.set_xlabel(r'$a2$')
    ax.set_ylabel(r'$a2_{err}$')
    #plt.savefig(f'../plots/errorsOfa2/N{Nran}.jpg')
# %%
