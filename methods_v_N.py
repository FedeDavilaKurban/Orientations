"""
Código para ver valores de a1 y eta con varios Nran
"""
#%%
import numpy as np
from functools import partial
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
            + np.pi*2*c*np.cos( 2.*np.pi*x ) \
            + np.pi*3*d*np.cos( 3.*np.pi*x ) \
            + np.pi*4*e*np.cos( 4.*np.pi*x )

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
    y = ecdf(cos1)-(cos1)#1+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 

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

def plot_a2N(Nrans,nseed,Nbs,cc,clrs,ax):

    for k in range(len(cc)):
        y=[]
        yerr=[]
        for Nran in Nrans:

            a2 = np.loadtxt(f'../data/zetaCos_v_zetaA2/a1_Nbs{Nbs}_nseed{nseed}_Nran{Nran}.gz')

            y.append( np.mean([a2[i:i+nseed] for i in range(0, len(a2), nseed)][k]) )
            yerr.append( np.std([a2[i:i+nseed] for i in range(0, len(a2), nseed)][k]) )
        x = Nrans

        y=np.array(y)
        ax.plot(x,y,'o-',ms=4,color=clrs[k],label=r'$e^2=$'+f'{1-cc[k]**2:.1f}')
        ax.fill_between(x,y-yerr,y+yerr,alpha=.3,color=clrs[k])
        ax.set_xscale('log')

def plot_cosN(Nrans,nseed,Nbs,cc,clrs,ax):

    for k in range(len(cc)):
        y=[]
        yerr=[]
        for Nran in Nrans:

            cos = np.loadtxt(f'../data/zetaCos_v_zetaA2/cos_Nbs{Nbs}_nseed{nseed}_Nran{Nran}.gz')

            y.append( np.mean([cos[i:i+nseed] for i in range(0, len(cos), nseed)][k]) )
            yerr.append( np.std([cos[i:i+nseed] for i in range(0, len(cos), nseed)][k]) )
        x = Nrans

        y=np.array(y)
        ax.plot(x,y,'o-',ms=4,color=clrs[k],label=r'$e^2=$'+f'{1-cc[k]**2:.1f}')
        ax.fill_between(x,y-yerr,y+yerr,alpha=.3,color=clrs[k])
        ax.set_xscale('log')

def plot_etaN(Nrans,nseed,Nbs,cc,clrs,ax):

    for k in range(len(cc)):
        y=[]
        yerr=[]
        for Nran in Nrans:

            eta = np.loadtxt(f'../data/methods_v_N/eta_Nbs{Nbs}_nseed{nseed}_Nran{Nran}.gz')

            y.append( np.mean([eta[i:i+nseed] for i in range(0, len(eta), nseed)][k]) )
            yerr.append( np.std([eta[i:i+nseed] for i in range(0, len(eta), nseed)][k]) )
        x = Nrans

        y=np.array(y)
        ax.plot(x,y,'o-',ms=4,color=clrs[k],label=r'$e^2=$'+f'{1-cc[k]**2:.1f}')
        ax.fill_between(x,y-yerr,y+yerr,alpha=.3,color=clrs[k])
        ax.set_xscale('log')
# %%
"""
Code parameters
"""
#Nran = 1000
Nbs = 100
nseed = 50

cc = np.array([.6,.8,1])
aa = np.array([1,1,1])
bb = aa
#%%
for Nran in [100,500,1000,5000,10000]:
    print(Nran)
    cos = []
    newcos = []
    ecdf = []
    residues = []
    residues_fit = []
    a2 = []
    a2_std = []
    a2_bs=[]

    beta = []
    eta = []
    eta_ran = 1/(np.sqrt(2)-1)
    eta_ran_std = np.sqrt(28.1421/Nran) #beta = perp/prll

    for k, (a, b, c) in enumerate(zip(aa,bb,cc)):
        f = partial(ellipsoid, a=a, b=b, c=c)
        #print(k)
        for rseed in range(nseed):
            np.random.seed(rseed)
            x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
            xs = x[0,:]
            ys = x[1,:]
            zs = x[2,:]

            cos.append( get_cos(xs,ys,zs) )
            newcos_, ecdf_, y_ = ecdf_residues(cos[-1])
            newcos.append(newcos_)
            ecdf.append(ecdf_)
            residues.append(y_)
            
            yfit_, d_yfit, a2_ = fits(newcos_,y_)
            residues_fit.append(yfit_)
            
            beta.append( get_beta(xs,ys,zs) )
            eta_m, eta_std, eta_array = get_eta_bs(beta[-1],Nbs=Nbs) 
            eta.append(eta_m)

            for i in range(Nbs):
                c_bs = np.random.choice(cos[-1],size=len(cos[-1]))
                nc_, e_, y_ = ecdf_residues(c_bs)
                yf_, d_yf, a2_ = fits(newcos_,y_)
                a2_bs.append(a2_)
            
            a2.append(np.mean([a2_bs[-Nbs:][i:i+Nbs] for i in range(0, len(a2_bs[-Nbs:]), Nbs)],axis=1))
            a2_std.append(np.std([a2_bs[-Nbs:][i:i+Nbs] for i in range(0, len(a2_bs[-Nbs:]), Nbs)],axis=1,ddof=1))

    np.savetxt('../data/zetaCos_v_zetaA2/a1_Nbs{}_nseed{}_Nran{}.gz'\
        .format(Nbs,nseed,Nran),a2)
    np.savetxt('../data/zetaCos_v_zetaA2/a1std_Nbs{}_nseed{}_Nran{}.gz'\
        .format(Nbs,nseed,Nran),a2_std)
    np.savetxt('../data/zetaCos_v_zetaA2/cos_Nbs{}_nseed{}_Nran{}.gz'\
        .format(Nbs,nseed,Nran),newcos)
    np.savetxt('../data/methods_v_N/eta_Nbs{}_nseed{}_Nran{}.gz'\
        .format(Nbs,nseed,Nran),eta)

#%%
from matplotlib.ticker import FormatStrFormatter

fig = plt.figure(figsize=(7,13))
#fig.subplots_adjust(hspace=0.2, wspace=0.4, bottom=.2)
plt.rcParams['font.size'] = 16
fs = 16
clrs = ['#361e15',\
        #'#402218', \
        '#865439', \
        #'#C68B59',\
        '#D7B19D']
#ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(211)
ax1 = fig.add_subplot(212,sharex=ax3)

Nrans = [300,400,500,1000,3000,4000,5000,10000]
plot_a2N(Nrans,nseed,Nbs,cc,clrs,ax1)
ax1.set_ylabel('a1')
ax1.set_xlabel(r'$\mathrm{N_{ran}}$')
ax1.axhline(0,ls=':',color='slategrey',label='Isotropy')
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

# plot_cosN(Nrans,nseed,Nbs,cc,clrs,ax2)
# ax2.set_ylabel(r'$\mathrm{cos}(\lambda)$')
# ax2.set_xlabel(r'$\mathrm{N_{ran}}$')
# ax2.axhline(.5,ls=':',color='slategrey')

plot_etaN(Nrans,nseed,Nbs,cc,clrs,ax3)
ax3.axhline(1/(np.sqrt(2)-1),ls=':',color='slategrey',label='Isotropy')
ax3.set_ylabel(r'$\eta$')
#ax3.set_xlabel(r'$\mathrm{N_{ran}}$')

ax3.legend(loc='upper right',ncol=2,framealpha=.6)

plt.savefig(f'../plots/methods_v_N/plot_Nbs{Nbs}_nseed{nseed}.png',bbox_inches = 'tight')
plt.savefig(f'../plots/methods_v_N/plot_Nbs{Nbs}_nseed{nseed}.pdf',bbox_inches = 'tight')
plt.tight_layout()
plt.show()
# %%
