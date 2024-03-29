"""
Quiero hacer un plot zeta_a2 vs zeta_cos
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
            + b*np.sin( np.pi*(x+1.)/2. ) \
            + c*np.sin( 2.*np.pi*(x+1.)/2. ) \
            + d*np.sin( 3.*np.pi*(x+1.)/2. ) \
            + e*np.sin( 4.*np.pi*(x+1.)/2. )
    def dfunc(x,b,c,d,e):
        return np.pi/2.*b*np.cos( np.pi*(x+1.)/2. ) \
            + np.pi*c*np.cos( 2.*np.pi*(x+1.)/2. ) \
            + np.pi*3./2.*d*np.cos( 3.*np.pi*(x+1.)/2. ) \
            + np.pi*2.*e*np.cos( 4.*np.pi*(x+1.)/2. )

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
    return yfit,d_yfit,a2

def ecdf_residues(cos_list,verbose=False):
    """
    Calculates ECDF and residues. Fits the residues.
    Returns: ecdf, fit, derivative of fit, and coefficient a2
    """
    #import numpy as np
    from statsmodels.distributions.empirical_distribution import ECDF

    cos_neg=-1*np.array(cos_list)
    cos1 = np.sort(np.concatenate([cos_list,cos_neg]))
    ecdf = ECDF(cos1) # Empirical cumulated distribution function
    y = ecdf(cos1)-(cos1+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 

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
# %%
"""
Code parameters
"""
Nran = 1000
Nbs = 100
nseed = 50

cc = np.array([.6,.8,1])
aa = np.array([1,1,1])
bb = aa
#%%
clrs = ['#361e15',\
        #'#402218', \
        '#865439', \
        #'#C68B59',\
        '#D7B19D']

cos = []
newcos = []
ecdf = []
residues = []
residues_fit = []
a2 = []

for k, (a, b, c) in enumerate(zip(aa,bb,cc)):
    f = partial(ellipsoid, a=a, b=b, c=c)

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
        a2.append(a2_)

#%%
k=0
for i in range(nseed):     
    plt.plot(newcos[i],residues[i],color=clrs[k],alpha=.2)
k=1
for i in range(nseed):
    plt.plot(newcos[i+nseed],residues[i+nseed],color=clrs[k],alpha=.2)
k=2
for i in range(nseed):
    plt.plot(newcos[i+nseed*2],residues[i+nseed*2],color=clrs[k],alpha=.2)
#plt.show()
#%%
def plot_fitcurves(newcos,residues_fit,ax):
    xm = np.linspace(-1,1,100)
    k = 0
    ys = []
    for i in range(nseed):
        xx = newcos[i].flatten()
        yy = residues_fit[i].flatten()
        ys.append(np.interp(xm,xx,yy))
    ym = np.mean(ys,axis=0)
    y_s1 = np.std(ys,axis=0,ddof=1)
    ax.fill_between(xm,ym-y_s1*3,ym+y_s1*3,color=clrs[k],alpha=.5)
    ax.plot(xm,ym,color=clrs[k],lw=3)

    k=1
    ys = []
    for i in range(nseed):
        xx = newcos[i+nseed]
        yy = residues_fit[i+nseed]
        ys.append(np.interp(xm,xx,yy))
    ym = np.mean(ys,axis=0)
    y_s1 = np.std(ys,axis=0,ddof=1)
    ax.fill_between(xm,ym-y_s1*3,ym+y_s1*3,color=clrs[k],alpha=.5)
    ax.plot(xm,ym,color=clrs[k],lw=3)

    k=2
    ys = []
    for i in range(nseed):
        xx = newcos[i+nseed*2]
        yy = residues_fit[i+nseed*2]
        ys.append(np.interp(xm,xx,yy))
    ym = np.mean(ys,axis=0)
    y_s1 = np.std(ys,axis=0,ddof=1)
    ax.fill_between(xm,ym-y_s1*3,ym+y_s1*3,color=clrs[k],alpha=.5)
    ax.plot(xm,ym,color=clrs[k],lw=3)

def plot_allcurves(newcos,residues,clrs,nseed,ax):
    k=0
    for i in range(nseed):     
        ax.plot(newcos[i],residues[i],color=clrs[k],alpha=.2)
    k=1
    for i in range(nseed):
        ax.plot(newcos[i+nseed],residues[i+nseed],color=clrs[k],alpha=.2)
    k=2
    for i in range(nseed):
        ax.plot(newcos[i+nseed*2],residues[i+nseed*2],color=clrs[k],alpha=.2)

def plot_eye(aa,bb,cc,clrs,ax):
    #EYE OF SAURON
    def ellipsp(t, a, c):
        e2 = 1 - (c/a)**2
        r = np.sqrt(c**2 / (1-e2*np.cos(t)**2))
        return r
    for k, (a, b, c) in enumerate(zip(aa,bb,cc)):
        t_ = np.linspace(0, 2*np.pi, 100)
        r_ = ellipsp(t_, c, a)
        if c<=a: e2=1-(c/a)**2
        if c>a: e2=1-(a/c)**2 
        ax0.plot(t_, r_, clrs[k],label=r'$e^2=$'+f'{e2:.1f}')
    ax.set_theta_zero_location("N")
    ax.legend(bbox_to_anchor=(-.2, 1.))

def plot_zetas(a2,cos,nseed,cc,clrs,ax):
    # def get_a2mean(a2,nseed):
    #     a2_m=[]
    #     parsed_a2=[a2[i:i+nseed] for i in range(0, len(a2), nseed)]
    #     for a2_chunk in parsed_a2:
    #         a2_m.append( np.mean(a2_chunk) )
    #     return np.array(a2_m)
    def get_a2stdran(a2,nseed):
        a2_ran=[a2[i:i+nseed] for i in range(0, len(a2), nseed)][-1]
        return np.array(np.std(a2_ran))
    # def get_cosmean(cos,nseed):
    #     cos_m = []
    #     parsed_cos=[cos[i:i+nseed] for i in range(0, len(cos), nseed)]
    #     for cos_chunk in parsed_cos:
    #         cos_m.append( np.mean(cos_chunk) )
    #     return np.array(cos_m)
    for k in range(len(cc)):
        y = np.array([a2[i:i+nseed] for i in range(0, len(a2), nseed)][k])
        x = np.mean([cos[i:i+nseed] for i in range(0, len(cos), nseed)][k],axis=1)
        #plt.scatter(x,y,color=clrs[k])
        zy = y/get_a2stdran(a2,nseed)
        zx = (x-0.5)/np.sqrt(12)
        ax.scatter(zx,zy,color=clrs[k])

fig = plt.figure(figsize=(13,17))
fig.subplots_adjust(hspace=0.6, wspace=0.4, bottom=.2)
plt.rcParams['font.size'] = 16
fs = 16

ncol = 2
nrow = 4

ax0 = fig.add_subplot(nrow,ncol,2,projection='polar')
ax1 = fig.add_subplot(nrow,ncol,1) #cos
ax2 = fig.add_subplot(nrow,ncol,3) #ECDF
ax3 = fig.add_subplot(nrow,ncol,4) #Residues 
ax4 = fig.add_subplot(nrow,ncol,5) #eta vs cos
ax5 = fig.add_subplot(nrow,ncol,6) #zetas

#AXIS LABELS
ax1.set_xlabel(r'$\cos(\lambda)$')
ax2.set_xlabel(r'$\cos(\lambda)$')
ax2.set_ylabel('ECDF')
ax3.set_xlabel(r'$\cos(\lambda)$'+' [mirrored]')
ax3.set_ylabel('Residues')
ax3.set_xlabel(r'$\cos(\lambda)$'+' [mirrored]')
ax4.set_ylabel('Residues Fits')
ax5.set_xlabel(r'$\zeta_{cos}=(\langle cos\rangle -cos_{ran})/\sigma_{cos}}$')
ax5.set_ylabel(r'$\zeta_{a2}=(a2-a2_{ran})/\sigma_{a2,ran}}$')

#VERTICAL/HORIZONTAL LINES
ax1.axhline(1,ls='--',lw=1,color='slategrey')
ax5.axvline(0,ls=':',color='slategrey')
ax5.axhline(0,ls=':',color='slategrey')

bins=8
clrs = ['#361e15',\
        #'#402218', \
        '#865439', \
        #'#C68B59',\
        '#D7B19D']

#from pylab import *
#clrs = cm.get_cmap('seismic', len(cc))  # matplotlib color palette name, n colors

plot_eye(aa,bb,cc,clrs,ax0)
plot_allcurves(newcos,residues,clrs,nseed,ax3)
plot_fitcurves(newcos,residues_fit,ax4)
for i in range(len(cc)):
    x = np.sort(cos[i*nseed])
    ax1.hist(x, bins=bins, histtype='step', \
        color=clrs[i], density=True, lw=2)
    ax2.step(x, ecdf[i*nseed](x), color=clrs[i], lw=2)
    
    ytext=np.arange(len(cc))/10+.1
    s = r'$\langle cos(\lambda)\rangle=$'+f'{cos[i*nseed].mean():.2f}'+\
        r'$\pm$'+f'{cos[i*nseed].std():.2f}'
    ax1.text(.1,ytext[i],s,\
        color=clrs[i],transform = ax1.transAxes)
plot_zetas(a2,cos,nseed,cc,clrs,ax5)

#TEXT
ax1.text(.6,.9,'Nran={}'.format(Nran),\
    fontsize=fs, alpha=.7, color='k', transform = ax1.transAxes)
ax3.text(.1,.7,'{} realizations \n for each {}'.format(nseed,r'$e^2$'),\
    fontsize=fs,alpha=.7,color='k',transform = ax3.transAxes)
ax4.text(.4,.1,'a2{}amplitude of the fit'.format(r'$\simeq$'),\
    fontsize=fs-1,alpha=.7,color='k',transform = ax4.transAxes)

plt.tight_layout()
# %%
