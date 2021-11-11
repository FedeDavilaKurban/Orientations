#%%
import matplotlib.pyplot as plt
import scipy.stats
import numpy as np
from astropy.io import ascii

def get_random_vec(R, N):

    phi = 2*np.pi*np.random.random(N)
    costheta = 1-2*np.random.random(N)
    u = np.random.random(N)

    theta = np.arccos( costheta )
    r = R * np.cbrt( u )

    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    
    return x,y,z

def get_beta_random(N_beta):

    sx, sy, sz = get_random_vec(1.,N_beta)

    s_randvec = np.column_stack((sx,sy,sz))

    s_norms = np.linalg.norm(s_randvec, axis=1)

    sx /= s_norms
    sy /= s_norms
    sz /= s_norms

    sp = np.sqrt(sx**2+sy**2)

    b = sp/np.abs(sz)

    return b

def get_eta_random(N_eta,N_beta):
    
    eta_ran = np.zeros(N_eta)

    for i in range(N_eta):

        b = get_beta_random(N_beta)
        eta_ran[i] = len(b[b>1]) / len(b[b<1])

    return eta_ran

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

rmin, rmax = 1.1,1.2

# %%
"""
B vs Vrad & Vtra
"""

minradV = 7.
maxradV = 0.
sec = 3
rmin, rmax = 1.1,1.2

plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 2, constrained_layout=True, sharex=False, sharey=True)

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))

    for v, ax_i, in zip([vel['vrad'],vel['vtra']],ax):

        p50 = np.percentile(v,50)
        v_lo = v[v<p50]
        v_hi = v[v>p50]
        b_lo = beta[v<p50]
        b_hi = beta[v>p50]

        m,b,rvalue,pvalue,std = scipy.stats.linregress(v,np.log10(beta)) 

        ax_i.scatter(v_lo,np.log10(b_lo),s=.5)
        ax_i.scatter(v_hi,np.log10(b_hi),s=.5)

        if vtype=='a': vtitle='All Voids'
        elif vtype=='r': vtitle='R Voids'
        elif vtype=='s': vtitle='S Voids'
        ax_i.set_title('{}'.format(vtitle))

        ax_i.plot(v,m*v+b,ls=':',color='k',label='r_lr={:.2f}'.format(rvalue))
        #ax_i.legend()

    ax[0].set_ylabel(r'$log_{10}(\beta)$')

axs[2,0].set_xlabel(r'$V_\mathrm{rad}$')
axs[2,1].set_xlabel(r'$V_\mathrm{tra}$')

plt.savefig('../plots/vel/BvsVel.jpg')
# %%
"""
Eta vs Low/High Vrad and Vtra
"""
minradV=7.
maxradV=0.
sec=3
rmin, rmax = 1.1,1.2


plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 2, constrained_layout=True, sharex=True, sharey=True)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for vtype, ax in zip(['a','r','s'],axs):
    beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

    vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))

    for v, color, fmt, ax_i, in zip([vel['vrad'],vel['vtra']],colors[:3],['o','x'],ax):

        p30 = np.percentile(v,30)
        p50 = np.percentile(v,50)
        p70 = np.percentile(v,70)

        lomask = v<p50
        himask = v>p50

        v1 = v[lomask]
        v2 = v[himask]
        #print(np.mean(v1),np.mean(v2))

        beta1 = beta[lomask]
        beta2 = beta[himask]

        logbeta1 = np.log10(beta1)
        logbeta2 = np.log10(beta2)

        eta1, eta1_std = get_eta_bs(logbeta1)
        eta2, eta2_std = get_eta_bs(logbeta2)
        print(eta1,eta2)
        ax_i.errorbar(1,eta1,eta1_std,fmt=fmt, capsize=5)
        ax_i.errorbar(2,eta2,eta2_std,fmt=fmt, capsize=5)
        ax_i.hlines(2.41,.5,2.5,ls='-.',color='k')

        if vtype=='a': vtitle='All Voids'
        elif vtype=='r': vtitle='R Voids'
        elif vtype=='s': vtitle='S Voids'
        ax_i.set_title(vtitle)
    ax[0].set_ylabel(r'$\eta$')

axs[2,1].set_xlim([0.,3.])
axs[2,0].set_xlabel(r'$V_\mathrm{rad}$')
axs[2,1].set_xlabel(r'$V_\mathrm{tra}$')
axs[2,1].set_xticks([1,2])
axs[2,1].set_xticklabels(['Low','High'])

plt.savefig('../plots/vel/eta_vs_lohiVel.jpg')
# %%
"""
Eta vs R for LowVrad
"""
plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15

minradV=7.
maxradV=0.

rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

for sec in [123]:
    print(sec)
    eta = []
    eta_std = []

    eta_random_mean = []
    eta_random_var = []
    eta_random_std = []
    
    n_gal = []
    for vtype in ['a','r','s']:
        print(vtype)
        
        for rmin,rmax in zip(r1,r2):
            
            beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

            vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

            p50 = np.percentile(vrad,50)
            beta = beta[vrad<p50]

            x = np.log10(beta.data)

            eta_, eta_std_ = get_eta_bs(x)
 
            eta.append( eta_ )
            eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

            # Obtain mean and var of control samples
            N = len(beta)
        
            eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
            eta_random_mean.append( np.mean(eta_random) )
            eta_random_var.append( np.var(eta_random,ddof=1) )
            eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )



    nbins = len(r1)
    x = r1


    ya = eta[:nbins]
    yr = eta[nbins:2*nbins]
    ys = eta[-nbins:]

    ya_err = eta_std[:nbins]
    yr_err = eta_std[nbins:2*nbins]
    ys_err = eta_std[-nbins:]

    ####################
    eta_random_mean = np.array(eta_random_mean)
    eta_random_var = np.array(eta_random_var)
    eta_random_std = np.array(eta_random_std)

    yran_mean_a = eta_random_mean[:nbins]
    yran_mean_r = eta_random_mean[nbins:2*nbins]
    yran_mean_s = eta_random_mean[-nbins:]

    yran_var_a = eta_random_var[:nbins]
    yran_var_r = eta_random_var[nbins:2*nbins]
    yran_var_s = eta_random_var[-nbins:]

    yran_std_a = eta_random_std[:nbins]
    yran_std_r = eta_random_std[nbins:2*nbins]
    yran_std_s = eta_random_std[-nbins:]
    ####################

    fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=False)

    axs[0].set_title('section {}'.format(sec))


    for ax, y, yerr, yran_mean, yran_err, label, in zip(axs,\
                                    (ya,yr,ys),\
                                    (ya_err,yr_err,ys_err),\
                                    (yran_mean_a,yran_mean_r,yran_mean_s),\
                                    (yran_std_a,yran_std_r,yran_std_s),('All','Rising','Shell')):

        # Theoretical n1/n2 value for random spins
        ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles='-.')
        ax.plot(x,yran_mean,c='k',alpha=.7)

        ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.3, color='k')
        ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.3, color='k')

        ax.errorbar(x,y,yerr=yerr,label=label)
        #print(y)
        ax.legend()

        # plt.ylabel(r'$n_1/n_2$')
        # plt.xlabel('R/Rv')
        # plt.legend()

    x_ticks_labels = []
    for i in range(nbins):
        x_ticks_labels.append( '{:.1f}-{:.1f}'.format(r1[i],r2[i]) )
    axs[2].set_xticks(x)
    axs[2].set_xticklabels(x_ticks_labels)

    #plt.savefig('../plots/BetaEta_vs_MassVel/Eta_LoVrad_allvtypes_sec{}_.jpg'.format(sec))
# %%
"""
Eta vs R for LoHiVrad
"""
plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15

minradV=7.
maxradV=0.

rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

for vfilter in ['hi','lo']:
    print(vfilter)

    for sec in [3,123]:
        print(sec)
        eta = []
        eta_std = []

        eta_random_mean = []
        eta_random_var = []
        eta_random_std = []
        
        n_gal = []
        for vtype in ['a','r','s']:
            print(vtype)
            
            for rmin,rmax in zip(r1,r2):
                
                beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

                vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

                p50 = np.percentile(vrad,50)
                
                if vfilter=='hi': beta = beta[vrad>p50]
                elif vfilter=='lo': beta = beta[vrad<p50]

                x = np.log10(beta.data)

                eta_, eta_std_ = get_eta_bs(x)
    
                eta.append( eta_ )
                eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

                # Obtain mean and var of control samples
                N = len(beta)
            
                eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
                eta_random_mean.append( np.mean(eta_random) )
                eta_random_var.append( np.var(eta_random,ddof=1) )
                eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )



        nbins = len(r1)
        x = (r1+r2)/2


        ya = eta[:nbins]
        yr = eta[nbins:2*nbins]
        ys = eta[-nbins:]

        ya_err = eta_std[:nbins]
        yr_err = eta_std[nbins:2*nbins]
        ys_err = eta_std[-nbins:]

        ####################
        eta_random_mean = np.array(eta_random_mean)
        eta_random_var = np.array(eta_random_var)
        eta_random_std = np.array(eta_random_std)

        yran_mean_a = eta_random_mean[:nbins]
        yran_mean_r = eta_random_mean[nbins:2*nbins]
        yran_mean_s = eta_random_mean[-nbins:]

        yran_var_a = eta_random_var[:nbins]
        yran_var_r = eta_random_var[nbins:2*nbins]
        yran_var_s = eta_random_var[-nbins:]

        yran_std_a = eta_random_std[:nbins]
        yran_std_r = eta_random_std[nbins:2*nbins]
        yran_std_s = eta_random_std[-nbins:]
        ####################

        fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

        if sec==3: 
            if vfilter=='hi': title='H-H Galaxies - High {}'.format(r'$V_\mathrm{rad}$')
            if vfilter=='lo': title='H-H Galaxies - Low {}'.format(r'$V_\mathrm{rad}$')

        elif sec==123: 
            if vfilter=='hi': title='High Spin Galaxies - High {}'.format(r'$V_\mathrm{rad}$')
            if vfilter=='lo': title='High Spin Galaxies - Low {}'.format(r'$V_\mathrm{rad}$')

        for ax, y, yerr, yran_mean, yran_err, label, in zip(axs,\
                                        (ya,yr,ys),\
                                        (ya_err,yr_err,ys_err),\
                                        (yran_mean_a,yran_mean_r,yran_mean_s),\
                                        (yran_std_a,yran_std_r,yran_std_s),('All','Rising','Shell')):

            # Theoretical n1/n2 value for random spins
            ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
            ax.plot(x,yran_mean,c='k',alpha=.7)

            ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.3, color='k')
            ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.3, color='k')
            ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.3, color='k')

            ax.errorbar(x,y,yerr=yerr,label=label)
            #print(y)
            #ax.legend()
            ax.set_ylabel(r'$\eta$')
            if sec==3:  
                ax.set_ylim([1.5,3.5])
            else: ax.set_ylim([2.,2.8])

        axs[2].set_xlabel('R/Rv')

        if sec==3: ytext=3.15
        elif sec==123: ytext=2.67
        axs[0].text(1.25,ytext,'All Voids')
        axs[1].text(1.25,ytext,'Rising Voids')
        axs[2].text(1.25,ytext,'Shell Voids')

        axs[0].set_title(title)
        axs[0].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[1].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[2].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

        # x_ticks_labels = []
        # for i in range(nbins):
        #     x_ticks_labels.append( '{:.1f}-{:.1f}'.format(r1[i],r2[i]) )
        #axs[2].set_xticklabels(x_ticks_labels)

        plt.savefig('../plots/vel/etavsr_{}Vrad_sec{}.jpg'.format(vfilter,sec))

#%%
"""
Eta vs R for LoHiVrad -- 2nd version (no "all voids")
"""
plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15

minradV=7.
maxradV=0.

rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

for vfilter in ['hi','lo']:
    print(vfilter)

    for sec in [3,123]:
        print(sec)
        eta = []
        eta_std = []

        eta_random_mean = []
        eta_random_var = []
        eta_random_std = []
        
        n_gal = []
        for vtype in ['r','s']:
            print(vtype)
            
            for rmin,rmax in zip(r1,r2):
                
                beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

                vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

                p50 = np.percentile(vrad,50)
                
                if vfilter=='hi': beta = beta[vrad>p50]
                elif vfilter=='lo': beta = beta[vrad<p50]

                x = np.log10(beta.data)

                eta_, eta_std_ = get_eta_bs(x)
    
                eta.append( eta_ )
                eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

                # Obtain mean and var of control samples
                N = len(beta)
            
                eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
                eta_random_mean.append( np.mean(eta_random) )
                eta_random_var.append( np.var(eta_random,ddof=1) )
                eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )



        nbins = len(r1)
        x = (r1+r2)/2


        yr = eta[:nbins]
        #yr = eta[nbins:2*nbins]
        ys = eta[-nbins:]

        yr_err = eta_std[:nbins]
        #yr_err = eta_std[nbins:2*nbins]
        ys_err = eta_std[-nbins:]

        ####################
        eta_random_mean = np.array(eta_random_mean)
        eta_random_var = np.array(eta_random_var)
        eta_random_std = np.array(eta_random_std)

        yran_mean_r = eta_random_mean[:nbins]
        #yran_mean_r = eta_random_mean[nbins:2*nbins]
        yran_mean_s = eta_random_mean[-nbins:]

        yran_var_r = eta_random_var[:nbins]
        #yran_var_r = eta_random_var[nbins:2*nbins]
        yran_var_s = eta_random_var[-nbins:]

        yran_std_r = eta_random_std[:nbins]
        #yran_std_r = eta_random_std[nbins:2*nbins]
        yran_std_s = eta_random_std[-nbins:]
        ####################

        fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=True)

        if sec==3: 
            if vfilter=='hi': title='H-H Galaxies - High {}'.format(r'$V_\mathrm{rad}$')
            if vfilter=='lo': title='H-H Galaxies - Low {}'.format(r'$V_\mathrm{rad}$')

        elif sec==123: 
            if vfilter=='hi': title='High Spin Galaxies - High {}'.format(r'$V_\mathrm{rad}$')
            if vfilter=='lo': title='High Spin Galaxies - Low {}'.format(r'$V_\mathrm{rad}$')

        for ax, y, yerr, yran_mean, yran_err, label, in zip(axs,\
                                        (yr,ys),\
                                        (yr_err,ys_err),\
                                        (yran_mean_r,yran_mean_s),\
                                        (yran_std_r,yran_std_s),('Rising','Shell')):

            # Theoretical n1/n2 value for random spins
            ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
            ax.plot(x,yran_mean,c='k',alpha=.7)

            ax.fill_between(x, yran_mean-yran_err, yran_mean+yran_err, alpha=.15, color='k')
            ax.fill_between(x, yran_mean-2*yran_err, yran_mean+2*yran_err, alpha=.15, color='k')
            ax.fill_between(x, yran_mean-3*yran_err, yran_mean+3*yran_err, alpha=.15, color='k')

            ax.errorbar(x,y,yerr=yerr,label=label)
            #print(y)
            #ax.legend()
            ax.set_ylabel(r'$\eta$')
            if sec==3:  
                ax.set_ylim([1.5,3.5])
            else: ax.set_ylim([2.,2.8])

        axs[1].set_xlabel('R/Rv')

        if sec==3: ytext=3.15
        elif sec==123: ytext=2.67
        #axs[0].text(1.25,ytext,'All Voids')
        axs[0].text(1.25,ytext,'Rising Voids')
        axs[1].text(1.25,ytext,'Shell Voids')

        axs[0].set_title(title)
        axs[0].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[1].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        #axs[2].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

        # x_ticks_labels = []
        # for i in range(nbins):
        #     x_ticks_labels.append( '{:.1f}-{:.1f}'.format(r1[i],r2[i]) )
        #axs[2].set_xticklabels(x_ticks_labels)

        plt.savefig('../plots/vel/etavsr_{}Vrad_sec{}_v2.jpg'.format(vfilter,sec))

#%%
"""
Eta vs R for LoHiVrad -- 3rd version (eta/sigma(eta))
"""

plt.rcParams['figure.figsize'] = (6, 8)
plt.rcParams['font.size'] = 15

minradV=7.
maxradV=0.

rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)


for sec in [3]:
    print(sec)

    fig, axs = plt.subplots(2, 1, constrained_layout=True, sharex=True, sharey=True)

    for vfilter in ['hi','lo']:
        print(vfilter)
        eta = []
        eta_std = []

        eta_random_mean = []
        eta_random_var = []
        eta_random_std = []
        
        for vtype in ['r','s']:
            print(vtype)

            n_gal = []
  
            for rmin,rmax in zip(r1,r2):
                
                beta = ascii.read('../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['beta']

                vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                    .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data

                p50 = np.percentile(vrad,50)
                print(p50)
                
                if vfilter=='hi': beta = beta[vrad>p50]
                elif vfilter=='lo': beta = beta[vrad<p50]

                x = np.log10(beta.data)

                eta_, eta_std_ = get_eta_bs(x)
    
                eta.append( eta_ )
                eta_std.append( eta_std_ )#/np.sqrt(len(bs_eta)) )

                # Obtain mean and var of control samples
                N = len(beta)
                n_gal.append(N)

                eta_random = get_eta_random(1000,N) # 1st parameter will be len(eta_random)
                eta_random_mean.append( np.mean(eta_random) )
                eta_random_var.append( np.var(eta_random,ddof=1) )
                eta_random_std.append( np.std(eta_random,ddof=1))#/np.sqrt(len(eta_random)) )

            print('N =', n_gal)
        nbins = len(r1)
        x = (r1+r2)/2

        ####################
        eta_random_mean = np.array(eta_random_mean)
        eta_random_std = np.array(eta_random_std)

        yran_mean_r = eta_random_mean[:nbins]
        yran_mean_s = eta_random_mean[-nbins:]

        yran_std_r = eta_random_std[:nbins]
        yran_std_s = eta_random_std[-nbins:]
        ####################

        yr = (eta[:nbins]-yran_mean_r)/yran_std_r
        ys = (eta[-nbins:]-yran_mean_s)/yran_std_s

        yr_err = eta_std[:nbins]/yran_std_r
        ys_err = eta_std[-nbins:]/yran_std_s

        ####################

        if sec==3: title = 'H-H Galaxies'
        elif sec==123: title = 'High Spin Galaxies'
        elif sec==0: title = 'All Galaxies'

        for ax, y, yerr, in zip(axs,(yr,ys),(yr_err,ys_err)):

            ax.fill_between(x, -1, 1, alpha=.025, color='k')
            ax.fill_between(x, -1, 3, alpha=.03, color='k')

            ax.hlines(0,x[0],x[-1],linestyles=':')

            if vfilter=='hi': label='High '+r'$V_\mathrm{rad}$'
            if vfilter=='lo': label='Low '+r'$V_\mathrm{rad}$'
            ax.errorbar(x,y,yerr=yerr,capsize=3,fmt='o-',ms=5,label=label)

            ax.set_ylabel(r'$(\eta-\eta_0)/\sigma(\eta)$')


        axs[1].set_xlabel('R/Rv')

        if sec==3: ytext=5.45
        elif sec==123: ytext=5.
        else: ytext=4.15

        axs[0].text(1.25,ytext,'Rising Voids')
        axs[1].text(1.25,ytext,'Shell Voids')

        axs[0].legend(loc='upper left',framealpha=.6)

        axs[0].set_title(title)
        axs[0].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        axs[1].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
        #axs[2].set_xticks([.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

        # x_ticks_labels = []
        # for i in range(nbins):
        #     x_ticks_labels.append( '{:.1f}-{:.1f}'.format(r1[i],r2[i]) )
        #axs[2].set_xticklabels(x_ticks_labels)

    plt.savefig('../plots/vel/etavsr_Vrad_sec{}_v3.jpg'.format(sec))
# %%
"""
Vrad vs R
"""
plt.rcParams['figure.figsize'] = (9, 5)
plt.rcParams['font.size'] = 18
import seaborn as sns

#vtype = 's'
sec = 0
colors = sns.color_palette()

minradV=7.
maxradV=0.
rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

for vtype, c in zip(['a','r','s'],colors[:3]):
    vrad_mean=[]
    vtra_mean=[]
    vrad_err=[]
    vtra_err=[]
    r_mean=[]

    for rmin,rmax in zip(r1,r2):
        vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                        .format(minradV,maxradV,rmin,rmax,sec,vtype))
        vrad_mean.append(np.mean(vel['vrad'].data))
        vtra_mean.append(np.mean(vel['vtra'].data))
        vrad_err.append(np.std(vel['vrad'].data)/(np.sqrt(len(vel['vrad']))))
        vtra_err.append(np.std(vel['vtra'].data)/(np.sqrt(len(vel['vtra']))))
        r_mean.append( (rmin+rmax)/2)

    if vtype=='a': label= r'$\mathrm{V_{rad}}$'+', All voids'
    if vtype=='r': label= r'$\mathrm{V_{rad}}$'+', R-type'
    if vtype=='s': label= r'$\mathrm{V_{rad}}$'+', S-type'

    #plt.plot(r_mean,vrad_mean,marker='o',color=c,label=label)
    #plt.plot(r_mean,vtra_mean,ls='--',marker='x',color=c,label=label)
    plt.errorbar(r_mean,vrad_mean,yerr=vrad_err,marker='o',color=c,label=label)
    plt.errorbar(r_mean,vtra_mean,yerr=vtra_err,ls='--',marker='x',markersize=10,color=c,label=label)

#plt.hlines(115,.8,1.5,linestyles=':',color='k')

plt.xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
plt.xlim([.8,1.5])
plt.legend(ncol=2)
plt.xlabel('r/Rv')
plt.ylabel(r'$\mathrm{V_{rad}, V_{tra}}$ (km/s)')
plt.tight_layout()
plt.savefig('../plots/vel/VradVtra_vs_R_allvtypes.pdf')
    # %%
