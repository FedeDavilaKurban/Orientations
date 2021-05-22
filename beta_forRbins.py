#%%
"""
Calculate Beta (J_perp/J_prll) for different bins in Radius
"""
import sys
import numpy as np
from scipy import spatial
import scipy.stats
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import random
from config import writePath, units
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
#%%
"""
Modified versions of some functions in orientationsTools
"""

def JvsM_(sec,gals,b,m,plot=True):
    """
    Performs linear regression of the Mass-Spin relation
    Determines galaxies of interest with 'sec'
    Returns galaxies of interest in 'gals_h'
    """
    import numpy as np
    import matplotlib.pyplot as plt

    #g_m,g_b,g_rvalue,g_pvalue,g_std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # Este ajuste es de las galaxias en total (ajuste global "g_...")
    
    # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
   #  m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) 
    m1 = -1./m
    b1 = -.7*(m-m1)+b
    b2 = -.3*(m-m1)+b

    M = np.log10(gals['mass'])
    S = np.log10(gals['sp_n'])

    #Alto Spin
    if sec == 1: gals_h = gals[(M < -.5)&(S > M*m+b)]  
    if sec == 3: gals_h = gals[(M > -.5)&(S > M*m+b)]  
    if sec == 13: gals_h = gals[(S > M*m+b)]  

    #Bajo Spin
    if sec == 4: gals_h = gals[(M < -.5)&(S < M*m+b)]  
    if sec == 6: gals_h = gals[(M > -.5)&(S < M*m+b)]   
    if sec == 46: gals_h = gals[(S < M*m+b)] 

    #Solo Masa
    if sec == 14: gals_h = gals[(M < -.5)] 
    if sec == 36: gals_h = gals[(M > -.5)]

    #Todas
    if sec == 0: gals_h = gals


    if plot:
        plt.scatter(M,S,c='gray', label='Galaxies Not Being Analyzed')
        plt.scatter(np.log10(gals_h['mass']),np.log10(gals_h['sp_n']),c='blue', label='Galaxies Being Analyzed')
        plt.plot(M,M*m+b,ls=':', label = 'Void Galaxies Spin-Mass Linear Regression')
        #plt.plot(M,M*g_m+g_b,ls='--', label = 'Box Galaxies Low/High Spin Linear Regression')
        plt.xlabel('log10( M/Msun )')
        plt.ylabel('Spin')
        plt.legend(loc=4)
        plt.show()
        #plt.savefig('../plots/JvsM_void{}.png'.format(nv),dpi=300)

    del M, S
    return gals_h

def orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5,jvsm_b,jvsm_m):
    """
    Determines galaxies of interest in the shell of the void 
    Returns cosine of angles of disks and void centric direction    
    """
    import numpy as np
    from astropy.table import Table

    global vx,vy,vz,vr,cos

    vx,vy,vz,vr = voids['x','y','z','r'][nv]
    if units=='kpc':
        vx*=1000.
        vy*=1000.
        vz*=1000.
        vr*=1000.

    idx1 = tree.query_ball_point([vx,vy,vz],vr*rmax)
    idx2 = tree.query_ball_point([vx,vy,vz],vr*rmin)
    shell = [g for g in idx1 if g not in idx2]
    gals = gxs[shell]
    #print('N of galaxies in shell:',len(gals))

    xv_mean = np.mean(gals['xv'])
    yv_mean = np.mean(gals['yv'])
    zv_mean = np.mean(gals['zv'])
    vshell_mean = [xv_mean, yv_mean, zv_mean]

    gals = JvsM_(sec,gals,jvsm_b,jvsm_m,plot=False) #Determine galaxies of interest (involves rmin,rmax,sec)

    #cos = cosCalc(gals_h,units,tree,xv,yv,zv,rv,s5) #Calculate the cosines of angle of J and void direction

    return gals,vshell_mean 

def readTNG_():
    """
    This function reads subhalos in the TNG300-1 simulation and returns 

    gxs: an ascii Table with the fields and filters I usually need for this: Position, Mass, Spin

    """
    import sys
    from config import basePath, illustrisPath
    sys.path.append(illustrisPath)
    import illustris_python as il
    import numpy as np
    from astropy.table import Table

    mass = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMass'])                                                                                                                      
    ids = np.where((np.log10(mass)>-1.5)&(np.log10(mass)<3.))
    mass = mass[ids]

    pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])[ids]

    spin = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSpin'])[ids]
    sp_n = np.sum(np.abs(spin)**2,axis=-1)**(1./2)
    
    vel = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloVel'])[ids]

    # sfr = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSFR'])[ids]
    
    # masstype = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMassType'])[ids]                                                                                                                      
    # gasmass = masstype[:,0]
    # starmass = masstype[:,4]

    gxs = Table(np.column_stack([pos[:,0],pos[:,1],pos[:,2],mass,spin,sp_n,\
        vel[:,0],vel[:,1],vel[:,2]]),\
        names=['x','y','z','mass','spx','spy','spz','sp_n','xv','yv','zv']) 

    del mass,pos,spin,sp_n,vel

    return gxs
#%%
"""
Reading Galaxies, Creating Tree
"""
gxs = readTNG_()
if units=='Mpc': 
    lbox=205.
    for col in ['x','y','z']:
        gxs[col]/=1000.
if units=='kpc': 
    lbox=205000.
#Applying some filters
gxs.remove_row(np.where(gxs['y']==lbox)[0][0])
gxs.remove_row(np.where(gxs['x']<0.)[0][0])

jvsm_m,jvsm_b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) 


tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

#%%
"""
Define Parameters 
"""

#exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(sys.argv[1])
minradV, maxradV, s5 = 7., 0., 0
#print('Codename of experiment:', exp)
print('minradVoid = {}Mpc'.format(minradV))
print('maxradVoid = {}Mpc'.format(maxradV))
#print('rmin = {}Rvoid'.format(rmin))
#print('rmax = {}Rvoid'.format(rmax))
#print('sec =',sec)
print('s5 =',s5)
#print('vtype =',vtype)

#voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)
#print('Num of voids:',len(voids))

# %%
"""
Calculating
"""

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
        eta_ran [i] = len(b[b>1]) / len(b[b<1])

    return eta_ran
    
def get_beta():
    nvs = range(len(voids))
    prll = []
    perp = []
    ratio = []
    r = []
    rrv = []

    for nv in nvs:

        vx, vy, vz, vr = voids[nv]['x'], voids[nv]['y'], voids[nv]['z'], voids[nv]['r']
        gals, vshell = orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5,jvsm_b,jvsm_m)

        for g in gals:

            gx, gy, gz = g['x'], g['y'], g['z']

            # A la velocidad de la galaxia le resto la "velocidad del shell/void" 
            #gxv, gyv, gzv = g['xv']-vshell[0], g['yv']-vshell[1], g['zv']-vshell[2]

            sx, sy, sz, sn = g['spx'], g['spy'], g['spz'], g['sp_n']

            s = [sx,sy,sz] #spin vector

            gr = [gx-vx,gy-vy,gz-vz] #galaxy position from void center
            for axis in range(len(gr)):
                if gr[axis]<=-(rmax+.5)*vr: gr[axis]+=lbox
                if gr[axis]>= (rmax+.5)*vr: gr[axis]-=lbox

            r.append( np.linalg.norm(gr) )
            rrv.append(r[-1]/vr)
            r_versor = gr/r[-1] #radial direction from void center

            prll.append( abs( np.dot(r_versor,s) ) )
            perp.append( np.sqrt(sn**2 - prll[-1]**2) ) # S**2 = S_perp**2 + S_paral**2


    r = np.array(r)
    rrv = np.array(rrv)

    perp =  np.array(perp)
    prll = np.array(prll)
    beta = perp/prll

    return beta,perp,prll,r,rrv
########################################################

plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

r1 = [0.4,0.5,0.6,0.7,0.8,0.9,1.,1.1,1.2]
r2 = [0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3]

for sec in [0]:

    eta = []
    eta_random_mean = []
    eta_random_var = []

    for vtype in ['a','r','s']:
        print(vtype)

        for rmin,rmax in zip(r1,r2):

            voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)

            print(rmin,rmax)
            beta, perp, prll, r, rrv = get_beta()

            x = np.log10(beta)
            n_perp = len(np.where(x>0.)[0])
            n_prll = len(np.where(x<0.)[0])
            eta.append( n_perp / n_prll )

            N = len(beta)
            eta_random = get_eta_random(500,N)
            eta_random_mean.append( np.mean(eta_random) )
            eta_random_var.append( np.var(eta_random,ddof=1) )


            plt.hist(x,bins=30,density=True)
            plt.vlines(0,0,1.,linestyles=':')

            plt.xlabel(r'$log_{10}\beta$')

            plt.text(-2.5,.6,r'$n_2(\beta<1)={}$'.format(n_prll))
            plt.text(1.,.6,r'$n_1(\beta>1)={}$'.format(n_perp))
            plt.text(3.,.2,r'$n_1/n_2={:.3f}$'.format(eta[-1]))

            plt.xlim([-3.,5.])
            plt.ylim([0.,1.])

            #plt.show()
            #plt.savefig('../plots/beta_forRbins/voids_{}/rmin{}_rmax{}.jpg'.\
            #            format(vtype, rmin, rmax))
            plt.close()

# %%

nbins = len(r1)
x = r1
ya = eta[:nbins]
yr = eta[nbins:2*nbins]
ys = eta[-nbins:]

# Theoretical n1/n2 value for random spins
plt.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
y = np.array(eta_random_mean[-nbins:])
plt.plot(x,y,c='k')
y1 = y - np.array(eta_random_var[-nbins:])
y2 = y + np.array(eta_random_var[-nbins:])
plt.fill_between(x, y1, y2, alpha=.5, color='k')

plt.plot(x,ya,label='All')
plt.plot(x,yr,label='R')
plt.plot(x,ys,label='S')


plt.ylabel(r'$n_1/n_2$')
plt.xlabel('R/Rv')
plt.legend()

x_ticks_labels = []
for i in range(nbins):
    x_ticks_labels.append( '{}-{}'.format(r1[i],r2[i]) )
plt.xticks(x, x_ticks_labels)

#plt.savefig('../plots/beta_forRbins/factorforallvoids.jpg')
#%%
nbins = len(r1)
x = r1

ya = eta[:nbins]
yr = eta[nbins:2*nbins]
ys = eta[-nbins:]

eta_random_mean = np.array(eta_random_mean)
eta_random_var = np.array(eta_random_var)

yran_mean_a = eta_random_mean[:nbins]
yran_mean_r = eta_random_mean[nbins:2*nbins]
yran_mean_s = eta_random_mean[-nbins:]

yran_var_a = eta_random_var[:nbins]
yran_var_r = eta_random_var[nbins:2*nbins]
yran_var_s = eta_random_var[-nbins:]

fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=False)

for ax, y, yran, yvar, label, in zip(axs,\
                                (ya,yr,ys),\
                                (yran_mean_a,yran_mean_r,yran_mean_s),\
                                (yran_var_a,yran_var_r,yran_var_s),('All','Rising','Shell')):

    # Theoretical n1/n2 value for random spins
    ax.hlines(1/(np.sqrt(2)-1),x[0],x[-1],linestyles=':')
    ax.plot(x,yran,c='k',alpha=.7)

    ax.fill_between(x, yran-yvar, yran+yvar, alpha=.3, color='k')
    ax.fill_between(x, yran-2*yvar, yran+2*yvar, alpha=.3, color='k')
    ax.fill_between(x, yran-3*yvar, yran+3*yvar, alpha=.3, color='k')

    ax.plot(x,y,label=label)

    ax.legend()

    # plt.ylabel(r'$n_1/n_2$')
    # plt.xlabel('R/Rv')
    # plt.legend()

x_ticks_labels = []
for i in range(nbins):
    x_ticks_labels.append( '{}-{}'.format(r1[i],r2[i]) )
axs[2].set_xticks(x)
axs[2].set_xticklabels(x_ticks_labels)

plt.savefig('../plots/beta_forRbins/factorforallvoids_err.jpg')
##########################################################
# %%
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

b = np.logspace(-4,4,100)
b_dist = np.sin(np.arctan(b))/(1+b**2)
b_dist2 = np.tan(np.arccos(np.random.random(10000)))
plt.plot(np.log10(b),b_dist)
#plt.hist(np.log10(b_dist2),density=True,bins=100)
#plt.vlines(np.log10(1/(np.sqrt(2)-1)),0,1,linestyles=':')
# %%

u=np.random.random(10000)

b_dist = np.log10(np.tan(np.arccos(u)))

plt.hist(b_dist,bins=100,density=True)
plt.vlines(np.log10(1/(np.sqrt(2)-1)),0,1,linestyles=':')
# %%
x = np.random.uniform(-10,10,1000)
y = np.cos(np.arctan(x))
plt.scatter(x,y)

x = np.random.uniform(-10,10,1000)
y = np.sin(np.arctan(x))
plt.scatter(x,y)

# %%

x = np.random.uniform(0,np.pi,1000)
y = np.abs(np.cos(x))
plt.scatter(x,y)

y = np.sin(x)
plt.scatter(x,y)
#%%

x = np.random.uniform(0,np.pi,1000)

y = np.tan(x)
plt.scatter(x,y)

y = np.sin(x)/np.abs(np.cos(x))
plt.scatter(x,y)

y = np.abs(np.tan(x))
plt.scatter(x,y)
# %%
