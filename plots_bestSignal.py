#%%
# import sys
# import numpy as np
# from scipy import spatial
# import scipy.stats
# from astropy.table import Table
# from astropy.io import ascii
# from orientationsTools import *
# import random
# from config import writePath, units
# import matplotlib.pyplot as plt

#%%
import sys
from config import illustrisPath
sys.path.append(illustrisPath)
import illustris_python as il
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import numpy as np
from scipy import spatial
import seaborn
import scipy.stats
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from orientationsTools import *
from config import writePath, units

# %%
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
    ids = np.where((np.log10(mass)>-1))#&(np.log10(mass)<3.))
    mass = mass[ids]

    pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])[ids]

    gxs = Table(np.column_stack([pos[:,0],pos[:,1],pos[:,2]\
                                ,mass\
                                #,spin,sp_n\
                                #,vel[:,0],vel[:,1],vel[:,2]\
                                #,gasmass,starmass,\
                                ]),
                names=['x','y','z','mass'])#,'spx','spy','spz','sp_n'])#,'xv','yv','zv'])#,'gasmass','starmass']) 

    del mass,pos

    return gxs

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

for col in ['x','y','z']:
    gxs[col]/=1000.
lbox /= 1000.

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

# %%
"""
Plot High Spin, Low Velocity, High/Low Mass
"""
import seaborn as sns
#colors = sns.color_palette()
#blue, oran = colors[0], colors[1]

plt.rcParams['figure.figsize'] = (9, 8)
plt.rcParams['font.size'] = 20
fig, axs = plt.subplots(2, 1, constrained_layout=True,\
    sharex=True, sharey=False, gridspec_kw={'height_ratios': [3, 1]})


minradV = 7.
maxradV = 0.
#axs[0].xscale


r1=[.8,.9,1.,1.1,1.2,1.3,1.4]
r2=[.9,1.,1.1,1.2,1.3,1.4,1.5]
x = (np.array(r1)+np.array(r2))/2

axs[0].fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
axs[0].fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')
axs[0].hlines(0,.8,1.5,linewidth=.6,color='k',alpha=.7)


sec=3
for vtype, label, fmt, c in zip(['a','r','s'],\
                            ['All Voids', 'R-Type', 'S-Type'],\
                            ['x:','o-','^-.'],\
                            ['k','C03','C00']):

    filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])

    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    

    axs[0].errorbar(x,y,yerr=yerr,label=label,fmt=fmt,\
                    capsize=3,ms=7,color=c)

axs[0].legend(loc='upper left')

axs[0].set_ylabel(r'$\zeta$')
#axs[0].set_xlabel(r'$\mathrm{r/R_v}$')
axs[0].set_xlim([.8,1.5])
axs[0].set_ylim([-4,8])


########################################

r_array = np.linspace(8,16.,30)

Ntotal = len(gxs)
rhototal = len(gxs)/(lbox**3)

r_array = np.linspace(6,10.6,30) #Rv
r_array = np.linspace(.8,1.5,30) #Rv

for vtype, c in zip(['s','r'],['C00','C03']):
    print('vtype=',vtype)

    voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)

    rho_arrays= []

    for nv in range(len(voids)):
        #print(nv)
        x,y,z,radius = voids['x','y','z','r'][nv].as_void()
        
        rho_array = []
        for r in r_array*radius:
            vol = 4.*np.pi*(r)**3./3
            rho_array.append( len(tree.query_ball_point([x,y,z],r))/vol)

        rho_arrays.append(rho_array)

    rho_mean = np.mean(rho_arrays,axis=0)
    rho_meanerr = np.std(rho_arrays,axis=0)/np.sqrt(len(voids))
    rho_std = np.std(rho_arrays,axis=0)

    if vtype=='r': label='R-type'
    elif vtype=='s': label='S-type'

    y = rho_mean/rhototal-1

    axs[1].fill_between(r_array,y-rho_std/rhototal,y+rho_std/rhototal,\
        color=c,alpha=.5)
    axs[1].errorbar(r_array,rho_mean/rhototal-1,yerr=rho_meanerr/rhototal,\
        marker='o',ms=5,lw=1, color=c,label=label)
    axs[1].set_ylabel(r'$\Delta(r/Rv)$')

axs[1].legend()
axs[1].set_yticks([-1.,-.75,-.5,-.25,0.])
axs[1].set_xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

axs[1].set_xlim([.8,1.5])
axs[1].set_ylim([-1.,0])

axs[1].set_xlabel(r'$\mathrm{r/R_v}$')

axs[0].text(0.825,-3.75,'High Mass, High Spin, Low {} galaxies'.format(r'$\mathrm{v_{rad}}$'))

#plt.tight_layout()
#plt.show()
plt.savefig('../plots/bestSignalwithProfiles.pdf')
# %%
"""
With Vrad
"""
import seaborn as sns
#colors = sns.color_palette()
#blue, oran = colors[0], colors[1]

plt.rcParams['figure.figsize'] = (9, 10)
plt.rcParams['font.size'] = 20
fig, axs = plt.subplots(3, 1, constrained_layout=True,\
    sharex=True, sharey=False, gridspec_kw={'height_ratios': [3, 1, 1]})


minradV = 7.
maxradV = 0.
#axs[0].xscale


r1=[.8,.9,1.,1.1,1.2,1.3,1.4]
r2=[.9,1.,1.1,1.2,1.3,1.4,1.5]
x = (np.array(r1)+np.array(r2))/2

axs[0].fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -1, 1, alpha=.1, color='k')
axs[0].fill_between([.8,.9,1.,1.1,1.2,1.3,1.4,1.5], -3, 3, alpha=.1, color='k')
axs[0].hlines(0,.8,1.5,linewidth=.6,color='k',alpha=.7)


sec=3
for vtype, label, fmt, c in zip(['a','r','s'],\
                            ['All Voids', 'R-Type', 'S-Type'],\
                            ['x:','o-','^-.'],\
                            ['k','C03','C00']):

    filename = '../data/eta/eta_vfilterlo_minradV{}_maxradV{}_sec{}_vtype{}.txt'\
            .format(minradV,maxradV,sec,vtype)
    etaTable = ascii.read(filename,\
            names=['eta','eta_std','rmin','rmax','N'])

    yran_mean = 1./(np.sqrt(2)-1)
    yran_err = np.sqrt(28.1421/etaTable['N'])

    y = (etaTable['eta'].data-yran_mean)/yran_err
    yerr = etaTable['eta_std'].data/yran_err

    

    axs[0].errorbar(x,y,yerr=yerr,label=label,fmt=fmt,\
                    capsize=3,ms=7,color=c)

axs[0].legend(loc='upper left')

axs[0].set_ylabel(r'$\zeta$')
#axs[0].set_xlabel(r'$\mathrm{r/R_v}$')
axs[0].set_xlim([.8,1.5])
axs[0].set_ylim([-4,8])


########################################

r_array = np.linspace(8,16.,30)

Ntotal = len(gxs)
rhototal = len(gxs)/(lbox**3)

r_array = np.linspace(6,10.6,30) #Rv
r_array = np.linspace(.8,1.5,30) #Rv

for vtype, c in zip(['s','r'],['C00','C03']):
    print('vtype=',vtype)

    voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)

    rho_arrays= []

    for nv in range(len(voids)):
        #print(nv)
        x,y,z,radius = voids['x','y','z','r'][nv].as_void()
        
        rho_array = []
        for r in r_array*radius:
            vol = 4.*np.pi*(r)**3./3
            rho_array.append( len(tree.query_ball_point([x,y,z],r))/vol)

        rho_arrays.append(rho_array)

    rho_mean = np.mean(rho_arrays,axis=0)
    rho_meanerr = np.std(rho_arrays,axis=0)/np.sqrt(len(voids))
    rho_std = np.std(rho_arrays,axis=0)

    if vtype=='r': label='R-type'
    elif vtype=='s': label='S-type'

    y = rho_mean/rhototal-1

    axs[1].fill_between(r_array,y-rho_std/rhototal,y+rho_std/rhototal,\
        color=c,alpha=.5)
    axs[1].errorbar(r_array,rho_mean/rhototal-1,yerr=rho_meanerr/rhototal,\
        marker='o',ms=5,lw=1, color=c,label=label)

axs[1].set_ylabel(r'$\Delta(r/Rv)$')

axs[1].legend()
axs[1].set_yticks([-1.,-.75,-.5,-.25,0.])
axs[1].set_xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

axs[1].set_xlim([.8,1.5])
axs[1].set_ylim([-1.,0])


axs[0].text(0.825,-3.75,'High Mass, High Spin, Low {} galaxies'.format(r'$\mathrm{v_{rad}}$'))

#####################################################

#vtype = 's'
sec = 3

minradV=7.
maxradV=0.
rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

for vtype, c in zip(['r','s'],['C00','C03']):
    vrad_mean=[]
    vtra_mean=[]
    vrad_err=[]
    vtra_err=[]
    vrad_std=[]
    r_mean=[]

    for rmin,rmax in zip(r1,r2):
        vel = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                        .format(minradV,maxradV,rmin,rmax,sec,vtype))

        vel = vel[vel['vrad']<np.median(vel['vrad'])]
        vrad_mean.append(np.mean(vel['vrad'].data))
        vrad_err.append(np.std(vel['vrad'].data)/(np.sqrt(len(vel['vrad']))))
        r_mean.append( (rmin+rmax)/2)
        vrad_std.append(np.std(vel['vrad']))

    if vtype=='r': label= 'R-type'
    if vtype=='s': label= 'S-type'

    vrad_mean=np.array(vrad_mean)
    vrad_std=np.array(vrad_std)

    axs[2].fill_between(r_mean,vrad_mean-vrad_std,vrad_mean+vrad_std,\
        alpha=.5,color=c)
    axs[2].errorbar(r_mean,vrad_mean,yerr=vrad_err,marker='o',color=c,\
        label=label)

axs[2].hlines(0,.8,1.5,linestyles=':',color='k')
#plt.hlines(15,.8,1.5,linestyles=':',color='k')

axs[2].set_xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
axs[2].set_yticks([-400,-200,0,200])
axs[2].set_xlim([.8,1.5])
#axs[2].legend(ncol=1)
axs[2].set_xlabel(r'$\mathrm{r/R_v}$')
axs[2].set_ylabel(r'$\mathrm{V_{rad}}$ (km/s)')
#plt.show()
plt.savefig('../plots/bestSignalwithProfilesVel.pdf')
# %%
