#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import spatial
import scipy.stats
from astropy.io import ascii
from astropy.table import Table
from orientationsTools import readTNG, readVoids, readExp
#import random
from config import writePath, units
#%%
voidsfilename = '../data/tng300-1_voids.dat'
names = ['r','x','y','z','vx','vy','vz','deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter']
voids = ascii.read(voidsfilename,names=names)
voids = voids[voids['r']>=7.]

gxs = readTNG()
lbox=205.
for col in ['x','y','z']:
    gxs[col]/=1000.
gxs.remove_row(np.where(gxs['y']==lbox)[0][0])
gxs.remove_row(np.where(gxs['x']<0.)[0][0])
tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

N=[]
for i in range(len(voids)):
	if i%1000==0: print(i)
	x,y,z,r = voids['x','y','z','r'][i]
	N.append(len(tree.query_ball_point([x,y,z],r)))

#%%
plt.rcParams['figure.figsize'] = (8, 8)
plt.rcParams['font.size'] = 15

fig , ax = plt.subplots(nrows = 1, ncols = 1, sharex=False, sharey=False, figsize=(10,8))

minradV = 7.
filename = writePath+'Proyectos/Orientations/data/vtype/minradV{}_vtype.dat'.format(float(minradV))
vtypes = ascii.read(filename)
for nv in range(len(voids[voids['r']>=minradV])):
    filename = writePath+'Proyectos/Orientations/data/profiles/minradV{}_void{}.dat'.format(minradV,nv)
    profile = ascii.read(filename,names=['r','rho'])
    if nv==38: profile = profile[1:]

    if vtypes[nv]['type']=='r': c, ls, a, label = 'C0', '--', 0.4, 'R-voids'
    if vtypes[nv]['type']=='s': c, ls, a, label = 'C1', '-', 0.5, 'S-voids'

    if nv==1:
        ax.plot(profile['r'],profile['rho'],c=c,ls=ls,alpha=a,label=label)
    elif nv==2:
        ax.plot(profile['r'],profile['rho'],c=c,ls=ls,alpha=a,label=label)
    else:
        ax.plot(profile['r'],profile['rho'],c=c,ls=ls,alpha=a)

ax.legend(loc = 'lower right')


ax.set_ylabel(r'$\rho/\bar{\rho}-1$')
ax.set_xlabel('r (Mpc h{})'.format(r'$^{-1}$'))

#%%

plt.rcParams['figure.figsize'] = (8, 8)
plt.rcParams['font.size'] = 15

fig , ax = plt.subplots(nrows = 3, ncols = 1, sharex=False, sharey=False, figsize=(10,8))
rmin, rmax = 1.0,1.1

#xv,yv,zv,rv = voids[voids['r']>=8.]['x','y','z','r'][0]
shell=[]
for nv in range(len(voids)):
    xv,yv,zv,rv = voids['x','y','z','r'][nv]
    idx1 = tree.query_ball_point([xv,yv,zv],rv*rmax)
    idx2 = tree.query_ball_point([xv,yv,zv],rv*rmin)
    shell.append( [g for g in idx1 if g not in idx2] )
shell = np.concatenate(shell)
gals = gxs[shell]

m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gals['mass']),np.log10(gals['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
m1 = -1./m
b1 = -.7*(m-m1)+b
b2 = -.3*(m-m1)+b

from matplotlib.colors import LogNorm
color='k'

filter = np.where(np.log10(gals['mass'])<=3)
x = np.log10(gals['mass'])
y = np.log10(gals['sp_n'])

#ax[1].scatter(x[filter],y[filter],s=5)
h = ax[0].hist2d(x[filter], y[filter], bins=30, norm=LogNorm(), cmap=plt.cm.Blues, cmin=1)
#h = ax[1].hist2d(x, y, bins=60, norm=LogNorm(), cmap=plt.cm.Blues, cmin=3)
fig.colorbar(h[3], ax=ax[0], orientation='vertical')


ax[0].set_xlabel(r'$\mathrm{log_{10}}M/M_{\odot}$')
ax[0].set_ylabel(r'$\mathrm{log_{10}}|\vec{S}|$')
ax[0].plot(x, x*m+b,ls='-',color=color,label='Linear regression')
ax[0].vlines(-.7, np.min(np.log10(gals['sp_n'])), \
    np.max(y),colors=color,linestyles=':')
ax[0].vlines(-.3, np.min(np.log10(gals['sp_n'])), \
    np.max(y),colors=color,linestyles=':')

s5 = ascii.read('../data/sigma5/sigma5_minradV7.0_maxradV0.0_rmin{}_rmax{}_sec0_vtypea.txt'.format(rmin,rmax))['sigma5']
ax[1].hist(np.log10(s5),bins=25,density=True)
ax[1].set_xlabel(r'$\mathrm{log_{10}}\Sigma_5$')

v = ascii.read('../data/vel/-1/vel_minradV7.0_maxradV0.0_rmin{}_rmax{}_sec0_vtypea.txt'.format(rmin,rmax))['vrad']
ax[2].hist(v,bins=25,density=True)
ax[2].set_xlabel(r'$\mathrm{log_{10}}v_{rad}$')

plt.tight_layout()
#plt.savefig('../plots/data.png')



# %%
