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
voids = voids[voids['r']>=5.]

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
fig , ax = plt.subplots(nrows = 2, ncols = 2, sharex=False, sharey=False, figsize=(10,8))

#Ngal vs R
ax[0][0].scatter(voids['r'],np.log10(N))
ax[0][0].set_xlabel('Void Radius')
ax[0][0].set_ylabel(r'$\mathrm{log}_{10}N_{gal}$')
ax[0][0].set_xscale('log')

#Histogram Void Radius
ax[0][1].hist(voids['r'],bins=30,cumulative=-1)
ax[0][1].set_xlabel('Void Radius (Mpc)')
ax[0][1].set_yscale('log')


minradV = 7.
for nv in range(len(voids[voids['r']>=minradV])):
    filename = writePath+'Proyectos/Orientations/data/profiles/minradV{}_void{}.dat'.format(minradV,nv)
    profile = ascii.read(filename,names=['r','rho'])
    ax[1][0].plot(profile['r'],profile['rho'])


ax[1][0].set_ylabel(r'$\rho/\bar{\rho}-1$')
ax[1][0].set_xlabel('r (Mpc)')

xv,yv,zv,rv = voids[voids['r']>=7.]['x','y','z','r'][0]

idx1 = tree.query_ball_point([xv,yv,zv],rv*1.2)
idx2 = tree.query_ball_point([xv,yv,zv],rv*0.9)
shell = [g for g in idx1 if g not in idx2]
gals = gxs[shell]

m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gals['mass']),np.log10(gals['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
m1 = -1./m
b1 = -.7*(m-m1)+b
b2 = -.3*(m-m1)+b


ax[1][1].scatter(np.log10(gals['mass']),np.log10(gals['sp_n']))
ax[1][1].set_xlabel(r'$\mathrm{log_{10}}M/M_{\odot}$')
ax[1][1].set_ylabel(r'$\mathrm{log_{10}}|\vec{J}|$')
ax[1][1].plot(np.log10(gals['mass']),np.log10(gals['mass'])*m+b,ls='-',color='orange',label='Linear regression')
ax[1][1].vlines(-.7, np.min(np.log10(gals['sp_n'])), np.max(np.log10(gals['sp_n'])),colors='orange',linestyles=':')
ax[1][1].vlines(-.3, np.min(np.log10(gals['sp_n'])), np.max(np.log10(gals['sp_n'])),colors='orange',linestyles=':')

ax[1][1].text(-0.99, 1.6, 'H-L', fontsize=15, color='orange')
ax[1][1].text(-0.64, 1.75, 'H-I', fontsize=15, color='orange')
ax[1][1].text(-0.25, 1.9, 'H-H', fontsize=15, color='orange')
ax[1][1].text(-0.99, 1.05, 'L-L', fontsize=15, color='orange')
ax[1][1].text(-0.64, 1.25, 'L-I', fontsize=15, color='orange')
ax[1][1].text(-0.25, 1.45, 'L-H', fontsize=15, color='orange')
ax[1][1].legend(loc ='lower right')

plt.tight_layout()
plt.savefig('../plots/data.png')

# %%
