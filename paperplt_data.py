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

plt.rcParams['figure.figsize'] = (8, 12)
plt.rcParams['font.size'] = 18

fig , ax = plt.subplots(nrows = 2, ncols = 1, sharex=False, sharey=False, figsize=(10,8))

# #Ngal vs R
# ax[0][0].scatter(voids['r'],np.log10(N))
# ax[0][0].set_xlabel('Void Radius')
# ax[0][0].set_ylabel(r'$\mathrm{log}_{10}N_{gal}$')
# ax[0][0].set_xscale('log')

# #Histogram Void Radius
# ax[0][1].hist(voids['r'],bins=30,cumulative=-1)
# ax[0][1].set_xlabel('Void Radius (Mpc)')
# ax[0][1].set_yscale('log')


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
        ax[0].plot(profile['r'],profile['rho'],c=c,ls=ls,alpha=a,label=label)
    elif nv==2:
        ax[0].plot(profile['r'],profile['rho'],c=c,ls=ls,alpha=a,label=label)
    else:
        ax[0].plot(profile['r'],profile['rho'],c=c,ls=ls,alpha=a)

ax[0].legend(loc = 'lower right')

ax[0].set_ylabel(r'$\Delta(r)$')
ax[0].set_xlabel('r (Mpc/h)')

xv,yv,zv,rv = voids[voids['r']>=8.]['x','y','z','r'][0]

idx1 = tree.query_ball_point([xv,yv,zv],rv*1.5)
idx2 = tree.query_ball_point([xv,yv,zv],rv*0.8)
shell = [g for g in idx1 if g not in idx2]
gals = gxs[shell]

m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gals['mass']),np.log10(gals['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
m1 = -1./m
b1 = -.8*(m-m1)+b
b2 = -.4*(m-m1)+b

from matplotlib.colors import LogNorm
color='k'

#filter = np.where(np.log10(gals['mass'])<=2)
x = np.log10(gals['mass'])
y = np.log10(gals['sp_n'])

#ax[1].scatter(x[filter],y[filter],s=5)
h = ax[1].hist2d(x, y, bins=25, norm=LogNorm(), cmap=plt.cm.Blues, cmin=0)
fig.colorbar(h[3], ax=ax[1], orientation='vertical')

#ax[1][1].scatter(np.log10(gals['mass']),np.log10(gals['sp_n']),edgecolor='blue',facecolor='none')

# ax[1][1].scatter(np.log10(gals['mass']),np.log10(gals['sp_n']),s=5)
# h = ax[1][1].hexbin(np.log10(gals['mass']),np.log10(gals['sp_n']),\
#     bins=100,norm=LogNorm(),mincnt=2,cmap=plt.cm.Blues)
# fig.colorbar(h, ax=ax[1][1], orientation='horizontal')

ax[1].set_xlabel(r'$\mathrm{log_{10}}(M/M_{\odot})$')
ax[1].set_ylabel(r'$\mathrm{log_{10}}|\vec{S}|$')
ax[1].plot(x, x*m+b,ls='-',color=color,label='Linear regression')
ax[1].vlines(-.8, np.min(np.log10(gals['sp_n'])), \
    np.max(y),colors=color,linestyles=':')
ax[1].vlines(-.4, np.min(np.log10(gals['sp_n'])), \
    np.max(y),colors=color,linestyles=':')
ax[1].set_xlim([-1.,2.5])

# ax[1][1].text(-0.95, 2.25, 'H-L', fontsize=15, color=color)
# ax[1][1].text(-0.64, 2.4, 'H-I', fontsize=15, color=color)
# ax[1][1].text(-0.25, 2.5, 'H-H', fontsize=15, color=color)
# ax[1][1].text(-0.95, 0.5, 'L-L', fontsize=15, color=color)
# ax[1][1].text(-0.64, 0.65, 'L-I', fontsize=15, color=color)
# ax[1][1].text(-0.25, 0.8, 'L-H', fontsize=15, color=color)
# ax[1].legend(loc ='lower right')

plt.tight_layout()
#plt.savefig('../plots/data2.png')
#plt.savefig('../plots/data.pdf')
# %%

# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.stats import kde
 
# # create data
# x = np.log10(gals['mass'])
# y = np.log10(gals['sp_n'])
 
# # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
# nbins=100
# k = kde.gaussian_kde([x,y])
# xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
# zi = k(np.vstack([xi.flatten(), yi.flatten()]))
 
# # Make the plot
# plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto')
# plt.colorbar()
# plt.show()
 
# # Change color palette
# # plt.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='auto', cmap=plt.cm.Greens_r)
# # plt.show()
# %%
