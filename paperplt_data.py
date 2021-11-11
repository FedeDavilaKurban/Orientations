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

plt.rcParams['figure.figsize'] = (6, 24)
plt.rcParams['font.size'] = 18

fig , ax = plt.subplots(nrows = 3, ncols = 1, sharex=False, sharey=False, figsize=(10,8))

# #Ngal vs R
# ax[0][0].scatter(voids['r'],np.log10(N))
# ax[0][0].set_xlabel('Void Radius')
# ax[0][0].set_ylabel(r'$\mathrm{log}_{10}N_{gal}$')
# ax[0][0].set_xscale('log')

# #Histogram Void Radius
# ax[0][1].hist(voids['r'],bins=30,cumulative=-1)
# ax[0][1].set_xlabel('Void Radius (Mpc)')
# ax[0][1].set_yscale('log')

##################################################################
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

##################################################################


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

##################################################################

sec = 0
import seaborn as sns
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
    ax[2].errorbar(r_mean,vrad_mean,yerr=vrad_err,marker='o',color=c,label=label)
    ax[2].errorbar(r_mean,vtra_mean,yerr=vtra_err,ls='--',marker='x',markersize=10,color=c,label=label)

ax[2].hlines(115,.8,1.5,linestyles=':',color='k')
ax[2].set_xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
ax[2].set_xlim([.8,1.5])
ax[2].legend(ncol=3)
ax[2].set_xlabel('r/Rv')
ax[2].set_ylabel(r'$\mathrm{V_{rad}, V_{tra}}$ (km/s)')

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
