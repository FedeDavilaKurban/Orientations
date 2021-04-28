#%%
"""
Calculate ratio of perpendicular and parallel J components
"""
import sys
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *
import random
from config import writePath, units
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
#%%
"""
New definitions based on orientationsTools
"""
def orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5):
    """
    Determines galaxies of interest in the shell of the void 
    Returns cosine of angles of disks and void centric direction    
    """
    import numpy as np
    from astropy.table import Table

    global xv,yv,zv,rv,cos

    xv,yv,zv,rv = voids['x','y','z','r'][nv]
    if units=='kpc':
        xv*=1000.
        yv*=1000.
        zv*=1000.
        rv*=1000.

    idx1 = tree.query_ball_point([xv,yv,zv],rv*rmax)
    idx2 = tree.query_ball_point([xv,yv,zv],rv*rmin)
    shell = [g for g in idx1 if g not in idx2]
    gals = gxs[shell]
    #print('N of galaxies in shell:',len(gals))

    gals = JvsM(sec,gals,gxs,plot=False) #Determine galaxies of interest (involves rmin,rmax,sec)

    #cos = cosCalc(gals_h,units,tree,xv,yv,zv,rv,s5) #Calculate the cosines of angle of J and void direction

    return gals 

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
    ids = np.where((np.log10(mass)>-1.)&(np.log10(mass)<3.))
    mass = mass[ids]

    pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])[ids]

    spin = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSpin'])[ids]
    sp_n = np.sum(np.abs(spin)**2,axis=-1)**(1./2)
    
    sfr = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSFR'])[ids]
    
    masstype = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMassType'])[ids]                                                                                                                      
    gasmass = masstype[:,0]
    starmass = masstype[:,4]

    gxs = Table(np.column_stack([pos[:,0],pos[:,1],pos[:,2],mass,spin,sp_n,sfr,gasmass,starmass]),\
        names=['x','y','z','mass','spx','spy','spz','sp_n','sfr','gasmass','starmass']) 

    del mass,pos,spin,sp_n,sfr,masstype

    return gxs

#%%
"""
Define Parameters, Reading Galaxies, Creating Tree
"""

#exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(sys.argv[1])
minradV, maxradV, rmin, rmax, sec, s5, vtype = 8., 9., .9, 1.2, 3, 0, 'a'
#print('Codename of experiment:', exp)
print('minradVoid = {}Mpc'.format(minradV))
print('maxradVoid = {}Mpc'.format(maxradV))
print('rmin = {}Rvoid'.format(rmin))
print('rmax = {}Rvoid'.format(rmax))
print('sec =',sec)
print('s5 =',s5)
print('vtype =',vtype)

voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)
print('Num of voids:',len(voids))

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

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

#%%
# g = gals[-1]
# v = voids[0]
# pg = g['x','y','z']
# pv = v['x','y','z']
# s = g['spx','spy','spz','sp_n']
# gx, gy, gz = pg[0], pg[1], pg[2]
# vx, vy, vz = pv[0]*1000, pv[1]*1000, pv[2]*1000
# sx, sy, sz = s[0], s[1], s[2]
# sn = s['sp_n']

# color = plt.rcParams['axes.prop_cycle'].by_key()['color']

# %matplotlib
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.quiver(0,0,0, vx, vy, vz, length=1, arrow_length_ratio=.02,color=color[0],\
#     label='Void Center')
# ax.quiver(gx, gy, gz, sx, sy, sz, length=1, arrow_length_ratio=.2,color=color[1],\
#     label='Spin')
# ax.quiver(0,0,0, gx, gy, gz,length=1, arrow_length_ratio=.02,color=color[2],\
#     label='Galaxy')
# plt.legend()
# plt.show()
# # %%
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# #ax.quiver(0,0,0, vx, vy, vz, length=1, arrow_length_ratio=.02,color=color[0],\
# #    label='Void Center')
# ax.quiver(gx, gy, gz, sx, sy, sz, length=1, arrow_length_ratio=.02,color=color[1],\
#     label='Spin')
# ax.quiver(vx, vy, vz, gx-vx, gy-vy, gz-vz,length=1, arrow_length_ratio=.02,color=color[2],\
#     label='Galaxy from Void Center')
# plt.legend()
# plt.show()
# %%
"""
Read Galaxy Properties (Spin, SFR, Masses) in the shells of Voids
"""

#nvs = [0]
nvs = range(len(voids))
prll = []
perp = []
ratio = []
sfr_ = []
gmass_ = []
smass_ = []

for nv in nvs:

    vx, vy, vz = voids[nv]['x'], voids[nv]['y'], voids[nv]['z']
    gals = orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5)

    sfr_.append(gals['sfr'].data)
    gmass_.append(gals['gasmass'].data)
    smass_.append(gals['starmass'].data)

    for g in gals:

        gx, gy, gz = g['x'], g['y'], g['z']

        sx, sy, sz, sn = g['spx'], g['spy'], g['spz'], g['sp_n']

        u = [gx-vx,gy-vy,gz-vz]
        u_versor = u/np.linalg.norm(u)
        v = [sx,sy,sz]
        prll.append( abs( np.dot(u_versor,v) ) )
        perp.append( np.sqrt(sn**2 - prll[-1]**2) )

ratio = np.array(perp)/np.array(prll)
sfr = np.concatenate(sfr_)
gmass = np.concatenate(gmass_)
smass = np.concatenate(smass_)

del sfr_,gmass_,smass_
# %%
"""
Histogramas de las Componentes Perpendicular y Paralela
"""
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

plt.hist(np.log10(perp),label='Perpendicular',bins=50,alpha=.8)
plt.hist(np.log10(prll),label='Parallel',bins=50,alpha=.8)

plt.xlabel(r'$\mathrm{log}_{10}\mathrm{J}_\perp,\,\mathrm{log}_{10}\mathrm{J}_\parallel$')

plt.legend()
plt.savefig('../plots/components_dist.png')
plt.show()

#%%
"""
Histogramas de los ratios de las componentes
"""
plt.hist(np.log10(ratio),bins=30,log=False)
plt.xlabel(r'$log_{10}(\frac{J_\perp}{/J_\parallel})$')
#plt.legend()
#plt.savefig('../plots/components_dist.png')
plt.show()


#%%
# """
# Scatter SFR vs Componentes Spin
# """
# plt.scatter(np.log10(prll),sfr,label='Parallel')
# plt.scatter(np.log10(perp),sfr,label='Perpendicular')
# plt.xlabel(r'$\mathrm{log}_{10}\mathrm{J}_\perp,\,\mathrm{log}_{10}\mathrm{J}_\parallel$')
# plt.ylabel('SFR')
# plt.legend()
#%%
"""
Scatter + Hist SFR vs Ratio Componentes Spin
"""
a = .8
fig = plt.figure(constrained_layout=True)
widths = [1]
heights = [.5, .5]
spec = fig.add_gridspec(ncols=1, nrows=2, width_ratios=widths,
                          height_ratios=heights)
ax00 = fig.add_subplot(spec[0,0])
ax10 = fig.add_subplot(spec[1,0])

y1 = np.log10(ratio[ratio>=1.]) #Perpendicular
y2 = np.log10(ratio[ratio<=1.]) #Paralelo

x1 = sfr[ratio>=1.]
x2 = sfr[ratio<=1.]


ax10.scatter(x1,y1,label='Perpendicular',alpha=a)
ax10.scatter(x2,y2,label='Parallel',alpha=a)
ax10.set_ylabel(r'$log_{10}(\frac{J_\perp}{/J_\parallel})$')
ax10.set_xlabel('SFR')

#plt.vlines(np.log10(1.),np.min(y),np.max(y))
ax10.legend()

ax00.hist(x1, bins=30, log=True, alpha=a, label='Perpendicular',density=True)
ax00.hist(x2, bins=20, log=True, alpha=a, label='Parallel',density=True)
ax00.legend()
plt.savefig('plot1.jpg')
#%%
# """
# Scatter+Hist SFR y Componentes J
# """
# fig5 = plt.figure(constrained_layout=True)
# widths = [2, .5]
# heights = [.5, 2]
# spec5 = fig5.add_gridspec(ncols=2, nrows=2, width_ratios=widths,
#                           height_ratios=heights)

# ax00 = fig5.add_subplot(spec5[0,0])
# ax10 = fig5.add_subplot(spec5[1,0])
# ax11 = fig5.add_subplot(spec5[1,1])

# a = .8
# y = sfr
# x1 = np.log10(prll)
# x2 = np.log10(perp)
# ax10.scatter(x1,y,label='Parallel',alpha=a)
# ax10.scatter(x2,y,label='Perpendicular',alpha=a)

# ax00.hist(x1,bins=30,alpha=a)
# ax00.hist(x2,bins=30,alpha=a)

# ax11.hist(y,bins=30,orientation='horizontal',alpha=a,log=True)
# #ax11.hist(y,bins=30,orientation='horizontal',alpha=a,log=True)

# ax10.legend()

# ax00.tick_params(axis="x", labelbottom=False)
# ax11.tick_params(axis="y", labelleft=False)
#%%
"""
Scatter + Hist Gmass/SMass vs Ratio Componentes Spin
"""
a = .8
fig = plt.figure(constrained_layout=True)
widths = [2,2]
heights = [.5, .5]
spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=widths,
                          height_ratios=heights)
ax00 = fig.add_subplot(spec[0,0])
ax10 = fig.add_subplot(spec[1,0])
ax01 = fig.add_subplot(spec[0,1])
ax11 = fig.add_subplot(spec[1,1])


y1 = np.log10(ratio[ratio>=1.]) #Perpendicular
y2 = np.log10(ratio[ratio<=1.]) #Paralelo

for x, ax in zip([gmass,smass],[ax10,ax11]):
    x1 = x[ratio>=1.]
    x2 = x[ratio<=1.]

    ax.scatter(x1,y1,label='Perpendicular',alpha=a)
    ax.scatter(x2,y2,label='Parallel',alpha=a)

for x, ax in zip([gmass,smass],[ax00,ax01]):
    
    x1 = x[ratio>=1.]
    x2 = x[ratio<=1.]

    ax.hist(x1, bins=30, log=True, alpha=a, label='Perpendicular', density=True)
    ax.hist(x2, bins=20, log=True, alpha=a, label='Parallel', density=True)

ax00.tick_params(axis="x", labelbottom=False)
ax01.tick_params(axis="x", labelbottom=False)
ax11.tick_params(axis="y", labelleft=False)
ax01.tick_params(axis="y", labelleft=False)

ax10.set_ylabel(r'$log_{10}(\frac{J_\perp}{/J_\parallel})$')
ax10.set_xlabel('Gass Mass')
ax11.set_xlabel('Star Mass')

ax00.legend()
ax11.legend()

plt.savefig('plot2.jpg')
# %%
"""
Scatter GMass/Smass vs Perp/Prll
"""
plt.scatter(gmass/smass,ratio)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Gas Mass / Star Mass')
plt.ylabel('Perpendicular / Parallel')
#plt.savefig('../plots/compRatio_massRatio.png')

# %%
