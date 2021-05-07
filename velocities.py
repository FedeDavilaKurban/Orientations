#%%
"""
Velocities
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

    gals = JvsM(sec,gals,gxs,plot=False) #Determine galaxies of interest (involves rmin,rmax,sec)

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
    ids = np.where((np.log10(mass)>-1.)&(np.log10(mass)<3.))
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
Define Parameters, Reading Galaxies, Creating Tree
"""

#exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(sys.argv[1])
minradV, maxradV, rmin, rmax, sec, s5, vtype = 7., 0., .9, float(sys.argv[1]), 3, 0, 'r'
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

#print('writing')
#ascii.write(gxs,'gxs.dat')

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)
#%%
print('calculating')
"""
Read Galaxy Velocity and Spin Components in the shells of Voids
"""

#nvs = [0]
nvs = range(len(voids))
prll = []
perp = []
ratio = []

# xvel_ = []
# yvel_ = []
# zvel_ = []

vrad = []
vtra = []
mass = []

for nv in nvs:
    # Ahora que estoy trabajando con velocidades
    # esta nomenclatura es un poco confusa :P
    # Cuando la 'v' está primero, significa que es una cantidad del void
    # Ej.: vx es la posición x del void; gxv es la velocidad en el eje x de la galaxia
    vx, vy, vz = voids[nv]['x'], voids[nv]['y'], voids[nv]['z']

    # Devuelve 'gals', un subset de 'gxs'. 
    # Gals son las gxs que estan en la cascara y
    # de la seccion 'sec' en el espacio Spin-Masa.
    # Estos parametros estan definidos al principio
    gals, vshell = orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5) 
    # vshell son las componentes promedio de la velocidad 
    # de todas las galaxias del shell (previo a la selección en por Spin-Masa)

    for g in gals:

        gx, gy, gz = g['x'], g['y'], g['z']

        # A la velocidad de la galaxia le resto la "velocidad del shell" 
        gxv, gyv, gzv = g['xv']-vshell[0], g['yv']-vshell[1], g['zv']-vshell[2]
        #vel = np.sqrt(gxv**2+gyv**2,gzv**2)

        sx, sy, sz, sn = g['spx'], g['spy'], g['spz'], g['sp_n']

        # Spin Components
        r = [gx-vx,gy-vy,gz-vz] #galaxy position from void center
        r_versor = r/np.linalg.norm(r) #radial direction from void center
        s = [sx,sy,sz] #spin vector
        prll.append( abs( np.dot(r_versor,s) ) )
        perp.append( np.sqrt(sn**2 - prll[-1]**2) ) # S**2 = S_perp**2 + S_paral**2

        # Velocities
        vrad.append( np.dot([gxv,gyv,gzv],r_versor) )

        v_norm_squared = gxv**2 + gyv**2 + gzv**2
        vrad_norm_squared = vrad[-1]**2
        vtrans_squared = v_norm_squared-vrad_norm_squared
        #Este error no deberia saltar:
        if vtrans_squared<0.: 
            print('Negative Value for V_trans_squared')
        vtra.append( np.sqrt(vtrans_squared) )

        mass.append( g['mass'] )

perp = np.array(perp)
prll = np.array(prll)

ratio = perp/prll

vrad = np.array(vrad)
vtra = np.array(vtra)
mass = np.array(mass)
#%%

print('plotting')
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

p1, p2 = 10, 90
pinf = np.percentile(ratio,p1)
psup = np.percentile(ratio,p2)

#pinf, psup = 1.,1.

fig, axs = plt.subplots(3, 1, constrained_layout=True)

x = np.log10(ratio)
axs[0].hist(x,bins=30,density=True,alpha=1.)
axs[0].vlines(np.log10(pinf),0,.9)
axs[0].vlines(np.log10(psup),0,.9)
axs[0].set_xlabel(r'$log_{10}(\frac{J_\perp}{J_\parallel})$')
# axs[0].hist(x[x<=np.log10(pinf)],bins=15,density=False)
# axs[0].hist(x[x>=np.log10(psup)],bins=15,density=False)


axs[0].text(2.,.25,'Min-Max Void R: {}, {}\nMin-Max Shell R: {}, {}\nSpin-Mass Section: {}\ns5: {}\nVoid Type: {}'\
    .format(minradV, maxradV, rmin, rmax, sec, s5, vtype))

x1 = vrad[ratio>=psup]
x2 = vrad[ratio<=pinf]

axs[1].hist(x1,label='Perpendicular',bins=25,density=True,alpha=0.7)
axs[1].hist(x2,label='Parallel',bins=25,density=True,alpha=0.7)

axs[1].legend()
axs[1].set_xlabel(r'$v_{rad}$')

x1 = vtra[ratio>=psup]
x2 = vtra[ratio<=pinf]

axs[2].hist(x1,label='Perpendicular',bins=25,density=True,alpha=0.7)
axs[2].hist(x2,label='Parallel',bins=25,density=True,alpha=0.7)

axs[2].legend()
axs[2].set_xlabel(r'$v_{tra}$')

axs[0].set_xlim([-2,3])
axs[0].set_ylim([0,1])
axs[1].set_xlim([-800,800])
axs[1].set_ylim([0.,0.006])
axs[2].set_xlim([0,1400])
axs[2].set_ylim([0,0.008])

plt.savefig('../plots/velocities_{}_{}_{}_{}_{}_{}_{}.jpg'.format(minradV, maxradV, rmin, rmax, sec, s5, vtype))
# %%
import statsmodels.api as sm
x1 = vrad[ratio>=psup]
x2 = vrad[ratio<=pinf]
fig, ax = plt.subplots(1, 1, constrained_layout=True)
color = plt.rcParams['axes.prop_cycle'].by_key()['color']
sm.qqplot(x1, line =None, ax=ax, color=color[0], label='Perpendicular')
sm.qqplot(x2, line =None, ax=ax, color=color[1], label='Parallel')
plt.legend()
plt.savefig('../plots/qqplot_vrad.jpg')
#%%
# fig, axs = plt.subplots(3, 1, sharex=True, constrained_layout=True)

# a = .4

# spin = np.sqrt(perp**2+prll**2) 

# x1 = vrad 
# x2 = vtra
# y = np.log10(spin/mass )

# axs[0].scatter(x1, y, label='Vrad',alpha=a)
# axs[0].scatter(x2, y, label='Vtra',alpha=a)
# axs[0].legend()
# axs[0].set_ylabel(r'$log_{10}(J/M)$')

# y = np.log10(perp/mass)

# axs[1].scatter(x1, y, label='Vrad',alpha=a)
# axs[1].scatter(x2, y, label='Vtra',alpha=a)
# axs[1].set_ylabel(r'$log_{10}(J_\perp/M)$')

# y = np.log10(prll/mass)

# axs[2].scatter(x1, y, label='Vrad',alpha=a)
# axs[2].scatter(x2, y, label='Vtra',alpha=a)
# axs[2].set_xlabel('Velocity')
# axs[2].set_ylabel(r'$log_{10}(J_\parallel/M)$')

# plt.savefig('../plots/JvsV.jpg')

# %%
import matplotlib.colors as colors
fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, constrained_layout=True)

spin = np.sqrt(perp**2+prll**2) 
bins = 80


x = vtra
y = np.log10(spin/mass )

axs[0,1].hist2d(x,y,bins=bins,density=True,norm=colors.LogNorm())
#axs[0,1].set_ylabel(r'$log_{10}(J/M)$')

y = np.log10(perp/mass)

axs[1,1].hist2d(x,y,bins=bins,density=True,norm=colors.LogNorm())
#axs[1,1].set_ylabel(r'$log_{10}(J_\perp/M)$')

y = np.log10(prll/mass)

axs[2,1].hist2d(x,y,bins=bins,density=True,norm=colors.LogNorm())
axs[2,1].set_xlabel(r'$V_{tra}$')
#axs[2,1].set_ylabel(r'$log_{10}(J_\parallel/M)$')

#----------

x = vrad
y = np.log10(spin/mass )

axs[0,0].hist2d(x,y,bins=bins,density=True,norm=colors.LogNorm())
axs[0,0].set_ylabel(r'$log_{10}(J/M)$')

y = np.log10(perp/mass)

axs[1,0].hist2d(x,y,bins=bins,density=True,norm=colors.LogNorm())
axs[1,0].set_ylabel(r'$log_{10}(J_\perp/M)$')

y = np.log10(prll/mass)

axs[2,0].hist2d(x,y,bins=bins,density=True,norm=colors.LogNorm())
axs[2,0].set_xlabel(r'$V_{rad}$')
axs[2,0].set_ylabel(r'$log_{10}(J_\parallel/M)$')
plt.savefig('../plots/JvsV_hist._{}_{}_{}_{}_{}_{}_{}.jpg'.format(minradV, maxradV, rmin, rmax, sec, s5, vtype))
# %%
