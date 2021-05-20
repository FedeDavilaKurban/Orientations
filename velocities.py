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
minradV, maxradV, rmin, rmax, sec, s5, vtype = 7., 0., .25, 1.5, 3, 0, 'r'
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

r = []
rrv = []

for nv in nvs:
    # Ahora que estoy trabajando con velocidades
    # esta nomenclatura es un poco confusa :P
    # Cuando la 'v' está primero, significa que es una cantidad del void
    # Ej.: vx es la posición x del void; gxv es la velocidad en el eje x de la galaxia
    vx, vy, vz, vr = voids[nv]['x'], voids[nv]['y'], voids[nv]['z'], voids[nv]['r']

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

        sx, sy, sz, sn = g['spx'], g['spy'], g['spz'], g['sp_n']

        # Spin Components
        gr = [gx-vx,gy-vy,gz-vz] #galaxy position from void center
        for axis in range(len(gr)):
            if gr[axis]<=-(rmax+.5)*vr: gr[axis]+=lbox
            if gr[axis]>= (rmax+.5)*vr: gr[axis]-=lbox


        r.append( np.linalg.norm(gr) )
        rrv.append(r[-1]/vr)
        r_versor = gr/r[-1] #radial direction from void center
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
v = np.sqrt(vrad**2+vtra**2)

r = np.array(r)
rrv = np.array(rrv)

mass = np.array(mass)
#%%
"""
Velocidad Radial vs R
"""
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15
color = plt.rcParams['axes.prop_cycle'].by_key()['color']

x = np.linspace(rmin,rmax,20)

y = []
for i in range(len(x)-1):
    y.append(np.mean(vrad[np.where(np.logical_and(rrv>x[i],rrv<x[i+1]))]))

x_step = (x[0]+x[1])/2
plt.scatter(x[:-1]+x_step,y)
plt.ylabel('Vrad (km/s)')
plt.xlabel('R/Rv')
plt.savefig('../plots/')

# Radio "dinámico" = Radio correspondiente a Vmax
rdin = x[np.where(y==np.max(y[15:]))][0]+x_step
rdin = 1.
print(rdin)
# %%
"""
Beta vs Velocidad
"""
outer_rdin = np.where(rrv>=rdin)[0]
inner_rdin = np.where(rrv<rdin)[0]

beta_inner = perp[inner_rdin]/prll[inner_rdin]
beta_outer = perp[outer_rdin]/prll[outer_rdin]

v_inner = vrad[inner_rdin]
v_outer = vrad[outer_rdin]

# Inner
#x = np.linspace(np.min(v_inner),np.max(v_inner),20)
x = np.linspace(0,200,10)

y = []
yerr = []
for i in range(len(x)-1):
    filter = np.where(np.logical_and(v_inner>x[i],v_inner<x[i+1]))[0]
    n = len(filter)
    y.append(np.mean(beta_inner[filter]))
    yerr.append(np.std(beta_inner[filter])/np.sqrt(n))
#plt.scatter(v_inner,np.log10(beta_inner),s=.01,color='k')

x_step = (x[0]+x[1])/2
plt.plot((x[:-1])+x_step,np.log10(y),label='Inner')
plt.errorbar((x[:-1])+x_step,np.log10(y),yerr=np.log10(yerr),color=color[0])


# Outer
#x = np.linspace(np.min(v_outer),np.max(v_outer),20)
y = []
yerr = []
for i in range(len(x)-1):
    filter = np.where(np.logical_and(v_outer>x[i],v_outer<x[i+1]))[0]
    n = len(filter)
    y.append(np.mean(beta_outer[filter]))
    yerr.append(np.std(beta_outer[filter])/np.sqrt(n))
x_step = (x[0]+x[1])/2
plt.plot((x[:-1])+x_step, np.log10(y), label='Outer')
plt.errorbar((x[:-1])+x_step, np.log10(y), yerr=np.log10(yerr), color=color[1])


plt.legend()
plt.xlabel('V (kms/s)')
plt.ylabel(r'$log_{10}(\beta)$')
# %%
"""
Hist de Betas
"""

plt.hist(np.log10(beta_inner),bins=22,density=True)
plt.hist(np.log10(beta_outer),bins=20,density=True,alpha=.5)
# %%
