#%%
"""
Sigma5 vs R
"""
import sys
import numpy as np
from scipy import spatial
from scipy import stats
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
minradV, maxradV, rmin, rmax, sec, s5, vtype = 7., 0., .5, 1.5, 3, 0, 'a'
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

# Hago otro tree con galaxias con masa > 10**11
ids = np.where((np.log10(gxs['mass'])>1.))
gxs11 = gxs[ids]
tree11 = spatial.cKDTree(data=np.column_stack((gxs11['x'],gxs11['y'],gxs11['z'])),boxsize=lbox)

#%%
print('calculating')
"""
Obtain a linear regression for Sigma5 vs R
"""

#nvs = [0]
nvs = range(len(voids))

sig5_list=[]
r_list=[]

for nv in nvs:
    print(nv)
    vx, vy, vz, vr = voids[nv]['x'], voids[nv]['y'], voids[nv]['z'], voids[nv]['r']

    gals, vshell = orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5) 

    data = np.column_stack((gals['x'],gals['y'],gals['z']))
    sig5_, inds5 = tree11.query(data, k=[6])  # k1 = identity
    sig5_list.append(sig5_)

    rx = gals['x']-vx
    ry = gals['y']-vy
    rz = gals['z']-vz


    for ri in [rx,ry,rz]: 
        ri[np.where(ri>2*vr)]-=lbox
        ri[np.where(ri<-2*vr)]+=lbox

    #plt.scatter(rx,ry,s=.5,c='blue')

    rs = [rx, ry, rz]
    #Normalized distance from void center:
    r_list.append( np.linalg.norm(rs, axis=0)/vr  )  

r = np.concatenate(r_list)
sig5 = np.concatenate(sig5_list)
sig5 = np.concatenate(sig5)

m,b,rvalue,pvalue,std = stats.linregress(x=r,y=sig5)

#%%
plt.rcParams['figure.figsize'] = (10, 8)
plt.rcParams['font.size'] = 15

#plt.vlines(1,0,3000,linestyles=':')
plt.scatter(r,sig5,alpha=.5,s=1)
plt.xlabel(r'$R/R_{v}$')
plt.ylabel(r'$\Sigma_5\,(kpc)$')
#plt.xlim(.5, 1.5)

plt.plot(r,r*m+b)

# %%
"""
Ahora que tengo la regressión lineal Sig5-R 
Leo de vuelta las galaxias en los shells (gals)
Guardo las componentes del Spin diferenciando entre galaxias
de Sig5 alto y bajo
"""
# hs5: high Sigma5
# ls5: low SIgma5
perp_hs5 = [] 
prll_hs5 = []
perp_ls5 = []
prll_ls5 = []
nvs = range(len(voids))

for nv in nvs:
    #print(nv)
    vx, vy, vz, vr = voids[nv]['x'], voids[nv]['y'], voids[nv]['z'], voids[nv]['r']

    gals, vshell = orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5) 

    rx = gals['x']-vx
    ry = gals['y']-vy
    rz = gals['z']-vz

    for ri in [rx,ry,rz]: 
        ri[np.where(ri>1.5*vr)]-=lbox
        ri[np.where(ri<-1.5*vr)]+=lbox

    rs = np.column_stack((rx.data,ry.data,rz.data))
    # Distance from void center:
    r_norms = np.linalg.norm(rs, axis=1)
    r_versors = np.zeros((np.shape(rs)))
    for i in range(len(r_norms)):
        r_versors[i] = rs[i]/r_norms[i]    

    s_vectors = np.column_stack((gals['spx'].data,\
                                gals['spy'].data,\
                                gals['spz'].data))

    sig5, inds5 = tree11.query(np.column_stack((gals['x'],gals['y'],gals['z'])), k=[6])  # k1 = identity

    for i in range(len(gals)):

        if sig5[i] > m*r_norms[i]/vr+b:

            prll_hs5.append(abs( np.dot(r_versors[i],s_vectors[i]))) 
            # S**2 = S_perp**2 + S_paral**2
            perp_hs5.append(np.sqrt(gals['sp_n'][i]**2 - prll_hs5[-1]**2))  
        
        if sig5[i] < m*r_norms[i]/vr+b:

            prll_ls5.append(abs( np.dot(r_versors[i],s_vectors[i]))) 
            perp_ls5.append(np.sqrt(gals['sp_n'][i]**2 - prll_ls5[-1]**2))  


perp_hs5 = np.array(perp_hs5)
prll_hs5 = np.array(prll_hs5)
perp_ls5 = np.array(perp_ls5)
prll_ls5 = np.array(prll_ls5)

beta_hs5 = perp_hs5/prll_hs5
beta_ls5 = perp_ls5/prll_ls5

plt.hist(np.log10(beta_hs5),bins=80,density=True)
plt.hist(np.log10(beta_ls5),bins=80,density=True,alpha=.7)

#plt.hist(np.log10(perp),label='Perpendicular',bins=50,density=True,alpha=.8)
#plt.hist(np.log10(prll),label='Parallel',bins=50,density=True,alpha=.8)
#plt.legend()
plt.show()


#%%


vrad = np.array(vrad)
vtra = np.array(vtra)
mass = np.array(mass)
