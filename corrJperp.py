#%%
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

def get_gals(gxs,tree,units,voids,nv,rmin,rmax,sec,s5,jvsm_b,jvsm_m,vshell=False):
    """
    Determines galaxies of interest in the shell of the void 
    Returns cosine of angles of disks and void centric direction    
    """
    # Esta funcion es una modificacion del orientations() que esta en el orientationTols
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

    if vshell:
        xv_mean = np.mean(gals['xv'])
        yv_mean = np.mean(gals['yv'])
        zv_mean = np.mean(gals['zv'])
        vshell_mean = [xv_mean, yv_mean, zv_mean]

    if sec!=0:
        gals = JvsM_(sec,gals,jvsm_b,jvsm_m,plot=False) #Determine galaxies of interest (involves rmin,rmax,sec)

    #cos = cosCalc(gals_h,units,tree,xv,yv,zv,rv,s5) #Calculate the cosines of angle of J and void direction

    if vshell:
        return gals,vshell_mean 
    else: return gals

def JvsM_(sec,gals,b,m,plot=True):
    """
    Performs linear regression of the Mass-Spin relation
    Determines galaxies of interest with 'sec'
    Returns galaxies of interest in 'gals_h'
    """
    import numpy as np
    import matplotlib.pyplot as plt

    m1 = -.8
    m2 = -.4

    M = np.log10(gals['mass'])
    S = np.log10(gals['sp_n'])

    #Alto Spin
    if sec == 1: gals_h = gals[(M < m1)&(S > M*m+b)]  
    if sec == 2: gals_h = gals[(M > m1)&(M < m2)&(S > M*m+b)]  
    if sec == 3: gals_h = gals[(M > m2)&(S > M*m+b)]  
    if sec == 12: gals_h = gals[(M < m2)&(S > M*m+b)]  
    if sec == 23: gals_h = gals[(M > m1)&(S > M*m+b)]  
    if sec == 123: gals_h = gals[(S > M*m+b)]  

    #Bajo Spin
    if sec == 4: gals_h = gals[(M < m1)&(S < M*m+b)]  
    if sec == 5: gals_h = gals[(M > m1)&(M < m2)&(S < M*m+b)]  
    if sec == 6: gals_h = gals[(M > m2)&(S < M*m+b)]  
    if sec == 45: gals_h = gals[(M < m2)&(S < M*m+b)]  
    if sec == 56: gals_h = gals[(M > m1)&(S < M*m+b)]  
    if sec == 456: gals_h = gals[(S < M*m+b)] 

    #Solo Masa
    if sec == 14: gals_h = gals[(M < m1)] 
    if sec == 25: gals_h = gals[(M > m1)&(M < m2)] 
    if sec == 36: gals_h = gals[(M > m2)]

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
    ids = np.where((np.log10(mass)>-1)&(np.log10(mass)<3.))
    mass = mass[ids]

    pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])[ids]

    spin = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSpin'])[ids]
    sp_n = np.sum(np.abs(spin)**2,axis=-1)**(1./2)


    gxs = Table(np.column_stack([pos[:,0],pos[:,1],pos[:,2]\
                                ,mass\
                                ,spin,sp_n\
                                #,vel[:,0],vel[:,1],vel[:,2]\
                                #,gasmass,starmass,\
                                ]),
                names=['x','y','z','mass','spx','spy','spz','sp_n'])#,'xv','yv','zv'])#,'gasmass','starmass']) 

    del mass,pos,spin,sp_n#,vel#,masstype,gasmass,starmass

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

# Linear regression on Spin-Mass
jvsm_m,jvsm_b,rvalue,pvalue,std = scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) 


tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)
#%%

minradV=7.
maxradV=0.
s5=0
sec = 0

r1 = [.9]
r2 = [1.4]

prll = []
perp = []

r = []
rrv = []

r_versor = []
phi = []
tita = []


for vtype in ['r']:
    print('vtype=',vtype)
    for rmin,rmax in zip(r1,r2):
        print('rmin, rmax =',rmin,rmax)

        voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)

        nvs = range(len(voids))
        #nvs=[0]

        for nv in nvs:

            vx, vy, vz, vr = voids[nv]['x'], voids[nv]['y'], voids[nv]['z'], voids[nv]['r']
            gals = get_gals(gxs,tree,units,voids,nv,rmin,rmax,sec,s5,jvsm_b,jvsm_m,vshell=False)

            # for g in gals:

            #     gx, gy, gz = g['x'], g['y'], g['z']

            #     sx, sy, sz, sn = g['spx'], g['spy'], g['spz'], g['sp_n']

            #     s = [sx,sy,sz] #spin vector

            #     gr = [gx-vx,gy-vy,gz-vz] #galaxy position from void center
            #     for axis in range(len(gr)):
            #         if gr[axis]<=-(rmax+.5)*vr: gr[axis]+=lbox
            #         if gr[axis]>= (rmax+.5)*vr: gr[axis]-=lbox

            #     r.append( np.linalg.norm(gr) )
            #     rrv.append(r[-1]/vr)
            #     r_versor.append( gr/r[-1] ) #radial direction from void center

            #     prll.append( np.dot(r_versor[-1],s)*r_versor[-1]  )
            #     perp.append( s - prll[-1] ) 

            #     tita.append( np.pi/2-np.arccos(gr[2]/np.sqrt(gr[0]**2+gr[1]**2+gr[2]**2)) )
            #     phi.append( np.arctan2(gr[1], gr[0]) )

            gx, gy, gz = gals['x'], gals['y'], gals['z']

            sx, sy, sz, sn = gals['spx'], gals['spy'], gals['spz'], gals['sp_n']

            s = [sx.data,sy.data,sz.data] #spin vector

            gr = [gx.data-vx.data,gy.data-vy.data,gz.data-vz.data] #galaxy position from void center
            for i in range(len(gr)):
                for j in range(len(gr[i])):
                    if gr[i][j]<=-(rmax+.5)*vr: gr[i][j]+=lbox
                    if gr[i][j]>= (rmax+.5)*vr: gr[i][j]-=lbox

            r.append( np.linalg.norm(gr,axis=0) )
            rrv.append(r[-1]/vr)
            r_versor.append( gr/r[-1] ) #radial direction from void center

            prll.append( np.multiply(r_versor[-1], s).sum(0)*r_versor[-1]  )
            perp.append( s - prll[-1] ) 

            tita.append( np.pi/2-np.arccos(gr[2]/np.sqrt(gr[0]**2+gr[1]**2+gr[2]**2)) )
            phi.append( np.arctan2(gr[1], gr[0]) )


        # r = np.array(r)
        # rrv = np.array(rrv)
        # r_versor = np.array(r_versor)

        # perp =  np.array(perp)
        # prll = np.array(prll)

        # tita = np.array(tita)
        # phi = np.array(phi)

#%%
#Para ver cuantas gxs hay en c/void
#s=0
for p in phi:
    s=len(p)
    print(s)

#%%
# # Build new basis
# ez = r_versor[0]

# ex = np.random.randn(3)  # take a random vector
# ex -= ex.dot(ez) * ez / np.linalg.norm(ez)**2    # make it orthogonal to ez
# ex /= np.linalg.norm(ex) # normalize it

# ey = np.cross(ez, ex)

# # print(np.linalg.norm(ex), np.linalg.norm(ey)) #Deberian dar 1
# # print(np.cross(ex, ey)-ez) #Deberian dar 0
# # print(np.dot(ex, ey),np.dot(ex,ez),np.dot(ez,ey)) #Deberian dar 0

# newbasis = [ey,ex,ez]

# print(perp[0])
# vec_new = np.linalg.solve(np.linalg.inv(newbasis),perp[0])
# print(vec_new) #Estoy buscando que la tercera componente me de 0

#%%
nv=27
print(len(phi[nv]))

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='mollweide')
ax.grid(True)
ax.scatter(phi[nv],tita[nv],s=4)
#ax.quiver(phi,tita,perp[0][0],perp[0][1],color='k',width=.005)

plt.show()
# %%

vrad = ascii.read('../data/vel/-1/vel_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
    .format(minradV,maxradV,rmin,rmax,sec,vtype))['vrad'].data
vrad = vrad[len(phi[nv-1]):len(phi[nv-1])+len(phi[nv])]
vmed = np.median(vrad)
print(vmed)
#%%

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='mollweide')
ax.grid(True)
sc=ax.scatter(phi[nv],tita[nv],s=2,c=vrad,cmap='plasma_r')
plt.colorbar(sc)
plt.show()

#%%
phi_v = phi[nv][vrad<vmed]
tita_v = tita[nv][vrad<vmed]

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='mollweide')
ax.grid(True)
sc=ax.scatter(phi_v,tita_v,s=5,c=vrad[vrad<vmed],cmap='plasma_r',alpha=.6)
plt.colorbar(sc)
plt.show()

    # %%
phi_1 = phi[nv][vrad<vmed]
tita_1 = tita[nv][vrad<vmed]

phi_2 = phi[nv][vrad>vmed]
tita_2 = tita[nv][vrad>vmed]

fig = plt.figure(figsize=(10,6))
ax = fig.add_subplot(111, projection='mollweide')
ax.grid(True)
ax.scatter(phi_1,tita_1,s=8,color='r',alpha=.7)
ax.scatter(phi_2,tita_2,s=8,color='b',alpha=.7)
plt.show()

# %%
