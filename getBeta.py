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
    
    # vel = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloVel'])[ids]

    # sfr = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSFR'])[ids]
    
    # masstype = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMassType'])[ids]                                                                                                                      
    # gasmass = masstype[:,0]
    # starmass = masstype[:,4]


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

# rinner_i = np.float64(0.6)
# rinner_f = np.float64(1.5)
# rstep = np.float64(0.1)
# r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
# r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)

r1 = np.array([0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4])
r2 = np.array([0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5])
r1 = [.9]
r2 = [1.4]


for sec in [0,123,456]:
    print('sec=',sec)
    for vtype in ['r','s','a']:
        print('vtype=',vtype)
        for rmin,rmax in zip(r1,r2):
            print('rmin, rmax =',rmin,rmax)

            voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)

            nvs = range(len(voids))
            prll = []
            perp = []

            r = []
            rrv = []

            for nv in nvs:

                vx, vy, vz, vr = voids[nv]['x'], voids[nv]['y'], voids[nv]['z'], voids[nv]['r']
                gals = get_gals(gxs,tree,units,voids,nv,rmin,rmax,sec,s5,jvsm_b,jvsm_m,vshell=False)

                for g in gals:

                    gx, gy, gz = g['x'], g['y'], g['z']

                    sx, sy, sz, sn = g['spx'], g['spy'], g['spz'], g['sp_n']

                    s = [sx,sy,sz] #spin vector

                    gr = [gx-vx,gy-vy,gz-vz] #galaxy position from void center
                    for axis in range(len(gr)):
                        if gr[axis]<=-(rmax+.5)*vr: gr[axis]+=lbox
                        if gr[axis]>= (rmax+.5)*vr: gr[axis]-=lbox

                    r.append( np.linalg.norm(gr) )
                    rrv.append(r[-1]/vr)
                    r_versor = gr/r[-1] #radial direction from void center

                    prll.append( abs( np.dot(r_versor,s) ) )
                    perp.append( np.sqrt(sn**2 - prll[-1]**2) ) # S**2 = S_perp**2 + S_paral**2


            r = np.array(r)
            rrv = np.array(rrv)

            perp =  np.array(perp)
            prll = np.array(prll)
            beta = perp/prll

            #print(rrv)
            
            ascii.write(np.column_stack([beta]),\
                        '../data/beta/-1/beta_minradV{}_maxradV{}_rmin{:.1f}_rmax{:.1f}_sec{}_vtype{}.txt'\
                        .format(minradV,maxradV,rmin,rmax,sec,vtype),\
                        names=['beta'],\
                        overwrite=True)
# %%
