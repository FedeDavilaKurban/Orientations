#%%
"""
Calculate Beta (J_perp/J_prll) for different bins in Radius
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
Reading Galaxies, Creating Tree
"""
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
"""
Define Parameters 
"""

#exp, minradV, maxradV, rmin, rmax, sec, s5, vtype = readExp(sys.argv[1])
minradV, maxradV, rmin, rmax, sec, s5, vtype = 7., 0., .4, 1.6, 3, 0, 'a'
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

# %%
"""
Calculating
"""

def get_beta():
    nvs = range(len(voids))
    prll = []
    perp = []
    ratio = []
    r = []
    rrv = []

    for nv in nvs:

        vx, vy, vz, vr = voids[nv]['x'], voids[nv]['y'], voids[nv]['z'], voids[nv]['r']
        gals, vshell = orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5)

        for g in gals:

            gx, gy, gz = g['x'], g['y'], g['z']

            # A la velocidad de la galaxia le resto la "velocidad del shell/void" 
            #gxv, gyv, gzv = g['xv']-vshell[0], g['yv']-vshell[1], g['zv']-vshell[2]

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

    return beta,perp,prll,r,rrv

for vtype in ['a']:
    print(vtype)
    for rmin,rmax in zip([0.6],[0.8]):

        voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)

        print(rmin,rmax)
        beta, perp, prll, r, rrv = get_beta()

        plt.rcParams['figure.figsize'] = (10, 8)
        plt.rcParams['font.size'] = 15


        x = np.log10(beta)
        n_perp = len(np.where(x>0.)[0])
        n_prll = len(np.where(x<0.)[0])
        b_factor = n_perp / n_prll 

        plt.hist(x,bins=30,density=True)
        plt.vlines(0,0,1.,linestyles=':')

        plt.xlabel(r'$log_{10}\beta$')

        plt.text(-2.5,.6,r'$n_2(\beta<0)={}$'.format(n_prll))
        plt.text(1.,.6,r'$n_1(\beta>0)={}$'.format(n_perp))
        plt.text(3.,.2,r'$n_1/n_2={:.3f}$'.format(b_factor))

        plt.xlim([-3.,5.])
        plt.ylim([0.,1.])

        #plt.savefig('../plots/beta_forRbins/voids_{}/rmin{}_rmax{}.jpg'.\
        #            format(vtype, rmin, rmax))
        #plt.close()
# %%
x = [1,2,3,4,5,6]
ya = [2.462, 2.606, 2.517, 2.429, 2.410, 2.379]
yr = [2.489, 2.576, 2.550, 2.439, 2.414, 2.366]
ys = [2.286, 2.783, 2.365, 2.394, 2.398, 2.423]
plt.plot(x,ya,label='All')
plt.plot(x,yr,label='R')
plt.plot(x,ys,label='S')

plt.ylabel(r'$n_1/n_2$')
plt.xlabel('R/Rv')
plt.legend()

x_ticks_labels = ['0.4-0.6','0.6-0.8','0.8-1.0','1.0-1.2','1.2-1.4','1.4-1.6']
plt.xticks(x, x_ticks_labels)

plt.savefig('../plots/beta_forRbins/factorforallvoids.jpg')
# %%
