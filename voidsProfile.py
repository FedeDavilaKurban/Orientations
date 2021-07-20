#%%
import sys
sys.path.append('/home/fede/')
import illustris_python as il
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import numpy as np
from scipy import spatial
import seaborn
import scipy.stats
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit
from orientationsTools import *
from config import writePath, units

# # %%

# for minradV in [7.]:
#     voids = readVoids(minradV)
#     gxs = readTNG()
#     lbox = 205000.

#     gxs.remove_row(np.where(gxs['y']==lbox)[0][0])
#     gxs.remove_row(np.where(gxs['x']<0.)[0][0])

#     for col in ['x','y','z']:
#         gxs[col]/=1000.
#     lbox /= 1000.

#     Ntotal = len(gxs)
#     tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)


#     r_array = np.linspace(0.5,50.,50) #Mpc

#     for nv in range(len(voids)):
#         #print(nv)
#         x,y,z,radius = voids['x','y','z','r'][nv].as_void()
#         rho_array = []
#         for r in r_array:
#             vol = 4.*np.pi*r**3./3
#             rho_array.append( len(tree.query_ball_point([x,y,z],r))/vol)

#         rho_total=Ntotal/lbox**3.
#         plt.plot(r_array,np.array(rho_array)/rho_total-1)
#         ascii.write(Table(np.column_stack([r_array,np.array(rho_array)/rho_total-1])),writePath+'Proyectos/Orientations/data/profiles/minradV{}_void{}.dat'.format(minradV,nv),names=['r','rho'],overwrite=True)
#     plt.savefig('../plots/profiles_minradV{}.png'.format(minradV),dpi=400)

# %%
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
    ids = np.where((np.log10(mass)>-1))#&(np.log10(mass)<3.))
    mass = mass[ids]

    pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])[ids]

    # spin = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSpin'])[ids]
    # sp_n = np.sum(np.abs(spin)**2,axis=-1)**(1./2)
    
    # vel = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloVel'])[ids]

    # sfr = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSFR'])[ids]
    
    # masstype = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMassType'])[ids]                                                                                                                      
    # gasmass = masstype[:,0]
    # starmass = masstype[:,4]


    gxs = Table(np.column_stack([pos[:,0],pos[:,1],pos[:,2]\
                                ,mass\
                                #,spin,sp_n\
                                #,vel[:,0],vel[:,1],vel[:,2]\
                                #,gasmass,starmass,\
                                ]),
                names=['x','y','z','mass'])#,'spx','spy','spz','sp_n'])#,'xv','yv','zv'])#,'gasmass','starmass']) 

    del mass,pos#,spin,sp_n#,vel#,masstype,gasmass,starmass

    return gxs

# def get_gals_forprofile(gxs,tree,units,voids,nv,rmin,rmax,sec,s5,vshell=False):
#     """
#     Determines galaxies of interest in the shell of the void 
#     Returns cosine of angles of disks and void centric direction    
#     """
#     # Esta funcion es una modificacion del orientations() que esta en el orientationTols
#     import numpy as np
#     from astropy.table import Table

#     global vx,vy,vz,vr,cos

#     vx,vy,vz,vr = voids['x','y','z','r'][nv]
#     if units=='kpc':
#         vx*=1000.
#         vy*=1000.
#         vz*=1000.
#         vr*=1000.

#     idx1 = tree.query_ball_point([vx,vy,vz],vr*rmax)
#     idx2 = tree.query_ball_point([vx,vy,vz],vr*rmin)
#     shell = [g for g in idx1 if g not in idx2]
#     gals = gxs[shell]
#     #print('N of galaxies in shell:',len(gals))

#     if vshell:
#         xv_mean = np.mean(gals['xv'])
#         yv_mean = np.mean(gals['yv'])
#         zv_mean = np.mean(gals['zv'])
#         vshell_mean = [xv_mean, yv_mean, zv_mean]

#     # if sec!=0:
#     #     gals = JvsM_(sec,gals,jvsm_b,jvsm_m,plot=False) #Determine galaxies of interest (involves rmin,rmax,sec)

#     #cos = cosCalc(gals_h,units,tree,xv,yv,zv,rv,s5) #Calculate the cosines of angle of J and void direction

#     if vshell:
#         return gals,vshell_mean 
#     else: return gals

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

for col in ['x','y','z']:
    gxs[col]/=1000.
lbox /= 1000.

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)
#%%
import seaborn as sns

plt.rcParams['figure.figsize'] = (8, 8)
plt.rcParams['font.size'] = 15
fig, axs = plt.subplots(3, 1, constrained_layout=True, sharex=True, sharey=True)

color = sns.color_palette()

minradV=7.
maxradV=0.

rinner_i = np.float64(0.8)
rinner_f = np.float64(1.5)
rstep = np.float64(0.1)
r1 = np.arange(rinner_i,rinner_f,rstep,dtype=np.float64)
r2 = np.arange(rinner_i+rstep,rinner_f+rstep,rstep,dtype=np.float64)
r_array = (r1+r2)/2

r_array = np.linspace(8,16.,30)

Ntotal = len(gxs)
rhototal = len(gxs)/(lbox**3)

r_array = np.linspace(6,10.6,30) #Rv

for vtype, c, ax in zip(['a','r','s'],color[:3],axs):
    print('vtype=',vtype)

    voids = readVoids(minrad=minradV,maxrad=maxradV,vtype=vtype)


    for nv in range(len(voids)):
        #print(nv)
        x,y,z,radius = voids['x','y','z','r'][nv].as_void()
        
        rho_array = []
        for r in r_array:
            vol = 4.*np.pi*(r)**3./3
            rho_array.append( len(tree.query_ball_point([x,y,z],r))/vol)

        ax.plot(r_array/radius,np.array(rho_array)/rhototal-1,marker='o',ms=1,lw=.5,color=c)
        ax.set_ylabel(r'$\rho(r)/\bar{\rho}-1$')

axs[2].set_xlim([.8,1.5])
axs[0].set_xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
axs[1].set_xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])
axs[2].set_xticks([0.8,.9,1.,1.1,1.2,1.3,1.4,1.5])

axs[0].text(0.825,-0.1,'All Voids')
axs[1].text(0.825,-0.1,'Rising Voids')
axs[2].text(0.825,-0.1,'Shell Voids')

axs[2].set_xlabel('R/Rv')

plt.savefig('../plots/vel/voidsprofiles.jpg')
# %%
