#%%
import sys
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import cosCalc, readVoids, JvsM
import random
from config import writePath, units
import matplotlib.pyplot as plt

def readTNG_():
    """
    This function reads subhalos in the TNG300-1 simulation and returns 

    gxs: an ascii Table with the fields and filters I usually need for this: Position, Mass, Spin

    """
    import sys
    sys.path.append('/home/fede/')
    import illustris_python as il
    import numpy as np
    from astropy.table import Table

    basePath = '/home/fede/TNG300-1/output/'

    mass = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMass'])                                                                                                                      
    ids = np.where((np.log10(mass)>-1.)&(np.log10(mass)<3.))
    mass=mass[ids]

    pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])
    pos = pos[ids]

    spin = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSpin'])
    spin = spin[ids]
    sp_n = np.sum(np.abs(spin)**2,axis=-1)**(1./2)

    sfr = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSFR'])
    sfr = sfr[ids]

    gxs = Table(np.column_stack([pos[:,0],pos[:,1],pos[:,2],mass,spin,sp_n,sfr]),names=['x','y','z','mass','spx','spy','spz','sp_n','sfr'])    
    del mass,pos,spin,sp_n,sfr

    return gxs

def orientations_(gxs,tree,units,voids,nv,rmin,rmax,sec,s5):
    """
    Determines galaxies of interest in the shell of the void 
    Returns cosine of angles of disks and void centric direction    
    """
    import numpy as np
    from astropy.table import Table

    #global xv,yv,zv,rv,cos

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

    gals_h = JvsM(sec,gals,gxs,plot=False) #Determine galaxies of interest (involves rmin,rmax,sec,s5)
    sfr_h = gals_h['sfr']

    cos = cosCalc(gals_h,units,tree,xv,yv,zv,rv,s5) #Calculate the cosines of angle of J and void direction

    return cos,sfr_h 

#··············································
#%%
print('Reading TNG')
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
print('Reading Voids')
voids = readVoids(7.,0.,'a')
print('Num of voids:',len(voids))

#%%

cos_list=[]
sfr_list=[]
for nv in range(len(voids)):
    #print(nv,'/',len(voids))
    cos, sfr = orientations_(gxs,tree,units,voids,nv,rmin=0.9,rmax=1.7,sec=0,s5=0)
    cos_list.append(cos)
    sfr_list.append(sfr)

    #ax.scatter(cos[np.where(sfr!=0)],sfr[np.where(sfr!=0)],color='b')

cos_flatlist = [i for j in cos_list for i in j]
sfr_flatlist = [i for j in sfr_list for i in j]

cos_array = np.array(cos_flatlist)
sfr_array = np.array(sfr_flatlist)

cos = cos_array[sfr_array>0.]
sfr = sfr_array[sfr_array>0.]

fig, ax = plt.subplots() 
plt.hist2d(cos,np.log10(sfr),bins=[10,10])

plt.xlabel(r'$\mathrm{cos}\theta$')
plt.ylabel('SFR')
plt.colorbar() 
plt.savefig('../plots/sfr.png')

# %%
