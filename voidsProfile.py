import sys
sys.path.append('/home/fede/')
import illustris_python as il
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import illustris_python as il
import numpy as np
from scipy import spatial
import seaborn
import scipy.stats
import matplotlib.gridspec as gridspec
from statsmodels.distributions.empirical_distribution import ECDF
from scipy.optimize import curve_fit

#%%
def readVoids(minrad=None,maxrad=None):
    """
    This function reads the file tng300-1_voids.dat
    'minrad': filter by a minimum radius (Mpc). Optional.
    'maxrad': filter by maximum radius (Mpc). Optional
    """
    voids=ascii.read('../data/tng300-1_voids.dat',names=['r','x','y','z','vx','vy','vz','deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter'])#,usecols=['r','x','y','z'])

    if minrad != None: voids=voids[np.where(voids['r']>=minrad)]
    if maxrad != None: voids=voids[np.where(voids['r']<=maxrad)]
    
    return voids
# %%
def readTNG(units,fields=['SubhaloPos','SubhaloMass']):
    """
    This function reads subhalos in the TNG300-1 simulation and returns 
    lbox: size of box in specified units
    Ntotal: num of galaxies
    gxs: an ascii Table with the fields and filters I usually need

    units=['Mpc','kpc']

    """
    basePath = '/home/fede/TNG300-1/output/'

    data = []
    if 'SubhaloPos' in fields:
        gxs = Table(il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos']),names=['x','y','z'] )
    if 'SubhaloMass' in fields: 
        gxs.add_column(Table.Column(il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMass'])),name='mass')
    if 'SubhaloSpin' in fields: 
        spin = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSpin'])
        sp_n = np.sum(np.abs(spin)**2,axis=-1)**(1./2)
        gxs.add_column(Table.Column(sp_n), name='sp_n')
        del spin, sp_n
    if 'SubhaloSFR' in fields: 
        gxs.add_column( il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSFR']),name='SFR' )
    if units=='Mpc': 
        lbox=205.
        for col in ['x','y','z']:
            gxs[col]/=1000.
    if units=='kpc': 
        lbox=205000.

    #Applying some filters
    gxs = gxs[(np.log10(gxs['mass'])>-1.)&(np.log10(gxs['mass'])<3)]
    gxs.remove_row(np.where(gxs['y']==lbox)[0][0])
    gxs.remove_row(np.where(gxs['x']<0.)[0][0])

    return lbox,gxs
# %%
voids=readVoids(7.)
lbox,Ntotal,gxs=readTNG('Mpc')
tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

# %%
r_array = np.linspace(0.5,50.,50) #Mpc

for nv in range(len(voids)):
    print(nv)
    x,y,z,radius = voids['x','y','z','r'][nv].as_void()
    rho_array = []
    for r in r_array:
        vol = 4.*np.pi*r**3./3
        rho_array.append( len(tree.query_ball_point([x,y,z],r))/vol)

    rho_total=Ntotal/lbox**3.
    plt.plot(r_array,np.array(rho_array)/rho_total-1)
plt.savefig('profiles.png',dpi=400)
# %%
