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
from config import writePath

# %%

for minradV in [7.]:
    voids = readVoids(minradV)
    gxs = readTNG()
    lbox = 205000.

    gxs.remove_row(np.where(gxs['y']==lbox)[0][0])
    gxs.remove_row(np.where(gxs['x']<0.)[0][0])

    for col in ['x','y','z']:
        gxs[col]/=1000.
    lbox /= 1000.

    Ntotal = len(gxs)
    tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)


    # %%
    r_array = np.linspace(0.5,50.,50) #Mpc

    for nv in range(len(voids)):
        #print(nv)
        x,y,z,radius = voids['x','y','z','r'][nv].as_void()
        rho_array = []
        for r in r_array:
            vol = 4.*np.pi*r**3./3
            rho_array.append( len(tree.query_ball_point([x,y,z],r))/vol)

        rho_total=Ntotal/lbox**3.
        plt.plot(r_array,np.array(rho_array)/rho_total-1)
        ascii.write(Table(np.column_stack([r_array,np.array(rho_array)/rho_total-1])),writePath+'Proyectos/Orientations/data/profiles/minradV{}_void{}.dat'.format(minradV,nv),names=['r','rho'],overwrite=True)
    plt.savefig('../plots/profiles_minradV{}.png'.format(minradV),dpi=400)
# %%
