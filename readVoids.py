import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
import illustris_python as il
import numpy as np
from scipy import spatial
import seaborn

voids=ascii.read('../data/tng300-1_voids.dat',names=['r','x','y','z','vx','vy','vz','deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter'])
tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)
