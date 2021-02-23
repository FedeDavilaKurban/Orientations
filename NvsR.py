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

lbox=205.
basePath = '/home/fede/TNG300-1/output/'
fields = ['SubhaloPos','SubhaloMass']
subgroups = il.groupcat.loadSubhalos(basePath,99,fields=fields)

gxs = Table(subgroups['SubhaloPos'],names=['x','y','z'],dtype=['float64','float64','float64'])

col_m = Table.Column(subgroups['SubhaloMass'],name='mass',dtype='float64')
gxs.add_column(col_m)

gxs = gxs[(np.log10(gxs['mass'])>-1.)&(np.log10(gxs['mass'])<3.)]

for col in ['x','y','z']:
	gxs[col]/=1000.
gxs.remove_row(np.where(gxs['y']==205.0)[0][0])
gxs.remove_row(np.where(gxs['x']<0.)[0][0])

tree = spatial.cKDTree(data=np.column_stack((gxs['x'],gxs['y'],gxs['z'])),boxsize=lbox)

voids = ascii.read('../data/tng300-1_voids.dat',names=['r','x','y','z','vx','vy','vz','deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter'])

N=[]
for i in range(len(voids)):
	if i%1000==0: print(i)
	x,y,z,r = voids['x','y','z','r'][i]
	N.append(len(tree.query_ball_point([x,y,z],r)))
	
plt.loglog(voids['r'],N,'.')
plt.xlabel('R')
plt.ylabel('N_gal')
plt.savefig('NvsR1.png')
plt.close()

plt.hist(voids['r'],bins=50)
plt.xlabel('R')
plt.yscale('log')
plt.savefig('Rvoids_hist1.png')
