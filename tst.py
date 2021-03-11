import numpy as np
from matplotlib import pyplot as plt
from astropy import coordinates
from mpl_toolkits.mplot3d import axes3d

from orientationsTools import ecdf_residues, fits

import localdev as tools

# Generar vectores de "spin" randoms dentro de una esfera de radio 3
# En la region 1<r<2, hay un exceso hacia el centro

#---------------
N = 500
Rmax = 3

a = -1
b = 1

Rex = 1
#---------------


# random locations (N points in a sphere of radius Rmax)
theta = np.arccos(np.random.rand(N)*2-1) - np.pi/2
phi = np.random.rand(N)*2*np.pi
r = np.sqrt(np.random.rand(N)*Rmax**2)
x, y, z = coordinates.spherical_to_cartesian(r, theta, phi)

# random spin vectors with radial excess

# ..random component
rtheta = np.arccos(np.random.rand(N)*2-1) - np.pi/2
rphi = np.random.rand(N)*2*np.pi
rr = np.ones(N)
xr, yr, zr = coordinates.spherical_to_cartesian(rr, rtheta, rphi)

# ..radial excess component

shell = np.logical_and(r > 1, r < 2)
xe, ye, ze = coordinates.spherical_to_cartesian(rr, theta, phi)
xe, ye, ze = -xe, -ye, -ze


xe, ye, ze = xe*Rex, ye*Rex, ze*Rex

xe[np.logical_not(shell)] = 0
ye[np.logical_not(shell)] = 0
ze[np.logical_not(shell)] = 0

Sx, Sy, Sz = xr + xe, yr + ye, zr + ze 

# ------------------------------------------------------------------


# VIS: 3D arrows

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.quiver(x, y, z, Sx, Sy, Sz, length=0.2, arrow_length_ratio=0.2)


fig.savefig('plot1.png')
plt.close('all')


# ------------------------------------------------------------------

# Tradicional: histograma de los cosenos
 
# ..1 calcular el coseno del ángulo

cosenos = []
for i in range(N):
    coso = np.dot([x[i], y[i], z[i]], [Sx[i], Sy[i], Sz[i]])
    coso = coso / np.sqrt(np.dot([x[i], y[i], z[i]], [x[i], y[i], z[i]]))
    coso = coso / np.sqrt(np.dot([Sx[i], Sy[i], Sz[i]], [Sx[i], Sy[i], Sz[i]]))
    cosenos.append(coso)

#cosenos = np.absolute(cosenos)
#cosenos = np.concatenate(cosenos, -cosenos)

 
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

Nbins = int(np.sqrt(300))
ax.hist(cosenos, bins=Nbins)
ax.set_xlim(a, b)

fig.savefig('plot2.png')
plt.close('all')


# ------------------------------------------------------------------

# Ahora aplicamos el procedimiento


# ..2 


cos, ecdf, res = tools.ecdf_residues(cosenos, a, b)

yfit, d_yfit, coefs = tools.fits(cos, res)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

cos_grid = np.linspace(a, b, len(yfit))

cos_grid_cdf = tools.ref_cdf(cos_grid, a, b)

ax.plot(cos, ecdf)
ax.plot(cos_grid, cos_grid_cdf)
ax.set_xlim(a, b)
ax.grid()

fig.savefig('plot3.png')
plt.close('all')     

 

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(cos, res)
ax.set_xlim(a, b)
ax.grid()

fig.savefig('plot4.png')
plt.close('all')     
                             


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(cos, yfit)
ax.scatter(cos, res)
ax.set_xlim(a, b)
ax.grid()

fig.savefig('plot5.png')
plt.close('all')     

        

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

ax.plot(cos_grid, d_yfit)
ax.set_xlim(a, b)
ax.grid()

fig.savefig('plot6.png')
plt.close('all')     
 

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
 
ax.plot(cos, ecdf)
ax.plot(cos_grid, cos_grid_cdf)
ax.set_xlim(a, b)
ax.grid()                                 


refcdf = tools.ref_cdf(cos, a, b)

y = yfit + refcdf
ax.plot(cos, y)


fig.savefig('plot7.png')
plt.close('all')     

 
                         
