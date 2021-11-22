#%%
"""
This is an experiment to interpret the meaning of the eta parameter
in terms of the angles
"""
import sys
from config import illustrisPath
sys.path.append(illustrisPath+'param_tools')
# use param_tools: https://github.com/maxkapur/param_tools
from param_tools import r_surface, surface_area

import numpy as np
from matplotlib import pyplot as plt
from functools import partial

#%%

# GENERATE POINTS UNIFORMLY DISTRIBUTED ON THE SURFACE OF AN ELLIPSE
# + VISUALIZE

def ellipsoid(t, u, a=1, b=1, c=1.1):
    return np.array([a*np.sin(u)*np.cos(t), b*np.sin(u)*np.sin(t), c*np.cos(u)])

domain_t = [0, 2*np.pi]
domain_u = [0, np.pi]

# Get random points
x, t, u, St, Su = r_surface(2000, ellipsoid, *domain_t, *domain_u, 20, 20)


# Set up plot
fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(hspace=0, wspace=0.1, bottom=.2)
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-1.3, 1.3)
ax.set_zlim(-1.3, 1.3)

# Plot random points
ax.scatter(*x, marker='o', alpha=0.3, color='cadetblue')

# Plot function
shape_t, shape_u = np.meshgrid(np.linspace(*domain_t, 25), np.linspace(*domain_u, 25))
shape_x = ellipsoid(shape_t, shape_u)
ax.plot_wireframe(*shape_x, color='black', ls="--", alpha = 0.3)

xs = x[0,:]
ys = x[1,:]
zs = x[2,:]

ax.set_box_aspect((np.ptp(xs), np.ptp(ys), np.ptp(zs))) # aspect ratio is 1:1:1 in data space

# Surface area
ellipsoid_A = surface_area(ellipsoid, *domain_t, *domain_u, 100, 100)
ax.set_title('Total surface area: {}'.format(np.round(ellipsoid_A, 3)), size=15, y=-0.05)

plt.show()
#%%

# COMPUTE THE BETA ALIGNMENT STATISTICS FROM THE POINTS
# + VISUALIZE

beta = []
for i in range(len(xs)):
    para = zs[i]
    perp = np.sqrt(xs[i]**2+ys[i]**2)
    beta.append(perp / abs(para))

fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(hspace=0, wspace=0.1, bottom=.2)
ax = fig.add_subplot(111)
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-1.3, 1.3)

ax.hist(np.log10(beta))
plt.show()
#%%

# COMPUTE THE ETA ALIGNMENT STATISTICS FROM THE POINTS
# + VISUALIZE
 
fig = plt.figure(figsize=(10,10))
fig.subplots_adjust(hspace=0, wspace=0.1, bottom=.2)
ax = fig.add_subplot(111)
bins = np.linspace(1.9, 4.8, 50)

clrs = ['#402218', '#865439', '#C68B59', '#D7B19D']

Nran = 200*5
Netas = 300*5

for k, c in enumerate(np.arange(0.7, 1, 0.1)):
    print(c)

    f = partial(ellipsoid, c=c)

    eta = np.zeros(Netas)
    for i in range(Netas):
        x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
        xs = x[0,:]
        ys = x[1,:]
        zs = x[2,:]
         
        beta = np.zeros(Nran)
        for j in range(len(xs)):
            para = zs[j]
            perp = np.sqrt(xs[j]**2+ys[j]**2)
            beta[j] = perp / abs(para)
        eta[i] = sum(beta>1)/sum(beta<1)

    ax.hist(eta, bins=bins, histtype='step', linewidth=2, 
            color=clrs[k], density=True, label=f'{c:.1f}')
    ax.axvline(eta.mean(), color=clrs[k], linestyle=':')

    print(eta.mean(), eta.mean()-1/(np.sqrt(2)-1) )

ax.axvline(2.41, color='slategrey', linestyle='--')
ax.legend()

#fig.savefig('signal_etas.pdf')
plt.show()

#%%
###########


# POLAR PLOT OF RANDOMS

Nran = 100000
Ntheta = 12
bins = np.arccos(np.linspace(0, 1, Ntheta))
bins = np.concatenate((-bins, bins[(Ntheta-2)::-1]))
bins = np.concatenate((bins, bins+np.pi))


fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

for c in [0.7, 0.8, 0.9, 1.]:

    f = partial(ellipsoid, c=c)

    x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
    xs = x[0,:]
    ys = x[1,:]
    zs = x[2,:]   

    beta = np.zeros(Nran)
    for j in range(len(xs)):
        para = zs[j]
        perp = np.sqrt(xs[j]**2 + ys[j]**2)
        beta[j] = perp / abs(para)

    beta = np.concatenate((-beta, beta))
    lmbda = np.arctan(beta)
    lmbda = np.concatenate((lmbda, lmbda+np.pi))

    h, t = np.histogram(lmbda, bins=bins)
    t = (t[1:]+t[:(len(t)-1)])/2

    ax.plot(t, h)

ax.grid(True)
plt.show()
#%%

# POLAR PLOT THEORETICAL

def ellips(t, a, b):
    return a*np.cos(t), b*np.sin(t)

def ellipsp(t, a, b):
    e2 = 1 - (b/a)**2
    r = np.sqrt(b**2 / (1-e2*np.cos(t)**2))
    return r

t = np.linspace(0, 2*np.pi, 100)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})

for k, c in enumerate([0.7, 0.8, 0.9, 1.]):

    r = ellipsp(t, 1, c)
    ax.plot(t, r, clrs[k])

ax.grid(True)
#fig.savefig('ellipses.pdf')           
plt.show()
#%%

# COS PLOT WITH RANDOMS

Nran = 100000
Ntheta = 12
bins = np.arccos(np.linspace(0, 1, Ntheta))
bins = np.concatenate((-bins, bins[(Ntheta-2)::-1]))
bins = np.concatenate((bins, bins+np.pi))
t = np.linspace(0, 2*np.pi, 100)

fig = plt.figure(figsize=(8, 3))
ax = fig.add_subplot(111)

for k, c in enumerate([0.7, 0.8, 0.9, 1.]):

    f = partial(ellipsoid, c=c)
    x, t, u, St, Su = r_surface(Nran, f, *domain_t, *domain_u, 20, 20)
    xs = x[0,:]
    ys = x[1,:]
    zs = x[2,:]   

    beta = np.zeros(Nran)
    for j in range(len(xs)):
        para = zs[j]
        perp = np.sqrt(xs[j]**2 + ys[j]**2)
        beta[j] = perp / abs(para)

    beta = np.concatenate((-beta, beta))
    lmbda = np.arctan(beta)
    lmbda = np.concatenate((lmbda, lmbda+np.pi))

    ax.hist(np.cos(lmbda), histtype='step', color=clrs[k])
#fig.savefig('cosenos.pdf')
plt.show()

# %%
