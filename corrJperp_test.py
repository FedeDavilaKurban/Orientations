"""
Objetivo final: hacer un grafico con delta tita (Jperp) vs delta tita (r)
"""
#%%
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

#%%
#Testing
np.random.seed(0)
ndim = 3

vc = np.random.rand(1,ndim)[0] #Void center

r1 = np.random.rand(1,ndim)[0] #Galaxia 1
r2 = r1+0.2*np.random.rand(1,ndim)[0] #Galaxia 2

s1 = 0.2*np.random.rand(1,ndim)[0] #Spin 1
s2 = 0.2*np.random.rand(1,ndim)[0] #Spin 2

gr1 = r1-vc #Gx1 from void center
gr2 = r2-vc #Gx2 from void center

###################################

n1 = gr1/np.linalg.norm(gr1)
n2 = gr2/np.linalg.norm(gr2)

s1_prll = np.dot(s1,n1)*n1 
s2_prll = np.dot(s2,n2)*n2 

s1_perp = s1 - np.dot(s1,n1)*n1 
s2_perp = s2 - np.dot(s2,n2)*n2 



# %%

tita_r = 180*np.arccos(np.dot(n1,n2))/np.pi
tita_s = 180*np.arccos((np.dot(s1_perp,s2_perp))/(np.linalg.norm(s1_perp)*(np.linalg.norm(s2_perp))))/np.pi

print(tita_r,tita_s)

# %%

def get_random_vec(R, N):

    phi = 2*np.pi*np.random.random(N)
    costheta = 1-2*np.random.random(N)
    u = np.random.random(N)

    theta = np.arccos( costheta )
    r = R * np.cbrt( u )

    x = r * np.sin( theta ) * np.cos( phi )
    y = r * np.sin( theta ) * np.sin( phi )
    z = r * np.cos( theta )
    
    return x,y,z

x,y,z = get_random_vec(1,1000)

tita = np.pi/2-np.arccos(z/np.sqrt(x**2+y**2+z**2))
phi = np.arctan2(y,x)

plt.scatter(180*phi/np.pi,180*tita/np.pi)
#%%
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='mollweide')
ax.scatter(phi,tita)
plt.show()
# %%
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='mollweide')
ax.quiver(phi,tita,.0,0.03,scale=1)
plt.show()
# %%
