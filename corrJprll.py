"""
Objetivo final: hacer un grafico con delta tita (Jprll) vs delta tita (r)
"""
#%%
import matplotlib.pyplot as plt
import numpy as np

#%%
#Testing
np.random.seed(0)
ndim = 2

vc = np.random.rand(1,ndim)[0] #Void center

r1 = np.random.rand(1,ndim)[0] #Galaxia 1
r2 = r1+0.2*np.random.rand(1,ndim)[0] #Galaxia 2

s1 = 0.2*np.random.rand(1,ndim)[0] #Spin 1
s2 = 0.2*np.random.rand(1,ndim)[0] #Spin 2

gr1 = r1-vc #Gx1 from void center
gr2 = r2-vc #Gx2 from void center


# %%

fig, ax = plt.subplots(1)

ax.arrow(0,0,vc[0],vc[1],head_width=.02)

ax.arrow(0,0,r1[0],r1[1],head_width=.02,color='g')
ax.arrow(0,0,r2[0],r2[1],head_width=.02,color='r')

ax.arrow(vc[0],vc[1],gr1[0],gr1[1],head_width=.02,color='g')
ax.arrow(vc[0],vc[1],gr2[0],gr2[1],head_width=.02,color='r')

ax.arrow(r1[0],r1[1],s1[0],s1[1],head_width=.02,color='b')
ax.arrow(r2[0],r2[1],s2[0],s2[1],head_width=.02,color='b')

ax.set_ylim((0,2))
ax.set_xlim((0,2))

# %%
