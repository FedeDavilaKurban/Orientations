"""
Objetivo final: hacer un grafico con delta tita (Jperp) vs delta tita (r)
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

n1 = gr1/np.linalg.norm(gr1)
n2 = gr2/np.linalg.norm(gr2)

s1_prll = np.dot(s1,n1)*n1 
s2_prll = np.dot(s2,n2)*n2 

s1_perp = s1 - np.dot(s1,n1)*n1 
s2_perp = s2 - np.dot(s2,n2)*n2 


# %%

fig, ax = plt.subplots(1)

ax.arrow(0,0,vc[0],vc[1],head_width=.02,alpha=.6)

ax.arrow(0,0,r1[0],r1[1],head_width=.02,color='g',alpha=.6)
ax.arrow(0,0,r2[0],r2[1],head_width=.02,color='g',alpha=.6)

ax.arrow(vc[0],vc[1],gr1[0],gr1[1],head_width=.02,color='g',alpha=.6)
ax.arrow(vc[0],vc[1],gr2[0],gr2[1],head_width=.02,color='g',alpha=.6)

ax.arrow(r1[0],r1[1],s1[0],s1[1],head_width=.02,color='b')
ax.arrow(r2[0],r2[1],s2[0],s2[1],head_width=.02,color='r')

ax.arrow(r1[0],r1[1],s1_prll[0],s1_prll[1],head_width=.02,color='k')
ax.arrow(r2[0],r2[1],s2_prll[0],s2_prll[1],head_width=.02,color='k')

ax.arrow(r1[0],r1[1],s1_perp[0],s1_perp[1],head_width=.02,color='k')
ax.arrow(r2[0],r2[1],s2_perp[0],s2_perp[1],head_width=.02,color='k')

ax.arrow(0,0,s1_perp[0],s1_perp[1],head_width=.02,color='k')
ax.arrow(0,0,s2_perp[0],s2_perp[1],head_width=.02,color='k')


ax.set_ylim((0,1))
ax.set_xlim((0,1))


# %%

tita_r = 180*np.arccos(np.dot(n1,n2))/np.pi
tita_s = 180*np.arccos(np.dot(s1_perp,s1_perp))/np.pi

print(tita_r,tita_s)

# %%
