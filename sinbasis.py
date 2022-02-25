#%%
import numpy as np
from matplotlib import pyplot as plt


t = np.linspace(0, np.pi/2, 200)

def f(t, n):
    return np.sin(2*n*t)


plt.close('all')
fig = plt.figure(figsize=(10, 5))
ax1 = fig.add_subplot(1,2,1)

ls = ['-','--','-.']
label = [r'$\phi_0$',r'$\phi_2$',r'$\phi_4$']
for i,n in enumerate(range(0,6,2)):
    ax1.plot(t, f(t,n), ls=ls[i], label=label[i])

ax2 = fig.add_subplot(1,2,2)

label = [r'$\phi_1$',r'$\phi_3$',r'$\phi_5$']
for i,n in enumerate(range(1,6,2)):
    ax2.plot(t, f(t,n), ls=ls[i], label=label[i])

ax1.legend(fontsize=18,framealpha=.5)
ax2.legend(fontsize=18,framealpha=.5)
ax1.set_xlabel('x', fontsize=16)
ax1.set_ylabel(r'$\phi_k(x)$', fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.tick_params(axis='both', which='minor', labelsize=16)
ax2.set_xlabel('x', fontsize=16)
ax2.set_ylabel(r'$\phi_k(x)$', fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.tick_params(axis='both', which='minor', labelsize=16)

plt.tight_layout()

fig.savefig('../plots/sinbasis.pdf')
plt.show()
# %%
