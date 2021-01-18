
#%%
"""
Quiero plotear los residuos de una sola seccion (generalmente va a ser la 3)
con los errores jackknife

Eventualmente añadire los errores sigma 1 y 2.
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
#%%
sec = 3 
rmin, rmax = .7, 1.4

dataList=[]
for jk in range(81):
    dataList.append( ascii.read('../data/fits_sec{}_rmin{}_rmax{}_jk{}'.format(sec,rmin,rmax,jk),names=['cos','ecdf','y','yfit','d_yfit']) )
    plt.plot(dataList[-1]['cos'],dataList[-1]['y'],alpha=.01,color='k')

#plt.plot(dataList[-1]['cos'],dataList[-1]['y'],alpha=1.)

#I need to create an array of x values where I will 
#evaluate the variance of the JK resampling
x_Var = np.linspace(0,1,100)
y_Var = []
y_Mean = []
for j in range(len(x_Var)):
    yList=[]
    for i in range(len(dataList)):
        x,y = dataList[i]['cos'], dataList[i]['y']
        yList.append( np.interp(x_Var[j],x,y) )

    y_Var.append( np.var(yList,ddof=1) )
    y_Mean.append( np.mean(yList) )

#Get Standard Deviation from Variance
y_sd = np.sqrt(y_Var)

#Actual plot
plt.fill_between(x_Var,y_Mean-y_sd,y_Mean+y_sd,alpha=.8,color='r')
plt.fill_between(-1*np.array(x_Var),-1*np.array(y_Mean)-y_sd,-1*np.array(y_Mean)+y_sd,alpha=.8,color='r')

plt.show()

# ax.set_title('Sec={}'.format(sec))
# #ax.step(cos,y)
# ax.plot(cos,yfit,'r--')
# ax.fill_between(xmean,ymean-ystd,ymean+ystd,color='k',alpha=.6,label=r'$1\sigma$')
# ax.fill_between(xmean,ymean-2*ystd,ymean+2*ystd,color='k',alpha=.4,label=r'$2\sigma$')
# #ax.legend()

# axes[0,0].legend(loc=2)
# plt.tight_layout()
# plt.savefig('../plots/cosines_allsections_rmin{}_rmax{}.png'.format(rmin,rmax),dpi=300)
# plt.show()
# %%
