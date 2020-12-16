#%%
import numpy as np
from orientationsTools import *
from astropy.io import ascii
from astropy.table import Table
import matplotlib.pyplot as plt

N=1000
np.random.seed(123)

cos = np.array([np.random.normal(.0,.3) for i in range(10000)])
cos0 = [c for c in cos if c>=0. and c<=1.]
cos = np.array([np.random.normal(1.,.3) for i in range(10000)])
cos1 = [c for c in cos if c>=0. and c<=1.]


plt.hist(cos0,color='b',alpha=.8)
plt.hist(cos1,color='g',alpha=.8)
plt.show()

cos,y,ecdf,yfit,d_yfit,a2 = fits(cos0)
print('a2=',a2)
ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_ones',names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)
plt.step(cos,y,'b-')
plt.plot(cos,yfit,'k--')

cos,y,ecdf,yfit,d_yfit,a2 = fits(cos1)
print('a2=',a2)
ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_zeros',names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)
plt.step(cos,y,'g-')
plt.plot(cos,yfit,'k--')
plt.show()
# %%

