import numpy as np
from orientationsTools import *
from astropy.io import ascii
from astropy.table import Table

N=1000

ones = np.ones((1,N))
zeros = np.zeros((1,N))

cos,y,ecdf,yfit,d_yfit,a2 = fits(ones.tolist())
print(a2)
ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_ones',names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)


cos,y,ecdf,yfit,d_yfit,a2 = fits(zeros.tolist())
print(a2)
ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_zeros',names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)
