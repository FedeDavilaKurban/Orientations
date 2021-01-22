#%%
import numpy as np
from scipy import spatial
from astropy.table import Table
from astropy.io import ascii
from orientationsTools import *

voids = readVoids(7.)

nvs = range(len(voids))
rmaxs = [1.4]
rmins = [0.7]

#Choose one of the following as parameter in orientations()
secs = [3]
sections = [4,5,45,56,456]
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

ran_iter = 100
s5 = 0 #Sigma5; considera sólo galaxias con sigma5 <= s5. Poner s5=0 para desactivar

"""
ECDF and residue fits
"""
#Read cosines of galaxies in shells of voids of interest

for jk in nvs:
    print('jk: {}/{}'.format(jk,len(nvs)-1))
    cos = []
    for nv in nvs:
        if jk==nv: continue
        for rmin in rmins:
            for rmax in rmaxs:
                for sec in secs:
                    cosTable = ascii.read('../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'])
                    cos.append( cosTable['cos'].data )

    cos_flattened = np.concatenate(cos).ravel()
    #print(cos_flattened)
    N = len(cos_flattened) #N/2 is the number of galaxies of interest

    #ECDF, fits
    cos,ecdf,y = ecdf_residues(cos_flattened)

    #ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_sec{}_rmin{}_rmax{}_jk{}'.format(sec,rmin,rmax,jk),names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)
    ascii.write(Table(np.column_stack([cos,ecdf(cos),y])),'../data/ecdf_sec{}_rmin{}_rmax{}_jk{}'.format(sec,rmin,rmax,jk),names=['cos','ecdf','y'],overwrite=True)

for rmin in rmins:
    for rmax in rmaxs:
        for sec in secs:
            dataList, x_ran, y_var, y_mean, y_sd = jk_mean_sd(1000,sec,rmin,rmax) #El argumento es el numero de valores random del eje x para evaluar mean y SD de las curvas JK

#Fits the mean values of the JK curves corresponding to the x_ran
yfit,d_yfit,a2 = fits(x_ran,y_mean)



#IN PROGRESS
#Random Sample

n_forsum = []
for nv in nvs:
    for rmin in rmins:
        for rmax in rmaxs:
            for sec in secs:
                cosTable = ascii.read('../data/cos_nv{}_sec{}_rmin{}_rmax{}_s5:{}'.format(nv,sec,rmin,rmax,s5),names=['cos'])
                n_forsum.append( len(cosTable['cos'].data) )

N_forRandom = np.sum(n_forsum) #Equals to the number of gxs in all the voids considered 
xmean, ymean, ystd = ranOrientations(ran_iter,N_forRandom)
ascii.write(Table(np.column_stack([xmean,ymean,ystd])),'../data/randomfits_sec{}_rmin{}_rmax{}'.format(sec,rmin,rmax),names=['xmean','ymean','ystd'],overwrite=True)

#PLOTS
import matplotlib.pyplot as plt
for i in range(len(dataList)):
    plt.plot(dataList[i]['cos'],dataList[i]['y'],alpha=.01,color='k')

plt.fill_between(x_ran,y_mean-y_sd,y_mean+y_sd,alpha=.8,color='r')
plt.plot(x_ran,yfit)

plt.fill_between(xmean,ymean-ystd,ymean+ystd,color='k',alpha=.6,label=r'$1\sigma$')
plt.fill_between(xmean,ymean-2*ystd,ymean+2*ystd,color='k',alpha=.4,label=r'$2\sigma$')
plt.legend()
plt.show()
# %%
