def orientations(rmins,rmaxs,sections,s5,ran_iter,write=True,ranfits=True):
    """
    Calculates orientations of disks, their ecdf, and fits.
    
    Writes files "fits_sec{}_rmin{}_rmax{}".format(sec,rmin,rmax)
    and
    "randomfits_sec{}_rmin{}_rmax{}".format(sec,rmin,rmax)

    """
    import numpy as np
    from astropy.table import Table

    global xv,yv,zv,rv, N

    for rmin in rmins:
        print(rmin)
        for rmax in rmaxs:
            print(rmax)
            for sec in sections:
                cos = []
                print('sec =',sec)
                for nv in range(len(voids)):
                    #if nv%10==0: print('nv=',nv)
                    xv,yv,zv,rv = voids['x','y','z','r'][nv]
                    if units=='kpc':
                        xv*=1000.
                        yv*=1000.
                        zv*=1000.
                        rv*=1000.

                    idx1 = tree.query_ball_point([xv,yv,zv],rv*rmax)
                    idx2 = tree.query_ball_point([xv,yv,zv],rv*rmin)
                    shell = [g for g in idx1 if g not in idx2]
                    gals = gxs[shell]
                    #print('N of galaxies in shell:',len(gals))

                    """
                    Spin-Mass linear regression 
                    Determine section of interest with 'sec'
                    """
                    gals_h = JvsM(sec,gals,plot=False)
                    #N_gals.append(len(gals_h))

                    #print('N of galaxies of interest:',len(gals_h))
                    """
                    Cosines
                    """
                    cos.append( cosines(gals_h,units,s5) )
                    #print(np.shape(cos))

                cos_flattened = [y for x in cos for y in x]
                N = len(cos_flattened)
                #print(N) #N/2 seria el numero de galaxias estudiadas

                #ECDF, fits
                cos,y,ecdf,yfit,d_yfit,a2 = fits(cos_flattened)
                print(a2)
                a2string="{:.3f}".format(a2)
                if write: ascii.write(Table(np.column_stack([cos,ecdf(cos),y,yfit,d_yfit])),'../data/fits_sec{}_rmin{}_rmax{}_s5:{}'.format(sec,rmin,rmax,s5),names=['cos','ecdf','y','yfit','d_yfit'],overwrite=True)

                #Random Sample
                if ranfits:
                    xmean, ymean, ystd = ranOrientations(ran_iter)
                    if write: ascii.write(Table(np.column_stack([xmean,ymean,ystd])),'../data/randomfits_sec{}_rmin{}_rmax{}_s5:{}'.format(sec,rmin,rmax,s5),names=['xmean','ymean','ystd'],overwrite=True)


rmaxs = [1.2]
rmins = [0.7]

sec = [6]
sections = [4,5,45,56,456]
allsections = [1,2,3,12,23,123,4,5,6,45,56,456]

ran_iter = 100
s5 = 0 #Sigma5; considera sólo galaxias con sigma5 <= s5. Poner s5=0 para desactivar

for s5 in []:
    orientations(rmins,rmaxs,sections,s5,ran_iter,write=True,ranfits=True)