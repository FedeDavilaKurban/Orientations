
#%%
def readVoids(minrad=None,maxrad=None):
    """
    This function reads the file tng300-1_voids.dat
    'minrad': filter by a minimum radius (Mpc). Optional.
    'maxrad': filter by maximum radius (Mpc). Optional
    """
    from astropy.io import ascii
    import numpy as np

    voids=ascii.read('../data/tng300-1_voids.dat',names=['r','x','y','z','vx','vy','vz','deltaint_1r','maxdeltaint_2-3r','log10Poisson','Nrecenter'])#,usecols=['r','x','y','z'])

    if minrad != None: voids=voids[np.where(voids['r']>=minrad)]
    if maxrad != None: voids=voids[np.where(voids['r']<=maxrad)]
    
    return voids

def readTNG():
    """
    This function reads subhalos in the TNG300-1 simulation and returns 

    gxs: an ascii Table with the fields and filters I usually need for this: Position, Mass, Spin

    """
    import sys
    sys.path.append('/home/fede/')
    import illustris_python as il
    import numpy as np
    from astropy.table import Table

    basePath = '/home/fede/TNG300-1/output/'

    mass = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloMass'])                                                                                                                      
    ids = np.where((np.log10(mass)>-1.)&(np.log10(mass)<3.))
    mass=mass[ids]

    pos = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloPos'])
    pos = pos[ids]

    spin = il.groupcat.loadSubhalos(basePath,99,fields=['SubhaloSpin'])
    spin = spin[ids]
    sp_n = np.sum(np.abs(spin)**2,axis=-1)**(1./2)
    gxs = Table(np.column_stack([pos[:,0],pos[:,1],pos[:,2],mass,spin,sp_n]),names=['x','y','z','mass','spx','spy','spz','sp_n'])    
    del mass,pos,spin,sp_n

    return gxs

def JvsM(sec,gals,plot=True):
    """
    Performs linear regression of the Mass-Spin relation
    Determines galaxies of interest with 'sec'
    Returns galaxies of interest in 'gals_h'
    """
    import numpy as np
    import scipy.stats
    import matplotlib.pyplot as plt

    #g_m,g_b,g_rvalue,g_pvalue,g_std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # Este ajuste es de las galaxias en total (ajuste global "g_...")
    m,b,rvalue,pvalue,std=scipy.stats.linregress(np.log10(gxs['mass']),np.log10(gxs['sp_n'])) # El ajuste tiene que ser con las 'gxs' (no con las 'gals')
    m1 = -1./m
    b1 = -.7*(m-m1)+b
    b2 = -.3*(m-m1)+b

    M = np.log10(gals['mass'])
    S = np.log10(gals['sp_n'])

    #Alto Spin
    if sec == 1: gals_h = gals[(M < -.7)&(S > M*m+b)] # Solo limite superior
    if sec == 2: gals_h = gals[(M > -.7)&(M < -.3)&(S > M*m+b)] # Limite superior e inferior en el espacio Spin-Masa
    if sec == 3: gals_h = gals[(M > -.3)&(S > M*m+b)] # Solo limite inferior
    if sec == 12: gals_h = gals[(M < -.3)&(S > M*m+b)] # Solo limite superior
    if sec == 23: gals_h = gals[(M < -.7)&(S > M*m+b)] # Solo limite inferior
    if sec == 123: gals_h = gals[(S > M*m+b)] # Solo limite en Spin

    #Bajo Spin
    if sec == 4: gals_h = gals[(M < -.7)&(S < M*m+b)] # Solo limite superior
    if sec == 5: gals_h = gals[(M > -.7)&(M < -.3)&(S < M*m+b)] # Limite superior e inferior en el espacio Spin-Masa
    if sec == 6: gals_h = gals[(M > -.3)&(S < M*m+b)] # Solo limite inferior
    if sec == 45: gals_h = gals[(M < -.3)&(S < M*m+b)] # Solo limite superior
    if sec == 56: gals_h = gals[(M < -.7)&(S < M*m+b)] # Solo limite inferior
    if sec == 456: gals_h = gals[(S < M*m+b)] # Solo limite en Spin

    if plot:
        plt.scatter(M,S,c='gray', label='Galaxies Not Being Analyzed')
        plt.scatter(np.log10(gals_h['mass']),np.log10(gals_h['sp_n']),c='blue', label='Galaxies Being Analyzed')
        plt.plot(M,M*m+b,ls=':', label = 'Void Galaxies Spin-Mass Linear Regression')
        #plt.plot(M,M*g_m+g_b,ls='--', label = 'Box Galaxies Low/High Spin Linear Regression')
        plt.xlabel('log10( M/Msun )')
        plt.ylabel('Spin')
        plt.legend(loc=4)
        plt.savefig('../plots/JvsM_void{}.png'.format(nv),dpi=300)

    return gals_h

def cosines(gals_h,units,s5=0):
    """
    Calculates absolute value of cosine of the angle between spin vector J 
    and void-centric direction.
    
    Returns: array of the cosines (will be double the amount of gals_h because it's duplicated with negative signs)
    """
    import numpy as np

    cos_list = []
    #sigma5 = []

    for i in range(len(gals_h)):
        if s5!=0:
            #Cuento cantidad de vecinos a un radio de 5Mpc, si es mayor a s5, rompo el loop
            if units=='kpc': rad = 5000.
            if units=='Mpc': rad = 5.

            sigma5 = len( tree.query_ball_point([gals_h[i]['x'],gals_h[i]['y'],gals_h[i]['z']],rad))
            #print(sigma5)
            if sigma5>=s5: continue
            #print("sigue")

        #Calculo el coseno
        u = [gals_h[i]['x']-xv,gals_h[i]['y']-yv,gals_h[i]['z']-zv]
        v = [gals_h[i]['spx'],gals_h[i]['spy'],gals_h[i]['spz']]
        c = abs( np.dot(u,v)/(np.linalg.norm(u)*np.linalg.norm(v)) )
        cos_list.append( c )
        #cos_list.append( -c )
                
    return cos_list

def fits(cos_list,verbose=False):
    """
    Calculates ECDF and residues. Fits the residues.
    Returns: ecdf, fit, derivative of fit, and coefficient a2
    """
    import numpy as np
    from statsmodels.distributions.empirical_distribution import ECDF
    from scipy.optimize import curve_fit

    cos_neg=-1*np.array(cos_list)
    cos1 = np.sort(np.concatenate([cos_list,cos_neg]))
    ecdf = ECDF(cos1) # Empirical cumulated distribution function
    y = ecdf(cos1)-(cos1+1.)/2. # Para cuando tomamos el valor verdadero de cos (no el absoluto) 

    def func(x,a,b,c,d,e):
        return a + b*np.sin( np.pi*(x+1.)/2. ) + c*np.sin( 2.*np.pi*(x+1.)/2. ) + d*np.sin( 3.*np.pi*(x+1.)/2. ) + e*np.sin( 4.*np.pi*(x+1.)/2. )
    def dfunc(x,b,c,d,e):
        return np.pi/2.*b*np.cos( np.pi*(x+1.)/2. ) + np.pi*c*np.cos( 2.*np.pi*(x+1.)/2. ) + np.pi*3./2.*d*np.cos( 3.*np.pi*(x+1.)/2. ) + np.pi*2.*e*np.cos( 4.*np.pi*(x+1.)/2. )


    x = np.array(cos1)

    coeffs, cov = curve_fit(func, x, y)
    yfit = func(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
    d_yfit = dfunc(x,coeffs[1],coeffs[2],coeffs[3],coeffs[4])

    a1 = coeffs[1]
    a2 = coeffs[2]
    a3 = coeffs[3]
    a4 = coeffs[4]
    if verbose==True: print('a1=',a1,'; a2=',a2,'; a3=',a3,'; a4=',a4)

    return cos1,y,ecdf,yfit,d_yfit,a2

def ranOrientations(n_iter):
    """
    Generate random orientations, calculate mean and standard deviations
    Returns: xmean (the random "cosines") , ymean (mean value of the fits), ystd (stand dev of fits)

    iter: integer number of iterations
    """
    import numpy as np

    #a2_ran = []
    yran_fit = []
    xran_fit = []
    for _ in range(n_iter):
        np.random.seed(_)
        rancos_pos = np.random.uniform(0.,1.,int(N/2))
        rancos_neg = rancos_pos*-1
        rancos = np.concatenate((rancos_pos,rancos_neg))

        cos_,y_,ecdf_,yfit_,d_yfit_,a2_ = fits(rancos)
        #a2_ran.append(a2_)
        yran_fit.append(yfit_)
        xran_fit.append(cos_)

    xmean = np.mean(xran_fit,axis=0)
    ymean = np.mean(yran_fit,axis=0)
    ystd = np.std(yran_fit,axis=0)

    return xmean,ymean,ystd

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
