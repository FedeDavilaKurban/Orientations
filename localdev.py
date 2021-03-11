
# lowpass filter density estimator
# 

import numpy as np

def ref_cdf(x, a, b):
    ref_cdf = (x-a)/(b-a)
    return ref_cdf



def ecdf_residues(x, a, b, verbose=False):
    """
    Calculates ECDF and residues. Fits the residues.
    Returns: ecdf, fit, derivative of fit, and coefficient a2
    """
    from statsmodels.distributions.empirical_distribution import ECDF

    x = np.array(x)
    x.sort()
    x_ecdf = ECDF(x)(x)

    refcdf = ref_cdf(x, a, b)

    residues = x_ecdf - refcdf

    return x, x_ecdf, residues


def func(x,a,b,c,d,e):
    f = a + \
        b*np.sin( np.pi*(x+1.)/2. ) + \
        c*np.sin( 2.*np.pi*(x+1.)/2. ) + \
        d*np.sin( 3.*np.pi*(x+1.)/2. ) + \
        e*np.sin( 4.*np.pi*(x+1.)/2. )
    return f

def dfunc(x,b,c,d,e):
    f = np.pi/2.*b*np.cos( np.pi*(x+1.)/2. ) + \
        np.pi*c*np.cos( 2.*np.pi*(x+1.)/2. ) + \
        np.pi*3./2.*d*np.cos( 3.*np.pi*(x+1.)/2. ) + \
        np.pi*2.*e*np.cos( 4.*np.pi*(x+1.)/2. )
    return f

 
def fits(x, y, verbose=False):
    "Perform fits on the residues of the ecdf"
    
    import numpy as np
    from scipy.optimize import curve_fit

    coeffs, cov = curve_fit(func, x, y)
    yfit = func(x,coeffs[0],coeffs[1],coeffs[2],coeffs[3],coeffs[4])
    d_yfit = dfunc(x,coeffs[1],coeffs[2],coeffs[3],coeffs[4])

    a1 = coeffs[1]
    a2 = coeffs[2]
    a3 = coeffs[3]
    a4 = coeffs[4]
    if verbose==True: print('a1=',a1,'; a2=',a2,'; a3=',a3,'; a4=',a4)

    return yfit, d_yfit, coeffs
 

