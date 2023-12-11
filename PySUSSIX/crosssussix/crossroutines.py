import numpy as np
import pandas as pd

import PySUSSIX.f90sussix.f90sussix as f90sussix
import PySUSSIX.crosssussix.crossref as crossref





def spectrum(x,px,number_of_harmonics = 5,method = 'hanning',return_values = False):
    """
    COMPUTE THE MAIN FREQUENCY OF A CANONICAL PAIR 
    (WITHOUT ORTHOGONALIZATION OF GRAM-SCHMIDT)
    """
    # Initialize arrays
    size_tbt    = 100000
    size_mterms = 300 

    _x    = np.zeros(size_tbt)
    _px   = np.zeros(size_tbt)
    tune  = np.zeros(size_mterms)
    zpesi = np.zeros(size_mterms) + 1j*np.zeros(size_mterms)

    assert number_of_harmonics >=2, 'number_of_harmonics needs to be > 2'
    narm  = number_of_harmonics
    meth  = {'hanning':1,'rectangular':2}[method]

    # Fill arrays with data
    maxn = len(x)
    _x[:maxn]  = x
    _px[:maxn] = px

    # Run fortran code
    f90sussix.spectrum(_x,_px,maxn,tune,zpesi,narm,meth)

    if return_values:
        return pd.DataFrame({'tune':tune[:narm],'zpesi':zpesi[:narm]})



def tunenewt(x,px):
    """
    COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
    IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
    """
    # Initialize arrays
    size_tbt    = 100000

    # zw: changed with funr and calcr
    _x    = np.zeros(size_tbt)
    _px   = np.zeros(size_tbt)

    # Fill arrays with data
    maxn = len(x)
    _x[:maxn]  = x
    _px[:maxn] = px

    # Run fortran code
    zw   = 0*1j
    tune = f90sussix.tunenewt(_x,_px,maxn,zw)
    zw   = f90sussix.data.zw_out[0]
    return tune,zw

def tunelasr(x,px):
    """
    SAME AS TUNENEWT BUT NO HANNING FILTER
    """
    # Initialize arrays
    size_tbt    = 100000

    # zw: ? not understood...
    _x    = np.zeros(size_tbt)
    _px   = np.zeros(size_tbt)

    # Fill arrays with data
    maxn = len(x)
    _x[:maxn]  = x
    _px[:maxn] = px

    # Run fortran code
    zw   = 0*1j
    tune = f90sussix.tunelasr(_x,_px,maxn,zw)
    zw   = f90sussix.data.zw_out[0]
    return tune,zw



def zfunr(z,tune0):
    size_tbt    = 100000

    _z     = np.zeros(size_tbt) + 1j*np.zeros(size_tbt)
    maxn = len(z)
    tune   = 0
    zw     = complex(0,0)
    tunea1 = tune0
    deltat = 1.0 / maxn
    _z[:maxn] = z

    npoint = 2**int(np.log2(maxn))
    deltat = 1.0 /npoint

    
    # Run fortran code
    f90sussix.zfunr(tune,zw,_z,maxn,tunea1,deltat)
    tune = f90sussix.zfunr_out.tune_out[0]
    zw   = f90sussix.zfunr_out.zw_out[0]
    return tune,zw


def calcr(z,tune):
    
    ztune1 = np.exp(-1j * 2 * np.pi * tune)
    zf     = z[-1]*0 + 1
    maxn   = len(z)


    f90sussix.calcr(ztune1,zf,z,maxn)

    return f90sussix.calcr_out.zpp_out[0]


def cfft(z):
    size_tbt = 100000
    maxn     = len(z)
    mft      = int(np.log2(maxn))

    _z       = np.zeros(size_tbt) + 1j*np.zeros(size_tbt)
    _z[:maxn]= z
    
    f90sussix.cfft(_z,-mft)
    
    return f90sussix.cfft_out.a_out