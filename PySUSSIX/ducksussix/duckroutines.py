import numpy as np
import pandas as pd

import PySUSSIX.f90sussix.f90sussix as f90sussix
import PySUSSIX.crossref as crossref





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
    zpesi = np.zeros(size_mterms)

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



def tunenewt(x,px,zw = 0*1j):
    """
    COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
    IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY.
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
    tune = f90sussix.tunenewt(_x,_px,maxn,zw)
    
    return tune

def tunelasr(x,px,zw = 0*1j):
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
    tune = f90sussix.tunelasr(_x,_px,maxn,zw)
    
    return tune