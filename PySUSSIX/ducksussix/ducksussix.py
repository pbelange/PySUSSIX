



import numpy as np
import pandas as pd

import PySUSSIX.crosssussix.crossroutines as crossroutines
import PySUSSIX.crosssussix.crosssussix as crosssussix



def calcr(tune_phasor, z):
    return np.polyval(z[::-1], tune_phasor)

def spectrum(x,px,number_of_harmonics = 5,method = 'hanning'):
    """
    COMPUTE THE MAIN FREQUENCY OF A CANONICAL PAIR 
    (WITHOUT ORTHOGONALIZATION OF GRAM-SCHMIDT)
    """

    assert number_of_harmonics >=1, 'number_of_harmonics needs to be > 1'
    
    # initialisation
    z  = x + 1j*px
    N  = np.arange(len(x))
    
    
    frequencies = []
    amplitudes  = [] 
    for _ in range(number_of_harmonics):

        # Computing frequency and amplitude
        freq,zw = crossroutines.tunenewt(x,px)
        zpesi   = zw / max(N+1)

        # Saving results
        frequencies.append(freq)
        amplitudes.append(zpesi)

        # Substraction procedure
        zgs  = zpesi * np.exp(2 * np.pi * 1j * freq * N)
        z   -= zgs
        x,px = np.real(z), np.imag(z)

    
    return pd.DataFrame({'amplitude':amplitudes,'frequency':frequencies})

    


def datspe(x = None,px = None,y = None,py = None,zeta = None,pzeta = None,number_of_harmonics = 5,method = 'hanning'):
    """
    THIS PROGRAM CALCULATES THE SPECTRUM OF A TRACKING DATA SET,
    THE TUNE AND THE LINES ARE CALCULATED WITH THE ROUTINE SPECTRUM
    """

    results = {}
    for pair in [(x,px,'x'),(y,py,'y'),(zeta,pzeta,'zeta')]:
        z,pz,plane = pair
        
        if z is not None:
            if pz is None:
                pz = np.zeros(len(z))

            # Computing spectrum
            # df = crossroutines.spectrum(z,pz,   number_of_harmonics = number_of_harmonics,
            #                                     method              = method,
            #                                     return_values       = True)
            # df.rename(columns={'tune':'frequency','zpesi':'amplitude'},inplace=True)
            # df = df[['amplitude','frequency']]
            df = spectrum(z,pz, number_of_harmonics = number_of_harmonics,
                                method              = method)
            
            results[plane] = df

        else:
            results[plane] = None


    return results


    