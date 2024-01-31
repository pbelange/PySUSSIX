



import numpy as np
import pandas as pd
import math
import numba

import sys
sys.path.append('/Users/pbelanger/ABPLocal/PySUSSIX')
import PySUSSIX.ducksussix.f90newton as f90newton



#===========================================
# Verbatim implementation of the SUSSIX algorithm
# Note: direct fortran translation found at the bottom of this file (will be deleted in the future)
#===========================================
profile = lambda x: x

def Hann(N,Nt = None,p=1):
    """
    Hann's window, centered at the middle of the dataset.
    ----------------------------------------------------
        N : index of the turn
        Nt: total number of turns
    ----------------------------------------------------
    """
    if Nt is None:
        Nt = np.max(N)
    if np.mod(Nt,2) == 0:
        center = Nt//2 - 1
    else:
        center = Nt//2
    return (2**p)*math.factorial(p)**2/(math.factorial(2*p)) * (1+np.cos(2*np.pi*(N-center)/Nt))**(p)



@numba.jit(nopython=True)
def raise2power(my_exp,N):
    out    = np.zeros(len(N)) + 1j*np.zeros(len(N))
    out[0] = 1
    out[1] = my_exp
    for i in range(1,len(out)):
        out[i] = out[i-1]*out[1]
    return out


@profile
def Laskar_DFFT(freq,N,z):
    """
    Discrete fourier transform of z, as defined in A. Wolski, Sec. 11.5.
    In a typical DFFT , freq = m/Nt where m is an integer. Here m could take any value.
    Note: this will differ from sussix.f.calcr by a factor 1/Nt (but does not slow the convergence of the Newton method, so not a problem)
    ----------------------------------------------------
        freq: discrete frequency to evaluate the DFFT at
        N   : turn numbers of the signal
        z   : complex array of the signal
    ----------------------------------------------------
    """
    Nt = len(z)
    deriv_factor = 1j*N

    # to_sum = 1/Nt*np.exp(-2*np.pi*1j*freq*N)*z
    to_sum = 1/Nt*raise2power(np.exp(-2*np.pi*1j*freq),N)*z

    _DFFT            = np.sum(to_sum)
    _DFFT_derivative = np.sum(deriv_factor*to_sum)
    return _DFFT,_DFFT_derivative

@profile
def Laskar_DFFT_derivative(freq,N,z):
    """
    Derivative of laskar_DFFT w.r.t freq
    Note: differs from exact derivative by a factor -2pi (but allows newton method to converge)
    """
    Nt = len(z)
    deriv_factor = 1j*N
    return np.sum(1/Nt*deriv_factor*1*np.exp(-2*np.pi*1j*freq*N)*z)




def FFT_tune_estimate(z,n_forced = None):
    """
    Estimate the tune using an FFT. The signal is cropped to the closest power of 2 for improved accuracy.
    ----------------------------------------------------
        z        : complex array of the signal
        n_forced : number of turns to use for the FFT. if > len(z), the signal is padded with zeros
    ----------------------------------------------------
    """
    # Cropping signal to closest power of 2
    if n_forced is None:
        n_forced = 2**int(np.log2(len(z)))

    # Search for maximum in Fourier spectrum
    z_spectrum = np.fft.fft(z,n=n_forced)
    idx_max    = np.argmax(np.abs(z_spectrum))
    
    # Estimation of Tune with FFT
    tune_est   = idx_max/n_forced
    resolution = 1/n_forced

    return tune_est,resolution

@profile
def fundamental_frequency(x,px,Hann_order = 1,optimization = 'fortran'):
    """
    Subroutine of the NAFF algorithm. 
    1. Applies a Hann window
    2. Estimates the fundamental frequency with an FFT 
    3. Refines the estimate with a complex Newton method
    4. Returns the main frequency and the amplitude.
    """

    # Windowing of the signal
    N   = np.arange(len(x))
    z   = np.array(x) - 1j*np.array(px)
    z_w = z * Hann(N, Nt=len(z),p=Hann_order)
    
    # Estimation of the tune with FFT
    tune_est,resolution = FFT_tune_estimate(z_w)

    # Preparing the estimate for the Newton refinement method
    if tune_est >= 0.5:
        tune_est = -(1.0 - tune_est)
    tune_est = tune_est - resolution

    # Refinement of the tune calulation
    if optimization == 'fortran':
        tune,amplitude = f90newton_method(z_w,tune_est,resolution)
    else:
        tune,amplitude = newton_method(z_w,N,tune_est,resolution)

    return tune,amplitude

@profile
def NAFF(x,px,number_of_harmonics = 5,Hann_order = 1,optimization = 'fortran'):
    """
    Applies the NAFF algorithm to find the spectral lines of a signal.
    """

    assert number_of_harmonics >=1, 'number_of_harmonics needs to be > 1'
    
    # Converting to numpy arrays
    x,px = np.array(x),np.array(px)

    # initialisation
    z  = x - 1j*px
    N  = np.arange(len(x))
    
    
    frequencies = []
    amplitudes  = [] 
    for _ in range(number_of_harmonics):

        # Computing frequency and amplitude
        freq,amp  = fundamental_frequency(x,px,Hann_order=Hann_order,optimization=optimization)

        # Saving results
        frequencies.append(freq)
        amplitudes.append(amp)

        # Substraction procedure
        zgs  = amp * np.exp(2 * np.pi * 1j * freq * N)
        z   -= zgs
        x,px = np.real(z), -np.imag(z)

    
    return pd.DataFrame({'amplitude':amplitudes,'frequency':frequencies})


@profile
def get_harmonics(x = None,px = None,y = None,py = None,zeta = None,pzeta = None,number_of_harmonics = 5,Hann_order = 1,optimization = 'fortran'):
    """
    Computes the spectrum of a tracking data set for all canonical pairs provided
    """

    results = {}
    for pair in [(x,px,'x'),(y,py,'y'),(zeta,pzeta,'zeta')]:
        z,pz,plane = pair
        
        if z is not None:
            if pz is None:
                pz = np.zeros(len(z))

            # Computing spectral lines
            df = NAFF(z,pz, number_of_harmonics = number_of_harmonics,
                                Hann_order      = Hann_order,
                                optimization    = optimization)
            
            results[plane] = df

        else:
            results[plane] = None


    return results


def find_linear_combinations(frequencies,fundamental_tunes = [],max_jklm = 10):
    """
    Categorisation of resonances. Returns the linear combinations of the fundamental tunes that are closest to the provided frequencies.
    This should be called after get_harmonics to have a list of frequencies.
    """

    assert len(fundamental_tunes) == 3, "Need 3 fundamental tunes"

    # Create a 3D array of all possible combinations of j, k, l
    j,k,l,m = np.mgrid[-max_jklm:max_jklm+1, -max_jklm:max_jklm+1,-max_jklm:max_jklm+1,-max_jklm:max_jklm+1]

    # Q_eff = j*Q_x + k*Q_y + l*Q_z + m
    all_combinations = j * fundamental_tunes[0] + k * fundamental_tunes[1] + l * fundamental_tunes[2] + m
    
    # Find the closest combination for each frequency
    jklm = []
    err = []
    for freq in frequencies:

        # Find the index of the closest combination
        closest_idx = np.unravel_index(np.argmin(np.abs(freq - all_combinations)), all_combinations.shape)

        # Get the corresponding values for l, j, k
        closest_combination = (j[closest_idx], k[closest_idx], l[closest_idx],m[closest_idx])
        closest_value = all_combinations[closest_idx]

        jklm.append(closest_combination)
        err.append(np.abs(closest_value-freq))

    return pd.DataFrame({'jklm':jklm,'err':err,'freq':frequencies})


def f90newton_method(z,tune_est,resolution):


    # Prepare fortran arguments
    maxn    = len(z)
    tunea1  = tune_est
    deltat  = resolution
    maxiter = len(z)
    
    # Run fortran code
    f90newton.zfunr(z,maxn,tunea1,deltat,maxiter)
    tune      = f90newton.zfunr_out.tune_out[0]
    amplitude = f90newton.zfunr_out.zw_out[0]
    return tune,amplitude/maxn

@profile
def newton_method(z,N,tune_est,resolution):
    # Increase resolution by factor 5
    resolution /= 5

    tune_test = np.zeros(10)
    tune_val  = np.zeros(10)


    # DFFT,DFFT_d = Laskar_DFFT(tune_est,N,z)
    # # DFFT_d = Laskar_DFFT_derivative(tune_est,N,z)
    # dtunea1 = DFFT.real*DFFT_d.real + DFFT.imag*DFFT_d.imag


    # old name
    tunea1 = tune_est
    deltat = resolution
    err    = 1e-10
    num    = 0
    

    DFFT,DFFT_d = Laskar_DFFT(tunea1,N,z)
    # DFFT_d = Laskar_DFFT_derivative(tunea1,N,z)

    
    dtunea1 = DFFT.real*DFFT_d.real + DFFT.imag*DFFT_d.imag



    tunea2 = 0
    dtunea2 = 0
    
    for ntest in range(1, 11):
        tunea2 = tunea1+deltat

        DFFT ,DFFT_d  = Laskar_DFFT(tunea2,N,z)
        # DFFT_d = Laskar_DFFT_derivative(tunea2,N,z)
        dtunea2 = DFFT.real*DFFT_d.real + DFFT.imag*DFFT_d.imag

        
        if (dtunea1 <= 0) and (dtunea2 >= 0):
            tune1, tune2, dtune1, dtune2 = tunea1, tunea2, dtunea1, dtunea2

            
            for ncont in range(1, 101):
                ratio = -dtune1 / dtune2 if abs(dtune2) > 0 else 0.0

                tune3 = (tune1 + ratio * tune2) / (1.0 + ratio)

                
                DFFT,DFFT_d = Laskar_DFFT(tune3,N,z)
                # DFFT_d = Laskar_DFFT_derivative(tune3,N,z)
                dtune3 = DFFT.real*DFFT_d.real + DFFT.imag*DFFT_d.imag


                if dtune3 <= 0.0:
                    if tune1 == tune3:
                        break
                    tune1, dtune1 = tune3, dtune3
                else:
                    if tune2 == tune3:
                        break
                    tune2, dtune2 = tune3, dtune3

                if abs(tune2 - tune1) <= err:
                    break

            
            tune_test[num] = tune3
            tune_val[num]  = np.abs(DFFT)
            num += 1
            


        tunea1, dtunea1 = tunea2, dtunea2

    idx_max = np.argmax(tune_val[:num])
    tune      = tune_test[idx_max]
    amplitude,_ = Laskar_DFFT(tune,N,z)
    return tune,amplitude


#===========================================
# 1-TO-1 PYTHON IMPLEMENTATION OF SUSSIX
#===========================================

def calcr(tune_phasor, z):
    return np.polyval(z[::-1], tune_phasor)


def tunenewt(x,px,hanning_order = 1):
    """COMPUTES THE TUNE USING A DISCRETE VERSION OF LASKAR METHOD.
        IT INCLUDES A NEWTON METHOD FOR THE SEARCH OF THE FREQUENCY."""
    
    import PySUSSIX.crosssussix.crossroutines as crossroutines
    import PySUSSIX.crosssussix.crosssussix as crosssussix

    # Estimation of Tune with FFT
    maxn = len(x)
    mft = int(np.log2(maxn))
    npoint = 2**mft
    maxn2 = maxn // 2
    step = 2 * np.pi / maxn

    mf_values = np.arange(1, maxn +1)


    z = (np.array(x) + 1j *np.array(px)) * (1.0 + np.cos(step * (mf_values - maxn2)))


    # Search for maximum in Fourier spectrum
    # zsing = crossroutines.cfft(z)
    zsing = np.fft.fft(z[:npoint])
    ftmax = np.max(np.abs(zsing[:npoint]))
    nftmax = np.argmax(np.abs(zsing[:npoint]))
    
    tunefou = float(nftmax ) / float(npoint)
    
    if tunefou >= 0.5:
        tunefou = -(1.0 - tunefou)
    deltat = 1.0 / float(npoint)

    
    tune1 = tunefou - deltat
    
    
    # Call zfun
    tune,zw = crossroutines.zfunr(z,tune1)

    return tune,zw


def spectrum(x,px,number_of_harmonics = 5,Hann_order = 1):
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

        # # Computing frequency and amplitude
        # freq,zw = tunenewt(x,px)
        # # freq,zw = crossroutines.tunenewt(x,px)
        # zpesi   = zw / max(N+1)

        freq,amp  = fundamental_frequency(x,px,Hann_order=Hann_order)

        # Saving results
        frequencies.append(freq)
        amplitudes.append(amp)

        # Substraction procedure
        zgs  = amp * np.exp(2 * np.pi * 1j * freq * N)
        z   -= zgs
        x,px = np.real(z), np.imag(z)

    
    return pd.DataFrame({'amplitude':amplitudes,'frequency':frequencies})

    


def datspe(x = None,px = None,y = None,py = None,zeta = None,pzeta = None,number_of_harmonics = 5,Hann_order = 1):
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
                                Hann_order          = Hann_order )
            
            results[plane] = df

        else:
            results[plane] = None


    return results


    