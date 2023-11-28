import numpy as np

def spectrum(x_values, xp_values, max_points, narm, method):
    """
    Compute the main frequency.

    Parameters:
    - x_values: Array of x-coordinate values.
    - xp_values: Array of xp-coordinate values.
    - max_points: Length of the orbit.
    - narm: Number of harmonics to compute.
    - method: Method to select the window (1: Hanning, 2: Rectangular).

    Without orthogonalization of Gram-Schmidt.

    Returns:
    - tune_results: Array of computed frequencies.
    - zpesi_results: Array of intermediate complex values.
    """

    duepi = 8.0 * np.arctan(1.0)
    max_iterations = 100000
    small_value, zero, one, two, pi, piph = 1e-17, 0.0, 1.0, 2.0, np.pi, duepi

    z_values = np.array([complex(x, xp) for x, xp in zip(x_values, xp_values)], dtype=np.complex128)
    original_z_values = np.copy(z_values)

    if max_points > max_iterations:
        print('ERROR IN SPECTRUM: MAX_POINTS TOO LARGE')
        return None, None

    if narm < 2:
        print('ERROR IN SPECTRUM: NARM SMALLER THAN 2')
        return None, None

    tune_results = np.zeros(narm, dtype=np.float64)
    zpesi_results = np.zeros(narm, dtype=np.complex128)

    for harmonic_index in range(narm):
        if method == 1:
            tune_results[harmonic_index] = tunenewt(x_values, xp_values, max_points)
        elif method == 2:
            tune_results[harmonic_index] = tunelasr(x_values, xp_values, max_points)

        frequency = tune_results[harmonic_index]
        zpesi_results[harmonic_index] = zfunr(complex(0.0, frequency * piph), z_values, max_points) / float(max_points)
        z_values -= np.power(zpesi_results[harmonic_index], np.arange(1, max_points + 1))

    x_values, xp_values = original_z_values.real, original_z_values.imag

    return tune_results, zpesi_results



import numpy as np
from scipy.fft import fft

def tunenewt(x_values, xp_values, max_points):
    """
    Compute the tune using a discrete version of Laskar method.
    It includes a Newton method for the search of the frequency.

    Parameters:
    - x_values: Array of x-coordinate values.
    - xp_values: Array of xp-coordinate values.
    - max_points: Length of the orbit.

    Returns:
    - tune_result: Computed tune.
    """

    max_iterations = 100000
    small_value, zero, one, two, pi, piph = 1e-17, 0.0, 1.0, 2.0, np.pi, 8.0 * np.arctan(1.0)

    maxn = len(x_values)
    tune_result = 0.0

    # Estimation of tune with FFT
    mft = int(np.log2(maxn))
    npoint = 2 ** mft
    maxn2 = maxn // 2
    step = piph / maxn

    z = np.array([complex(x, xp) * (1.0 + np.cos(step * (mf - maxn2)))
                  for mf, (x, xp) in enumerate(zip(x_values, xp_values))], dtype=np.complex128)

    zsing = np.copy(z)
    zsing = fft(zsing)

    # Search for maximum of Fourier spectrum
    nftmax = np.argmax(np.abs(zsing))
    tunefou = (nftmax - 1) / npoint if nftmax > 0 else 0.0
    tune_result = zfunr(complex(0.0, tunefou * piph), z, maxn) / float(maxn)

    return tune_result


import numpy as np
from scipy.fft import fft

def tunelasr(x_values, xp_values, maxn):
    """
    Same as tunenewt but no Hanning filter.

    Parameters:
    - x_values: Array of floating-point values.
    - xp_values: Array of floating-point values.
    - maxn: Integer value.

    Returns:
    - Floating-point value.
    """
    duepi = np.pi * 8.0
    tune = 0.0
    zw = complex(0.0, 0.0)

    # Estimation of tune with FFT
    mft = int(np.log2(maxn))
    npoint = 2 ** mft
    maxn2 = maxn // 2
    step = duepi / maxn

    z = np.array([complex(xi, xpi) for xi, xpi in zip(x_values, xp_values)], dtype=np.complex128)
    zsing = np.fft.fft(z)

    # Search for the maximum of the Fourier spectrum
    nftmax = np.argmax(np.abs(zsing))
    tunefou = (nftmax - 1) / float(npoint) if nftmax > 0 else 0.0
    tune1 = tunefou - 1.0 / float(npoint)

    # Call zfunr to compute the tune
    tune = zfunr(complex(0.0, tune1 * duepi), z, maxn) / float(maxn)

    return tune


import numpy as np

def zfunr(ztune, z_values, maxn):
    """
    Auxiliary routine used by tunenewt.

    Parameters:
    - ztune: Complex value.
    - z_values: Array of complex values.
    - maxn: Maximum value for the loop.

    Returns:
    - tune: Tune value.
    - zw: Complex value.
    """

    small_value, zero, one, two, pi, piph = 1e-17, 0.0, 1.0, 2.0, np.pi, 8.0 * np.arctan(1.0)

    # Initialization
    err = 1e-10
    zu = 0.0 + 1.0j

    # We divide deltat in 5 parts
    deltat = 1.0 / float(maxn)
    deltat /= 5.0

    zd = zu * np.arange(1, maxn + 1) * z_values

    tunetest = np.zeros(10)
    tuneval = np.zeros(10)

    # Calculate ztune1
    zw = calcr(ztune, z_values, maxn)
    calcr(ztune, zd, maxn)
    dtunea1 = np.real(zw) * np.real(zd) + np.imag(zw) * np.imag(zd)
    num = 1

    for ntest in range(10):
        ztune2 = np.exp(-zu * piph * (ztune - 1j * deltat))
        zw = calcr(ztune2, z_values, maxn)
        calcr(ztune2, zd, maxn)
        dtunea2 = np.real(zw) * np.real(zd) + np.imag(zw) * np.imag(zd)

        if dtunea1 <= 0.0 and dtunea2 >= 0.0:
            tune1, tune2 = ztune, ztune2
            dtune1, dtune2 = dtunea1, dtunea2

            for ncont in range(100):
                if abs(dtune2) > 0:
                    ratio = -dtune1 / dtune2
                else:
                    ratio = 0.0

                tune3 = (tune1 + ratio * tune2) / (1.0 + ratio)
                ztune3 = np.exp(-zu * piph * tune3)
                zw = calcr(ztune3, z_values, maxn)
                calcr(ztune3, zd, maxn)
                dtune3 = np.real(zw) * np.real(zd) + np.imag(zw) * np.imag(zd)

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

            tunetest[num] = tune3
            tuneval[num] = abs(zw)
            num += 1

        ztune, dtunea1 = ztune2, dtunea2

    tune = tunetest[0]
    tunevmax = tuneval[0]

    for nc in range(1, num - 1):
        if tunevmax <= tuneval[nc]:
            tunevmax = tuneval[nc]
            tune = tunetest[nc]

    zw = calcr(np.exp(-zu * piph * tune), z_values, maxn)

    return tune, zw



import numpy as np

def calcr(zv, zp_values, maxd):
    """
    Auxiliary routine used by tunenewt.

    Parameters:
    - zv: Complex value.
    - zp_values: Array of complex values.
    - maxd: Maximum value for the loop.

    Returns:
    - zpp: Complex value.
    """
    zpp = zp_values[maxd - 1]

    for np in range(maxd - 2, -1, -1):
        zpp = zpp * zv + zp_values[np]

    return zpp
