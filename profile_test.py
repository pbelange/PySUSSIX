import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import scipy.stats as sciStat
import time

import PySUSSIX.ducksussix.ducksussix as ducksussix


Q0 = 0.31025793875089835
Qs = 0.002
dQ = Qs/12
Jx = (0.5*(10**2))

n_bands_Qs = 1
n_bands_dQ = 5
N   = np.arange(int(1e5))
j   = np.arange(-n_bands_Qs,n_bands_Qs+1)
i   = np.arange(-n_bands_dQ,n_bands_dQ+1)


Ai  = sciStat.cauchy.pdf(i/np.max(i),0,0.05)
Aj  = sciStat.cauchy.pdf(j/np.max(j),0,0.05)
Ai  = Ai/np.max(Ai)
Aj  = Aj/np.max(Aj)
np.random.seed(0)
phii = np.zeros(len(i))#np.random.uniform(-np.pi/2,np.pi/2,len(i))
phij = np.zeros(len(j))#np.random.uniform(-np.pi/2,np.pi/2,len(j)) 
# x = sum([sum([ _Ai*_Aj * np.sin(2*np.pi*(Q0+ _j*Qs + _i*dQ)*N) for _i,_Ai in zip(i,Ai) ]) for _j,_Aj in zip(j,Aj) ])
hx_full = sum([sum([ np.sqrt(2*Jx)*_Ai*_Aj * np.exp(1j*2*np.pi*(Q0+ _j*Qs + _i*dQ)*N + _phii + _phij) for _i,_Ai,_phii in zip(i,Ai,phii) ]) for _j,_Aj,_phij in zip(j,Aj,phij) ])


#=================================================

x,px = hx_full.real,-hx_full.imag
N_vec = np.logspace(3,5,50).astype(int)
# N_vec = list(set(2**np.linspace(np.log2(10**3),np.log2(10**5)+1,100).astype(int)))
freq_df = []
amp_df = []
exc_df = []
N = int(1e5)

t1 = time.perf_counter()
results = ducksussix.get_harmonics( x       = x[:N], 
                                    px      = px[:N],
                                    y       = None,
                                    py      = None,
                                    zeta    = None,
                                    pzeta   = None,
                                    number_of_harmonics = 32,Hann_order = 1,
                                    optimization='python')

t2 = time.perf_counter()
print(results['x'])
print(f'elapsed time: {t2-t1:.4f} s')