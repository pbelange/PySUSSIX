
import numpy as np
import matplotlib.pyplot as plt
import time
import numba


Q0 = 0.31025793875089835

Nt = int(1e5)
N  = np.arange(Nt)

my_exp = np.exp(1j*2*np.pi*Q0)
target = np.exp(1j*2*np.pi*Q0*N)

# trick  = my_exp**N

#=========================
def reference(freq,N,z):
    Nt = len(z)
    to_sum = 1/Nt*np.exp(-2*np.pi*1j*freq*N)*z
    _DFFT  = np.sum(to_sum)
    return _DFFT
#=========================


#=========================
@numba.jit(nopython=True)
def custom_exp(N,my_exp):
    out    = np.zeros(len(N)) + 1j*np.zeros(len(N))
    out[0] = 1
    out[1] = my_exp
    for i in range(1,len(out)):
        out[i] = out[i-1]*out[1]
    return out
#----
# @numba.jit(nopython=True)
def test1(freq,N,z):
    Nt     = len(z)
    to_sum = 1/Nt*custom_exp(N,np.exp(-2*np.pi*1j*freq))*z
    _DFFT  = np.sum(to_sum)
    return _DFFT

#=========================


#=========================
@numba.jit(nopython=True)
def custom_exp2(N,my_exp):
    return np.power(my_exp,N)
#----
def test2(freq,N,z):
    Nt     = len(z)
    to_sum = 1/Nt*custom_exp2(N,np.exp(-2*np.pi*1j*freq))*z
    _DFFT  = np.sum(to_sum)
    return _DFFT

#=========================


#=========================
@numba.jit(nopython=True)
def custom_exp3(N, my_exp):
    exp_array = np.repeat(my_exp,len(N))
    return np.cumprod(exp_array)
#----
# @numba.jit(nopython=True)
def test3(freq,N,z):
    Nt     = len(z)
    to_sum = 1/Nt*custom_exp3(N,np.exp(-2*np.pi*1j*freq))*z
    _DFFT  = np.sum(to_sum)
    return _DFFT

#=========================


#=========================
@numba.jit(nopython=True)
def custom_exp4(N,my_exp,out):
    
    out[0] = 1
    out[1] = my_exp
    for i in range(1,len(out)):
        out[i] = out[i-1]*out[1]
    return out
#----
# @numba.jit(nopython=True)
def test4(freq,N,z):
    Nt     = len(z)
    out    = np.zeros(len(N)) + 1j*np.zeros(len(N))
    out    = custom_exp4(N,np.exp(-2*np.pi*1j*freq),out)
    to_sum = 1/Nt*out*z
    _DFFT  = np.sum(to_sum)
    return _DFFT

#=========================

for fun in [reference,test1,test2,test3,test4]:
    # Compiling
    _DFFT = fun(Q0,N,target)
    # Timing
    t1 = time.perf_counter()
    for _ in range(100):
        _DFFT = fun(Q0,N,target)
    t2 = time.perf_counter()
    print('=====================')
    print('function: ',fun.__name__)
    print(f'elapsed time: {(t2-t1)/100/1e-3:.4f} ms')
    print('=====================')




#=========================
def test5(freq,N,z,tune_phasor):
    Nt = len(z)
    to_sum = 1/Nt*(tune_phasor**freq)*z
    _DFFT  = np.sum(to_sum)
    return _DFFT
#=========================

# New test
fun  = test5
# Compiling

# Timing
t1 = time.perf_counter()
tune_phasor = np.exp(-2*np.pi*1j*N)
for _ in range(500):
    _DFFT = fun(Q0,N,target,tune_phasor)
t2 = time.perf_counter()
print('=====================')
print('function: ','test5')
print(f'elapsed time: {(t2-t1)/500/1e-3:.4f} ms')
print('=====================')