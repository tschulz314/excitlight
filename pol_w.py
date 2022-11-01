# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:00:37 2022

@author: derto
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

path = os.getcwd()
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))
import Integration
import constants as C
import parameters as P
import misc as misc

import time
st = time.time()



### initial conditions

w = np.linspace(-100, 100, 2000)

#grid = Integration.grids.equi(P.k0, P.k1, P.nk)
grid = Integration.grids.gts2(P.k0, P.k1, P.nk)
k = grid[0]



### calculte psi_k(w) (rotating frame)


def psi_over_w():
    k_2D = k.reshape((1, len(k)))
    w_2D = w.reshape((len(w), 1))
    E_bracket = C.hbar*w_2D + C.hbar**2*k_2D**2/(2*P.mu) - P.Phi - 1j*P.damp
    psi = P.dipole * misc.E0w(w_2D) / E_bracket
    return psi
# psi = psi_over_w() 
 
# plt.plot(w, psi[:,0].real)
# plt.show()
# plt.close()


### calculte psi_k(t) (rotating frame)

def psi_over_t():
    t = np.linspace(P.t0, P.t1, 1000)
    #t = np.linspace(0, 10, 1000)
    psiw = psi_over_w() 
    psit = np.zeros((len(t), len(k)), dtype=np.complex_)
    for ii in range(len(k)):    
        t, psit[:, ii] = misc.fourier_trafo(w, psiw[:, ii], t[0], t[-1], len(t),
                                            inverse=True)
    np.savetxt("sol_w/time", t)
    np.savetxt("sol_w/momentum", k)
    np.savetxt("sol_w/psit", psit.view(float))    
    return t, psit
#t, psit = psi_over_t()


# plt.plot(t, psit[:, -1])
# plt.plot(k, psit[int(0), :].real)
# plt.show()
# plt.close()


### calculate chi(w) (rotating frame)

def chi_over_w():
    chiw = np.zeros(len(w), dtype=np.complex_)
    #alpha = np.zeros(len(w))
    k_2D = k.reshape((1, len(k)))
    w_2D = w.reshape((len(w), 1))
    epsk = C.hbar**2 * k_2D**2 / (2*P.mu)
    integrand = k_2D**2/(-C.hbar*w_2D + epsk + P.Eg - C.hbar*P.w0 -1j*P.damp)
    #integrand2 = epsk/((C.hbar*w_2D + epsk + P.Eg - C.hbar*P.w0)**2  + P.damp**2)
    for ii in range(len(w)):
        chiw[ii] = Integration.integrater.int_disc(integrand[ii, :], grid) 
        #alpha[ii] = Integration.integrater.int_disc(integrand2[ii, :], grid)     
    chiw *= np.abs(P.dipole)**2 / (2*np.pi**2)   
    #alpha *= np.abs(P.dipole)**2 * np.sqrt(2*P.mu**3) / (np.pi**2*C.hbar**3) * P.damp
    np.savetxt("sol_w/freq", w)
    np.savetxt("sol_w/chi", chiw.view(float))  
    return w, chiw     

w, chiw = chi_over_w()

# plt.plot(w, chiw.imag)    
# plt.plot(w, alpha)
# plt.show()
# plt.close()



et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
