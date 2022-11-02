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

w = np.linspace(-40, 50, 1000)

#grid = Integration.grids.equi(P.k0, P.k1, P.nk)
grid = Integration.grids.gts2(P.k0, P.k1, P.nk)
k = grid[0]


### calculte psi_k(w) (rotating frame)

def psi_of_w():
    k_2D = k.reshape((1, len(k)))
    w_2D = w.reshape((len(w), 1))
    E_bracket = - C.hbar*w_2D + C.hbar**2*k_2D**2/(2*P.mu) - P.Phi - 1j*P.damp
    psiw = P.dipole * misc.E0w(w_2D) / E_bracket
    np.savetxt("sol_w/momentum", k)
    np.savetxt("sol_w/psiw", psiw.view(float)) 
psi_of_w()    
 

### calculte psi_k(t) (rotating frame)

def psi_of_t():
    t = np.linspace(P.t0, P.t1, 1000)
    psiw = np.loadtxt("sol_w/psiw").view(complex)
    psit = np.zeros((len(t), len(k)), dtype=np.complex_)
    for ii in range(len(k)):    
        psit[:, ii] = misc.fourier_trafo(w, psiw[:, ii], t, inverse=False)
    np.savetxt("sol_w/time", t)
    np.savetxt("sol_w/psit", psit.view(float))    
psi_of_t()


### calculate chi(w) (rotating frame)

def chi_of_w():
    chiw = np.zeros(len(w), dtype=np.complex_)
    k_2D = k.reshape((1, len(k)))
    w_2D = w.reshape((len(w), 1))
    epsk = C.hbar**2 * k_2D**2 / (2*P.mu)
    integrand = k_2D**2/(- C.hbar*w_2D + epsk - P.Phi - 1j*P.damp)
    for ii in range(len(w)):
        chiw[ii] = Integration.integrater.int_disc(integrand[ii, :], grid)   
    chiw *= np.abs(P.dipole)**2 / (2*np.pi**2)   
    np.savetxt("sol_w/frequency", w)
    np.savetxt("sol_w/chi", chiw.view(float))     
chi_of_w()




et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
