# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:00:37 2022

@author: derto
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import sys
import os


# get current directory
path = os.getcwd()
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))
import ODE
import Integration
import constants as C
import parameters as P

import time
st = time.time()



    
### initial conditions

k0 = 0
k1 = 500
nk = 100

k = Integration.grids.gts1(k0, k1, nk)[0]

psik0 = np.zeros(len(k), dtype=np.complex_)

w = np.linspace(0, 100, 1000)
chiw = np.zeros(len(w), dtype=np.complex_)
alpha = np.zeros(len(w))

k_2D = k.reshape((1, len(k)))
w_2D = w.reshape((len(w), 1))


epsk = C.hbar**2 * k_2D**2 / (2*P.mu)
integrand = k_2D**2/(C.hbar*w_2D + C.hbar**2*k_2D**2/(2*P.mu) + P.Eg - C.hbar*P.w0 -1j*P.damp)
integrand2 = epsk/((C.hbar*w_2D + epsk + P.Eg - C.hbar*P.w0)**2  + P.damp**2)

for ii in range(len(w)):
    chiw[ii] = Integration.integrater.int_disc(integrand[ii, :], Integration.grids.gts1(k0, k1, nk))
    alpha[ii] = Integration.integrater.int_disc(integrand2[ii, :], Integration.grids.gts1(k0, k1, nk)) 
    
chiw *= np.abs(P.dipole)**2 / (2*np.pi**2)   
alpha *= np.abs(P.dipole)**2 * np.sqrt(2*P.mu**3) / (np.pi**2*C.hbar**3) * P.damp     

#plt.plot(w, chiw.imag)    
#plt.plot(w, alpha)
#print(chiw.imag)
    
    
#pol[ii] = Integration.integrater.int_disc(
#            integrand[ii, :], Integration.grids.gts1(P.k0, P.k1, P.nk))    

    
#w = np.loadtxt("solutions/frequency")
#Ew = np.loadtxt("solutions/Ew").view(complex)


 

#psi_k_w = P.dipole * 

et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
