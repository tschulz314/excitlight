# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:00:37 2022

@author: derto
"""

import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

import sys
import os


# get current directory
path = os.getcwd()
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))
import Integration
import constants as C
import parameters as P
import misc as misc

import time
st = time.time()




def TMDC_pot(q):
    chi =   7.112 * 0.1 
    kappa = 2 
    r0 = 2 * np.pi * chi / kappa
    V = C.e**2 / (2 * C.eps0 * q) 
    W = V * 1 / kappa * 1 / (1 + r0 * q)
    return  W
#TMDC_pot(k)


def int_grid_eff():
    double = np.zeros(len(k))
    nkjj = 200
    nphi = 200
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], 0, 10, nkjj)
        kjj = kjj_grid[0]
        kjj_integrand = np.zeros(len(kjj))
        for jj in range(len(kjj)):
            sqrt = np.sqrt(k[ii]**2 + kjj[jj]**2 - 2*k[ii]*kjj[jj]*np.cos(phi)) 
            phi_integrand = kjj[jj] * TMDC_pot(sqrt) * 4 * k[ii]**4 / (k[ii]**2+kjj[jj]**2)**2
            kjj_integrand[jj] = Integration.integrater.int_disc(phi_integrand, phi_grid)  
        double[ii] = np.sum(kjj_grid[1]*kjj_integrand)
    
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                sqrt = np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi))
                phi_integrand = k[jj] * TMDC_pot(sqrt)
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)  
    return w_ij


grid = Integration.grids.gts2(0, 5, 50)
k = grid[0]

def eig():
    eps = C.hbar**2 * k**2 / (2*P.mu) - P.Phi #-1j*P.damp
    w_ij = int_grid_eff()
    mat = np.diag(eps) - 1/(2*np.pi)**2*w_ij 
    energy, psi = linalg.eigh(mat)
    np.savetxt("sol_eig/momentum", k)
    np.savetxt("sol_eig/omega", energy/C.hbar)
    np.savetxt("sol_eig/psi", psi.view(float)) 
    #return energy, psi
eig()



#plt.plot(energy, energy-energy, 'x')
# ind = 10
# print(energy[ind])
# plt.plot(k, psi[:, ind])
#plt.xlim(0, 0.3)




# w = energy / C.hbar
# def chi_of_w():
#     chiw = np.zeros(len(w), dtype=np.complex_)
#     #integrand = 
#     for ii in range(len(w)):
#         chiw[ii] = Integration.integrater.int_disc(psi, grid) 
#     k_2D = k.reshape((1, len(k)))
#     w_2D = w.reshape((len(w), 1))
#     epsk = C.hbar**2 * k_2D**2 / (2*P.mu)
#     integrand = k_2D**2/(- C.hbar*w_2D + epsk - P.Phi - 1j*P.damp)
#     for ii in range(len(w)):
#         chiw[ii] = Integration.integrater.int_disc(integrand[ii, :], grid)   
#     chiw *= np.abs(P.dipole)**2 / (2*np.pi**2)   
#     np.savetxt("sol_w/frequency", w)
#     np.savetxt("sol_w/chi", chiw.view(float))     
# chi_of_w()





et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
