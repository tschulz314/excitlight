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
    return  W #V / P.eps
#TMDC_pot(k)


grid = Integration.grids.gts2(0, 5, 100)
k = grid[0]


def int_grid_eff():
    double = np.zeros(len(k))
    nkjj = 1000
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


def norm(psi):
    psinew = np.zeros(psi.shape) # psi(k, w)
    for ii in range(psi.shape[1]):    
        norm = 1/(2*np.pi)**2 * np.sum(np.abs(psi[:, ii])**2*grid[1])
        psinew[:, ii] = psi[:, ii] / np.sqrt(norm)
        #print(norm)
    return psinew    


def eig():
    eps = C.hbar**2 * k**2 / (2 * P.mu) - P.Phi #-1j*P.damp
    w_ij = int_grid_eff()
    mat = np.diag(eps) - 1/(2*np.pi)**2*w_ij 
    energy, psi = linalg.eig(mat) # energy[i], psi[k, i]
    energy, psi = misc.sort(energy, psi)
    psi = norm(psi)
    np.savetxt("sol_eig/momentum", k)
    np.savetxt("sol_eig/frequency", energy.real/C.hbar)
    np.savetxt("sol_eig/psi", psi)#.view(float))
    #return energy, psi
eig()


#w = np.loadtxt("sol_eig/frequency")
#psi = np.loadtxt("sol_eig/psi")#.view(complex)
#print(psi.shape)

#psi2 = norm(psi)



#print(w.shape)
#nw = 20
#w = w[:nw]
#psi = psi[:, :nw]
#w = np.sort(w)
#print(w)

#plt.plot(w, w-w, 'x')
#plt.axvline(4*P.ryd_frq)
#plt.xlim(-600, 50)
#ind = 0
#print(k.shape, psi[:, :].shape)
#print(C.hbar * w[ind])
#plt.plot(k, psi[:, ind].imag)
#plt.xlim(0, 0.3)








et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
