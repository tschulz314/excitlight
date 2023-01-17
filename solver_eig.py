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
import interaction as inter

import time
st = time.time()



### k grid

k, kgrid = Integration.grids.gts2(0, 4, 200) # 4, 30
#k = grid[0]


def norm(psi):
    psinew = np.zeros(psi.shape) # psi(k, w)
    for ii in range(psi.shape[1]):    
        norm = 1/(2*np.pi) * np.sum(np.abs(psi[:, ii])**2*kgrid*k)
        psinew[:, ii] = psi[:, ii] / np.sqrt(norm)
    #print(norm)
    return psinew    


def eig():
    eps = C.hbar**2 * k**2 / (2 * P.mu) - P.Phi #-1j*P.damp
    w_ij = inter.int_grid_eff(k, kgrid) #int_grid_eff()
    mat = np.diag(eps) - 1/(2*np.pi)**2*w_ij 
    energy, psi = linalg.eig(mat) # energy[i], psi[k, i]
    #print(psi[0, :])
    #np.savetxt("sol_eig/psi_complex", psi.view(float))
    #print(psi.imag)
    #psi=psi.real
    energy, psi = misc.sort(energy, psi)
    psi = norm(psi)
    #psi = norm(psi)
    np.savetxt("sol_eig/momentum", k)
    np.savetxt("sol_eig/momentum_weights", kgrid)
    np.savetxt("sol_eig/frequency", energy.real/C.hbar)
    np.savetxt("sol_eig/psi", psi)#.view(float))
    #return energy, psi
eig()


def chi():
    w = np.loadtxt("sol_eig/frequency")
    psi = np.loadtxt("sol_eig/psi")
    ### calc psi(r=0)
    psi0 = np.zeros(len(w))
    for ii in range(len(w)):
        psi0[ii] = 1/(2*np.pi) * np.sum(kgrid*psi[:, ii]*k)
    ### calc chi      
    #print(len(w))
    nw = 5000
    #w = w.reshape((1, len(w)))
    #psi0 = psi0.reshape(w.shape)
    w_for_chi = np.linspace(-600, 50, nw).reshape((nw, 1))
    sumover = np.abs(psi0)**2 * 1 / (C.hbar*(w_for_chi - w) + 1j*P.damp) #-1/(C.hbar*(w_for_chi + w + 1j*P.damp)))
    chi = - np.abs(P.dipole)**2 * np.sum(sumover, axis=1) / C.eps0
    np.savetxt("sol_eig/w", w_for_chi)
    np.savetxt("sol_eig/chiw", chi.imag)
    #return w_for_chi, chi
chi()

#w_for_chi, chi = chi()

#print(w2.shape, chi.shape, chi)

#plt.plot(w_for_chi[:, 0], chi.imag)


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
