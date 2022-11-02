# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:00:37 2022

@author: derto
"""

import numpy as np
#from scipy import interpolate
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
import misc as misc

import time
st = time.time()



    
### initial conditions

#grid = Integration.grids.equi(P.k0, P.k1, P.nk)
grid = Integration.grids.gts2(P.k0, P.k1, P.nk)
k = grid[0]
epsk = C.hbar**2 * k**2 / (2*P.mu)
psik0 = np.zeros(len(k), dtype=np.complex_)


### calculate psi(t) (rotating frame)
    
def rhs(k, coulomb=True):
    def ode(t, psik):
        Dpsik = (epsk - P.Phi - 1j*P.damp)*psik - P.dipole*misc.E0(t)
        if coulomb is True:
            Dpsik -= P.cI * int_integral(psik)
        Dpsik  *= (-1j) / C.hbar 
        #print(t)
        return Dpsik
    return ode


def int_integral(psik):
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    np.seterr(divide='ignore')
    v_ij = np.where(np.abs(ki-kj) > 0,
                    grid[1] * kj / ki * (np.log(np.abs(ki+kj))-np.log(np.abs(ki-kj))), 0) 
    np.seterr(divide='warn')
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+np.pi*ki)
    I = np.sum(w_ij*psik, axis=1)
    #print(I.shape)                
    return I 

   
def psi_of_t(coulomb=True):
    method = ODE.solvers.RK4(rhs(k, coulomb)) 
    t, sol = ODE.main.solve(method, P.t0, P.t1, P.nt, psik0)
    #sol = sol*np.exp( -1j * P.w0*t.reshape((len(t), 1)) )
    #t = t[::2]
    #sol = sol[::2, :]
    np.savetxt("sol_t/time", t)
    np.savetxt("sol_t/momentum", k)
    np.savetxt("sol_t/sol", sol.view(float))
psi_of_t()
#psi_of_t(False)



### calculate P(t) (rotating frame)

def pol_of_t():
    t = np.loadtxt("sol_t/time")
    sol = np.loadtxt("sol_t/sol").view(complex)
    pol = np.zeros(len(t), dtype=np.complex_)
    integrand = 1/(2*np.pi**2)*np.conjugate(P.dipole)*sol*k**2 #* np.exp(-1j*P.w0*t)
    #print(integrand.shape)
    for ii in range(len(t)):
        pol[ii] = Integration.integrater.int_disc(
            integrand[ii, :], grid)
    np.savetxt("sol_t/polt", pol.view(float))
pol_of_t()



### calculate P(w) (rotating frame)
    
def pol_of_w(): # t-> w: exp(iwt), i.e. inverse
    polt = np.loadtxt("sol_t/polt").view(complex)
    t = np.loadtxt("sol_t/time")
    w = np.linspace(-40, 50, 1000)
    polw = misc.fourier_trafo(t, polt, w, inverse=True)
    np.savetxt("sol_t/frequency", w)    
    np.savetxt("sol_t/polw", polw.view(float)) 
pol_of_w()    





   



et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
