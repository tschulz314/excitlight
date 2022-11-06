# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:00:37 2022

@author: derto
"""

import numpy as np
from scipy import interpolate
from scipy import integrate
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
    

def int_grid_3D():
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    np.seterr(divide='ignore')
    v_ij = np.where(np.abs(ki-kj) > 0,
                    grid[1] * kj / ki * (np.log(np.abs(ki+kj))-np.log(np.abs(ki-kj))), 0) 
    np.seterr(divide='warn')
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+np.pi*ki)              
    return w_ij 


def int_grid_2D():
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nkjj = 5000
    nphi = 500
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], k[0], 10*k[-1], nkjj)
        kjj = kjj_grid[0]
        kjj_integrand = np.zeros(len(kjj))
        for jj in range(len(kjj)):
            sqrt = np.sqrt(k[ii]**2 + kjj[jj]**2 - 2*k[ii]*kjj[jj]*np.cos(phi)) 
            phi_integrand = kjj[jj] / sqrt * 4 * k[ii]**4 / (k[ii]**2+kjj[jj]**2)**2
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
                phi_integrand = k[jj]/sqrt
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)  
    return w_ij


def int_grid_eff(D):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nkjj = 5000
    nphi = 500
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], k[0], 5*k[-1], nkjj)
        kjj = kjj_grid[0]
        kjj_integrand = np.zeros(len(kjj))
        for jj in range(len(kjj)):
            sqrt = np.sqrt(k[ii]**2 + kjj[jj]**2 - 2*k[ii]*kjj[jj]*np.cos(phi)) 
            phi_integrand = kjj[jj] / sqrt**D * 4 * k[ii]**4 / (k[ii]**2+kjj[jj]**2)**2
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
                phi_integrand = k[jj] / sqrt**D
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)  
    return w_ij


def rhs(k, dim=3, coulomb=True):
    if coulomb is True:
        if dim == 3:
            w_ij = int_grid_3D()
            cI = P.cI_3D
        if dim == 2:
            w_ij = int_grid_2D()
            cI = P.cI_2D
        else:
            w_ij = int_grid_eff(1.2)
            cI = P.cI_2D
    else:
        w_ij = np.zeros((len(k), len(k)))
        cI = 0
        
    def ode(t, psik):
        Dpsik = (epsk - P.Phi - 1j*P.damp)*psik - P.dipole*misc.E0(t)
        Dpsik -= cI * np.sum(w_ij*psik, axis=1)
        Dpsik  *= (-1j) / C.hbar 
        #print(t)
        return Dpsik
    return ode


def psi_of_t(dim=3, coulomb=True):
    method = ODE.solvers.RK4(rhs(k, dim, coulomb)) 
    t, sol = ODE.main.solve(method, P.t0, P.t1, P.nt, psik0)
    #sol = sol*np.exp( -1j * P.w0*t.reshape((len(t), 1)) )
    #t = t[::2]
    #sol = sol[::2, :]
    np.savetxt("sol_t/time", t)
    np.savetxt("sol_t/momentum", k)
    np.savetxt("sol_t/sol", sol.view(float))
#psi_of_t(dim=2, coulomb=True)
psi_of_t(dim=2.5, coulomb=True)



### calculate P(t) (rotating frame)

def pol_of_t(dim=3):
    t = np.loadtxt("sol_t/time")
    sol = np.loadtxt("sol_t/sol").view(complex)
    pol = np.zeros(len(t), dtype=np.complex_)
    if dim == 3:
        integrand = 1/(2*np.pi**2)*np.conjugate(P.dipole)*sol*k**2 #* np.exp(-1j*P.w0*t)
    if dim == 2:
        integrand = 1/(2*np.pi)*np.conjugate(P.dipole)*sol*k    
    #print(integrand.shape)
    for ii in range(len(t)):
        pol[ii] = Integration.integrater.int_disc(
            integrand[ii, :], grid)
    np.savetxt("sol_t/polt", pol.view(float))
pol_of_t(dim=2)




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
