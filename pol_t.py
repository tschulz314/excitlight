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
import interaction as inter

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
    nkjj = 1000
    nphi = 1000
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], k[0], 4*k[-1], nkjj)
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



def TMDC_pot(q):
    # def eps_tilde
    # eps_inf =  2
    # chi =   7.112 * 0.1 #0.66  
    # eps_up = eps_inf
    # eps_down = eps_inf
    # eps_tilde = (eps_up+eps_down)/2 + 2*np.pi*chi*q   # chi/(2*C.eps0)*q*C.eps0*4*np.pi
    # V = C.e**2 / (2 * C.eps0 * q) 
    # eps_HS_inv = 1 / eps_tilde
    # W1 = V * eps_HS_inv 

    #d = 6.18
    #eps = 15.46
    #chi = d*(eps-1)/(4*np.pi) * 0.1
    #print(chi)
    
    chi =   7.112 * 0.1 # 7.112 * e-10
    kappa = 2 
    r0 = 2 * np.pi * chi / kappa
    V = C.e**2 / (2 * C.eps0 * q) 
    W = V * 1 / kappa * 1 / (1 + r0 * q)
    
    # plt.plot(q, W)
    # #plt.plot(q, W1, '--')
    # plt.plot(q, V, '--')
    # plt.ylim(0, 1e5)
    # plt.show()
    
    return  W
    #return  V / P.eps
#TMDC_pot(k)


def int_grid_eff():
    ### calucalte the double integral 
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
    #print(w_ij)
    return w_ij


def rhs(k, dim=3, coulomb=True):
    if coulomb is True:
        if dim == 3:
            w_ij = int_grid_3D()
            cI = 1 / (2*np.pi)**2 * C.e**2 / (C.eps0*P.eps)   
        if dim == 2:
            w_ij = int_grid_2D()
            cI = 1 / (2*np.pi)**2 * C.e**2 / (C.eps0*P.eps*2) 
        else:
            w_ij = inter.int_grid_eff(k, grid)
            cI = 1 / (2*np.pi)**2 
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
    #print(sol.shape)
    np.savetxt("sol_t/time", t)
    np.savetxt("sol_t/momentum", k)
    np.savetxt("sol_t/sol", sol.view(float))


### calculate P(t) (rotating frame)

def pol_of_t(dim=3):
    t = np.loadtxt("sol_t/time")
    sol = np.loadtxt("sol_t/sol").view(complex)
    #print(sol.shape)
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


### calculate P(w) (rotating frame)
    
def pol_of_w(): # t-> w: exp(iwt), i.e. inverse
    polt = np.loadtxt("sol_t/polt").view(complex)
    t = np.loadtxt("sol_t/time")
    #w = np.linspace(-40, 50, 2000)
    w = Integration.grids.gts2(-600, 50, 10000)[0]
    polw = misc.fourier_trafo(t, polt, w, inverse=True)
    np.savetxt("sol_t/frequency", w)    
    np.savetxt("sol_t/polw", polw.view(float)) 


### execute programm 

#psi_of_t(dim=2, coulomb=True)
psi_of_t(dim='eff', coulomb=True)
pol_of_t(dim=2)
pol_of_w()
     


### calculate the wavefunctions
    
#def psi_of_w(): # t-> w: exp(iwt), i.e. inverse
#    polt = np.loadtxt("sol_t/polt").view(complex)
#    t = np.loadtxt("sol_t/time")
#    #w = np.linspace(-40, 50, 2000)
#    w = Integration.grids.gts2(-600, 50, 10000)[0]
#    polw = misc.fourier_trafo(t, polt, w, inverse=True)
#    np.savetxt("sol_t/frequency", w)    
#    np.savetxt("sol_t/polw", polw.view(float))     
    

   



et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
