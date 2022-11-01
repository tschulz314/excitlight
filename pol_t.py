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
    
def rhs(k):
    def ode(t, psik):
        Dpsik = (epsk - P.Phi - 1j*P.damp)*psik - P.dipole*misc.E0(t) - P.cI*int_integral3(psik)
        Dpsik  *= (-1j) / C.hbar 
        #print(t)
        return Dpsik
    return ode


def int_integral(psik):
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    v_ij = np.where(np.abs(ki-kj) > 0,
                    grid[1] * kj  / ki * np.log(np.abs(ki+kj)/np.abs(ki-kj)), 0) 
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=0)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+np.pi*ki)
    I = np.sum(w_ij*psik, axis=0)
    #print(I.shape)                
    return I 


def int_integral2(psik):
    I = np.zeros(len(k), dtype=np.complex_)
    ### calculate the grid v_(i, j)
    v_ij = np.zeros((len(k), len(k)))
    for ii in range(len(k)):  # loop over outer k grid
        for jj in range(len(k)):
            if ii == jj:
                v_ij[ii, jj] = 0
            else:
                v_ij[ii, jj] = grid[1][jj] * k[jj] / k[ii] * np.log(np.abs(k[jj]+k[ii])/np.abs(k[jj]-k[ii])) # v_(i,j)
    ### calculate the modified grid w_(i, j)            
    w_ij = np.zeros((len(k), len(k)))
    for ii in range(len(k)):  # loop over outer k grid
        for jj in range(len(k)):
            if ii == jj:
                #v_ij[ii, jj] = 0
                A = 0
                for mm in range(len(k)):
                    A += v_ij[ii, mm]
                w_ij[ii, jj] = -A + np.pi*k[ii]    
            else:
                w_ij[ii, jj] = v_ij[ii, jj]
    
    for ii in range(len(k)):            
        I[ii] = np.sum(w_ij[ii, :]*psik) # sum over j       
    #print(I) 
    #print(v_ij)
    return I     


def int_integral3(psik):
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    integ = kj / ki * psik * np.log(np.abs(ki+kj)/np.abs(ki-kj+1e-10))
    I = np.zeros(len(k), dtype=np.complex_)
    for ii in range(len(k)):
    #print(ii)
        I[ii] = Integration.integrater.int_disc(integ[ii, :], grid)
    return I


#psi = np.ones(len(k))
#Ik1 = int_integral(psi)
#Ik2 = int_integral2(psi)
#Ik3 = int_integral3(psi)
#
#
#psi = np.ones(len(k))
#plt.plot(k, Ik1)
#plt.plot(k, Ik2)
#plt.plot(k, Ik3)

    
#ki = k.reshape((len(k), 1))
#kj = k.reshape((1, len(k)))
#v_ij = np.where(np.abs(ki-kj) > 0,
#                grid[1] * kj / ki * np.log(np.abs(ki+kj)/np.abs(ki-kj)), 0) 
#v_im = v_ij * 4 * ki**2 / (ki**2+kj**2)**2
#sum_i = np.sum(v_im, axis=0)

#x = np.linspace(0, 5, 6)
#x1 = x.reshape((len(x), 1))
#x2 = x.reshape((1, len(x)))
##u = np.where(np.abs(x1-x2) == 0, 0, 1/np.abs((x1+x2)))
#u = np.where(np.abs(x1-x2) > 0, (x1-x2), 0)
#v = np.sum(u, axis=0)
#print(u)
#print(v)

def psiovertime():
    method = ODE.solvers.RK4(rhs(k)) 
    t, sol = ODE.main.solve(method, P.t0, P.t1, P.nt, psik0)
    #sol = sol*np.exp( -1j * P.w0*t.reshape((len(t), 1)) )
    #print(len(t), sol.shape)
    #t = t[::10]
    #sol = sol[::10, :] #int(P.nt/P.ntprint)
    #print(sol)
    #print(len(t), sol.shape)
    np.savetxt("sol_t/time", t)
    np.savetxt("sol_t/momentum", k)
    np.savetxt("sol_t/sol", sol.view(float))
psiovertime()



### calculate P(t) (rotating frame)

def polovertime():
    t = np.loadtxt("sol_t/time")
    sol = np.loadtxt("sol_t/sol").view(complex)
    pol = np.zeros(len(t), dtype=np.complex_)
    integrand = 1/(2*np.pi**2)*np.conjugate(P.dipole)*sol*k**2 #* np.exp(-1j*P.w0*t)
    #print(integrand.shape)
    for ii in range(len(t)):
        pol[ii] = Integration.integrater.int_disc(
            integrand[ii, :], grid)
    np.savetxt("sol_t/polt", pol.view(float))
polovertime()



### calculate P(w) (rotating frame)
    
def poloverfreq():
    polt = np.loadtxt("sol_t/polt").view(complex)
    t = np.loadtxt("sol_t/time")
    nw = 2000
    w0 = -30
    w1 = 50
    #w, Ew = misc.fourier_trafo(t, misc.E0(t), w0, w1, nw)#, inverse=True)
    w, polw = misc.fourier_trafo(t, polt, w0, w1, nw, inverse=True)
    np.savetxt("sol_t/frequency", w)    
    #np.savetxt("solutions/Ew", Ew.view(float))
    np.savetxt("sol_t/polw", polw.view(float)) 
poloverfreq()    





   



et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
