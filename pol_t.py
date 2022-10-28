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

k = Integration.grids.gts2(P.k0, P.k1, P.nk)[0]
epsk = C.hbar**2 * k**2 / (2*P.mu)
psik0 = np.zeros(len(k), dtype=np.complex_)



### calculate psi(t) (rotating frame)
    
def rhs(k):
    def ode(t, psik):
        Dpsik = (epsk - P.Phi - 1j*P.damp)*psik - P.dipole*misc.E0(t)
        Dpsik  *= (-1j) / C.hbar 
        return Dpsik
    return ode


def psiovertime():
    method = ODE.solvers.RK4(rhs(k)) 
    t, sol = ODE.main.solve(method, P.t0, P.t1, P.nt, psik0)
    #sol = sol*np.exp( -1j * P.w0*t.reshape((len(t), 1)) )
    np.savetxt("solutions/time", t)
    np.savetxt("solutions/momentum", k)
    np.savetxt("solutions/sol", sol.view(float))
psiovertime()



### calculate P(t) (rotating frame)

def polovertime():
    t = np.loadtxt("solutions/time")
    sol = np.loadtxt("solutions/sol").view(complex)
    pol = np.zeros(len(t), dtype=np.complex_)
    integrand = 1/(2*np.pi**2)*np.conjugate(P.dipole)*sol*k**2 #* np.exp(-1j*P.w0*t)
    #print(integrand.shape)
    for ii in range(len(t)):
        pol[ii] = Integration.integrater.int_disc(
            integrand[ii, :], Integration.grids.gts2(P.k0, P.k1, P.nk))
    np.savetxt("solutions/polt", pol.view(float))
polovertime()



### calculate P(w) (rotating frame)
    
def poloverfreq():
    polt = np.loadtxt("solutions/polt").view(complex)
    t = np.loadtxt("solutions/time")
    nw = 8000
    w0 = -2000
    w1 = 2000
    #w, Ew = misc.fourier_trafo(t, misc.E0(t), w0, w1, nw)#, inverse=True)
    w, polw = misc.fourier_trafo(t, polt, w0, w1, nw) #, inverse=True)
    np.savetxt("solutions/frequency", w)    
    #np.savetxt("solutions/Ew", Ew.view(float))
    np.savetxt("solutions/polw", polw.view(float)) 
poloverfreq()    


   



et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
