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

def Efield(t):
    return P.E0 * np.exp(-1/2*((t-P.texp) / P.FWHM)**2) * 4 * np.log(2) 
    

#def interaction_integral(fx, grid, N):
#    x, g = grid
#    I=0
#    for ii in range(N):
#            I += fx[ii] * g[ii]
#    return I   

def int_grid(k, N):
    a = P.k0
    b = k
    c = P.k1
    ii = np.arange(0, N)
    # Häufungspunkt am Ende
    x1 = ( (ii+N+1)/(2*N+1) - 1 / ( 2 * np.pi ) * np.sin ( 2 * np.pi * (ii+N+1)/(2*N+1) ) ) * 2 * ( b-a ) + a - ( b-a )
    v1 = ( 1 - np.cos( 2 * np.pi * (ii+N+1)/(2*N+1) ) ) * ( 2 * ( b-a )/(2*N+1) )
    # Häufungspunkt am Anfang
    x2 = ( (ii+1)/(2*N+1) - 1 / ( 2 * np.pi ) * np.sin ( 2 * np.pi * (ii+1)/(2*N+1) ) ) * 2 * ( c-b ) + b
    v2 = ( 1 - np.cos( 2 * np.pi * (ii+1)/(2*N+1) ) ) * ( 2 * ( c-b )/(2*N+1) )
    # Insgesamt 
    x = np.concatenate((x1, x2))
    v = np.concatenate((v1, v2))
    #w = 
    return x, v
#x, w = int_grid(0.5, 20)
#plt.plot(x, x-x, 'x')

#        Ik = np.pi*k*psik 
#        #Ik = Ik + 
#        for ii in range(len(k)):
#            integrand = k/k[ii] * np.log(abs(k+k[ii])/abs(k-k[ii]))*(psik -
#                           psik[ii]*4*k/((k**2 + k[ii]**2)**2))
#            N = 200
#            Re = Integration.integrater.integrater(
#                    interpolate.interp1d(k, integrand.real),
#                                Integration.grids.gauß(P.k0, P.k1, N), N)
#            Im = Integration.integrater.integrater(
#                    interpolate.interp1d(k, integrand.imag),
#                                Integration.grids.gauß(P.k0, P.k1, N), N)
#            Ik[ii] += Re + 1j*Im 


def rhs(k):
    def ode(t, psik):
#        Ik = np.pi*k*psik 
#        #Ik = Ik + 
#        for ii in range(len(k)):
#            integrand = k/k[ii] * np.log(abs(k+k[ii])/abs(k-k[ii]))*(psik -
#                           psik[ii]*4*k/((k**2 + k[ii]**2)**2))
#            N = 200
#            Re = Integration.integrater.integrater(
#                    interpolate.interp1d(k, integrand.real),
#                                int_grid(k[ii], N), 2*N)
#            Im = Integration.integrater.integrater(
#                    interpolate.interp1d(k, integrand.imag),
#                                int_grid(k[ii], N), 2*N)
#            Ik[ii] += Re + 1j*Im 
        Dpsik = (P.Eg - P.w0 + epsk -1j*P.damp)*psik - P.dipole * Efield(t) #- P.c * Ik 
        Dpsik  = Dpsik * (-1j) / C.hbar 
        return Dpsik
    return ode


### calculate psi(t)


#k = np.linspace(P.k0, P.k1, P.nk)
k = Integration.grids.gts1(P.k0, P.k1, P.nk)[0]
epsk = C.hbar**2 * k**2 / (2 * P.mu)
psik0 = np.zeros(len(k), dtype=np.complex_)

def psiovertime():
    method = ODE.solvers.RK4(rhs(k)) 
    t, sol = ODE.main.solve(method, P.t0, P.t1, P.nt, psik0)
    np.savetxt("solutions/time", t)
    np.savetxt("solutions/momentum", k)
    np.savetxt("solutions/sol", sol.view(float))
psiovertime()


### calculate P(t)

#t = np.loadtxt("solutions/time")
#plt.plot(t, t-t, 'x')



def polovertime():
    t = np.loadtxt("solutions/time")
    sol = np.loadtxt("solutions/sol").view(complex)
    pol = np.zeros(P.nt+1, dtype=np.complex_)
    t = t.reshape((len(t), 1))
    integrand = 1/(2*np.pi**2)*(np.conjugate(P.dipole)*sol*np.exp(-1j*P.w0*t) +
                                    P.dipole * np.conjugate(sol*np.exp(-1j*P.w0*t)))
    #print(t.shape)
    #t = t.reshape(())
    #integrand = np.transpose(integrand)
    #print(integrand.shape)
    for ii in range(P.nt+1):
        pol[ii] = Integration.integrater.int_disc(
            integrand[ii], Integration.grids.gts1(P.k0, P.k1, P.nk))
    np.savetxt("solutions/polt", pol.view(float))
polovertime()


### calculate P(w)


def poloverfreq():
    polt = np.loadtxt("solutions/polt").view(complex)
    t = np.loadtxt("solutions/time")
    #print(polt.shape)
    nw = 1000
    w0 = 2000
    w1 = 3000
    w = Integration.grids.gts1(w0, w1, nw)[0]
    polw = np.zeros(nw, dtype=np.complex_)
    Ew = np.zeros(nw, dtype=np.complex_)
    for ii in range(nw):
        fpolt = interpolate.interp1d(t, polt*np.exp(-1j*w[ii]*t))
        fEt = lambda time: Efield(time) * np.cos(P.w0*time)  * np.exp(-1j*w[ii]*time) #  *np.exp(-1j*P.w0*time)
        #print(fpolt(0.1))
        polw[ii] = Integration.integrater.int_func(
            fpolt, Integration.grids.gts1(P.t0, P.t1, P.nt))
        Ew[ii] = Integration.integrater.int_func(
            fEt, Integration.grids.gts1(P.t0, P.t1, P.nt))
    np.savetxt("solutions/frequency", w)     
    np.savetxt("solutions/polw", polw.view(float))
    np.savetxt("solutions/Ew", Ew.view(float))
    
    
    #print(polw)
    
poloverfreq()    





et = time.time()
elapsed_time = et - st
print('Execution time:', elapsed_time, 'seconds')
