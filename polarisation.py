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
        Ik = np.pi*k*psik 
        #Ik = Ik + 
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
        Dpsik = (P.Eg + epsk -1j*P.damp)*psik - P.dipole * Efield(t) #- P.c * Ik 
        Dpsik  = Dpsik * (-1j) / C.hbar 
        return Dpsik
    return ode



k = np.linspace(P.k0, P.k1, P.nk)
epsk = C.hbar**2 * k**2 / (2 * P.mu)
psik0 = np.zeros(len(k), dtype=np.complex_)

method = ODE.solvers.RK4(rhs(k)) 
t, sol = ODE.main.solve(method, P.t0, P.t1, P.nt, psik0)

np.savetxt("solutions/time", t)
np.savetxt("solutions/momentum", k)
np.savetxt("solutions/sol", sol.view(float))

#print(sol)


### calculate P(t)
pol = np.zeros(P.nt+1, dtype=np.complex_)
for ii in range(P.nt+1):
    integrand = 1/(2*np.pi**2)*(np.conjugate(P.dipole)*sol[ii, :] + 
                   P.dipole * np.conjugate(sol[ii, :]))
    N = 200
    Re = Integration.integrater.integrater(
        interpolate.interp1d(k, integrand.real),
                                Integration.grids.gauß(P.k0, P.k1, N), N)
    Im = Integration.integrater.integrater(
        interpolate.interp1d(k, integrand.imag),
                                Integration.grids.gauß(P.k0, P.k1, N), N)
    
    pol[ii] = Re + 1j*Im
np.savetxt("solutions/pol", pol.view(float))


#A= Integration.integrater.integrater(lambda x: np.sin(x), Integration.grids.gauß(0, np.pi, 200), 200)
#A = Integration.grids.gauß(P.k0, P.k1, 200)



