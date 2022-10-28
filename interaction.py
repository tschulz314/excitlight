#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 12:42:44 2022

@author: tschulz
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
    

def int_integral(k):
    N = 1000
    x, g = int_grid(k, round(N/2))
    U = x, g
    func = lambda x2: x2/k * np.log(np.abs(k + x2)/np.abs(k - x2))    
    v = g * x/k * np.log(np.abs(x + k)/np.abs(x - k))
    #w = np.zeros((N, N))
    w = v
    delta = 0.001
    for ii in range(N):
        if abs(x[ii]-k) < delta: 
            print(x[ii])
            w[ii] = -np.sum(v) + v[ii] + Integration.integrater.integrater2(func, U, N)        
    #return x, w
    #print(func(x))
    return np.sum(w)

x = int_integral(1)
#np.integrate()    
    
    
    
#    Integration.integrater.integrater(
#                    interpolate.interp1d(k, integrand.real),
#                                int_grid(k[ii], N), 2*N)
    

    
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



#def rhs(k):
#    def ode(t, psik):
##        Ik = np.pi*k*psik 
##        #Ik = Ik + 
##        for ii in range(len(k)):
##            integrand = k/k[ii] * np.log(abs(k+k[ii])/abs(k-k[ii]))*(psik -
##                           psik[ii]*4*k/((k**2 + k[ii]**2)**2))
##            N = 200
##            Re = Integration.integrater.integrater(
##                    interpolate.interp1d(k, integrand.real),
##                                int_grid(k[ii], N), 2*N)
##            Im = Integration.integrater.integrater(
##                    interpolate.interp1d(k, integrand.imag),
##                                int_grid(k[ii], N), 2*N)
##            Ik[ii] += Re + 1j*Im 
#        Dpsik = (P.Eg - P.w0 + epsk -1j*P.damp)*psik - P.dipole * Efield(t) #- P.c * Ik 
#        Dpsik  = Dpsik * (-1j) / C.hbar 
#        return Dpsik
#    return ode    
    