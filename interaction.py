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
    




def int_integral(psik):
    I = np.zeros(len(k), dtype=np.complex_)
    # calculate the modified grid w_(i, j)
    w_ij = np.zeros((len(k), len(k)))
    for ii in range(len(k)):  # loop over outer k grid
        v_ij = grid[1] * k / k[ii] * np.log(np.abs(k+k[ii])/np.abs(k-k[ii])) # v_(i,j)
        v_im = v_ij * 4 * k[ii]**4 / (k[ii]**2+k**2)**2 # v_(i, m) in sum for i = j
        vsum = np.sum(v_im)
        for jj in range(len(k)):
            if ii == jj:
                w_ij[ii, jj] = -vsum + v_im[ii] + np.pi*k[ii]  
            else:
                w_ij[ii, jj] = v_ij[jj]
        I[ii] = np.sum(w_ij[ii, :]*psik) # sum over j       
    #print(I)    
    return I            
#x = int_integral(1)







def int_integral2(psik):
    I = np.zeros(len(k), dtype=np.complex_)
    #print(psik)
    ### calculate the grid v_(i, j)
    v_ij = np.zeros((len(k), len(k)))
    for ii in range(len(k)):  # loop over outer k grid
        for jj in range(len(k)):
            if ii == jj:
                v_ij[ii, jj] = None
            else:
                v_ij[ii, jj] = grid[1][jj] * k[jj] / k[ii] * np.log(np.abs(k[jj]+k[ii])/np.abs(k[jj]-k[ii])) # v_(i,j)
                #v_ij[ii, jj] *= 4*k[ii]**4/(k[ii]**2 + k[jj]**2)**2   
    #print(v_ij)         
    ### calculate the modified grid w_(i, j)            
    w_ij = np.zeros((len(k), len(k)))
    for ii in range(len(k)):  # loop over outer k grid
        for jj in range(len(k)):
            if ii == jj:
                #print(ii, jj)
                #v_ij[ii, jj] = 0
                A = 0
                for mm in range(len(k)):
                    if mm != jj:
                        A += v_ij[ii, mm]*4*k[ii]**4/(k[ii]**2 + k[mm]**2)**2 
                    
                w_ij[ii, jj] = -A + np.pi*k[ii]    
            else:
                #print(ii, jj)
                w_ij[ii, jj] = v_ij[ii, jj] #*4*k[ii]**4/(k[ii]**2 + k[mm]**2)**2 
    #I = np.sum(w_ij*psik, axis=1)
    #for ii in range(len(k)):   
        #I[ii] = Integration.integrater.int_disc(psik[:], (k, w_ij[ii, :]))         
        #I[ii] = np.sum(w_ij[ii, :]*psik[:]) # sum over j       
    #print(I) 
    #print(w_ij)
    #print(A)
    return I     




def int_integral3(psik):
    #print(psik)
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    integ = kj / ki * psik * np.log(np.abs(ki+kj)/np.abs(ki-kj+1e-12))
    I = np.zeros(len(k), dtype=np.complex_)
    for ii in range(len(k)):
    #print(ii)
        I[ii] = Integration.integrater.int_disc(integ[ii, :], grid)
    return I