#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 16:09:17 2022

@author: tschulz
"""

import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt

import sys
import os
path = os.getcwd()
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))
import Integration
import constants as C
import parameters as P




def k_int_grid(k, k0, k1, N):
    a = k0
    b = k
    c = k1
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


### find value in array

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


### elecrtic fields

def E0(t):
    return P.E0 * np.exp(-1/2*((t-P.texp) / P.FWHM)**2) * 2*np.log(2)


def E(t):
    return E0(t) * np.exp(-1j*P.w0*t)


def E0w(w): # t-> w: exp(iwt), i.e. inverse
    return P.E0 * 2*np.log(2) * np.exp(-1/2*(w*P.FWHM)**2) * P.FWHM * np.exp(1j*w*P.texp)


### fourier transform

def fourier_trafo(x, fx, y, inverse=False):
    f = interpolate.interp1d(x, fx)
    grid = Integration.grids.equi(x[0], x[-1], len(x))
    x = grid[0]
    gy = np.zeros(len(y), dtype=np.complex_)
    for ii in range(len(y)):
        if inverse is True:
            func = lambda x: 1/np.sqrt(2*np.pi) * f(x) * np.exp(1j*y[ii]*x) 
        else:
            func = lambda x: 1/np.sqrt(2*np.pi) * f(x) * np.exp(-1j*y[ii]*x) 
        gy[ii] = Integration.integrater.int_disc(func(x), grid)   
    return gy


### test 1

#x = np.linspace(-10, 10, 10000)
#y, gy = fourier_trafo(x, 1/(x**2+1), -10, 10, 10000)
#plt.plot(y, gy.real)
#plt.plot(y, np.pi*np.exp(-np.abs(y))*1/np.sqrt(2*np.pi), '--')
#x2, fx2 = fourier_trafo(y, gy, -10, 10, 1000, inverse=True)
#plt.plot(x, 1/(x**2+1))
#plt.plot(x2, fx2.real, '--')


### test 2

#def gauss(t, texp, FWHM):
#    return  np.exp(-1/2*((t-texp) / FWHM)**2) * 4 * np.log(2) * 1/2 # * 1 / np.sqrt(2*np.pi*P.FWHM**2) #P.E0 * P.E0 *

#t0 = 0
#t1 = 5
#nt = 1000
#
#texp = 1
#FWHM = 0.06

#w0 = -1000
#w1 = 1000
#
#t = np.linspace(t0, t1, nt)
#w, gw = fourier_trafo(t, gauss(t, texp, FWHM), w0, w1, nt) #, inverse=False)
#gw = gw * np.exp(1j*w*texp)
#plt.plot(t, gauss(t, texp, FWHM))
#plt.show()
#plt.close()
#plt.plot(w, gw.real)
#plt.plot(w, gw.imag)
#plt.plot(w, gauss(w, 0, 1/FWHM)*FWHM, '--')
##plt.plot(w, E0w(w), '--')
#plt.show()


###test 3

# w0 = -100
# w1 = 100

# t = np.linspace(P.t0, P.t1, P.nt)
# w, gw = fourier_trafo(t, E0(t), w0, w1, P.nt, inverse=True)
# #gw = gw * np.exp(1j*w*texp)
# plt.plot(t, E0(t))
# plt.show()
# plt.close()
# plt.plot(w, gw.real)
# #plt.plot(w, gw.imag)
# plt.plot(w, E0w(w).real, '--')
# #plt.plot(w, E0w(w).imag, '--')
# plt.show()







