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

### elecrtic field



def E0(t):
    return P.E0 * np.exp(-1/2*((t-P.texp) / P.FWHM)**2) * 2*np.log(2)


def E(t):
    return E0(t) * np.exp(-1j*P.w0*t)


def E0w(w):
    return P.E0 * 2*np.log(2) * np.exp(-1/2*(w*P.FWHM)**2) * P.FWHM * np.exp(-1j*w*P.texp)



def fourier_trafo(x, fx, y0, y1, ny, inverse=False):
    f = interpolate.interp1d(x, fx)
    #plt.plot(x, f(x))
    #plt.show()
    #plt.close()
    y = np.linspace(y0, y1, ny)
    gy = np.zeros(len(y), dtype=np.complex_)
    grid = Integration.grids.equi(x[0], x[-1], len(x))
    x = grid[0]
    for ii in range(ny):
        func = lambda xx: 1/np.sqrt(2*np.pi) * f(xx) * np.exp(-1j*y[ii]*xx) 
        if inverse is True:
            func = lambda x: 1/np.sqrt(2*np.pi) * f(x) * np.exp(1j*y[ii]*x) 
        gy[ii] = Integration.integrater.int_disc(func(x), grid) 
    #print(x[0], x[-1], len(x))    
    return y, gy


def fourier_trafo2(x, fx, y, inverse=False):
    f = interpolate.interp1d(x, fx)
    #plt.plot(x, f(x))
    #plt.show()
    #plt.close()
    #y = np.linspace(y0, y1, ny)
    gy = np.zeros(len(y), dtype=np.complex_)
    grid = Integration.grids.equi(x[0], x[-1], len(x))
    x = grid[0]
    for ii in range(len(y)):
        func = lambda xx: 1/np.sqrt(2*np.pi) * f(xx) * np.exp(-1j*y[ii]*xx) 
        if inverse is True:
            func = lambda x: 1/np.sqrt(2*np.pi) * f(x) * np.exp(1j*y[ii]*x) 
        gy[ii] = Integration.integrater.int_disc(func(x), grid) 
    #print(x[0], x[-1], len(x))    
    return y, gy

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

#w0 = -100
#w1 = 100
#
#t = np.linspace(P.t0, P.t1, P.nt)
#w, gw = fourier_trafo(t, E0(t), w0, w1, P.nt) #, inverse=False)
##gw = gw * np.exp(1j*w*texp)
#plt.plot(t, E0(t))
#plt.show()
#plt.close()
#plt.plot(w, gw.real)
##plt.plot(w, gw.imag)
#plt.plot(w, E0w(w).real, '--')
##plt.plot(w, E0w(w).imag, '--')
#plt.show()




