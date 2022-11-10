# -*- coding: utf-8 -*-
"""
Parameters for electron-phonon.scattering
"""
import constants as C
import numpy as np


# material spicific
#me = 0.067 * C.me 
#mh = 0.15 * C.me
me = 0.47 * C.me # MoS2
mh = 0.54 * C.me # MoS2

mu = 1/(1/me + 1/mh)

eps = 12  # dielectric constant
dipole = 0.5 * C.e 
Eg = 2500 # gap energy
damp = C.hbar*1 # damping of the polaization


# eletric field 
E0 = 1e-07  # field strength
w0 = 2500 / C.hbar
texp = 1
FWHM = 0.001 # 0.005
Phi = C.hbar*w0 - Eg


# grids
k0 = 0.
k1 = 10.
nk = 100 


# ODE conditions
t0          = 0.
t1          = 10
nt          = 20000
#ntprint     = 100
tres = nt/(t1-t0) # 1 Stützstelle per ps


# fourier space
#w0 = -40
#w1 = 50
#nw = 1000

# general
debug = True

ryd = - mu * C.e**4 / (8*C.eps0**2*eps**2*C.hbar**2*4*np.pi**2)
ryd_frq = ryd / C.hbar

ryd_TMDC = - mu * C.e**4 / (8*C.eps0**2*C.hbar**2*4*np.pi**2)
ryd_frq_TMDC = ryd_TMDC / C.hbar



