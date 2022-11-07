# -*- coding: utf-8 -*-
"""
Parameters for electron-phonon.scattering
"""
import constants as C
import numpy as np


# material spicific
me = 0.067 * C.me #0.067
mh = 0.15 * C.me
mu = 1/(1/me + 1/mh)
#mu = C.me
eps = 12  # dielectric constant
#eps = 1
dipole = 0.5 * C.e
Eg = 2500 # gap energy
damp = C.hbar*1 # damping of the polaization


# interraction 
cI_3D = 1 / (2*np.pi)**2 * C.e**2 / (C.eps0*eps) 
cI_2D = 1 / (2*np.pi)**2 * C.e**2 / (C.eps0*eps*2) 


# eletric field 
E0 = 1e-07  # field strength
w0 = 2500 / C.hbar
texp = 1
FWHM = 0.05 #0.06
Phi = C.hbar*w0 - Eg


# grids
k0 = 0.
k1 = 2.
nk = 100 #2000


# ODE conditions
t0          = 0.
t1          = 10
nt          = 20000
ntprint     = 100
tres = nt/(t1-t0) # 1 Stützstelle per ps


# fourier space
#w0 = -40
#w1 = 50
#nw = 1000

# general
debug = True

ryd = - mu * C.e**4 / (8*C.eps0**2*eps**2*C.hbar**2*4*np.pi**2)
ryd_frq = ryd / C.hbar



