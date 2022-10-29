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
eps = 12  # dielectric constant
dipole = 0.5 * C.e
Eg = 2500 # gap energy
damp = C.hbar*1 # damping of the polaization


# interraction 
c = C.e**2/((2*np.pi)**3 * 2 * C.eps0 * eps) 


# eletric field 
E0 = 1e-07  # field strength
w0 = 2500 / C.hbar
texp = 1
FWHM = 0.05 #0.06
Phi = C.hbar*w0 - Eg


# grids
k0 = 0.
k1 = 1.
nk = 2000


# ODE conditions
t0          = 0.
t1          = 10
nt          = 5000
tres = nt/(t1-t0) # 1 Stützstelle per ps

# general
debug = True




