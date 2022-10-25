# -*- coding: utf-8 -*-
"""
Parameters for electron-phonon.scattering
"""
import constants as C
import numpy as np


# material spicific
me = 0.067 * C.me 
mh = 0.015 * C.me
mu = 1/(1/me + 1/mh)
eps = 12  # dielectric constant
dipole = 0.5 * C.e
Eg = 0 # gap energy
damp = 1 # damping of the polaization

# interraction constant
c = C.e**2/((2*np.pi)**3 * 2 * C.eps0 * eps) 

# eletric pulse
E0 = 1e-07  # field strength
texp = 0.1
FWHM = 0.015


# grids
k0 = 0.
k1 = 2.
nk = 100



# ODE conditions
t0          = 0.
t1          = .2
nt          = 4000
#dt          = (t1-t0)/nt
#dens_start  = 1e12 # cm^-2
#n = dens_start * 1E-14 # nm^-2
#temp        = 300.
#beta        = (C.kb * temp)**-1

# general
debug = True




