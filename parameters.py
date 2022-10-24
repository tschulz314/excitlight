# -*- coding: utf-8 -*-
"""
Parameters for electron-phonon.scattering
"""
import constants as C
import numpy as np
pi = np.pi

# grids
k0 = 0.
k1 = 3.
nk = 400

qopt0 = 0.
qopt1 = 3.
nqopt = 400

qac0 = 0.
qac1 = 3.
nqac = 400

phi0 = 0.
phi1 = 2. * pi
nphi = 400

# phonons
# MoS2:
c_acoustic      = 6.6       # speed of longitudinal acoustic phonons in nm*ps^-1 (MoS2, taken from PRB 90, 045422)
rho_ion         = 1.92e4    # ion mass density in meV*ps^2*nm^-4 (MoS2)
D_1_acoustic_e  = 4.5e3     # deformation potential for longitudinal acoustic phonons in meV (MoS2)
D_0_optical_e   = 5.8e4     # deformation potential for longitudinal optical  phonons in meV/nm (MoS2)
#D_1_acoustic_h  = 2.5e3     # deformation potential for longitudinal acoustic phonons in meV (MoS2)
#D_0_optical_h   = 4.6e4     # deformation potential for longitudinal optical  phonons in meV/nm (MoS2)
E_LO            =  48.9 #49.       # LO phonon energy in meV

# electron-phonon
D_1_acoustic    = D_1_acoustic_e
D_0_optical     = D_0_optical_e

# hole-phonon
#D_1_acoustic    = D_1_acoustic_h
#D_0_optical     = D_0_optical_h

# electrons
meff = 0.51 * C.me #0.41 #0.51 (MoS2)
K = C.hbar**2  / (2 * meff)

# holes
# meff = 0.46 * C.me

# ODE conditions
t0          = 0.
t1          = 10. 
nt          = 800
dt          = (t1-t0)/nt
dens_start  = 1e12 # cm^-2
n = dens_start * 1E-14 # nm^-2
temp        = 300.
beta        = (C.kb * temp)**-1

# general
debug = True


# Weitere
a = E_LO # = hbar * Omega_opt = const
b = C.hbar * c_acoustic  # = hbar * omega_acc / q = hbar * c_ac
tau_opt = 100
tau_ac = 1
q_LO = np.sqrt(a / K)
q_LA = b / K


