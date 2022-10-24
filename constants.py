# -*- coding: utf-8 -*-
"""
Constants
"""
import numpy as np
pi = np.pi

# Physikalische Konstanten
# Einheiten in [ nm, ps, pA, meV ]
c0        =   299792.458                                    # nm / ps
e         =   1.602176565e5                                 # pA * ps
kb        =   0.08617343                                    # meV / K
me        =   5.686e-3                                      # meV * ps^2 / nm^2
hbar      =   0.6582119514                                  # meV * ps
eps0      =   1.41844e6                                     # pA^2 * ps^2 / meV / nm
coul_pref =   e**2 / 2. / eps0                              # meV * nm
Ry        =   me * e**4 / 8. / eps0**2 / (hbar * 2. * pi)**2# 1 Rydberg in meV
a0        =   4. * pi * eps0 * hbar**2 / me / e**2          # Bohrradius in nm                  
mu0       =   1. / ( c0**2 * eps0 )                         # pV * ps / pA / pm
