# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:00:37 2022

@author: derto
"""

import numpy as np
from scipy import interpolate

import sys
import os

# get current directory
path = os.getcwd()
#print("Current Directory", path)
# prints parent directory
#print(os.path.abspath(os.path.join(path, os.pardir)))
#sys.path.append('C:\\Users\\derto\\OneDrive\\Dokumente\\Masterarbeit\\Numerik\\')
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))
import ODE
import Integration
import constants as C
import parameters as P


def Efield(t):
    return P.E0 * np.exp(-1/2*((t-P.texp) / P.FWHM)**2) * 4 * np.log(2) 
    

def rhs(k):
    def ode(t, psik):
        Dpsik = (P.Eg + epsk -1j*P.damp)*psik #- P.dipole * Efield(t)
        #Dpsik = epsk*psik
        Dpsik  = Dpsik * (-1j) / C.hbar 
        #print(t)
        #Dpsik = np.cos(t)
        return Dpsik
    return ode



k = np.linspace(P.k0, P.k1, P.nk)
epsk = C.hbar**2 * k**2 / (2 * P.mu)
psik0 = np.zeros(len(k), dtype=np.complex_)

method = ODE.solvers.RK4(rhs(k)) 
t, sol = ODE.main.solve(method, P.t0, P.t1, P.nt, psik0)

np.savetxt("solutions/time", t)
np.savetxt("solutions/momentum", k)
#np.savetxt("solutions/sol", sol)
np.savetxt("solutions/sol", sol.view(float))

#print(sol)

### calculate P(t)
pol = np.zeros(P.nt+1, dtype=np.complex_)
#print(len(pol))
for ii in range(P.nt+1):
    integrand = 1/(2*np.pi**2)*(np.conjugate(P.dipole)*sol[ii, :].real +
                           np.conjugate(sol[ii, :].imag)*P.dipole)
    N = 200
    Re = Integration.integrater.integrater(
        interpolate.interp1d(k, integrand.real),
                                Integration.grids.gauß(P.k0, P.k1, N), N)
    Im = Integration.integrater.integrater(
        interpolate.interp1d(k, integrand.imag),
                                Integration.grids.gauß(P.k0, P.k1, N), N)
    
    pol[ii] = Re + 1j*Im
np.savetxt("solutions/pol", pol.view(float))


#A= Integration.integrater.integrater(lambda x: np.sin(x), Integration.grids.gauß(0, np.pi, 200), 200)
#A = Integration.grids.gauß(P.k0, P.k1, 200)



