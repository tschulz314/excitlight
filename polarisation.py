# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 19:00:37 2022

@author: derto
"""

import numpy as np

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
import constants as C

#hbar=1
#mu = 1

damp = 0
mu = 0.2*C.me
dipole = 1 


k = np.linspace(0, 2)
epsk = C.hbar**2 * k**2 / (2 * mu)
Eg = 0


def gauss(t, texp, sigma):
    return np.exp(-1/2*((t-texp) / sigma)**2) 



def rhs(k):
    def ode(t, psik):
        Dpsik = (Eg + epsk -1j*damp)*psik - dipole * gauss(t, 0.05, 0.015)
        #Dpsik = epsk*psik
        Dpsik  = Dpsik * (-1j) / C.hbar 
        #print(t)
        #Dpsik = np.cos(t)
        return Dpsik
    return ode




psik0 = np.zeros(len(k), dtype=np.complex_)

method = ODE.solvers.RK4(rhs(k)) 
t, sol = ODE.main.solve(method, 0, 0.4, 2000, psik0)

#print(sol[4, 0])

np.savetxt("solutions/time", t)
np.savetxt("solutions/momentum", k)
#np.savetxt("solutions/sol", sol)
np.savetxt("solutions/sol", sol.view(float))

#print(sol)







