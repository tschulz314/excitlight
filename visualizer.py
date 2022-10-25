# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:40:27 2021

@author: derto
"""



import numpy as np
import matplotlib.pyplot as plt

import parameters as P


def read_data():
    t = np.loadtxt("solutions/time")
    k = np.loadtxt("solutions/momentum")
    sol = np.loadtxt("solutions/sol").view(complex) #, dtype=np.complex_
    pol = np.loadtxt("solutions/pol").view(complex)
    return t, k, sol, pol
t, k, sol, pol = read_data()


def time_evo():
    T = len(t)-1
    t_ind = [0, int(T/100), int(T/10), int(T/5), int(T/2), int(T)]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True) 
    #plt.figure() #figsize=(7, 7) dpi=300
    for index in t_ind:
        ax1.plot(k, sol[index, :].real, '-', label=r"$t={}\,$ps".format(t[index]))
        ax2.plot(k, sol[index, :].imag, '-', label=r"$t={}\,$ps".format(t[index]))
    ax1.set_xlabel(r"$k$ in nm$^{-1}$")
    ax2.set_xlabel(r"$k$ in nm$^{-1}$")
    ax1.set_ylabel(r"$Re(\psi_k)$")
    ax2.set_ylabel(r"$Im(\psi_k)$")
    #fig.title('Anfangsverteilung')
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    ax1.grid()
    ax2.grid()
    #plt.savefig("fig/fk_" + branch + '.png', dpi=300)
    #plt.savefig("fig/fk_0.png", dpi=300)
    plt.show()
#time_evo() 


def Efield(t):
    return P.E0 * np.exp(-1/2*((t-P.texp) / P.FWHM)**2) * 4 * np.log(2) 


def psit():
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, dpi=300) 
    K = len(k) - 1
    #print(K)
    ktoplot = [round(K/2), round(K*3/4)]
    for index in ktoplot:
        ax1.plot(t, sol[:, index].real, '-')
        ax2.plot(t, sol[:, index].imag, '-', label=r"$k={}\,$nm$^{{-1}}$".format(round(k[index], 1)))
    ax3.plot(t, Efield(t))    
    ax3.set_xlabel(r"$t$ in ps")
    #ax3.set_xlim(0, 0.1)
    ax1.set_ylabel(r"$Re(\psi)$")
    ax2.set_ylabel(r"$Im(\psi)$")
    ax3.set_ylabel(r"$E$")
    #fig.title('Anfangsverteilung')
    ax2.legend(loc="right") #bbox_to_anchor=(1.1, 1.05)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()
psit() 


def pol_evo():
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, dpi=300) 
    ax1.plot(t, pol[:].real, '-')
    ax2.plot(t, pol[:].imag, '-')
    ax3.plot(t, Efield(t))  
    #ax3.set_ylim(0, 100)
    ax3.set_xlabel(r"$t$ in ps")
    #ax3.set_xlim(0, 0.1)
    ax1.set_ylabel(r"$Re(P)$")
    ax2.set_ylabel(r"$Im(P)$")
    ax3.set_ylabel(r"$E$")
    ax2.legend(loc="right") 
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()
pol_evo() 
