# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:40:27 2021

@author: derto
"""

import numpy as np
import matplotlib.pyplot as plt


import constants as C
import parameters as P
import misc as misc


### loading data 
def load_all():
    t = np.loadtxt("sol_t/time")
    k = np.loadtxt("sol_t/momentum")
    sol = np.loadtxt("sol_t/sol").view(complex) 
    polt = np.loadtxt("sol_t/polt").view(complex)
    w = np.loadtxt("sol_t/frequency")
    polw = np.loadtxt("sol_t/polw").view(complex)
    Ew = misc.E0w(w)
#load_all()
    

def psi_of_k():
    """
    Plot psi of k for various times t
    """
    t = np.loadtxt("sol_t/time")
    k = np.loadtxt("sol_t/momentum")
    sol = np.loadtxt("sol_t/sol").view(complex) 
    T = len(t)-1
    #t_ind = [0, int(T/100), int(T/10), int(T/5), int(T/2), int(T)]
    t_ind = [int(T/100), int(T/10)]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=300) 
    for index in t_ind:
        ax1.plot(k, sol[index, :].real, '-', label=r"$t={}\,$ps".format(t[index]))
        ax2.plot(k, sol[index, :].imag, '-', label=r"$t={}\,$ps".format(t[index]))
    ax2.set_xlabel(r"$k$ in nm$^{-1}$")
    ax1.set_ylabel(r"$Re(\tilde{\psi}_k)$")
    ax2.set_ylabel(r"$Im(\tilde{\psi}_k)$")
    plt.legend() #bbox_to_anchor=(1.1, 1.05)
    ax1.grid()
    ax2.grid()
    plt.show()
#psi_of_k() 


def spicific_psik(ktoplot):
    """
    Plot psi of t for a specific k 
    """
    t = np.loadtxt("sol_t/time")
    k = np.loadtxt("sol_t/momentum")
    sol = np.loadtxt("sol_t/sol").view(complex) 
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=300) 
    plt.tight_layout()
    kind = misc.find_nearest(k, ktoplot)
    ax1.plot(t, sol[:, kind].real, '-',
             label=r"Re, $k={}\,$nm$^{{-1}}$".format(round(k[kind], 2)))
    ax1.plot(t, sol[:, kind].imag, '--',
             label=r"Im, $k={}\,$nm$^{{-1}}$".format(round(k[kind], 2))) 
    ax2.plot(t, misc.E0(t)) 
    ax2.set_xlabel(r"$t$ in ps")
    ax1.set_ylabel(r"$\tilde{\psi}_k(t)$")
    #ax1.set_xlim(6.6, 6.7)
    ax2.set_ylabel(r"$E_0(t)$")
    ax1.legend()
    ax1.grid()
    ax2.grid()
    plt.show()
#spicific_psik(0) 
#spicific_psik(0.4) 
#spicific_psik(1) 

 
def pol_of_t():
    """
    Plot Polaization P over time t
    """
    t = np.loadtxt("sol_t/time")
    polt = np.loadtxt("sol_t/polt").view(complex)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, dpi=300) 
    ax1.plot(t, polt[:].real, '-')
    ax2.plot(t, polt[:].imag, '-')
    ax3.plot(t, misc.E0(t))  
    #ax3.plot(t, misc.E(t))  
    #ax3.set_ylim(0, 100)
    ax3.set_xlabel(r"$t$ in ps")
    #ax3.set_xlim(0, 2)
    ax1.set_ylabel(r"$Re(P)$")
    ax2.set_ylabel(r"$Im(P)$")
    ax3.set_ylabel(r"$E_0$")
    #ax2.legend(loc="right") 
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()
#pol_of_t() 


def all_over_w():
    """
    Plot Polaization P, electric Field E and the suszebility chi over
    frequency w
    """
    w = np.loadtxt("sol_t/frequency")
    polw = np.loadtxt("sol_t/polw").view(complex)
    Ew = misc.E0w(w)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, dpi=300) 
    #plt.tight_layout()
    chi = polw / (Ew*C.eps0)  
    ax1.plot(w, polw.real, '-')
    ax1.plot(w, polw.imag, '-')
    ax2.plot(w, Ew.real, '-')
    ax2.plot(w, Ew.imag, '-')
    ax3.plot(w, chi.real, '-')
    ax3.plot(w, chi.imag, '-')
    #ax3.set_xlim(0, 100)
    ax3.set_xlabel(r"$\omega$ in ps$^{-1}$")
    #ax3.set_xlim(-100, 100)
    #ax3.set_ylim(-1e2, 1e2)
    ax1.set_ylabel(r"$P(\omega)$")
    ax2.set_ylabel(r"$E(\omega)$")
    ax3.set_ylabel(r"$\chi(\omega)$")
    #ax2.set_ylabel(r"$Im(\chi)$")
    #ax2.legend(loc="right") 
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()
#all_over_w()    
    
    
def abs_of_w():
    """
    Plot the absorption obver the frequency w
    """
    w = np.loadtxt("sol_t/frequency")
    polw = np.loadtxt("sol_t/polw").view(complex)
    Ew = misc.E0w(w)
    plt.figure(dpi=300) 
    chi = polw / (Ew*C.eps0) 
    #plt.plot(w, chi.real, '-')
    plt.plot(w, chi.imag, '-')
    print(w[np.argmax(chi.imag)], C.hbar*w[np.argmax(chi.imag)])

    # plt.axvline(P.ryd_frq, color='green', linewidth=1,
    #             label=r'$\omega_{\mathrm{Ryd}}$')
    # plt.axvline(4*P.ryd_frq, color='red', linewidth=1,
    #             label=r'$4\omega_{\mathrm{Ryd}}$')
    # plt.axvline(4*P.ryd_frq_TMDC, color='red', linewidth=1,
    #             label=r'$4\omega_{\mathrm{Ryd}}$')
    # plt.legend(loc="right") 
    plt.xlabel(r"$\omega$ in ps$^{-1}$")
    plt.ylabel(r"$Im(\chi)$")
    #plt.yscale('log')
    plt.grid()
    plt.show()
abs_of_w()     

# w = 529.6166367969113, E = 348.6


# 1.5,  10k: -26.30594880788261
# 3.5,  60k: -26.532283227888374
# 5.0, 120k: -26.56443677149896
# 6.0, 160k: -26.58049667978786    
# 8.0, 300k: -26.586654772040177
# Lit: -26.600097478572224
def convergence_2D():
    kmax = [1.5, 3.5, 5., 6., 8.]
    w = [-26.30594880788261, -26.532283227888374, -26.56443677149896,
         -26.58049667978786, -26.586654772040177]
    plt.figure(dpi=300)
    plt.plot(kmax, w, 'x', label='numerical calculation')
    plt.axhline(4*P.ryd_frq, color='red', linewidth=1,
                label=r'$4\omega_{\mathrm{Ryd}}$')
    plt.xlabel(r"$k$ in nm")
    plt.ylabel(r"$\omega$ in ps$^{-1}$")
    plt.legend()
    plt.grid()
    plt.show()
#convergence_2D()


