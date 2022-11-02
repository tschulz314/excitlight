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



t = np.loadtxt("sol_t/time")
k = np.loadtxt("sol_t/momentum")
sol = np.loadtxt("sol_t/sol").view(complex) 
#sol = sol*np.exp(- 1j * P.w0*t.reshape((len(t), 1)) )
polt = np.loadtxt("sol_t/polt").view(complex)
w = np.loadtxt("sol_t/frequency")
polw = np.loadtxt("sol_t/polw").view(complex)
#Ew = np.loadtxt("sol_t/Ew").view(complex)
Ew = misc.E0w(-w)

def psi_over_k():
    T = len(t)-1
    #t_ind = [0, int(T/100), int(T/10), int(T/5), int(T/2), int(T)]
    t_ind = [int(T/100)] #int(T/10), int(T/5), 
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=300) 
    #plt.figure() #figsize=(7, 7) 
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
#psi_over_k() 


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def spicific_psik(ktoplot):
    plt.plot(t, sol[:, 0].real)
    plt.show()
    plt.close()
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=300) 
    plt.tight_layout()
    kind = find_nearest(k, ktoplot)
    ax1.plot(t, sol[:, kind].real, '-', label=r"$k={}\,$nm$^{{-1}}$".format(round(k[kind], 2)))
    #ax1.plot(t, sol[:, kind].imag, '-') 
    ax2.plot(t, misc.E0(t)) 
    #ax2.plot(t, misc.E(t).real) 
    ax2.set_xlabel(r"$t$ in ps")
    ax1.set_ylabel(r"$\tilde{\psi}_k(t)$")
    #ax1.set_xlim(6.6, 6.7)
    ax2.set_ylabel(r"$E(t)$")
    ax1.legend()
    ax1.grid()
    ax2.grid()
    plt.show()
#spicific_psik(0) 
#spicific_psik(0.4) 
#spicific_psik(1) 


def psit():
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, dpi=300) 
    K = len(k) - 1
    #print(K)
    #for ii in range(len(k)):
    #    if k[ii] <
    ktoplot = [round(K/2), round(K*7/10), round(K*3/10)]
    for index in ktoplot:
        ax1.plot(t, sol[:, index].real, '-')
        ax2.plot(t, sol[:, index].imag, '-', label=r"$k={}\,$nm$^{{-1}}$".format(round(k[index], 1)))
    ax3.plot(t, misc.E0(t))    
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
#psit() 


def psit2():
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, dpi=300) 
    plt.tight_layout()
    K = len(k) - 1
    #print(round(k[round(K*7/10)], 1))
    ax1.plot(t, sol[:, round(K*3/10)].real, '-')#, label=r"$k={}\,$nm$^{{-1}}$".format(round(k[round(K*3/10)])))
    #ax1.plot(t, sol[:, round(K*3/10)].imag, '-') #, label=r"$k={}\,$nm$^{{-1}}$".format(round(k[index], 1)))
    ax2.plot(t, sol[:, round(K/2)].real, '-')
    #ax2.plot(t, sol[:, round(K/2)].imag, '-')
    ax3.plot(t, sol[:, round(K*7/10)].real, '-')
    #ax3.plot(t, sol[:, round(K*7/10)].imag, '-')
    #ax1.set_xlim(0, 0.1)
    #ax4.plot(t, E0(t)) 
    ax4.plot(t, misc.E(t)) 
    ax3.set_xlabel(r"$t$ in ps")
    ax1.set_ylabel(r"$\psi_k, k={}\,$nm$^{{-1}}$".format(round(k[round(K*3/10)], 1)))
    ax2.set_ylabel(r"$\psi_k, k={}\,$nm$^{{-1}}$".format(round(k[round(K/2)], 1)))
    ax3.set_ylabel(r"$\psi_k, k={}\,$nm$^{{-1}}$".format(round(k[round(K*7/10)], 1)))
    ax4.set_ylabel(r"$E(t)$")
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()
#psit2() 
    
    
    



def polt_evo():
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
    ax3.set_ylabel(r"$E$")
    #ax2.legend(loc="right") 
    ax1.grid()
    ax2.grid()
    ax3.grid()
    plt.show()
#polt_evo() 


def polw_evo():
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=300) 
    ax1.plot(w, polw.real, '-')
    ax2.plot(w, polw.imag, '-')
    ax2.set_xlabel(r"$\omega$ in ps$^{-1}$")
    #ax3.set_xlim(0, 0.1)
    ax1.set_ylabel(r"$Re(P)$")
    ax2.set_ylabel(r"$Im(P)$")
    #ax2.legend(loc="right") 
    ax1.grid()
    ax2.grid()
    plt.show()
#polw_evo() 



def chi_over_w():
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=300) 
    chi = polw / (Ew*C.eps0) 
    ax1.plot(w, chi.real, '-')
    ax2.plot(w, chi.imag, '-')
    ax2.set_xlabel(r"$\omega$ in ps$^{-1}$")
    #ax2.set_xlim(-50, 50)
    #ax2.set_ylim(0, 100)
    ax1.set_ylabel(r"$Re(\chi)$")
    ax2.set_ylabel(r"$Im(\chi)$")
    #ax2.legend(loc="right") 
    ax1.grid()
    ax2.grid()
    plt.show()
#chi_over_w()     


def allw_evo():
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, dpi=300) 
    #plt.tight_layout()
    chi = polw / Ew 
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
#allw_evo()    
    
    
def abs_over_w():
    plt.figure(dpi=300) #
    chi = polw / (Ew*C.eps0) 
    #plt.plot(w, chi.real, '-')
    plt.plot(w, chi.imag, '-')
    plt.xlabel(r"$\omega$ in ps$^{-1}$")
    #ax2.set_xlim(-50, 50)
    #plt.ylim(0, 0.025)
    plt.axvline(P.ryd_frq, color='grey')
    plt.ylabel(r"$Im(\chi)$")
    #ax2.legend(loc="right") 
    plt.grid()
    plt.show()
abs_over_w()     
