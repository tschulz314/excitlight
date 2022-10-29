# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:40:27 2021

@author: derto
"""



import numpy as np
import matplotlib.pyplot as plt

#import parameters as P
#import misc as misc






def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def comp_psi_over_k(tval):
    # sol from t space
    t1 = np.loadtxt("sol_t/time")
    k1 = np.loadtxt("sol_t/momentum")
    psit1 = np.loadtxt("sol_t/sol").view(complex) 
    index1 = find_nearest(t1, tval)
    # sol from w space
    t2 = np.loadtxt("sol_w/time")
    k2 = np.loadtxt("sol_w/momentum")
    psit2 = np.loadtxt("sol_w/psit").view(complex) 
    index2 = find_nearest(t2, tval)
    # plot
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, dpi=300) 
    #ax1.plot(k1, psit1[index1, :].real, '-', label=r"$t={}\,$ps".format(t1[index1]))
    #ax1.plot(k1, psit1[index1, :].real, 'x', label=r"$t={}\,$ps".format(t1[index1]))
    ax1.plot(k2, psit2[index2, :].real, '--', label=r"$t={}\,$ps".format(t2[index2]))
    #ax2.plot(k1, psit1[index1, :].imag, '-', label=r"$t={}\,$ps".format(t1[index1]))
    #ax2.plot(k1, psit1[index1, :].imag, 'x', label=r"$t={}\,$ps".format(t1[index1]))
    ax2.plot(k2, psit2[index2, :].imag, '--', label=r"$t={}\,$ps".format(t2[index2]))
    ax1.set_xlabel(r"$k$ in nm$^{-1}$")
    ax2.set_xlabel(r"$k$ in nm$^{-1}$")
    ax1.set_ylabel(r"$Re(\psi_k)$")
    ax2.set_ylabel(r"$Im(\psi_k)$")
    ax2.set_xlim(0.1, 0.2)
    #fig.title('Anfangsverteilung')
    plt.legend()
    ax1.grid()
    ax2.grid()
    #plt.savefig("fig/fk_" + branch + '.png', dpi=300)
    #plt.savefig("fig/fk_0.png", dpi=300)
    plt.show()
comp_psi_over_k(5) 


  