# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:40:27 2021

@author: derto
"""



import numpy as np
import matplotlib.pyplot as plt


def read_data():

    t = np.loadtxt("solutions/time")
    k = np.loadtxt("solutions/momentum")
    sol = np.loadtxt("solutions/sol").view(complex) #, dtype=np.complex_
    return t, k, sol


t, k, sol = read_data()


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
    
    #plt.close()
    #plt.plot(t, sol[:, 0])
    #plt.plot(t, np.sin(t))
    
#time_evo() 

def gauss(t, texp, sigma):
    return np.exp(-1/2*((t-texp) / sigma)**2) 

def psit():
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True) 
    ktoplot = [20, 40]
    for index in ktoplot:
        ax1.plot(t, sol[:, index].real, '-')
        ax2.plot(t, sol[:, index].imag, '-', label=r"$k={}\,$nm$^{{-1}}$".format(round(k[index])))
     
    ax3.plot(t, gauss(t, 0.05, 0.015))    
    #ax1.set_xlabel(r"$t$ in ps")
    ax3.set_xlabel(r"$t$ in ps")
    ax1.set_ylabel(r"$Re(\psi)$")
    ax2.set_ylabel(r"$Im(\psi)$")
    #fig.title('Anfangsverteilung')
    ax2.legend(bbox_to_anchor=(1.1, 1.05))
    ax1.grid()
    ax2.grid()
    ax3.grid()
    #plt.savefig("fig/fk_" + branch + '.png', dpi=300)
    #plt.savefig("fig/fk_0.png", dpi=300)
    plt.show()

    #plt.close()
    #plt.plot(t, sol[:, 0])
    #plt.plot(t, np.sin(t))
    
psit() 


