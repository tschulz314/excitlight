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
sys.path.append(os.path.abspath(os.path.join(path, os.pardir)))
import Integration
import constants as C
import parameters as P
import misc as misc


    
def int_grid_3D(k, grid):
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    np.seterr(divide='ignore')
    v_ij = np.where(np.abs(ki-kj) > 0,
                    grid[1] * kj / ki * (np.log(np.abs(ki+kj))-np.log(np.abs(ki-kj))), 0) 
    np.seterr(divide='warn')
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+np.pi*ki)              
    return w_ij 


def int_grid_2D(k, grid):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nkjj = 1000
    nphi = 1000
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], k[0], 4*k[-1], nkjj)
        kjj = kjj_grid[0]
        kjj_integrand = np.zeros(len(kjj))
        for jj in range(len(kjj)):
            sqrt = np.sqrt(k[ii]**2 + kjj[jj]**2 - 2*k[ii]*kjj[jj]*np.cos(phi)) 
            phi_integrand = kjj[jj] / sqrt * 4 * k[ii]**4 / (k[ii]**2+kjj[jj]**2)**2
            kjj_integrand[jj] = Integration.integrater.int_disc(phi_integrand, phi_grid)  
        double[ii] = np.sum(kjj_grid[1]*kjj_integrand)
    
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                sqrt = np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi))
                phi_integrand = k[jj]/sqrt
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)  
    return w_ij



def TMDC_pot(q):
    chi =   7.112 * 0.1 # 7.112 * e-10
    kappa = 2 
    r0 = 2 * np.pi * chi / kappa
    V = C.e**2 / (2 * C.eps0 * q) 
    W = V * 1 / kappa * 1 / (1 + r0 * q)
    return  W
    #return  V / P.eps


def k_int_grid2(k, k0, k1, N):
    a = k0
    b = k.reshape((len(k), 1))
    c = k1
    ii = np.arange(0, N)
    ii = ii.reshape((1, N))
    # Häufungspunkt am Ende
    x1 = ( (ii+N+1)/(2*N+1) - 1 / ( 2 * np.pi ) * np.sin ( 2 * np.pi * (ii+N+1)/(2*N+1) ) ) * 2 * ( b-a ) + a - ( b-a )
    v1 = ( 1 - np.cos( 2 * np.pi * (ii+N+1)/(2*N+1) ) ) * ( 2 * ( b-a )/(2*N+1) )
    # Häufungspunkt am Anfang
    x2 = ( (ii+1)/(2*N+1) - 1 / ( 2 * np.pi ) * np.sin ( 2 * np.pi * (ii+1)/(2*N+1) ) ) * 2 * ( c-b ) + b
    v2 = ( 1 - np.cos( 2 * np.pi * (ii+1)/(2*N+1) ) ) * ( 2 * ( c-b )/(2*N+1) )
    # Insgesamt 
    x = np.concatenate((x1, x2), axis=1)
    v = np.concatenate((v1, v2), axis=1)
    #w = 
    return x, v



def int_grid_eff(k, kgrid):
    # phi grid
    fpot = TMDC_pot
    nphi = 200
    phi, phigrid = Integration.grids.gts1(0, 2*np.pi, nphi)
    ### calculate the double integral over kj and phi 
    # calculate integrand toint(ki, kj, phi) for doubnle(ki)
    nkjj = 200
    kj, kjgrid = k_int_grid2(k, 0, 10, nkjj)
    k_3d = k.reshape((len(k), 1, 1))
    kj_3d = kj.reshape((len(k), 2*nkjj, 1))
    phi_3d = phi.reshape((1, 1, len(phi)))
    sqrt = np.sqrt( k_3d**2 + kj_3d**2 - 2*k_3d*kj_3d*np.cos(phi_3d) )
    toint = kj_3d * fpot(sqrt) * 4 * k_3d**4 / (k_3d**2+kj_3d**2)**2 
    toint = np.sum(toint*phigrid, axis=2) # phi integral
    double = np.sum(toint*kjgrid, axis=1) # kj integral
    
    ### calculate the weights
    ki_3d = k.reshape((len(k), 1, 1))
    kj_3d = k.reshape((1, len(k), 1))
    phi_3d = phi.reshape((1, 1, len(phi)))
    sqrt = np.sqrt( ki_3d**2 + kj_3d**2 - 2*ki_3d*kj_3d*np.cos(phi_3d) )
    toint = kj_3d * fpot(sqrt)
    phiint = np.sum(toint * phigrid, axis=2) 
    ki_2d = k.reshape((len(k), 1))
    kj_2d = k.reshape((1, len(k)))
    v_ij = np.where(np.abs(ki_2d - kj_2d)==0, 0, kgrid*phiint)
    v_im = v_ij * 4 * ki_2d**4 / (ki_2d**2+kj_2d**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki_2d-kj_2d) > 0, v_ij, -sum_i+double)  
    #print(w_ij)
    return w_ij




def int_grid_eff2(k, kgrid):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nkjj = 200
    nphi = 200
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], 0, 10, nkjj)
        kjj = kjj_grid[0]
        kjj_integrand = np.zeros(len(kjj))
        for jj in range(len(kjj)):
            sqrt = np.sqrt(k[ii]**2 + kjj[jj]**2 - 2*k[ii]*kjj[jj]*np.cos(phi)) 
            phi_integrand = kjj[jj] * TMDC_pot(sqrt) * 4 * k[ii]**4 / (k[ii]**2+kjj[jj]**2)**2
            kjj_integrand[jj] = Integration.integrater.int_disc(phi_integrand, phi_grid)  
        double[ii] = np.sum(kjj_grid[1]*kjj_integrand)
    
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                sqrt = np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi))
                phi_integrand = k[jj] * TMDC_pot(sqrt)
                v_ij[ii, jj] = kgrid[jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)  
    #print(w_ij)
    return w_ij


def int_grid_eff3(k, kgrid):
    # calucalte the double integral 
    double = np.zeros(len(k))
    nkjj = 200
    nphi = 200
    phi, phigrid = Integration.grids.gts1(0, 2*np.pi, nphi)
    # double(ki) # toint(ki, kj, phi)
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], 0, 10, nkjj)
        kjj = kjj_grid[0]
        kjj_integrand = np.zeros(len(kjj))
        for jj in range(len(kjj)):
            sqrt = np.sqrt(k[ii]**2 + kjj[jj]**2 - 2*k[ii]*kjj[jj]*np.cos(phi)) 
            phi_integrand = kjj[jj] * TMDC_pot(sqrt) * 4 * k[ii]**4 / (k[ii]**2+kjj[jj]**2)**2
            kjj_integrand[jj] = np.sum(phi_integrand * phigrid)  
        double[ii] = np.sum(kjj_grid[1]*kjj_integrand)
    
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                sqrt = np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi))
                phi_integrand = k[jj] * TMDC_pot(sqrt)
                v_ij[ii, jj] = kgrid[jj]*np.sum(phi_integrand * phigrid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)  
    #print(w_ij)
    return w_ij


def int_grid_eff4(k, kgrid):
    # calucalte the double integral 
    nkjj = 200
    nphi = 200
    phi, phigrid = Integration.grids.gts1(0, 2*np.pi, nphi)

    double = np.zeros(len(k)) # double(ki)
    toint = np.zeros((len(k), 2*nkjj, len(phi)))  # toint(ki, kj, phi) 
    kj = np.zeros((len(k), 2*nkjj)) # kj(ki, kj)    
    kjgrid = np.zeros((len(k), 2*nkjj)) # kjgrid(ki, kj)             
    for ii in range(len(k)):
        kj[ii, :], kjgrid[ii, :] = misc.k_int_grid(k[ii], 0, 10, nkjj)
    k_3d = k.reshape((len(k), 1, 1))
    kj_3d = kj.reshape((len(k), 2*nkjj, 1))
    phi_3d = phi.reshape((1, 1, len(phi)))
    sqrt = np.sqrt( k_3d**2 + kj_3d**2 - 2*k_3d*kj_3d*np.cos(phi_3d) )
    toint = kj_3d * TMDC_pot(sqrt) * 4 * k_3d**4 / (k_3d**2+kj_3d**2)**2
    toint = np.sum(toint*phigrid, axis=2)
    double = np.sum(toint*kjgrid, axis=1)
    
    
    
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                sqrt = np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi))
                phi_integrand = k[jj] * TMDC_pot(sqrt)
                v_ij[ii, jj] = kgrid[jj]*np.sum(phi_integrand * phigrid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)  
    #print(w_ij)
    return w_ij


def int_grid_eff5(k, kgrid):
    # phi grid
    nphi = 200
    phi, phigrid = Integration.grids.gts1(0, 2*np.pi, nphi)
    ### calculate the double integral over kj and phi 
    nkjj = 200
    double = np.zeros(len(k)) # double(ki)
    toint = np.zeros((len(k), 2*nkjj, len(phi)))  # toint(ki, kj, phi) 
    kj = np.zeros((len(k), 2*nkjj)) # kj(ki, kj)    
    kjgrid = np.zeros((len(k), 2*nkjj)) # kjgrid(ki, kj)             
    for ii in range(len(k)):
        kj[ii, :], kjgrid[ii, :] = misc.k_int_grid(k[ii], 0, 10, nkjj)
    k_3d = k.reshape((len(k), 1, 1))
    kj_3d = kj.reshape((len(k), 2*nkjj, 1))
    phi_3d = phi.reshape((1, 1, len(phi)))
    sqrt = np.sqrt( k_3d**2 + kj_3d**2 - 2*k_3d*kj_3d*np.cos(phi_3d) )
    toint = kj_3d * TMDC_pot(sqrt) * 4 * k_3d**4 / (k_3d**2+kj_3d**2)**2 
    toint = np.sum(toint*phigrid, axis=2) # phi integral
    double = np.sum(toint*kjgrid, axis=1) # kj integral
    
    ### calculate the weights
    ki_3d = k.reshape((len(k), 1, 1))
    kj_3d = k.reshape((1, len(k), 1))
    phi_3d = phi.reshape((1, 1, len(phi)))
    sqrt = np.sqrt( ki_3d**2 + kj_3d**2 - 2*ki_3d*kj_3d*np.cos(phi_3d) )
    toint = kj_3d * TMDC_pot(sqrt)
    phiint = np.sum(toint * phigrid, axis=2) 
    ki_2d = k.reshape((len(k), 1))
    kj_2d = k.reshape((1, len(k)))
    v_ij = np.where(np.abs(ki_2d - kj_2d)==0, 0, kgrid*phiint)
    v_im = v_ij * 4 * ki_2d**4 / (ki_2d**2+kj_2d**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki_2d-kj_2d) > 0, v_ij, -sum_i+double)  
    #print(w_ij)
    return w_ij




def int_grid_eff6(k, kgrid):
    # phi grid
    fpot = TMDC_pot
    nphi = 200
    phi, phigrid = Integration.grids.gts1(0, 2*np.pi, nphi)
    ### calculate the double integral over kj and phi 
    # calculate integrand toint(ki, kj, phi) for doubnle(ki)
    nkjj = 200
    kj, kjgrid = k_int_grid2(k, 0, 10, nkjj)
    k_3d = k.reshape((len(k), 1, 1))
    kj_3d = kj.reshape((len(k), 2*nkjj, 1))
    phi_3d = phi.reshape((1, 1, len(phi)))
    sqrt = np.sqrt( k_3d**2 + kj_3d**2 - 2*k_3d*kj_3d*np.cos(phi_3d) )
    toint = kj_3d * fpot(sqrt) * 4 * k_3d**4 / (k_3d**2+kj_3d**2)**2 
    toint = np.sum(toint*phigrid, axis=2) # phi integral
    double = np.sum(toint*kjgrid, axis=1) # kj integral
    
    ### calculate the weights
    ki_3d = k.reshape((len(k), 1, 1))
    kj_3d = k.reshape((1, len(k), 1))
    phi_3d = phi.reshape((1, 1, len(phi)))
    sqrt = np.sqrt( ki_3d**2 + kj_3d**2 - 2*ki_3d*kj_3d*np.cos(phi_3d) )
    toint = kj_3d * fpot(sqrt)
    phiint = np.sum(toint * phigrid, axis=2) 
    ki_2d = k.reshape((len(k), 1))
    kj_2d = k.reshape((1, len(k)))
    v_ij = np.where(np.abs(ki_2d - kj_2d)==0, 0, kgrid*phiint)
    v_im = v_ij * 4 * ki_2d**4 / (ki_2d**2+kj_2d**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki_2d-kj_2d) > 0, v_ij, -sum_i+double)  
    #print(w_ij)
    return w_ij