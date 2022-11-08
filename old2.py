# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 14:37:30 2022

@author: derto
"""

def int_integral_2D4(psik):
    nphi = 100
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    I = np.zeros(len(k))
    for ii in range(len(k)):
        k_integrand = np.zeros(len(k))
        for jj in range(len(k)):
            phi_integrand = k[jj]/np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi_grid[0])) 
            k_integrand[jj] = np.sum(phi_grid[1]*phi_integrand) * psik[jj]   
        I[ii] = np.sum(grid[1]*k_integrand) 
    return I



def int_integral_2D3(psik):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nphi = 1000
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    for ii in range(len(k)):
        k_integrand = np.zeros(len(k))
        for jj in range(len(k)):
            phi_integrand = k[jj]/np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi_grid[0])) 
            phi_integrand *= 4*k[ii]**4/(k[ii]**2 + k[jj]**2)**2
            k_integrand[jj] = np.sum(phi_grid[1]*phi_integrand)   
        double[ii] = np.sum(grid[1]*k_integrand)  
    #print(double)
    
    ### calculate the remidning terms
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                phi_integrand = k[jj]/(np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi_grid[0])))
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
            # phi_integrand = k[jj]/(np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi_grid[0])))
            # v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)
    I = np.sum(w_ij*psik, axis=1)    
    return I


def int_integral_2D5(psik):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nphi = 100
    phi_grid = Integration.grids.gts2(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        phi_integrand = np.zeros(len(phi))
        for ff in range(len(phi)):
            kjj_grid = k_int_grid(k[ii], 200)
            kjj = kjj_grid[0]
            k_integrand = kjj / np.sqrt(k[ii]**2 + kjj**2 - 2*k[ii]*kjj*np.cos(phi[ff])) *4*k[ii]**4/(k[ii]**2 + kjj**2)**2
            phi_integrand[ff] = Integration.integrater.int_disc(k_integrand, kjj_grid) 
            # k_integrand = lambda kjj: kjj / np.sqrt(k[ii]**2 + kjj**2 - 2*k[ii]*kjj*np.cos(phi[ff])) *4*k[ii]**4/(k[ii]**2 + kjj**2)**2
            # phi_integrand[ff] = integrate.quad(k_integrand, k[0], k[1])[0]
            #print(integrate.quad(k_integrand, k[0], k[1])[1])
        double[ii] = np.sum(phi_grid[1]*phi_integrand)
    print(double)    
    ### calculate the remidning terms
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                phi_integrand = k[jj]/(np.sqrt(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi)))
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    ### w
    w_ij = np.zeros((len(k), len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii == jj:
                w_ij[ii, jj] = - sum_i[ii] + double[ii]
            else:
                w_ij[ii, jj] = v_ij[ii, jj]  
    #w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)            
    I = np.sum(w_ij*psik, axis=1)   
    return I








def int_integral_3D3(psik):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nphi = 100
    phi_grid = Integration.grids.gts1(0, np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        k_integrand = np.zeros(len(k))
        for jj in range(len(k)):
            phi_integrand = k[jj]**2*np.sin(phi)/(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi)) 
            phi_integrand *= 4*k[ii]**4/(k[ii]**2 + k[jj]**2)**2
            k_integrand[jj] = np.sum(phi_grid[1]*phi_integrand)   
        plt.plot(k, k_integrand)    
        plt.axvline(k[ii])
        plt.show()
        plt.close()
        double[ii] = np.sum(grid[1]*k_integrand)  
        #print(double[ii]/k[ii])
    #print(double/k)
    
    ### calculate the remidning terms
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                phi_integrand = k[jj]**2*np.sin(phi)/(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi))
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)          
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+double)
    I = np.sum(w_ij*psik, axis=1)
    #print('M2', '\n', sum_i)   
    return I





def int_integral_3D4(psik):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nphi = 100
    phi_grid = Integration.grids.gts2(0, np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        phi_integrand = np.zeros(len(phi))
        for ff in range(len(phi)):
            kjj_grid = k_int_grid(k[ii], k[0], k[-1], 200)
            kjj = kjj_grid[0]
            k_integrand = kjj**2*np.sin(phi[ff]) / (k[ii]**2 + kjj**2 - 2*k[ii]*kjj*np.cos(phi[ff])) *4*k[ii]**4/(k[ii]**2 + kjj**2)**2
            phi_integrand[ff] = Integration.integrater.int_disc(k_integrand, kjj_grid) 
            # plt.plot(kjj, k_integrand)
            # plt.show()
            # plt.close()
            # k_integrand = lambda kjj: kjj**2*np.sin(phi[ff]) / (k[ii]**2 + kjj**2 - 2*k[ii]*kjj*np.cos(phi[ff])) *4*k[ii]**4/(k[ii]**2 + kjj**2)**2
            # phi_integrand[ff] = integrate.quad(k_integrand, k[0], k[1])[0]
        double[ii] = np.sum(phi_grid[1]*phi_integrand)

    ### calculate the remidning terms
    v_ij = np.zeros((len(k), len(k)))
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii ==  jj:
                v_ij[ii, jj] = 0
            else:
                phi_integrand = k[jj]**2*np.sin(phi)/(k[ii]**2 + k[jj]**2 - 2*k[ii]*k[jj]*np.cos(phi))
                v_ij[ii, jj] = grid[1][jj]*Integration.integrater.int_disc(phi_integrand, phi_grid)   
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    
    ### w
    w_ij = np.zeros((len(k), len(k)))
    for ii in range(len(k)):
        for jj in range(len(k)):
            if ii == jj:
                w_ij[ii, jj] = - sum_i[ii] + double[ii]
            else:
                w_ij[ii, jj] = v_ij[ii, jj]          
    I = np.sum(w_ij*psik, axis=1)    
    #print('M2', '\n', sum_i)   
    return I



##### gut

def int_integral_3D(psik):
    ki = k.reshape((len(k), 1))
    kj = k.reshape((1, len(k)))
    np.seterr(divide='ignore')
    v_ij = np.where(np.abs(ki-kj) > 0,
                    grid[1] * kj / ki * (np.log(np.abs(ki+kj))-np.log(np.abs(ki-kj))), 0) 
    np.seterr(divide='warn')
    v_im = v_ij * 4 * ki**4 / (ki**2+kj**2)**2
    sum_i = np.sum(v_im, axis=1)
    #print(-sum_i + np.pi)
    #print('M1', '\n', sum_i) 
    w_ij = np.where(np.abs(ki-kj) > 0, v_ij, -sum_i+np.pi*ki)
    I = np.sum(w_ij*psik, axis=1)
    #print(I.shape)                
    return I 





def int_integral_3D2(psik):
    psik_func = interpolate.interp1d(k, psik)
    I = np.zeros(len(k), dtype=np.complex_)
    nphi = 100
    nkjj = 500
    phi_grid = Integration.grids.gts2(0, np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        phi_integrand = np.zeros(len(phi))
        for ff in range(len(phi)):
            kjj_grid = misc.k_int_grid(k[ii], k[0], k[-1], nkjj)
            kjj = kjj_grid[0]
            denom = k[ii]**2 + kjj**2 - 2*k[ii]*kjj*np.cos(phi[ff])
            k_integrand = kjj**2*np.sin(phi[ff]) / denom  * psik_func(kjj)
            phi_integrand[ff] = Integration.integrater.int_disc(k_integrand, kjj_grid) 
        I[ii] = np.sum(phi_grid[1]*phi_integrand)
    return I


def int_integral_3D3(psik):
    psik_func = interpolate.interp1d(k, psik)
    I = np.zeros(len(k), dtype=np.complex_)
    nphi = 500
    nkjj = 1000
    phi_grid = Integration.grids.gts2(0, np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], k[0], k[-1], nkjj)
        kjj = kjj_grid[0]
        kjj_integrand = np.zeros(len(kjj))
        for jj in range(len(kjj)):
            denom = k[ii]**2 + kjj[jj]**2 - 2*k[ii]*kjj[jj]*np.cos(phi) 
            phi_integrand = kjj[jj]**2*np.sin(phi) / denom  * psik_func(kjj[jj])
            kjj_integrand[jj] = Integration.integrater.int_disc(phi_integrand, phi_grid)  
        I[ii] = np.sum(kjj_grid[1]*kjj_integrand)
    return I

# sol = np.loadtxt("sol_t/sol").view(complex)
# psi = sol[100, :]
# #psi = np.exp(-10*k)
# # psi = 1/(k)
# # plt.plot(k, psi)
# # plt.show()
# # plt.close()
# Ik = int_integral_3D(psi.real)
# Ik2 = int_integral_3D2(psi.real)
# plt.plot(k, Ik)
# plt.plot(k, Ik2,  '--')


def int_integral_2D(psik):
    psik_func = interpolate.interp1d(k, psik)
    I = np.zeros(len(k), dtype=np.complex_)
    nphi = 100
    nkjj = 500
    phi_grid = Integration.grids.gts2(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        phi_integrand = np.zeros(len(phi))
        for ff in range(len(phi)):
            kjj_grid = misc.k_int_grid(k[ii], k[0], k[-1], nkjj)
            kjj = kjj_grid[0]
            sqrt = np.sqrt(k[ii]**2 + kjj**2 - 2*k[ii]*kjj*np.cos(phi[ff]))
            k_integrand = kjj / sqrt * psik_func(kjj)
            phi_integrand[ff] = Integration.integrater.int_disc(k_integrand, kjj_grid) 
        I[ii] = np.sum(phi_grid[1]*phi_integrand)
    return I


def int_integral_2D2(psik):
    ### calucalte the double integral 
    double = np.zeros(len(k))
    nkjj = 100
    nphi = 100
    phi_grid = Integration.grids.gts1(0, 2*np.pi, nphi)
    phi = phi_grid[0]
    for ii in range(len(k)):
        kjj_grid = misc.k_int_grid(k[ii], k[0], k[-1], nkjj)
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
    I = np.sum(w_ij*psik, axis=1)    
    return I


# psi = np.exp(-10*k)
# Ik = int_integral_2D(psi)
# Ik2 = int_integral_2D2(psi)
# plt.plot(k, Ik)
# plt.plot(k, Ik2, '--')



#####################################################
    


def TMDC_pot(q):
    omega = 1
    gamma = P.damp
    
    # def eps_tilde
    eps_inf = 1
    chi = 1
    eps_tilde = (1+eps_inf)/2 + chi/(2*C.eps0)*q
    
    # def s²
    n = 1
    m = 1
    Omega_pl = np.sqrt(C.e**2 * n / (C.eps0 * eps_inf * m))
    s_sq = eps_inf * Omega_pl**2 / 2
    
    # def w_0,q
    v = C.hbar**2 * q**2 / (2 * m)
    kF = (3 * np.pi**2 * n)**1/3
    kappa = np.sqrt(C.e**2 / (np.pi**2 * C.eps0 * eps_inf * C.hbar**2) * kF * m)
    omega_0 = Omega_pl**2*q**2/kappa**2 + v**2
    
    # calculate W
    V = C.e**2 / (2 * C.eps0 * q) 
    bracket = omega_0**2 - omega**2 - 1j*gamma*omega
    eps_HS = (bracket)/(eps_tilde*bracket + s_sq)
    W = V * eps_HS
    return W



def TMDC_pot(q):
    # def eps_tilde
    eps_inf = 2
    chi = 0.66  #####
    #eps_tilde = (1+eps_inf)/2 + chi/(2*C.eps0)*q *C.eps0*4*np.pi
    eps_up = eps_inf
    eps_down = eps_inf
    eps_tilde = (eps_up+eps_down)/2 + chi/(2*C.eps0)*q *C.eps0*4*np.pi
    #print((1+eps_inf)/2, chi/(2*C.eps0)*q)
    
    # def s²
    n = 1e-3 # 10**18 cm**-3 
    m = 0.4 * C.me
    Omega_pl = np.sqrt(C.e**2 * n / (C.eps0 * eps_inf * m))
    s = np.sqrt(eps_inf * Omega_pl**2 / 2)
    s = 0
    # def w_0,q
    v = C.hbar**2 * q**2 / (2 * m)
    kF = (3 * np.pi**2 * n)**(1/3)
    kappa = np.sqrt(C.e**2 / (np.pi**2 * C.eps0 * eps_inf * C.hbar**2) * kF * m)
    omega_0 = np.sqrt(Omega_pl**2*q**2/kappa**2 + v**2)
    
    # calculate W
    V = C.e**2 / (2 * C.eps0 * q) 
    #V = 1 / q
    eps_HS_inv = omega_0**2/(eps_tilde*omega_0**2  + s**2)
    W = V * eps_HS_inv 
    
    #eps_HS = eps_tilde + s**2/omega_0**2
    #plt.figure(dpi=300, figsize=(3.5, 5))
    #plt.plot(q, 1/eps_HS_inv)
    #plt.plot(q, eps_HS, '-')
    #plt.ylim(0, 300)
    #plt.plot(q, W)
    #plt.plot(q, V, '-')
    #plt.ylim(0, 40000)
    return W
#TMDC_pot(k)
