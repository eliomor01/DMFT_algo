#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 09:04:19 2024

@author: cerradavelascoj
"""

import numpy as np


# Dyson equation

def solve_dyson(G_0, sigma):
    
    """
    Solves the Dyson equation given the impurity self-energy and a (1) non-interacting /
    (2) local Green's function. 
    
    In the first case, yields the impurity Green's function. In the second case, yields
    the non interacting Green's function.
    
    Input:
        G_0:   (1) Non-interacting Green's function in Matsubara axis
               (2) Local Green's function in Matsubara axis
             
        Sigma: The impurity self-energy
        
    Output:
        G_loc: (1) Impurity Green's function in Matsubara axis
               (2) Non-interacting Green's function in Matsubara axis
    """
    
    inverse_G0 = 1/G_0
    G_loc = 1/(inverse_G0 + sigma)
    
    return G_loc


# Self consistent test
    

def convergence_test(G_imp, G_loc, delta):
    
    """
    Given the newly computed impurity Green's function and the old local Green's
    function, performs the convergence test.
    
    Returns True if convergence has been reached with precision delta. Returns 
    False otherwise.
    
    Input:
        G_imp: New impurity Green's function in Matsubara axis
        G_loc: Old local Green's function in Matsubara axis
        delta: Precision level of the DMFT loop
        
    Output:
        bool: True if convergence has been achieved; false otherwise
    """

    difference = np.absolute(G_imp - G_loc)**2
    D = np.sqrt(np.sum(difference))
    
    return bool(np.where(D < delta, True, False))


# Pade approximation for analytical continuation

def real_axis_gf(G_loc, Momega, Romega, delta):
    
    '''
    Given the local Green's function in the Matsubara axis, computes the analytical
    continuation to the real axis with the Maximum Entropy method.
    
    Input:
        G_loc:  Local Green's function in Matsubara axis
        Momega: Matsubara axis grid
        Romega: Real frequency axis grid
        delta:  Precision level of the DMFT loop
        
    Output:
        G_loc_real: Local Green's function in real frequency axis
    '''
    
    import continuation as cont
    
    G_loc_real_double = np.zeros((Romega.shape[0],2))
    
    err = np.ones_like(Romega)*delta
    model = np.ones_like(Romega)
    model /= np.trapz(model, Romega)
    
    probl = cont.AnalyticContinuationProblem(im_axis=Momega, re_axis=Romega,
                                             im_data=G_loc, kernel_mode='freq_fermionic')
    
    sol, _ = probl.solve(method='maxent_svd', alpha_determination='chi2kink', optimizer='newton',
                         model=model, stdev=err, interactive=True, alpha_start=1e12, alpha_end=1e-2,
                         preblur=True, blur_width=0.5)
    
    G_loc_real_double[:,0] = sol.backtransform.real
    G_loc_real_double[:,1] = sol.backtransform.imag
    
    G_loc_real = G_loc_real_double[:,0] + 1j*G_loc_real_double[:,1]
    
    return G_loc_real
    
    

# Spectral function

def compute_spectral(G_loc_real):
    
    '''
    Takes the Green's function in real freauency axis and computes the spectral
    function.
    
    Input:
        G_loc_real: Local Green's function in real frequency axis
        
    Output:
        A_real: Spectral function in real frequency axis
        
    '''
    
    A_real = - G_loc_real.imag / np.pi
    
    return A_real
