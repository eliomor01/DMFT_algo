import numpy as np
import matplotlib.pyplot as plt

########################################################################
# INPUTS
########################################################################

#D0, u0, nloop, xmix, max_omega, T, tol  = np.loadtxt('input.bin')
D0, u0, nloop, xmix, max_omega, T, tol  = 1.0, 2.0, 40, 0.7, 4.0, 0.1, 1e-6

try:

    DOS = np.loadtxt('DOS.txt')

except:

    energies = np.linspace(-D0, D0, 1000)
    DOS = 2*np.sqrt(1-(energies/D0)**2)/(np.pi * D0)



########################################################################
# DMFT ROUTINES
########################################################################

def DMFT_iteration_first_half(u0,S_imp,Energies,Momega,tau,T):
    N_tau = len(tau)
    G_loc = Dyson_Green(S_imp, energies, DOS, Momega)
    G_0 = G_loc/(1+G_loc*S_imp)
    G_0_tau = Inv_Matsubara_Fourier_div (G_0,T,tau,Momega)
    new_S_imp = Self_Energy_fromG(G_0_tau,u0,tau,Momega)
    return G_loc, G_0, G_0_tau, new_S_imp

def Self_Energy_fromG (G,U,tau,omega):
    Ntau = len(tau)
    dtau = tau[1] - tau[0]
    esp = np.exp(1.0j*np.tensordot(tau,omega,0))
    integrand = G*G*G
    return U*U*dtau*np.matmul(integrand,esp)

def Matsubara_Fourier (v):
    return np.fft.fft(v)

def Inv_Matsubara_Fourier_div (G,T,tau,omega):
    N_tau = len(tau)
    Gp = G - 1/(omega)
    return 2*T*np.sum(np.real(Gp.reshape((-1,1)).repeat(N_tau,1)*np.exp(-np.tensordot(omega,tau,axes=0))),axis=0) - np.sign(tau)/2
    
    
def Dyson_Green (Self_Energy, Energy, Density_Energy, omega):
    N_freqs = len(omega)
    epsrep = Energy.reshape((1,-1)).repeat(N_freqs,0)
    desrep = Density_Energy.reshape((1,-1)).repeat(N_freqs,0)
    if len(epsrep) == len(desrep):
        N_e = epsrep.shape[1]
        summand =omega.reshape((-1,1)).repeat(N_e,1) - epsrep - Self_Energy.reshape((-1,1)).repeat(N_e,1)
        return np.sum(desrep/summand, axis=1)/(N_e)


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


########################################################################
# INITIALIZE FUNCTIONS
########################################################################

Momega = (2*np.arange(1000) + 1)*np.pi*T
Romega = Momega
tau = np.linspace(0, 1/T, 1000)
S_imp = np.zeros((len(Momega)))


########################################################################
# DMFT LOOP
########################################################################


for k in range(nloop):

    G_loc, G_0, G_0_tau, new_S_imp = DMFT_iteration_first_half(u0,S_imp, DOS, Momega, tau, T)
    G_imp = solve_dyson(G_0, new_S_imp)
    
    if convergence_test(G_imp, G_loc,tol):
        print('Convergence has been reached, iteration '+str(k))
        break
    else:
        S_imp = xmix * S_imp + (1-xmix) * new_S_imp

G_loc_real = real_axis_gf(G_loc, Momega, Romega, tol)
spectral   = compute_spectral(G_loc_real)

plt.plot(Romega, spectral)
plt.show()