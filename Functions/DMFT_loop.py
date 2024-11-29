import numpy as np
import matplotlib.pyplot as plt

########################################################################
# INPUTS
########################################################################

#D0, u0, nloop, xmix, max_omega, T, tol  = np.loadtxt('input.bin')
D0, u0, nloop, xmix, n_omega, T, tol  = 1.0, 1.5, 40, 0.7, 500, 0.1, 1e-6

try:

    data1 = np.loadtxt('.\\DOS.txt')
    energies = data1[:,0]
    DOS   = data1[:,1]

except:

    energies = np.linspace(-D0, D0, 1000)
    DOS = 2*np.sqrt(1-(energies/D0)**2)/(np.pi * D0)



########################################################################
# DMFT ROUTINES
########################################################################

def DMFT_first_iteration(u0,S_imp,Energies,dos,Momega,tau,T):
    N_tau = len(tau)
    G_loc = Dyson_Green(S_imp, Energies, dos, Momega)
    g_0 = G_loc/(1+G_loc*S_imp)
    G_0_tau = Inv_Matsubara_Fourier_div (g_0,T,tau,Momega)
    plt.plot(tau, G_0_tau.real)
    plt.show()
    new_S_imp = Self_Energy_fromG(G_0_tau,u0,tau,Momega)
    return G_loc, g_0, G_0_tau, new_S_imp

def DMFT_iteration_first_half(u0,S_imp,Energies,dos,Momega,tau,T):
    N_tau = len(tau)
    G_loc = Dyson_Green(S_imp, Energies, dos, Momega)
    g_0 = G_loc/(1+G_loc*S_imp)
    G_0_tau = Inv_Matsubara_Fourier_div (g_0,T,tau,Momega)
    new_S_imp = Self_Energy_fromG(G_0_tau,u0,tau,Momega)
    return G_loc, g_0, G_0_tau, new_S_imp

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
    Gp = G - 1/(1.0j*omega)
    return 2*T*np.sum(np.real(Gp.reshape((-1,1)).repeat(N_tau,1)*np.exp(-1.0j*np.tensordot(omega,tau,axes=0))),axis=0) - np.sign(tau)/2
    
    
def Dyson_Green (Self_Energy, Energy, Density_Energy, omega):
    N_freqs = len(omega)
    epsrep = Energy.reshape((1,-1)).repeat(N_freqs,0)
    desrep = Density_Energy.reshape((1,-1)).repeat(N_freqs,0)
    if len(epsrep) == len(desrep):
        N_e = epsrep.shape[1]
        N_states = desrep.sum()
        summand =omega.reshape((-1,1)).repeat(N_e,1) - epsrep - Self_Energy.reshape((-1,1)).repeat(N_e,1)
        return np.sum(desrep/summand, axis=1)/(N_states)


# Dyson equation

def solve_dyson(G_0, sigma):
    
    """
    Solves the Dyson equation given the impurity self-energy and a non-interacting Green's function. 
    
    Yields the impurity Green's function

    Input:
        G_0:   Non-interacting Green's function in Matsubara axis
             
        Sigma: The impurity self-energy
        
    Output:
        G_loc: Impurity Green's function in Matsubara axis 
    """
    
    inverse_G0 = 1/G_0
    G_loc = 1/(inverse_G0 - sigma)
    
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

Momega = (2*np.arange(n_omega) + 1)*np.pi*T
tau = np.linspace(0, 1/T, 10000)
S_imp = np.zeros((len(Momega)))


########################################################################
# DMFT LOOP
########################################################################

G_loc, G_0, G_0_tau, new_S_imp = DMFT_first_iteration(u0,S_imp, energies, DOS, Momega, tau, T)
G_imp = solve_dyson(G_0, new_S_imp)
S_imp = new_S_imp

for k_loop in range(nloop):

    G_loc, G_0, G_0_tau, new_S_imp = DMFT_iteration_first_half(u0,S_imp, energies, DOS, Momega, tau, T)
    G_imp = solve_dyson(G_0, new_S_imp)
    
    if convergence_test(G_imp, G_loc, tol):
        print('Convergence has been reached, iteration '+str(k_loop))
        break
    else:
        S_imp = xmix * new_S_imp + (1-xmix) * S_imp

path1 = '.\\g_loc.dat'
with open(path1, 'w') as writer:

    for i in range(len(Momega)):
        writer.write(str(Momega[i]) + ' ' + str(np.real(G_loc[i])) + ' ' + str(np.imag(G_loc[i])) +'\n')

path2 = '.\\sigma.dat'
with open(path2, 'w') as writer:

    for i in range(len(Momega)):
        writer.write(str(Momega[i]) + ' ' + str(np.real(S_imp[i])) + ' ' + str(np.imag(S_imp[i])) +'\n')


import pade as pade
pade

G_loc_real_table = np.loadtxt('.\\pade.dat')
G_loc_real = G_loc_real_table[:,1] + 1.0j*G_loc_real_table[:,2]

spectral = compute_spectral(G_loc_real)
Romega = G_loc_real_table[:,0]

plt.plot(Romega, spectral)
plt.show()
