import numpy as np

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
    
def Self_Energy_fromG (G,U,tau,omega):
    Ntau = len(tau)
    dtau = tau[1] - tau[0]
    esp = np.exp(1.0j*np.tensordot(tau,omega))
    integrand = G*G*G
    return U*U*dtau*np.matmul(integrand,esp)

def DMFT_iteration_first_half(S_imp,Energies,U,Momega,tau,T):
    N_tau = len(tau)
    G_loc = Dyson_Green(S_imp, Energies, Momega)
    G_0 = G_loc/(1+G_loc*S_imp)
    G_0_tau = Inv_Matsubara_Fourier_div (G_0,T,tau,Momega)
    S_new = Self_Energy_fromG(G_0_tau,U,tau,Momega)
    return G_0_tau

