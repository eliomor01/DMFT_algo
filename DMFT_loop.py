import numpy as np


def DMFT_iteration_first_half(S_imp,Energies,Momega,tau,T):
    N_freqs = len(Momega)
    N_tau = len(tau)
    epsrep = Energies.reshape((1,-1)).repeat(N_freqs,0)
    N_e = epsrep.shape[1]
    summand =Momega.reshape((-1,1)).repeat(N_e,1) - epsrep - S_imp.reshape((-1,1)).repeat(N_e,1)
    G_loc = np.sum(1/summand, axis=1)/(N_e)
    G_0 = G_loc/(1+G_loc*S_imp)
    Gp_0 = G_0 - 1/(Momega)
    G_0_tau = 2*T*np.sum(np.real(Gp_0.reshape((-1,1)).repeat(N_tau,1)*np.exp(-np.tensordot(Momega,tau,axes=0))),axis=0) - np.sign(tau)/2
    return G_0_tau
