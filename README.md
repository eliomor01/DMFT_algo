DMFT algorithm in python

Hyperparameters:
T = Temperature
N_tau = Number of discretization points
Energies = From the input

Variable names:

G_loc:dim 1 array $G_{loc}(i\omega_n)$
G_loc_tau: $G_{loc}(\tau)$

S_imp: $\Sigma (i\omega_n)$
Sp_imp: $\Sigma ' (i\omega_n)$

Momega: matsubara frequencies vector, already saved as imaginary, we save only the positive ones since the negatives and functions of them can be got by symmetry


##Guidelines:
Use a routine based approach

Make a routine for Fourier transform to go from imaginar time and viceversa

Make a Dyson solver routine

Make a routine to get self energy from the Green function

Check the dependance on chemical potential

For the moment do only the fermionic case

Check input structure on the virtualbox DMFT and retake it to do comparisons after



##v1:

Note that DMFT_solved.zip includes a complete implementation of the IPT solver. Check for reference

A numerical integration yields different (and maybe better?) results when computing the local Green's function. We should take care here

The FFT approach seems the most reasonable one. In DMFT_solved, check ./Documentation/report_git2.pdf; they discuss the choice of Fourier transform formula.

The pade.py code is either not working or very bad. When using a Maximum Entropy approach, the results are metallic, but quite off I'd say.
However, the Green's function in Matsubara axis in the metallic phase seems ok.

