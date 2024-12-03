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


#########################################################################################################
##    v1
#########################################################################################################
Note that DMFT_solved.zip includes a complete implementation of the IPT solver. Check for reference

A numerical integration yields different (and maybe better?) results when computing the local Green's function. We should take care here

The FFT approach seems the most reasonable one. In DMFT_solved, check ./Documentation/report_git2.pdf; they discuss the choice of Fourier transform formula.

The pade.py code is either not working or very bad. When using a Maximum Entropy approach, the results are metallic, but quite off I'd say.
However, the Green's function in Matsubara axis in the metallic phase seems ok.

Now the FFT in the first iteration goes smooth. The price to pay is a weird handling of the Matsubara frequencies, which I accomodated to the numpy sintax.
In the end, we can (and shall) only keeep positive freqs.

SOLVED ISSUES:
  - FFT seems to work
  - The DMFT loop is actually running
  - Im(GF) of the metallic phase in Matsubara axis goes to -2 at zero; good ! Sanity check

TASKS:
  - Either learn to use properly the ana_cont package, so that MaxEnt yields a better result, or ask Michele why pade.py is not working
  - Insulating phases are still a problem. We should address how to initialize insulating phases. My proposal is to change a bit the DMFT loop so that we start
from a GF, not a self-energy. That way, the typical insulating seed goes as 1/i*omega
  - We should carefully check if the implementation of Dyson_Green is ok. The direct numerical integration is surely fine, look at the GF.

