# Gadget-2-openmp
This is a modified version of Gadget-2 with a hybrid (openmp+MPI) paralellization scheme for tidal disruption simulations.  Gadget-2 is originally available from http://wwwmpa.mpa-garching.mpg.de/gadget/right.html .

It has been modified to include an openmp parallelization routine in the Gravity and Hydro solvers, in parts following MPJ Express Meets Gadget: Towards a Java Code for Cosmological Simulations (https://doi.org/10.1007/11846802_50 or http://mpj-express.org/docs/papers/mpj-gadget-parsim06.pdf) 

This code was used in some of the simulations of 
Turbovelocity Stars: Kicks Resulting from the Tidal Disruption of Solitary Stars
Manukian, H., Guillochon, J., Ramirez-Ruiz, E., & O'Leary, R.~M. 2013, ApJL, 771, L28 




Original License:

GADGET is free software, distributed under the GNU General Public License. This implies that you may freely distribute and copy the software. You may also modify it as you wish, and distribute these modified versions as long as you indicate prominently any changes you made in the original code, and as long as you leave the copyright notices, and the no-warranty notice intact. Please read the General Public License for more details. Note that the authors retain their copyright on the code.

If you use GADGET-2 for scientific work, we kindly ask you to reference the code paper(s) on GADGET, i.e.

Springel V., 2005, MNRAS, 364, 1105
Springel V., Yoshida N., White S. D. M., 2001, New Astronomy, 6, 51
