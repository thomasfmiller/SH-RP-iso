Example code for the paper:

"Path-integral isomorphic Hamiltonian for including nuclear quantum effects in non-adiabatic dynamics"
X. Tao, P. Shushkov, and T. F. Miller III, submitted, available at: arXiv:1709.06722.
(Citation to be updated upon journal acceptance.)

Any use of this code should cite the paper indicated above.

The enclosed fortran code illustrates the use of the SH-RP-iso method applied to a two-level 
scattering problem with a single nuclear degree of freedom, as employed in Fig. 4 of the paper above.

Description of included files:

[1] rpiso.f
This subroutine computes the rpmd version of the isomorphic potential for a general 
two-level system with one nuclear coordinate (Eqs. 16-19 in the paper), 
as well as its derivative with respect to ring polymer bead positions. 
Following evaluation of the isomorphic potential in the diabatic representation, 
the isomorphic potential is diagonalized to obtain the adiabatic representation and 
the first-derivative non-adiabatic coupling, as described in the paper.
INPUT:
the number of beads for the ring-polymer, nbead,
the imaginary time step, \beta_n=\beta/nbead=\frac{1}/{nbead k_B T},
the position array of the ring-polymer \vec{x}, x[nbead].
OUTPUT:
adiabatic ring-polymer surfaces, V_{ii}ia; their derivative with respect to the ring-polymer bead
position, dV_{ii}ia; first-derivative non-adiabatic coupling, D12i.


[2] model.f
A module that provides the diabatic physical potential matrix for the subsequent calculations. 
INPUT: 
the position variable x.
OUTPUT:  
the value of the diabatic potential energy surfaces, couplings, and their derivatives with respect
to position, V_ij and dV_ij.

[3] rpsh_iso_integrator.f
The mododule uses the subroutines in [1] and [2] to evolve a SH-RP-iso trajectory. 
Initially, thermalized ring-polymer internal modes are sampled with a given 
centroid position and velocity. The velocity-Verlet algorithm with a fourth-order 
Runge-Kutta integrator is applied, as described in the paper.
INPUT:
temp: temperature, dt: timestep, 
xconstr & vcin: initial centroid position and velocity for the ring-polymer 
OUTPUT:
A time-dependent trajectory \vec{x}(t), saved in the array x, which will be used for post-evolution calculations.

[4] main.f
A simple driver main routine.

Given a set of input parameters, such as temperature, timestep, initial adiabatic surface,
initial centroid position, initial centroid momentum, number of ring polymer beads, and random number generator seed,
this routine integrates a surface-hopping ring-polymer trajectory on the isomorphic Hamiltonian.
The driver routine requires a stopping criterion for trajectory termination and allows for 
the sampling of various properties to be included by the user.

