# SH-RP-iso
Example of ring polymer surface hopping (SH-RP-iso) for a two-level scattering problem

#### Any use of this code should cite the following paper:

#### "Path-integral isomorphic Hamiltonian for including nuclear quantum effects in non-adiabatic dynamics."  X. Tao, P. Shushkov, and T. F. Miller III, submitted, available at: arXiv:1709.06722.  (Submitted to J. Chem. Phys.  Citation to be updated upon journal acceptance.)

The enclosed fortran code illustrates the use of the SH-RP-iso method applied to a two-level 
scattering problem with a single nuclear degree of freedom, as employed in Fig. 4 of the paper above.

Description of the enclosed files:

## [1] rp_iso.f
 
This subroutine computes the rpmd version of the isomorphic potential for a general 
two-level system with one nuclear coordinate (Eqs. 16-19 in the paper), 
as well as its derivative with respect to ring polymer bead positions. 
Following evaluation of the isomorphic potential in the diabatic representation, 
the isomorphic potential is diagonalized to obtain the adiabatic representation and 
the first-derivative non-adiabatic coupling, as described in the paper.

###### INPUT:
- non-negative integer number log_{2}{nbead}, where n is the number of beads for the ring-polymer
- the imaginary time step, \beta_n=\beta/nbead=\frac{1}/{nbead k_B T}
- the position array of the ring-polymer \vec{x}, x[nbead]

###### OUTPUT:
- adiabatic ring-polymer surfaces, V11ia and V22ia
- their derivative with respect to the ring-polymer bead position, dV11ia and dV22ia
- first-derivative non-adiabatic couplings, D12i

## [2] model.f

A module that provides the diabatic physical potential matrix for a two-level scattering problem in one dimension.

###### INPUT: 
- the position variable, x

###### OUTPUT:
- the diabatic potential energy matrix elements, Vij
- their derivatives with respect to position, dVij

## [3] sh_rp_iso_integrator.f

A module that uses  subroutines in [1] and [2] to evolve a SH-RP-iso trajectory. 
Initially, thermalized ring-polymer internal modes are sampled with a given 
centroid position and velocity. Uses the velocity-Verlet algorithm with a fourth-order 
Runge-Kutta integrator, as described in the paper.

###### INPUT:
- temperature, temp
- timestep, dt
- initial ring-polymer centroid position, xcin
- initial ring-polymer centroid velocity, vcin

###### OUTPUT:
- The SH-RP-iso trajectory.

## [4] main.f

A simple driver main routine.

Given a set of input parameters, such as temperature, timestep, initial adiabatic surface,
initial centroid position, initial centroid momentum, number of ring polymer beads, and random number generator seed,
this routine integrates a surface-hopping ring-polymer trajectory on the isomorphic Hamiltonian.
The driver routine requires a stopping criterion for trajectory termination and allows for 
the sampling of various properties to be included by the user.
