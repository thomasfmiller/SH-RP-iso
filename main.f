c     simple surface-hopping ring-polymer isomorphic Hamiltonian driver routine
      program SHRPiso
c     this module is located in rpsh_iso_integrator.f
      use dynamics
      implicit none
      integer ISEED
      double precision xcin,vcin
      real rrr

c      read input parameters from the command line
c      initial temperature (K)
       call GETARG(1,buffer)
       read(buffer,*) temp
c      timestep (a.u.)
       call GETARG(2,buffer)
       read(buffer,*) dt
c      initial surface
       call GETARG(3,buffer)
       read(buffer,*) nsurf
c      initial centroid position (a.u.)
       call GETARG(4,buffer)
       read(buffer,*) xcin
c      initial centroid velocity (a.u.)
       call GETARG(5,buffer)
       read(buffer,*) vcin
c      number of beads
       call GETARG(6,buffer)
       read(buffer,*) nbead
c      initial SEED
       call GETARG(7,buffer)
       read(buffer,*) ISEED

c      initialize free ring-polymer quantities 
       call InitParams()

c      initialize random number generator
       rrr = ran2(ISEED)

c      initialize the free ring-polymer normal modes
       call Samplefree(xcin,vcin,ISEED)

c      initialize the electronic coefficients
       call InitEl(nsurf)

c      initialize the forces for the nuclear degrees of freedom propagation
       call Force_Long()

c      all other initialization over sampled quantities go here ...

c      loop over integration steps
       do

c        integrate for a single timestep
c        using fourth-order Runge-Kutta method for the electronic coefficients
c        and alternation between velocity Verlet and free ring-polymer propagation for the nuclear degrees of freedom
         call Integrate(ISEED) 

c        problem specific criterion to terminate integration goes here ...

       enddo


      stop
      end
