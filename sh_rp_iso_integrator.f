      module dynamics
c     requires module model with the physical mass of the particle
      use model, only:mass
      implicit none
      save

c     Dynamics parameters      
c     number of beads
      integer nbead
c     Boltzmann constant $ pi
      double precision kb,pi
      parameter(kb=3.16682E-06)
      parameter(pi=3.14159265359d0)
c     Dynamics arrays and parameters
c     Bead positions
      double precision, allocatable :: x(:)
c     Bead velocities
      double precision, allocatable :: v(:)
c     Forces from the potential energy surface (divided by mass)
      double precision, allocatable :: fl(:)
c     Bead normal mode positions
      double precision, allocatable :: xnm(:)
c     Bead normal mode velocities
      double precision, allocatable :: vnm(:)
c     Gradient of the potential energy surfaces
      double precision, allocatable :: grd(:,:)
c     Bead first-derivative non-adiabatic coupling 
      double precision, allocatable :: d12(:)
c     Normal mode-to-Cartesian position transformation matrix 
      double precision, allocatable :: nm2p(:,:)
c     Free ring-polymer normal mode propagator for dt/2
      double precision, allocatable :: propfr(:,:,:)
c     Free ring-polymer normal mode frequencies
      double precision, allocatable :: omfr(:)
c     Free ring-polymer normal mode thermal widths
      double precision, allocatable :: sigfr(:)
c     beta/Nbead
      double precision betan
c     (Nbead*temp)^2
      double precision omega2
c     half-timestep, timestep, 1/mass, fact=1/3
      double precision hdt,dt,m_1,fact
c     temperature (a.u.), thermal width, thermal centroid width
      double precision temp,sigma,sigmac
c     active surface
      integer nsurf
c     electronic coefficients
      double precision C(4)
c     off-diagonal Hamiltonian matrix element in adiabatic representation - H12 = -D12V 
      double precision H12
c     trajectory phase = \int_{0}^{t} (E_1-E_2) dt
      double precision Phase
c     D12V = \sum_i^{Nbead} d_{12,i} * v_{i} / Nbead
      double precision dc(2)
c     isomorphic adiabatic potential energy surfaces
      double precision epot(2)

      contains

      subroutine InitParams()
      implicit none
      integer i,j,k
      double precision ph,norm,omn,cs,sn

c      the program uses only even number of beads
       if(mod(nbead,2).ne.0)stop'Odd number of beads!'

c      set-up common parameters
c      temperature in kbT
       temp = kb*temp
       betan= 1.d0/(temp*dble(nbead))
       sigma = DSQRT(temp*dble(nbead)/mass)
       sigmac = DSQRT(temp/mass)
       omn = dble(nbead)*temp
       omega2 = omn*omn
       m_1 = 1.d0/mass
       fact = 1.d0/3.d0
       hdt = 0.5d0*dt
       
c      allocate main arrays
       allocate(x(nbead),v(nbead),fl(nbead),grd(nbead,2),
     .          d12(nbead),nm2p(nbead,nbead),propfr(2,2,nbead),
     .          omfr(nbead),sigfr(nbead),xnm(nbead),vnm(nbead))

       ph=2.d0*pi/dble(nbead)
       norm=sqrt(1.d0/dble(nbead)) 

c      normal mode-to-Cartesian position transformation matrix
       do j = 1, nbead
        nm2p(j,1)=norm
       enddo

       do j = 1, nbead
        nm2p(j,nbead/2+1)=norm*(-1.d0)**j
       enddo

       norm=norm*sqrt(2.d0) 

       do k = 2, nbead/2
        do j = 1, nbead
         nm2p(j,k)=norm*cos(ph*dble(j*(k-1)))
        enddo
       enddo

       do k = nbead/2+2,nbead
        do j = 1, nbead
         nm2p(j,k)=norm*sin(ph*dble(j*(k-1)))
        enddo
       enddo

c      normal mode frequencies
       ph=pi/dble(nbead)
       do i = 1, nbead
        omfr(i)=2.d0*omn*sin(dble(i-1)*ph)
       enddo

c      normal mode thermal widths
       sigfr(1)=0.d0
       do i = 2, nbead
        sigfr(i)=sqrt(dble(nbead)*temp/(mass*omfr(i)**2))
       enddo

c      free ring-polymer propagator for 1/2 dt
       propfr(1,1,1)=1.d0
       propfr(2,1,1)=hdt
       propfr(1,2,1)=0.d0
       propfr(2,2,1)=1.d0
       do i = 2, nbead
        ph = omfr(i)*hdt
        cs = cos(ph)
        sn = sin(ph)
        propfr(1,1,i)=cs
        propfr(2,1,i)=sn/omfr(i)
        propfr(1,2,i)=-omfr(i)*sn
        propfr(2,2,i)=cs
       enddo

      return
      end subroutine initparams

c     initialize the electronic states to nsurfin - initial surface
      subroutine InitEl(nsurfin)
      implicit none
      integer nsurfin

       if(nsurfin.eq.1)then
        C(1)=1.0d0
        C(2)=0.0d0
        C(3)=0.0d0
        C(4)=0.0d0
       elseif(nsurfin.eq.2)then
        C(1)=0.0d0
        C(2)=0.0d0
        C(3)=1.0d0
        C(4)=0.0d0
       endif

      return
      end subroutine initel

c     sample free ring-polymer modes from a thermal distribution
c     given centroid initial velocity and position returns a ring-polymer with the rest of the modes thermalized 
      subroutine Samplefree(xcin,vcin,ISEED)
      implicit none
      double precision xcin,vcin
      integer ISEED
      integer i,j,k
      double precision xc,vc

c      sample normal modes from thermal distribution
       xnm(1)=0.d0
       do i = 2, nbead
        xnm(i) = GASDEV(sigfr(i),ISEED)
       enddo

       x=0.d0
       do k = 1, nbead
        do j = 1, nbead
         x(j)=x(j)+nm2p(j,k)*xnm(k)
        enddo
       enddo

c      sample velocities from thermal distribution
       do i=1,nbead
        v(i) = GASDEV(sigma,ISEED)
       enddo

c      adjust centroid position to xcin
       xc = Centroid(x,nbead)
       do i=1,nbead
        x(i) = x(i)-xc+xcin
       enddo

c      adjust centroid velocity to vcin
       vc = Centroid(v,nbead)
       do i=1,nbead
        v(i) = v(i)-vc+vcin
       enddo

      return
      end subroutine samplefree

c     computes the centroid of a(n)
      double precision function Centroid(a,n)
      implicit none
      integer n,i
      double precision a(n),ac

      ac = 0.0d0
      do i=1,n
       ac = ac + a(i)
      enddo
      Centroid = ac/dble(n)

      return
      end function centroid

C     Note that the isomorphic adiabatic surfaces and the non-adiabatic coupling matrix element are divided by Nbead, 
c     according to the definition of the isomorphic Hamiltonian in the paper 
c     see also the appendix A in the paper for the equivalent forms of the ring-polymer equations of motion
c     we use eqn A9 for propagation of the ring polymer positions and velocities

c     integrator for a ring-polymer surface-hopping trajectory on the isomorphic potential energy surfaces
c     the electronic time-propagation is carried out in interaction representation 
c     uses fourth-order Runge-Kutta (RK4) integrator for the electronic coefficients
c     alternates between free ring-polymer and velocity Verlet steps with the isomorphic adiabatic surface forces
      subroutine Integrate(ISEED)
      implicit none
      integer ISEED 
      double precision k1(4),k2(4),k3(4),k4(4),Cc(4),prob,dV1,dV2
      integer I,J

c     first step of RK4
      call Hamiltonian(dV1)
      call IntStep(k1,C,hdt)

c     propagate velocities with the isomorphic adiabatic surface forces for a half timestep
      do I=1,nbead
       v(I) = v(I) + fl(I)*hdt
      enddo

c     propagate the free ring-polymer normal modes for a half timestep
      call IntegrateFree()

c     compute adiabatic potential energy surfaces, 
c     forces of the adiabatic surfaces, 
c     and first-derivative non-adiabatic coupling 
      call Force_Long()

c     second and third step of RK4 at t+1/2dt
      call Hamiltonian(dV2)
      Phase = 0.5d0*(dV2 + dV1)*hdt + Phase

      Cc = C + k1
      call IntStep(k2,Cc,hdt)

      Cc = C + k2
      call IntStep(k3,Cc,dt)

c     propagate the free ring-polymer normal modes for a half timestep
      call IntegrateFree()

c     compute adiabatic potential energy surfaces, 
c     forces of the adiabatic surfaces, 
c     and first-derivative non-adiabatic coupling
c     at t+dt
      call Force_Long()

c     propagate velocities with the isomorphic adiabatic surface forces for a half timestep
      do I=1,nbead
       v(I) = v(I) + fl(I)*hdt
      enddo

c     fourth step of RK4 at t+dt
      call Hamiltonian(dV1)
      Phase = 0.5d0*(dV2 + dV1)*hdt + Phase

      Cc = C + k3
      call IntStep(k4,Cc,hdt)

c     RK4 propagation of the electronic coefficients to t+dt
      do I=1,4
       C(I) = C(I) + (k1(I)+k4(I)+2.d0*k2(I)+k3(I))/3.d0               
      enddo

c     compute probability for surface hop
c     prob(k\rightarrow l)=\Delta{t}\frac{1}{a_{kk}}\left( \frac{1}{N}\sum_i^N D_{l'l,i}\text{v}_i \right)2\text{Re}(a_{l'l})
      prob = (C(1)*C(3)+C(2)*C(4))*DCOS(Phase)
      prob =  prob + (C(1)*C(4)-C(2)*C(3))*DSIN(Phase)
      prob = 2.d0*dc(nsurf)*dt*prob 
      prob = prob / (C(2*nsurf-1)*C(2*nsurf-1)+C(2*nsurf)*C(2*nsurf))
c     if the hopping probability is negative, set it to zero
      if(prob.lt.0.0d0) prob = 0.0d0

c     draw a random number and deside whether to hop
      if(ran2(ISEED).lt.prob) call SurfaceSwitch()

      return
      end subroutine integrate 

C     IntStep computes \dot{\widetilde{c}}(t)=-i\widetilde{H}(t)\widetilde{c}(t) and k=incr*\dot{\widetilde{c}}, where incr=dt or dt/2
      subroutine IntStep(k,Cc,incr)
      implicit none
      integer i
      double precision k(4),Cc(4),incr,Hi(2)

c      interaction representation Hamiltonian elements - \tilde{H_{ij}}=H_{ij}e^{i\int_0^{t} (H_{ii} - H_{jj}) dt'}
c      time-dependent Schrodinger equation in adiabatic interaction representation for a two-level system is:
c      \dot{\tilde{c_k}}(t)=-\left(\frac{1}{N}\sum_i^N D_{kl,i}\text{v}_i\right)e^{i\int_0^{t} (H_{kk} - H_{ll}) dt'}\tilde{c_l}
       Hi(1) = H12*DCOS(Phase)
       Hi(2) = -H12*DSIN(Phase)
 
       k(1) = Hi(1)*Cc(3) - Hi(2)*Cc(4)
       k(2) = Hi(1)*Cc(4) + Hi(2)*Cc(3)
       k(3) = - Hi(1)*Cc(1) - Hi(2)*Cc(2)
       k(4) = - Hi(1)*Cc(2) + Hi(2)*Cc(1)

       k = k * incr

      return
      end subroutine intstep

C     Hamiltonian computes the electronic Hamiltonian coupling element
      subroutine Hamiltonian(dE)
      implicit none
      double precision dE
      integer i

       dc(1) = 0.d0
       do i = 1, nbead
        dc(1) = dc(1) + d12(i)*v(i)
       enddo
       dc(1) = dc(1)/dble(nbead)
       dc(2) = -dc(1)
       H12 = -dc(1)
       dE = (epot(2)-epot(1))/dble(nbead)

      return
      end subroutine hamiltonian

c     Perform energy conserving surface hop
      subroutine SurfaceSwitch()
      implicit none
      integer surf_new,I
      double precision Pot(2),dE
      double precision dnorm,proj,proj_new

       if(nsurf.eq.1) surf_new = 2
       if(nsurf.eq.2) surf_new = 1

       pot=epot

c      unit vector along the first-derivative non-adiabatic coupling
       dnorm = 0.0d0
       do I=1,nbead
        dnorm = dnorm + d12(I)*d12(I)
       enddo
       dnorm = DSQRT(dnorm)

c      projection of the velocity along the first-derivative non-adiabatic coupling direction
       proj = 0.0d0
       do I=1,nbead
        proj = proj + v(I)*d12(I)
       enddo
       proj = proj / dnorm

c      adjust this projection of velocity to conserve total energy in the surface switch 
       proj_new = proj*proj - 2.d0*m_1*(Pot(surf_new)-Pot(nsurf))

c      do not hop if there is not enough projection of velocity to conserve total energy,
c      i.e. reject frustrated hops
c      does not reverse velocity (if velocity reversal is implemented, then it should be included here)
       if(proj_new.lt.0.0d0) return

c      take the square root and preserve the direction of projection of velocity 
       proj_new = DSQRT(proj_new)
       if(proj.lt.0.0d0) proj_new = -proj_new

c      adjust the bead velocities according to the new projection 
       proj = proj - proj_new
       do I=1,nbead
        v(I) = v(I) - proj*d12(I)/dnorm
       enddo

c      switch the active surface
       nsurf = surf_new

c      compute the forces to restart the nuclear degrees of freedom integration
       call Force_Long()

      return
      end subroutine surfaceswitch

c     integrator of the free ring-polymer modes
      subroutine IntegrateFree()
      implicit none
      integer j,k
      double precision t,t2

c      transform to free ring-polymer normal modes
       xnm=0.d0
       do k = 1, nbead
        do j = 1, nbead
         xnm(k)=xnm(k)+nm2p(j,k)*x(j)
        enddo
       enddo

c      transform the velocities, as well
       vnm=0.d0
       do k = 1, nbead
        do j = 1, nbead
         vnm(k)=vnm(k)+nm2p(j,k)*v(j)
        enddo
       enddo

c      propagate the normal modes for a half timestep
       do k = 1, nbead
        t=propfr(1,1,k)*vnm(k)+propfr(1,2,k)*xnm(k)
        t2=propfr(2,1,k)*vnm(k)+propfr(2,2,k)*xnm(k)
        vnm(k)=t
        xnm(k)=t2
       enddo

c      back-transform to Cartesian positions and velocities
       x=0.d0
       do k = 1, nbead
        do j = 1, nbead
         x(j)=x(j)+nm2p(j,k)*xnm(k)
        enddo
       enddo

       v=0.d0
       do k = 1, nbead
        do j = 1, nbead
         v(j)=v(j)+nm2p(j,k)*vnm(k)
        enddo
       enddo

      return
      end subroutine integratefree

c     this subroutine (from rp_iso.f) computes
c     forces from the isomorphic adiabatic potential energy surfaces
c     isomorphic adiabatic potential energy surfaces
c     first-derivative non-adiabatic coupling from the isomorphic Hamiltonian
      subroutine Force_Long()
      implicit none

        call PotDeriv_iso(x,nbead,betan,
     .       epot(1),epot(2),grd(:,1),grd(:,2),d12)

        fl(:)=-grd(:,nsurf)*m_1

      return
      end subroutine force_long

c     Numerical recipes random number generators
      DOUBLE PRECISION FUNCTION gasdev(sigma,idum)
      INTEGER idum
      DOUBLE PRECISION sigma
      INTEGER iset
      REAL fac,gset,rsq,v1,v2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1      v1=2.*ran2(idum)-1.
       v2=2.*ran2(idum)-1.
       rsq=v1**2+v2**2
       if(rsq.ge.1..or.rsq.eq.0.) goto 1
       fac=sqrt(-2.*log(rsq)/rsq)
       gset=v1*fac
       gasdev=v2*fac
       iset=1
      else
       gasdev=gset
       iset=0
      endif
       GASDEV = DBLE(GASDEV)*SIGMA
      return
      END FUNCTION gasdev

      REAL FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2
      INTEGER IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1)
      PARAMETER (IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211)
      PARAMETER (IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB)
      PARAMETER (EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do j=NTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if (idum.lt.0) idum=idum+IM1
           if (j.le.NTAB) iv(j)=idum
         enddo
         iy=iv(1)
       endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if(idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      RETURN
      END FUNCTION ran2

      end module dynamics


