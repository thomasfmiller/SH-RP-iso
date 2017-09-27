c     subroutine, which computes the adiabatic ring-polymer isomorphic potential energy surfaces
c     their derivatives and the first-derivative non-adiabatic coupling
      subroutine PotDeriv_iso(x,nbead,betan,
     .           V11ia,V22ia,dV11ia,dV22ia,D12i)
      implicit none
      integer nbead
      double precision x(nbead),betan
      double precision V11ia,V22ia,D12i(nbead)
      double precision dV11ia(nbead),dV22ia(nbead)
      double precision V11i,V22i,V12i
      double precision dV11i(nbead),dV22i(nbead),dV12i(nbead)
      double precision mu,dmu(nbead)
      double precision root,droot(nbead)
      double precision holel(2,2,nbead)
      double precision holer(2,2,nbead)
      double precision m(2,2,nbead)
      double precision dm(2,2,nbead)
      double precision V11,dV11,V22,dV22,V12,dV12
      double precision t(2,2),t2(2,2)
      double precision exps,beta,snh,dvs,ddvs,dv
      double precision expsmu,expsn
      double precision rt2,rt2tol,pf
      parameter(rt2tol=1.d-20)
      integer i

        V11i=0.d0
        V22i=0.d0
        do i = 1, nbead
c        computes M=[exp(-\betan*hat{V(x_i)})]_{\alpha\beta} and dM(x_i)_{\alpha\beta}/dx_i
         call m_dm(x(i),betan,m(:,:,i),dm(:,:,i),V11,dV11,V22,dV22)
         V11i=V11i+V11
         V22i=V22i+V22
         dV11i(i)=dV11
         dV22i(i)=dV22 
        enddo

c       recursive computation of the electronic ring-polymer partition function
c       left and right hole matrices are computed in linear Nbead effort
c       Follows Bell's algorithm, see the paper for a reference
        holel(:,:,1)=m(:,:,1)
        do i = 2, nbead-1
          holel(:,:,i)=matmul(holel(:,:,i-1),m(:,:,i))
        enddo

        holer(:,:,nbead)=m(:,:,nbead)
        do i = nbead-1, 2, -1
          holer(:,:,i)=matmul(m(:,:,i),holer(:,:,i+1))
        enddo

c       electronic ring-polymer partition function (\Mu from the paper)
        t=matmul(holel(:,:,nbead-1),m(:,:,nbead))
        mu=t(1,1)+t(2,2) 

c       dMu/dXi - gradients
        t=matmul(dm(:,:,1),holer(:,:,2))
        dmu(1)=t(1,1)+t(2,2)
        do i = 2, nbead-1
          t=matmul(dm(:,:,i),holer(:,:,i+1))
          t2=matmul(holel(:,:,i-1),t)
          dmu(i)=t2(1,1)+t2(2,2)
        enddo
        t=matmul(holel(:,:,nbead-1),dm(:,:,nbead))
        dmu(nbead)=t(1,1)+t(2,2)

c       implementation of the formula for the off-diagonal potential energy matrix element
        exps=exp(0.5d0*betan*(V11i+V22i))
        expsmu=exps*mu
        root=2.d0*acosh(0.5d0*expsmu)/betan
        
c       derivative prefactor
        snh=exps/(betan*sinh(0.5d0*betan*root))

c       isomorphic adiabatic potential energy surfaces
        dvs=V11i+V22i
        V11ia=0.5d0*(dvs-root)
        V22ia=0.5d0*(dvs+root)

c       gradient of the isomorphic adiabatic potential energy surfaces
        do i = 1, nbead
         ddvs=dV22i(i)+dV11i(i)
         droot(i)=(0.5d0*betan*ddvs*mu+dmu(i))*snh
         dV11ia(i)=0.5d0*(ddvs-droot(i))
         dV22ia(i)=0.5d0*(ddvs+droot(i))
        enddo

c       isomorphic off-diagonal potential energy matrix element (K_12^{iso} from the paper)
        dv=V22i-V11i
        rt2=root**2-dv**2

c       cut-off for very small or negative coupling
        if(rt2.lt.rt2tol)then
         V12i=0.d0
         D12i=0.d0
         dV12i=0.d0
         write(*,'(''WARNING: The off-diagonal coupling is'')')
         write(*,'(''         imaginary or very close to zero!'')')
         write(*,'(''Setting to zero and continuing'')')
        else
         V12i=0.5d0*sqrt(root**2-dv**2)
c        gradient and first derivative non-adiabatic coupling
         do i = 1, nbead
          dV12i(i)=0.25d0*(root*droot(i)-dv*(dV22i(i)-dV11i(i)))/V12i
          D12i(i)=(dV12i(i)*dv-V12i*(dV22i(i)-dV11i(i)))/root**2
         enddo
c        scale by the number of beads
         D12i=dble(nbead)*D12i
        endif

      return
      end 

c     subroutine, which computes the matrix M=exp(-betan*H) and its derivative
c     this subroutine requires the function potderiv in a module model, 
c     which computes the physical diabatic potentials and their derivatives for a given position
      subroutine m_dm(x,betan,m,dm,V11,dV11,V22,dV22)
      use model, only: potderiv
      implicit none
      double precision x, betan, m(2,2), dm(2,2)
      double precision V11,V22,V12,dV11,dV22,dV12 
      double precision dv,vs,ddv,dvs,dis,cs,sn,dcs,dsn
      double precision th,dth,root,droot,l(2),dl(2)
      double precision u(2,2),du(2,2),ut(2,2),dut(2,2)
      double precision expl(2,2),dexpl(2,2)
      double precision t1(2,2)

c      PotDeriv computes the diabatic potential energy matrix and its derivatives
c      It is located in module model.
       call PotDeriv(x,V11,V22,V12,dV11,dV22,dV12)

c      adiabatic potential energy
       vs=V22+V11
       dv=V22-V11
       dvs=dV22+dV11
       ddv=dV22-dV11

       dis = dv*dv+4.d0*V12*V12
       root = dsqrt(dis)
       l(1) = 0.5d0*(vs-root)
       l(2) = 0.5d0*(vs+root)
       droot=(dv*ddv+4d0*V12*dV12)/root
       dl(1)= 0.5d0*(dvs-droot)
       dl(2)= 0.5d0*(dvs+droot)

       th =0.5d0*DATAN2(2.d0*V12,dv)
       dth=(dV12*dv-V12*ddv)/dis

c      unitary transformation matrix for adiabatic to diabatic representation
       cs=dcos(th)
       sn=dsin(th)
       dcs=-sn*dth
       dsn=cs*dth

       u(1,1)=cs
       u(2,1)=-sn
       u(1,2)=sn
       u(2,2)=cs

       ut(1,1)=cs
       ut(2,1)=sn
       ut(1,2)=-sn
       ut(2,2)=cs

c      derivative of the unitary matrix
       du(1,1)=dcs
       du(2,1)=-dsn
       du(1,2)=dsn
       du(2,2)=dcs

       dut(1,1)=dcs
       dut(2,1)=dsn
       dut(1,2)=-dsn
       dut(2,2)=dcs

c      Boltzmann operator in adiabatic representation
       expl(1,1)=dexp(-betan*l(1))
       expl(2,1)=0.d0
       expl(1,2)=0.d0
       expl(2,2)=dexp(-betan*l(2))

c      derivative of the Boltzmann operator in adiabatic representation
       dexpl(1,1)=(-betan)*expl(1,1)*dl(1)
       dexpl(2,1)=0.d0
       dexpl(1,2)=0.d0
       dexpl(2,2)=(-betan)*expl(2,2)*dl(2)

c      transformation to diabatic representation
       t1=matmul(expl,ut)
       m =matmul(u,t1)

c      transformation to diabatic representation of the derivative
       dm=matmul(du,t1)
       
       t1=matmul(dexpl,ut)
       dm=dm+matmul(u,t1)

       t1=matmul(expl,dut)
       dm=dm+matmul(u,t1) 

      return
      end
