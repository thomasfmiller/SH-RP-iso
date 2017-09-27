      module model
      implicit none
      save

      integer, parameter :: nsurf = 2
      double precision, parameter :: A1 = -5.729577951D0 
      double precision, parameter :: B1 = 17.18873385D0
      double precision, parameter :: aw = 0.767495031D0
      double precision, parameter :: A2 = 7.D0
      double precision, parameter :: aw2 = 1.D0
      double precision, parameter :: x01 = -1.6D0
      double precision, parameter :: B2 = -0.75D0
      double precision, parameter :: D  = 0.25D0
      double precision, parameter :: aw3 = 0.25D0
      double precision, parameter :: x02 = -2.625D0
      double precision, parameter :: mass = 1.D0

      contains

      subroutine PotDeriv(x,V11,V22,V12,dV11,dV22,dV12)
      implicit none
      double precision, intent(in) :: x
      double precision, intent(out) :: V11,V22,V12,dV11,dV22,dV12
      double precision :: expa,cosha
       
       expa=dexp(-aw*x)
       cosha=dcosh(0.5d0*aw*x)
       V11 = A1/(1.d0+expa) + B1/(4.d0*cosha**2)
       dV11 = A1*expa*aw/(1.d0+expa)**2 - 
     .        0.25d0*B1*aw*sinh(0.5d0*aw*x)/cosha**3

       expa=dexp(-aw2*(x-x01))
       V22 = A2/(1.d0+expa) + B2
       dV22 = A2*expa*aw2/(1.d0+expa)**2

       V12 = D*dexp(-aw3*(x-x02)**2)
       dV12 = -aw3*2.d0*(x-x02)*V12

      return
      end subroutine potderiv

      end module model 
