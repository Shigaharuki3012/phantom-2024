!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module moddump
use prompting, only: prompt
!
! set solid body rotation for all gas particles about the z-axis given angular frequency
!
! :References: None
!
! :Owner: Mike Lau
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 implicit none

contains

subroutine modify_dump(npart,npartoftype,massoftype,xyzh,vxyzu)
 integer, intent(inout) :: npart
 integer, intent(inout) :: npartoftype(:)
 real,    intent(inout) :: massoftype(:)
 real,    intent(inout) :: xyzh(:,:),vxyzu(:,:)
 real                   :: omega,vphi,R
 integer                :: i
 real                   :: d0,k,omega_s,omega_s_max,d,R_star,M_star,alpha

 ! Assume rotation axis is the z-axis
!  omega = 7.92e-3  ! Omega in code units
 d0 = 0.02
 R_star = 1 
 M_star = 1
 alpha = 3.0
 call prompt('Core Relative Size:', d0, 0.)
 call prompt('Star Radius:', R_star, 0.)
 call prompt('Star Mass:', M_star, 0.)
 call prompt('Power law Index:', alpha, 0.)
 k = d0**(alpha-1.5)
 omega_s_max=sqrt(M_star/R_star**3.0)
 omega_s = k/sqrt(2.0)*omega_s_max
 print*,"surface rotation frequency:", omega_s

 do i = 1,npart
   R = sqrt(xyzh(1, i)**2 + xyzh(2, i)**2)
   d = R/R_star
   if (d > d0) then
      omega = omega_s*d**(-alpha)
   elseif (d <= d0) then
      omega = omega_s*d0**(-alpha)
   end if
   vphi = omega*R
   vxyzu(1, i) = -omega*xyzh(2, i)
   vxyzu(2, i) = omega*xyzh(1, i)
   vxyzu(3, i) = 0.
 enddo

 return
end subroutine modify_dump

end module moddump

