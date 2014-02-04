!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE OPTDEPTH_V2(CHI,NZ,IP,DO_P_RAY)
      USE MOD_SPACE_GRID_V2
      IMPLICIT NONE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Calculates dtau using a formula based on the Euler-Maclaurin
! summation formula (Auer 84, -> Knuth 68)
! Determine differences: dtau == dtau(d) == tau(d+1)-tau(d)
!
! written   11-94 DLM
! altered   12-95 DLM Added abs(z) so a positive tau and dtau would
!                     be found for negative z (ie s_m).
! Altered 2/21/96 DLM Updated to F90 standard
! altered 5/27/97 DLM Removed include "nd_parameters.f".  Changed to
!                     allocate statements.
! altered 6/20/97 DLM Corrected logic for inward and outward direction.
!                     Previously used absolute values, now check direction
!                     and use separated loops.
!
!--------------------------------------------------------------------
!
      integer nz
      integer ip
      real*8 chi(nz)
      logical do_p_ray
!
! Local variables
!
      integer iz
      real*8 dchidz(nz)
!
!--------------------------------------------------------------------
!
! Determine dtau for either inward or outward direction
!
      IF(DO_P_RAY)THEN
!
! For outward, determine the derivative of chi
!
        dchidz(1)=(chi(1)-chi(2))/(ray(ip)%s_p(1)-ray(ip)%s_p(2))
        do iz=2,nz-1
          dchidz(iz)=(chi(iz-1)-chi(iz+1))/(ray(ip)%s_p(iz-1)-ray(ip)%s_p(iz+1))
        enddo
        dchidz(nz)=(chi(nz-1)-chi(nz))/(ray(ip)%s_p(nz-1)-ray(ip)%s_p(nz))
	dchidz(1:nz)=0.0d0
!
! Determine dtau
!
        tau(1)=0.0D0
        do iz=1,nz-1
          dtau(iz)=0.5d0*(ray(ip)%s_p(iz)-ray(ip)%s_p(iz+1))
     *         *( chi(iz)+chi(iz+1)
     *         +(ray(ip)%s_p(iz)-ray(ip)%s_p(iz+1))
     *         *(dchidz(iz+1)-dchidz(iz))/6.0d0 )
          tau(iz+1)=tau(iz)+dtau(iz)
        enddo
!
      else
!
! For inward, determine the derivative of chi
!
        dchidz(1)=(chi(2)-chi(1))/(ray(ip)%s_m(2)-ray(ip)%s_m(1))
        do iz=2,nz-1
          dchidz(iz)=(chi(iz+1)-chi(iz-1))/(ray(ip)%s_m(iz+1)-ray(ip)%s_m(iz-1))
        enddo
        dchidz(nz)=(chi(nz)-chi(nz-1))/(ray(ip)%s_m(nz)-ray(ip)%s_m(nz-1))
	dchidz(1:nz)=0.0d0
!
! Determine dtau
!
        tau(1)=0.0D0
        do iz=1,nz-1
          dtau(iz)=0.5d0*(ray(ip)%s_m(iz+1)-ray(ip)%s_m(iz))
     *         *(chi(iz+1)+chi(iz)
     *         +(ray(ip)%s_m(iz+1)-ray(ip)%s_m(iz))
     *         *(dchidz(iz)-dchidz(iz+1))/6.0d0)
          tau(iz+1)=tau(iz)+dtau(iz)
        enddo
!
      endif
!
      return
      end
