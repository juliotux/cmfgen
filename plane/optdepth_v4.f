!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
      SUBROUTINE OPTDEPTH_V4(DTAU_LOC,CHI,NZ,IP,DO_P_RAY,METHOD)
      USE MOD_SPACE_GRID_V2
      IMPLICIT NONE
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Calculates dtau using a formula based on the Euler-Maclaurin
! summation formula (Auer 84, -> Knuth 68)
! Determine differences: dtau == dtau(d) == tau(d+1)-tau(d)
!
! written   11-94  DLM
! altered   12-95  DLM Added abs(z) so a positive tau and dtau would
!                      be found for negative z (ie s_m).
! Altered 2/21/96  DLM Updated to F90 standard
! altered 5/27/97  DLM Removed include "nd_parameters.f".  Changed to
!                     allocate statements.
! altered 6/20/97  DLM Corrected logic for inward and outward direction.
!                      Previously used absolute values, now check direction
!                      and use separated loops.
! Altered 05/09/10 DJH Changed to V4. Deleted TAU from call, and added METHOD to
!                      call. Code nolonger compute dCHIdZ when using ZERO option.
!
!--------------------------------------------------------------------
!
      INTEGER NZ
      INTEGER IP
      REAL*8 CHI(NZ)
      REAL*8 DTAU_LOC(NZ)
      LOGICAL DO_P_RAY
      CHARACTER(LEN=*) METHOD
!
! Local variables
!
      INTEGER IZ
      REAL*8 T1
      REAL*8 DCHIDZ(NZ)
!
!--------------------------------------------------------------------
!
! Determine dtau for either inward or outward direction
!
      IF(DO_P_RAY)THEN
!
	IF(METHOD .EQ. 'ZERO')THEN
	  DCHIDZ(1:NZ)=0.0D0
          DO IZ=1,NZ-1
            DTAU_LOC(IZ)=0.5D0*(RAY(IP)%S_P(IZ)-RAY(IP)%S_P(IZ+1))*(CHI(IZ)+CHI(IZ+1))
          END DO
	ELSE
!
! For outward direction, determine the derivative of chi
!
          DCHIDZ(1)=(CHI(1)-CHI(2))/(RAY(IP)%S_P(1)-ray(ip)%s_p(2))
          DO IZ=2,NZ-1
            DCHIDZ(IZ)=(CHI(IZ-1)-CHI(IZ+1))/(RAY(IP)%S_P(IZ-1)-RAY(IP)%S_P(IZ+1))
          ENDDO
          DCHIDZ(NZ)=(CHI(NZ-1)-CHI(NZ))/(RAY(IP)%S_P(NZ-1)-RAY(IP)%S_P(NZ))
!
! Determine dtau
!
          DO IZ=1,NZ-1
            DTAU_LOC(IZ)=0.5D0*(RAY(IP)%S_P(IZ)-RAY(IP)%S_P(IZ+1))*(CHI(IZ)+CHI(IZ+1))
	    T1=( (RAY(IP)%S_P(IZ)-RAY(IP)%S_P(IZ+1))**2 )*(DCHIDZ(IZ+1)-DCHIDZ(IZ))/12.0D0
            IF(T1 .LT. 0.0D0)THEN
	      DTAU_LOC(IZ)=MAX(0.5D0*DTAU_LOC(IZ),DTAU_LOC(IZ)+T1)
	    ELSE
	      DTAU_LOC(IZ)=MIN(1.5D0*DTAU_LOC(IZ),DTAU_LOC(IZ)+T1)
	    END IF
	  END DO
	END IF
!
      ELSE
	IF(METHOD .EQ. 'ZERO')THEN
	  DCHIDZ(1:NZ)=0.0D0
          DO IZ=1,NZ-1
            DTAU_LOC(IZ)=0.5D0*(RAY(IP)%S_M(IZ+1)-RAY(IP)%S_M(IZ))*(CHI(IZ+1)+CHI(IZ))
          END DO
	ELSE
!
! For inward, determine the derivative of chi
!
	  DCHIDZ(1)=(CHI(2)-CHI(1))/(RAY(IP)%S_M(2)-RAY(IP)%S_M(1))
          DO IZ=2,NZ-1
            DCHIDZ(IZ)=(CHI(IZ+1)-CHI(IZ-1))/(RAY(IP)%S_M(IZ+1)-RAY(IP)%S_M(IZ-1))
          END DO
          DCHIDZ(NZ)=(CHI(NZ)-CHI(NZ-1))/(RAY(IP)%S_M(NZ)-RAY(IP)%S_M(NZ-1))
!
! Determine dtau
!
          DO IZ=1,NZ-1
            DTAU_LOC(IZ)=0.5D0*(RAY(IP)%S_M(IZ+1)-RAY(IP)%S_M(IZ))*(CHI(IZ+1)+CHI(IZ))
     	    T1=((RAY(IP)%S_M(IZ+1)-RAY(IP)%S_M(IZ))**2)*(DCHIDZ(IZ)-DCHIDZ(IZ+1))/12.0D0
            IF(T1 .LT. 0.0D0)THEN
	      DTAU_LOC(IZ)=MAX(0.5D0*DTAU_LOC(IZ),DTAU_LOC(IZ)+T1)
	    ELSE
	      DTAU_LOC(IZ)=MIN(1.5D0*DTAU_LOC(IZ),DTAU_LOC(IZ)+T1)
	    END IF
          END DO
	END IF
!
      END IF
!
      RETURN
      END
