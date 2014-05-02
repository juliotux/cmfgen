	SUBROUTINE SET_DC_LTE_V2(XzVLTE,DXzV,EDGE,NXzV,T,TMIN,ND)
	IMPLICIT NONE
!
! Altered  9-Dec-2013:  Changed way code handles TMIN_LOC -- still needs development.
! Created 24-Feb-2013.  Returns Log (dep. coef.)
!                       Based on SET_DC_LTE
!
	INTEGER NXzV
	INTEGER ND
	REAL*8 XzVLTE(NXzV,ND)
	REAL*8 DXzV(ND)
	REAL*8 T(ND)
	REAL*8 EDGE(NXzV)
	REAL*8 TMIN
!
! Constants for opacity etc. These are set in CMFGEN.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 TMIN_LOC
	INTEGER I,J
!
	TMIN_LOC=TMIN*(1.0D0+(MAX(0.0D0,(EDGE(1)-3.5D0))/20.0D0))
	DO I=1,ND
	  DXzV(I)=1.0D0
	  IF(T(I) .GT. TMIN_LOC)THEN
	    XzVLTE(:,I)=0.0D0
	  ELSE
	    DO J=1,NXzV
	      XzVLTE(J,I)=HDKT*EDGE(J)*(1.0D0/TMIN_LOC-1.0D0/T(I))
	    END DO
	  END IF
	END DO
!
	RETURN
	END
