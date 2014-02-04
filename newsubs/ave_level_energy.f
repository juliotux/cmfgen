!
! Routine to compute the Average energy of a super-level. The temperature is
! not taken into account.
!
	SUBROUTINE AVE_LEVEL_ENERGY(AVE_ENERGY,EDGE,STAT_WT,F_TO_S,
	1                  EQHYD,N_S,N_F,NT,HYD_PRES)
	IMPLICIT NONE
!
! Input
	LOGICAL HYD_PRES
	INTEGER EQHYD,NT
	INTEGER N_S,N_F
	REAL*8 EDGE(N_F)
	REAL*8 STAT_WT(N_F)
	INTEGER F_TO_S(N_F)
!
! Output

	REAL*8 AVE_ENERGY(NT)
!
! Local variables
!
	REAL*8 G_SUM(N_S)
	REAL*8 EDGE_SUM(N_S)
	INTEGER I,J
!
	IF(HYD_PRES)THEN
	  G_SUM(:)=0.0D0
	  EDGE_SUM(:)=0.0D0
	  DO J=1,N_F
	    I=F_TO_S(J)
            G_SUM(I)=G_SUM(I)+STAT_WT(J)
	    EDGE_SUM(I)=EDGE_SUM(I)+EDGE(J)*STAT_WT(J)
	  END DO
	  DO I=1,N_S
	    AVE_ENERGY(EQHYD+I-1)=EDGE_SUM(I)/G_SUM(I)
	  END DO
	END IF
!
	RETURN
	END
