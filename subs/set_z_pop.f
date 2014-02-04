C
C Routines set the vector ZPOP to contain the charge of each atomic species.
C ZPOP should be initialized before the fist call.
C
	SUBROUTINE SET_Z_POP(ZPOP,ZHYD,EQHYD,NLEV,NT,HYD_PRES)
	IMPLICIT NONE
C
C Created 09-Aug-1994
C
	LOGICAL HYD_PRES
	INTEGER EQHYD,NLEV,NT,I
	REAL*8 ZPOP(NT),ZHYD
C
	IF(HYD_PRES)THEN
	  DO I=EQHYD,EQHYD+NLEV-1
	    ZPOP(I)=ZHYD-1
	  END DO
	  ZPOP(EQHYD+NLEV)=ZHYD
	END IF
C
	RETURN
	END
