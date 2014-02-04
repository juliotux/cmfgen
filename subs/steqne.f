C
C Routine to increment the Charge conservation equation, and the
C variation charge equation. This was originally done in
C steqheii.
C
	SUBROUTINE STEQNE(BA,STEQ,ED,EQNE,NT,NUM_BNDS,ND,DST,DEND)
	IMPLICIT NONE
C
C Altered 26-May-1996 :  ERROR_LU installed.
C Altered 05-Oct-1989 :  DST,DEND variable installed to avoid reading in
C                          enitre BA matrix for each ion.
C Created 15-Feb-1989
C
	INTEGER EQNE,NT,NUM_BNDS,ND,DST,DEND
C
	REAL*8 BA(NT,NT,NUM_BNDS,ND)
	REAL*8 STEQ(NT,ND),ED(ND)
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C Local varaiables.
C
	INTEGER K,M
C
C NB - all S. E. routines assume EQNE = NT-1
C
	IF(EQNE .NE. NT-1)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error - EQNE must equal NT-1'
	  STOP
	END IF
C
	DO K=DST,DEND
	  IF(ND .EQ. NUM_BNDS)THEN
	    M=K
	  ELSE
	    M=(NUM_BNDS/2)+1
	  END IF
	  BA(NT-1,NT-1,M,K)=BA(NT-1,NT-1,M,K)-1.0D0
	  STEQ(NT-1,K)=STEQ(NT-1,K)-ED(K)
	END DO
C
	RETURN
	END
