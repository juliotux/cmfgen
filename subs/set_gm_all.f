C
C Routine to store the statistical weights, atomic mass, and name for each
C species in 3 vectors, each of LENGTH NT.
C
C G_ALL,MASS_ALL should be zeroed before first call.
C
	SUBROUTINE SET_GM_ALL(G_ALL,MASS_ALL,LEVEL_ID,NT,
	1                 GCI,MASS_CARB,EQCI,CILEVNAME,NCI,CI_PRES)
	IMPLICIT NONE
C
C Altered 26-May-1996 : Variable format <> removed.
C                       Limit of 9999 levels.
C Created 17-May-1995
C
	INTEGER NCI,NT,EQCI
	REAL*8 G_ALL(NT)
	REAL*8 MASS_ALL(NT)
	CHARACTER*(*) LEVEL_ID(NT)
C
	REAL*8 GCI(NCI)
	REAL*8 MASS_CARB
	CHARACTER*(*) CILEVNAME(NCI)
	LOGICAL CI_PRES
C
	INTEGER ICHRLEN,ERROR_LU
	EXTERNAL ICHRLEN,ERROR_LU
C
	INTEGER I,J,K,LUER,ID_LENGTH
C
	IF(CI_PRES)THEN
	  DO I=1,NCI
	    G_ALL(EQCI+I-1)=GCI(I)
	    MASS_ALL(EQCI+I-1)=MASS_CARB
	  END DO
	  IF(CILEVNAME(1) .NE. ' ')THEN
	    DO I=1,NCI
	      LEVEL_ID(EQCI+I-1)=CILEVNAME(I)
	    END DO
	  ELSE
	    DO I=1,NCI
	      LEVEL_ID(EQCI+I-1)=' '
	    END DO
	  END IF
C
	  ID_LENGTH=LEN(LEVEL_ID(1))
	  DO I=1,NCI
     	     J=ICHRLEN(LEVEL_ID(EQCI+I-1))
	     K=LOG10(0.1D0+I)+1
	     IF( J+K+2 .GT. ID_LENGTH)THEN
	       LUER=ERROR_LU()
	       WRITE(LUER,*)'Error in SET_GM_ALL --- LEVEL_ID is too short'
	       WRITE(LUER,*)'Operating on:',CILEVNAME(I)
	     ELSE
	       IF(I .LT. 10)THEN
	         WRITE(LEVEL_ID(EQCI+I-1)(J+1:),'(A,I1,A)')'<',I,'>'
	       ELSE IF(I .LT. 100)THEN
	         WRITE(LEVEL_ID(EQCI+I-1)(J+1:),'(A,I2,A)')'<',I,'>'
	       ELSE IF(I .LT. 1000)THEN
	         WRITE(LEVEL_ID(EQCI+I-1)(J+1:),'(A,I3,A)')'<',I,'>'
	       ELSE IF(I .LT. 10000)THEN
	         WRITE(LEVEL_ID(EQCI+I-1)(J+1:),'(A,I4,A)')'<',I,'>'
	       ELSE
	         LUER=ERROR_LU()
	         WRITE(LUER,*)'Error in SET_GM_ALL --- Index is too big'
	       END IF
	     END IF
	  END DO
	END IF
C
	RETURN
	END
