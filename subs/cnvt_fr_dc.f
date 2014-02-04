C
C Routine to convert from Departure coefficients to POPULATIONS.
C In addition, the SPECIES population (==SUM) is incremented.
C The ION contribution is only added in if the next ionization
C state is not present.
C
C The ion population for the previous ionization stage is also set.
C
C Altered 14-Apr-1989
C
	SUBROUTINE CNVT_FR_DC(C2,C2LTE,DC2,NC2,DCI,SUM,ND,
	1                      FIRST,CIII_PRES)
	IMPLICIT NONE
C
	LOGICAL FIRST,CIII_PRES
	INTEGER NC2,ND
	REAL*8 C2(NC2,ND),C2LTE(NC2,ND),DC2(ND),DCI(ND),SUM(ND)
C
	INTEGER I,J
C
	IF(FIRST)THEN
	  DO J=1,ND
	    SUM(J)=0.0D0
	  END DO
	END IF
C
	DO J=1,ND
	  DO I=1,NC2
	    C2(I,J)=C2(I,J)*C2LTE(I,J)
	    SUM(J)=SUM(J)+C2(I,J)
	  END DO
	END DO
C
	IF(.NOT. CIII_PRES)THEN
	  DO J=1,ND
	    SUM(J)=SUM(J)+DC2(J)
	  END DO
	END IF
C
C Set ion population for next lower ionization stage,
C
	DO J=1,ND
	  DCI(J)=C2(1,J)
	END DO
C
	FIRST=.FALSE.
	RETURN
	END
