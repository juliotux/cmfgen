C
C Routine to check whether equataion numbers are consistent with
C the atomic species present, and the number of levels. Must be
C called for each ion.
C
C CV is the current equation numer.
C ATOM_PRES is set TRUE when ion is present. Used to indicate whether
C any ionic species of a particular element (eg carbon) is present.
C
	SUBROUTINE CHK_EQ_NUM(C2_PRES,EQC2,NC2,CV,ATOM_PRES,STRING)
	IMPLICIT NONE
C
C Altered 24-May-1996 : ERROR_LU, LUER inserted
C Created 01-Sep-1989.
C
	INTEGER EQC2,NC2,CV,L,ICHRLEN
	LOGICAL C2_PRES,ATOM_PRES
	CHARACTER*(*) STRING
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	IF(C2_PRES)THEN
	  IF(EQC2 .NE. CV)THEN
	   LUER=ERROR_LU()
	   L=ICHRLEN(STRING)
	   WRITE(LUER,100)string(1:L)
100	   FORMAT(1X,'Error - invalid equation number for ',A)
	   STOP
	  ELSE
	    ATOM_PRES=.TRUE.
	    CV=CV+NC2
	  END IF
	END IF
C
	RETURN
	END
