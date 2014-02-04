C Converted to PGPLOT from MONGO by Gregson Vaux, September 1996
C
C ALTERED 31-MAY-2000 - ENCODE removed.
C                       Negative exponents now handled.
C Altered 12-APR-1989 - MGO prefixed to MOMGHO calls.
C
	SUBROUTINE MON_NUM(X,Y,VAL,LOC,NDEC)
	IMPLICIT NONE
	INTEGER NDEC,NCHAR,LOC
	REAL*4 X,Y,VAL,T1,PGLOC
	CHARACTER LABEL*30,FMT*7
	CHARACTER*80 OUT_FRMAT
C
C Convert MONGO LOC to PGPLOT LOC
C
	IF ((LOC .EQ. 1) .OR. (LOC .EQ. 4) .OR. (LOC .EQ. 1)) THEN
	  PGLOC=1.0
	ELSE IF ((LOC .EQ. 2) .OR. (LOC .EQ. 5) .OR. (LOC .EQ. 8)) THEN
	  PGLOC=.5
	ELSE IF ((LOC .EQ. 3) .OR. (LOC .EQ. 6) .OR. (LOC .EQ. 9)) THEN
	  PGLOC=0.0
	END IF

C
C Compute number of digits to be encoded
C
	LABEL=' '
	IF( ABS(VAL) .GT. 1.0D+05)THEN
	  NCHAR=6+NDEC
	  IF(VAL .LT. 0)NCHAR=NCHAR+1
	  WRITE(OUT_FRMAT,'(A5,I2.2,A1,I2.2,A1)')'(1P,E',NCHAR,'.',NDEC,')'
	  WRITE(LABEL,OUT_FRMAT)VAL
	ELSE IF( ABS(VAL) .LT. 1.0D-04 .AND. VAL .NE. 0)THEN
	  NCHAR=6+NDEC
	  IF(VAL .LT. 0)NCHAR=NCHAR+1
	  WRITE(OUT_FRMAT,'(A5,I2.2,A1,I2.2,A1)')'(1P,E',NCHAR,'.',NDEC,')'
	  WRITE(LABEL,OUT_FRMAT)VAL
	ELSE
	  T1=ABS(VAL)
	  NCHAR=1
	  IF(T1 .GT. 1.0)THEN
	    NCHAR=NCHAR+LOG10(T1)
	  END IF
	  IF(VAL .LT. 0)NCHAR=NCHAR+1		!Minus sign
	  NCHAR=NCHAR+1				!+1 for decimal pont
	  IF(NDEC .NE. 0)NCHAR=NCHAR+NDEC
C
	  FMT(1:7)='(F10.1)'
	  WRITE(FMT(3:4),'(I2.2)')NCHAR
	  IF(NDEC .NE. 0)THEN
	    WRITE(FMT(6:6),'(I1)')NDEC
	  ELSE
	    FMT(6:6)='0'
	  END IF
	  WRITE(LABEL,FMT)VAL
C
C Output label centered on the X position, and with the current Y pointer
C at the top of the string.
C
	  IF(NDEC .EQ. 0)NCHAR=NCHAR-1		!Don't print out decimal point.
	END IF
C
	CALL PGPTXT(X,Y,0.0,PGLOC,LABEL(1:NCHAR))

C
	RETURN
	END

