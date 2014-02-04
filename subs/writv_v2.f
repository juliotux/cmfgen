C
	SUBROUTINE WRITV_V2(F,ND,NDEC,A,LU)
	IMPLICIT NONE
C
C Altered 28-Apr-2000 : Bug fix: Length of FORM (at 15) was one character 
C                                too short.
C Altered 28-May-1996 : IMPLICIT NONE installed.
C Altered 30-APR-1985 : Now writes 10 columns instead of five across a page.)
C
	INTEGER ND
	INTEGER LU
	INTEGER NDEC		!Number of decimal digits.
	REAL*8 F(ND)
	CHARACTER*(*) A
C
C Local variables.
C
	INTEGER I,J,L
	INTEGER N_PER_LINE
	INTEGER NX
	CHARACTER*16 FORM
C
	NX=NDEC+8
	N_PER_LINE=131/NX
C
	WRITE(FORM,'(A,I2.2,A,I2.2,A,I2.2,A)')
	1             '(1X,1P,',N_PER_LINE,'E',NX,'.',NDEC,')'
C
	L=1
	WRITE(LU,'(//,1X,A,/)')A
C
	IF(ND .GE. N_PER_LINE)THEN
	  DO I=1,ND+1-N_PER_LINE,N_PER_LINE
	    WRITE(LU,FORM)(F(J),J=I,I+N_PER_LINE-1,1)
	    L=I+N_PER_LINE
	  END DO
	END IF
	IF(L .LE. ND)WRITE(LU,FORM)(F(J),J=L,ND)
C
	RETURN
	END
