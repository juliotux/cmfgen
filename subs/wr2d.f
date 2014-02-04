C
C Routine to write a 2D Matrix out as a "MATRIX" . A maximum of ten
C numbers are written across the page.
C
C Altered 13-Dec-1989 - Implicit none installed. I index written out.
C
	SUBROUTINE WR2D(A,N,M,MES,LU)
	IMPLICIT NONE
	INTEGER N,M,LU
	REAL*8 A(N,M)
	CHARACTER*(*) MES
C
	INTEGER MS,MF,ML,I,J
C
	WRITE(LU,200)MES
200	FORMAT(/,1X,A)
C
	MS=1
	DO 10 ML=0,M-1,10
	  MF=ML+10
	  IF(MF .GT. M)MF=M
	  WRITE(LU,100)
	  DO 20 I=1,N
	    WRITE(LU,110)I,(A(I,J),J=MS,MF)
20	  CONTINUE
	  MS=MS+10
10	CONTINUE
C
100	FORMAT(/)
110	FORMAT(1X,I4,'* ',1P,10E12.4)
	RETURN
	END
