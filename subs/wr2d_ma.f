!
! Routine to write a 2D Matrix out as a "MATRIX" . A maximum of five
! numbers are written across the page. This routine gives greater
! output accuracy than WR2D.
!
! Created 06-Feb-2004
!
	SUBROUTINE WR2D_MA(A,N,M,MES,LU)
	IMPLICIT NONE
	INTEGER N,M,LU
	REAL*8 A(N,M)
	CHARACTER*(*) MES
!
	INTEGER MS,MF,ML,I,J
!
	WRITE(LU,200)MES
200	FORMAT(/,1X,A)
!
	MS=1
	DO 10 ML=0,M-1,5
	  MF=ML+5
	  IF(MF .GT. M)MF=M
	  WRITE(LU,100)
	  DO 20 I=1,N
	    WRITE(LU,110)I,(A(I,J),J=MS,MF)
20	  CONTINUE
	  MS=MS+5
10	CONTINUE
!
100	FORMAT(/)
110	FORMAT(1X,I4,'* ',1P,5E24.14)
	RETURN
	END
