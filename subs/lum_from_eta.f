C
C Routine to compute the luminosity arising from lines (non-blanketed mode),
C dIELECTRONIC lines (non-blankted mode) or from shock X_RAY emisison.
C On entry LINE_LUM should contain the EMISSIVITY at radius R corrected for the
C LOCAL escape probablity.
C We compute the amount of energy emitted between R(I) and R(I+1), and
C store it in LINE_LUM(I). This is computed using the Euler-McLaurin summation
C rule (i.e. the trapazoidal rule with a correction for the first derivatives).
C
	SUBROUTINE LUM_FROM_ETA(LINE_LUM,R,ND)
	IMPLICIT NONE
C
C Altered 30-Mar-2000 : Bug fixed with derivative computation. Derivatives were
C                         offset.
C Created 25-Jun-1998
C
	INTEGER ND
	REAL*8 LINE_LUM(ND)
	REAL*8 R(ND)
C
	REAL*8 DERIV(ND)
	INTEGER I
C
C We compute the amount of energy emitted between R(I) and R(I+1), and
C store it in LINE(I). This is computed using the Euler-McLaurin summation
C rule (i.e. the trapazoidal rule with a correction for the first derivatives).
C
	DERIV(1)=(LINE_LUM(1)-LINE_LUM(2))/(R(1)-R(2))
	DO I=2,ND-1
	  DERIV(I)=(LINE_LUM(I-1)-LINE_LUM(I+1))/(R(I-1)-R(I+1))
	END DO
	DERIV(ND)=(LINE_LUM(ND-1)-LINE_LUM(ND))/(R(ND-1)-R(ND))
C
	DO I=1,ND-1
	  LINE_LUM(I)=0.5D0*(R(I)-R(I+1))*( LINE_LUM(I)+LINE_LUM(I+1) 
	1            +(R(I)-R(I+1))*(DERIV(I+1)-DERIV(I))/6.0D0 )
	END DO
	LINE_LUM(ND)=0.0D0
C
	RETURN
	END
