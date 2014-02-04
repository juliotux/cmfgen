!
! Program to evaluate a multiple pseudo-Gaussian fit to a set of
! data points passed by the module GAUS_FIT_DATA.
!
! The function has the form
!
!   Y =P1 + P2*(X-X(1))+ P5*EXP( -((X-P3)/P4)^P6) + ... 
!
! The routine returns the fit (YFIT) and the SQUARED error.
! between the fit and data.
!
! Altered   -Sep-2007 : Use of alternative exponent to 2 installed.
! Created 21-Jul-2005
!
	REAL*8 FUNCTION GAUS_FIT_FUNC(PARAMS)
	USE GAUS_FIT_DATA
	IMPLICIT NONE
	REAL*8 PARAMS(NG_PAR)
!
	REAL*8 SUM
	REAL*8 T1
	INTEGER I,J,K
!
	GAUS_FIT_FUNC=0.0D0
	DO J=1,NG_DATA
	  SUM=PARAMS(1)+PARAMS(2)*(X_GAUS(J)-X_GAUS(1))
	  DO K=3,NG_PAR,4
	   T1=PARAMS(K+2)*EXP(-(ABS((X_GAUS(J)-PARAMS(K))/PARAMS(K+1)))**PARAMS(K+3))
	   IF(ABS(PARAMS(K+3)) .LT. 1.0D0)T1=100.0*T1
	   SUM=SUM+T1
	  END DO
	  YFIT(J)=SUM
	  GAUS_FIT_FUNC=GAUS_FIT_FUNC+(Y_GAUS(J)-SUM)**2
	END DO
!
	RETURN
	END
