!
! Program to evaluate the error assoiciated with a multiple pseudo-Gaussian fit 
! to a set of data points passed by the module GAUS_FIT_DATA.
!
! The function has the form
!
!   Y =P1 + P2*(X-X(1))+ P5*EXP( -((X-P3)/P4)^P6) + ... 
!
! This routine must be kept compatible with GAUS_FIT_FUNC
!
! Created 05-Oct-2007
!
	SUBROUTINE GAUS_FIT_ER(PARAMS)
	USE GAUS_FIT_DATA
	IMPLICIT NONE
	REAL*8 PARAMS(NG_PAR)
!
	REAL*8 T1
	REAL*8 SUM
	INTEGER I,J,K
!
! Make sure YFIT is up to date.
!
	DO J=1,NG_DATA
	  SUM=PARAMS(1)+PARAMS(2)*(X_GAUS(J)-X_GAUS(1))
	  DO K=3,NG_PAR,4
	   SUM=SUM+PARAMS(K+2)*EXP(-(ABS((X_GAUS(J)-PARAMS(K))/PARAMS(K+1)))**PARAMS(K+3))
	  END DO
	  YFIT(J)=SUM
	END DO
!
! Conservative error estimate based on fit quality. We use abs value, with 1.0-03 of
! line center. Possible influence of multiple lines is ignored.
!
	IF(ALLOCATED(EW_ERROR))DEALLOCATE(EW_ERROR,ALT_ERROR,MIN_ERROR)
	ALLOCATE (EW_ERROR(NUM_GAUS),ALT_ERROR(NUM_GAUS),MIN_ERROR(NUM_GAUS))
	EW_ERROR(:)=0.0D0; ALT_ERROR(:)=0.0D0; MIN_ERROR(:)=0.0D0
	DO J=1,NG_DATA
	  DO I=1,NUM_GAUS
	    K=3+4*(I-1)
	    T1=EXP(-(ABS((X_GAUS(J)-PARAMS(K))/PARAMS(K+1)))**PARAMS(K+3))
	    IF(T1 .GT. 1.0D-03)THEN
	      EW_ERROR(I)=EW_ERROR(I)+(X_GAUS(MIN(J+1,NG_DATA))-X_GAUS(MAX(1,J-1)))*
	1                ABS(Y_GAUS(J)-YFIT(J))
	    END IF
	    IF(ABS((X_GAUS(J)-PARAMS(K))/PARAMS(K+1)) .LT. 4)THEN
	      ALT_ERROR(I)=ALT_ERROR(I)+(X_GAUS(MIN(J+1,NG_DATA))-X_GAUS(MAX(1,J-1)))*
	1                ABS(Y_GAUS(J)-YFIT(J))
	      MIN_ERROR(I)=MIN_ERROR(I)+(X_GAUS(MIN(J+1,NG_DATA))-X_GAUS(MAX(1,J-1)))*
	1                (Y_GAUS(J)-YFIT(J))
	    END IF
	  END DO
	END DO
	EW_ERROR=0.5D0*EW_ERROR
	ALT_ERROR=0.5D0*ALT_ERROR
	MIN_ERROR=0.5D0*ABS(MIN_ERROR)
!
	RETURN
	END
