C
C Subroutine to perform integration using a fit by a cubic polynomial.
C Based on NAG routin D01GAF which is in turn based on the work of
C Gill and Miller.
C
	SUBROUTINE INTEGRATE(X,F,TA,IFAIL,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 - ERROR_LU installed.
C Altered 07-Jan-1991 - Main loop made vectozable by separating computation
C                     of total integrand.
C Altered 29-Jan-1987 - Error handling for non monotonic function inserted.
C Created 29-APR-1985
C
	INTEGER ND,I,IFAIL
	REAL*8 X(ND),F(ND),TA(ND)
	REAL*8 H0,H1,H2,F01,F12,F23,F012,F123,F0123
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	IF(ND .LT. 4)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in INTEGRATE - not enough points'
	  IFAIL=1
	  RETURN
	END IF
C
	IF(X(ND) .GT. X(1))THEN
	  DO I=1,ND-1
	    IF(X(I+1) .LT. X(I))THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in INTEGRATE - X function not monotonic'
	      IFAIL=2
	      RETURN
	    END IF
	  END DO
	ELSE
	  DO I=1,ND-1
	    IF(X(I+1) .GT. X(I))THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in INTEGRATE - X function not monotonic'
	      IFAIL=2
	      RETURN
	    END IF
	  END DO
	END IF
C
C Initialize values.
C
	IFAIL=0
	H0=X(2)-X(1)
	H1=X(3)-X(2)
	H2=X(4)-X(3)
C
	F01=(F(2)-F(1))/H0
	F12=(F(3)-F(2))/H1
	F23=(F(4)-F(3))/H2
C
	F012=(F12-F01)/(X(3)-X(1))
	F123=(F23-F12)/(X(4)-X(2))
C
	F0123=(F123-F012)/(X(4)-X(1))
C
C Compute integral for first step.
C
	TA(1)=H0*(F(1)+0.5D0*H0*(F01-H0*(F012/3.0D0-(H0*0.5D0
	1             +(H1-H0)/3.0D0)*F0123)))
C
C Compute integral for all subsequent steps. We first evaluate DINT
C (but store in TA) so that loop is vectorizable. The total integrand
C is subsequently formed in a non-vectorizable do loop.
C
	DO I=2,ND-2
C
	  IF( I .NE. 2)THEN
	    H0=H1
	    H1=H2
	    F01=F12
	    F12=F23
	    F012=F123
	  END IF
C
	  H2=(X(I+2)-X(I+1))
	  F23=(F(I+2)-F(I+1))/H2
	  F123=(F23-F12)/(X(I+2)-X(I))
	  F0123=(F123-F012)/(X(I+2)-X(I-1))
C
	  TA(I)=0.5D0*H1*(F(I)+F(I+1)-H1*H1/6.0D0*
	1            (F012+F123+(H0-H2)*F0123))
C
	END DO
C
	DO I=2,ND-2
	  TA(I)=TA(I)+TA(I-1)
	END DO
C
C Compute integral for last step.
C
	TA(ND-1)=TA(ND-2)+H2*(F(ND)-H2*(0.5D0*F23+H2*(F123/6.0D0
	1          +F0123*(0.25D0*H2+(H1-H2)/6.0D0))))
C
	IFAIL=0
C
	RETURN
	END
