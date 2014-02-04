C
C Created Jan-9189
C Revised 02-Feb-1989 - Variables A,B etc made REAL*8. Check for
C                       division by zero installed.
C
C Subroutine to perform an NG acceleration on estimates obtained
C using an operator which has linear convergence. The method
C is described in detail by Auer (p101, Numerical Radiative Transfer).
C
C IF the logical variable WEIGHT is true,  weighting inversely proportional
C to the value is used.
C
	SUBROUTINE NGACCEL(RJ,PREVRJ,ND,WEIGHT)
	IMPLICIT NONE
C
C Altered 10-Apr-2006 - Check if weight is zero (possible for sigma)
C Altered 24-May-1996 - ERROR_LU installed.
C Altered 23-Jun-1989 - Test for singularity (relative to machine
C                       accuracy) installed. Before just checked whether
C                       DIVISOR was zero. Computation of D1 and D2 improved
C                       to help preserve precsion.
C
	INTEGER ND
	REAL*8 RJ(ND),PREVRJ(4,ND)
	LOGICAL WEIGHT,EQUAL
C
	INTEGER I
	REAL*8 A,B,A1,B1,B2,C1,C2,D0,D1,D2,W,DIVISOR
	REAL*8 DIV1,DIV2,PRECIS,X02AJF
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C X02AJF returns the machine precision.
C
	PRECIS=100.0D0*X02AJF()
C
	A1=0.0D0
	B1=0.0D0
	B2=0.0D0
	C1=0.0D0
	C2=0.0D0
	LUER=ERROR_LU()
C
	DO I=1,ND
	  W=1.0D0
          IF(WEIGHT .AND. PREVRJ(1,I) .NE. 0.0D0)W=1.0D0/PREVRJ(1,I)
C
	  D0=PREVRJ(1,I)-PREVRJ(2,I)
	  D1=D0 + (PREVRJ(3,I)-PREVRJ(2,I))
	  D2=D0 + (PREVRJ(4,I)-PREVRJ(3,I))
	  A1=A1+W*D1*D1
	  B1=B1+W*D1*D2
	  C1=C1+W*D0*D1
	  B2=B2+W*D2*D2
	  C2=C2+W*D0*D2
	END DO
C
C Machine checks that DIVISOR is know with sufficient accuracy for us to
C be able to compute both constants A and B. If matrix is singular, one
C of the constants is set to zero. The special case of no solution can only
C happen in very special circumstance (i.e values all constant).
C
C The DIVISOR (i.e.DETERMINAT) may be zero (or nearly so) if succeeding
C iterations are related exactly by the SAME geometric series.
C In this case we can arbitrarly set one of the constants to zero,
C and compute the value of the other.
C
	DIV1=B2*A1
	DIV2=B1*B1
	IF( EQUAL(DIV1,DIV2,PRECIS) )THEN
	  WRITE(LUER,*)'Warning - Singular determinant in NGACCEL'
	  IF(A1 .NE. 0.0D0)THEN
	    A=C1/A1
	    B=0.0D0
	  ELSE IF(B2 .NE. 0.0D0)THEN
	    A=0.0D0
	    B=C2/B2
	  ELSE
	    WRITE(LUER,*)'Error - no solution for C1 and C2 in NGACCEL'
	    A=0.0D0
	    B=0.0D0
	  END IF
	ELSE
	  DIVISOR=DIV1-DIV2
	  A=(B2*C1-B1*C2)/DIVISOR
	  B=(A1*C2-B1*C1)/DIVISOR
	END IF
C
	DO I=1,ND
	  RJ(I)=PREVRJ(1,I)+A*( PREVRJ(2,I)-PREVRJ(1,I) )
	1      +B*( PREVRJ(3,I)-PREVRJ(1,I) )
	END DO
C
	RETURN
	END
