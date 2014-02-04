C
C Subroutine to define the cubic interpolation coefficent to interpolate a 
C vecto onto a new grid. The grid vector must be 
C either a monotonically decreasing or increasing function. A modified cubic
C polynomial is used to do the interpolation. Instead of using
C the excact cubic estiamtes for the first derivative at the two nodes,
C we use revised estimates which insure that the interpolating function
C is mononotonic in the interpolating interval.
C
C Both the interpolated function and its first derivatives are continuous.
C
C The techniques is somewhat similar to that suggested by Nordulund.
C
C Disadvantages: The interpolating coefificents can only be defined when the
C                function is known. In principal could use these modified
C                first derivatives to compute an accurate integration
C                formulae. However, the integration weights cannot be defined
C                independently of the function values, as desired in many
C                situations.
C
C Ref: Steffen. M, 1990, A/&A, 239, 443-450
C
	SUBROUTINE MON_INT_FUNS_V2(COEF,CHI,R,ND)
	IMPLICIT NONE
C
C Altered 16-Jan-2005 - COEF(ND,1:4) was previously zero. For ease of use
C                         elsewhere, COEF(ND,4) now contains CHI(ND) and
C                         COEF(ND,3) contains dCHIdR for the last point. The
C                         other two coefficients are zero. These two
C                         assignments are consistent with those at other depths.
C                         The other coefficients are left zero. CHI(ND,:) should
C                         not be used for function fitting.
C Altered 25-Sep-1997 - Rewritten for speed and to vectorize efficiently.
C Created 25-Mar-1996 - Based on MON_INTERP
C
	INTEGER ND
	REAL*8 COEF(ND,4)
	REAL*8 R(ND)
	REAL*8 CHI(ND)
C
	REAL*8 ONE
	PARAMETER (ONE=1.0D0)
	INTEGER I
C
	REAL*8 H(ND)			!Delta R [ R(I+1)-R(I) ]
	REAL*8 S(ND)			!Slope in interval (I to I+1)
	REAL*8 D(ND)			!First derivative at node I
C
C The array R may be either monotonically increasing, or decreasing.
C
	DO I=1,ND-1
	  H(I)=R(I+1)-R(I)
	  S(I)=(CHI(I+1)-CHI(I))/H(I)
	END DO
C
C Compute the first derivatives at node I.
C
        D(1)=S(1) +(S(1)-S(2))*H(1)/(H(1)+H(2))
	DO I=2,ND-1
          D(I)=(S(I-1)*H(I)+S(I)*H(I-1))/(H(I-1)+H(I))
	END DO
        D(ND)=S(ND-1)+(S(ND-1)-S(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
C
C Adjust first derivatives so that function is monotonic  in each interval.
C
	D(1)=( SIGN(ONE,S(1))+SIGN(ONE,D(1)) )*MIN(ABS(S(1)),0.5D0*ABS(D(1)))
	DO I=2,ND-1
	  D(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5D0*ABS(D(I)))
	END DO
	D(ND)=( SIGN(ONE,S(ND-1))+SIGN(ONE,D(ND)) )*
	1        MIN(ABS(S(ND-1)),0.5D0*ABS(D(ND)))
C
C Determine the coeffients of the monotonic cubic polynmial.
C
C If T1=X-R(I) then
C             Y=COEF(I,1)*T1^3 + COEF(I,2)*T1^3 +COEF(I,3)*T1 +COEF(I,4)
C
	DO I=1,ND-1
          COEF(I,1)=(D(I)+D(I+1)-2.0D0*S(I))/H(I)/H(I)
	  COEF(I,2)=(3.0D0*S(I)-2.0D0*D(I)-D(I+1))/H(I)
	  COEF(I,3)=D(I)
	  COEF(I,4)=CHI(I)
	END DO
	COEF(ND,1:2)=0.0D0
	COEF(ND,3)=D(ND)
	COEF(ND,4)=CHI(ND)
C
	RETURN
	END
