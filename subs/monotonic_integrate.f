!
! Subroutine to comput the integral of a numerically defined function using
! monotonic cubic integration. The grid vector must be 
! either a monotonically decreasing or increasing function. A modified cubic
! polynomial is used to do the interpolation. Instead of using
! the excact cubic estiamtes for the first derivative at the two nodes,
! we use revised estimates which insure that the interpolating function
! is mononotonic in the interpolating interval.
!
! Both the interpolated function and its first derivatives are continuous.
!
! The techniques is somewhat similar to that suggested by Nordulund.
!
! Disadvantages: The interpolating coefificents can only be defined when the
!                function is known. In principal could use these modified
!                first derivatives to compute an accurate integration
!                formulae. However, the integration weights cannot be defined
!                independently of the function values, as desired in many
!                situations.
!
! Ref: Steffen. M, 1990, A/&A, 239, 443-450
!
	SUBROUTINE MONOTONIC_INTEGRATE(INTEGRAL,CHI,R,ND)
	IMPLICIT NONE
!
! Created 23-Oct-2006 - Based on MON_INT_FUNS_V2.
!
	INTEGER ND
	REAL*8 INTEGRAL(ND)
	REAL*8 CHI(ND)
	REAL*8 R(ND)
!
	REAL*8 ONE
	PARAMETER (ONE=1.0D0)
	INTEGER I
!
	REAL*8 H(ND)			!Delta R [ R(I+1)-R(I) ]
	REAL*8 S(ND)			!Slope in interval (I to I+1)
	REAL*8 D(ND)			!First derivative at node I
!
! The array R may be either monotonically increasing, or decreasing.
!
	DO I=1,ND-1
	  H(I)=R(I+1)-R(I)
	  S(I)=(CHI(I+1)-CHI(I))/H(I)
	END DO
!
! Compute the first derivatives at node I.
!
        D(1)=S(1) +(S(1)-S(2))*H(1)/(H(1)+H(2))
	DO I=2,ND-1
          D(I)=(S(I-1)*H(I)+S(I)*H(I-1))/(H(I-1)+H(I))
	END DO
        D(ND)=S(ND-1)+(S(ND-1)-S(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
!
! Adjust first derivatives so that function is monotonic  in each interval.
!
	D(1)=( SIGN(ONE,S(1))+SIGN(ONE,D(1)) )*MIN(ABS(S(1)),0.5D0*ABS(D(1)))
	DO I=2,ND-1
	  D(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5D0*ABS(D(I)))
	END DO
	D(ND)=( SIGN(ONE,S(ND-1))+SIGN(ONE,D(ND)) )*
	1        MIN(ABS(S(ND-1)),0.5D0*ABS(D(ND)))
!
! Now do the integration.
!
	INTEGRAL(1)=0.5D0*H(1)*( CHI(1)+CHI(2)+H(1)*(D(1)-D(2))/6.0D0 )
	DO I=2,ND-1
	  INTEGRAL(I)=INTEGRAL(I-1)+0.5D0*H(I)*( CHI(I)+CHI(I+1)+H(I)*(D(I)-D(I+1))/6.0D0 )
	END DO
!
	RETURN
	END
