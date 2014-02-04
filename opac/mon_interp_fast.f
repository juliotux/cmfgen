C
C Subroutine to interpolate an array onto a new grid. The grid vector must be 
C either a monotonically decreasing or increasing function. A modified cubic
C polynomial is used to do the interpolation. Instead of using
C the exact cubic estimates for the first derivative at the two nodes,
C we use revised estimates which insure that the interpolating function
C is monotonic in the interpolating interval.
C
C Disadvantages: The interpolating weights can only be defined when the
C                function is known. In principal could use these modified
C                first derivatives to compute an accurate integration
C                formulae. However, the integration weights cannot be defined
C                independently of the function values, as desired in many
C                situations.
C
C Ref: Steffen. M, 1990, A/&A, 239, 443-450
C
	SUBROUTINE MON_INTERP_FAST(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND)
	IMPLICIT NONE
C
C Created 09-Dec-1998 : Based on MON_INTERP.
C                       Design to be fast when from creating a large 
C                           array from a much smaller array.
C
	INTEGER*4 NQ,LIN_END,NX,NV,ND
	REAL*8 QZ(NQ,LIN_END),QZR(NX)
	REAL*8 VARRAY(NV,LIN_END),R(ND)
C
	REAL*8 S(ND)		!Slopes
	REAL*8 H(ND)
	INTEGER*4 LST_INTERVAL
	INTEGER*4 ND_SM
	INTEGER*4 IVEC(NX)
C
	REAL*8 ONE
	PARAMETER (ONE=1.0D0)
	INTEGER*4 I,J,ML
	REAL*8 T1
	REAL*8 A(ND)
	REAL*8 B(ND)
	REAL*8 C(ND)
	REAL*8 D(ND)		!Used for derivative at I.
	REAL*8 E(ND)
	REAL*8 SGN
C
	INTEGER*4 ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C The array R may be either monotonically increasing, or decreasing.
C
	SGN=SIGN(ONE,R(ND)-R(1))
	IF( (SGN*QZR(1) .LT. SGN*R(1)) .OR.
	1   (SGN*QZR(NX) .GT. SGN*R(ND)) )THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in MON_INTERP - values outside range'
	  STOP
	END IF
C
C Determine intervals and slopes to minimize computational effort.
C
	DO I=1,ND-1
	  H(I)=R(I+1)-R(I)
	END DO
C
C Check that R and QZR vectors are monotonic.
C
	DO I=1,ND-2
	  IF(H(I)*H(I+1) .LE. 0)THEN
	    WRITE(LUER,*)'Error in MON_INTERP_FAST'
	    WRITE(LUER,*)'R values must be monotonic'
	    STOP
	  END IF
	END DO
	DO I=1,NX-2
	  IF( (QZR(I+2)-QZR(I+1))*(QZR(I+1)-QZR(I)) .LE. 0)THEN
	    WRITE(LUER,*)'Error in MON_INTERP_FAST'
	    WRITE(LUER,*)'QZR values must be monotonic'
	    STOP
	  END IF
	END DO
C
C Determine the interval (R(I) to R(I+1)) containing QZR(J).
C
	I=1
	DO J=1,NX
500	  IF( SGN*QZR(J) .LE. SGN*R(I+1))THEN
	    IVEC(J)=I
	  ELSE
	    I=I+1
	    GOTO 500
	  END IF
	END DO
	LST_INTERVAL=MAXVAL(IVEC)
	ND_SM=LST_INTERVAL+1
C
C Loop over frequency space.
C
	DO ML=1,LIN_END
C
C Compute the slopes.
C
	  DO I=1,MIN(ND-1,ND_SM)
	    S(I)=(VARRAY(I+1,ML)-VARRAY(I,ML))/H(I)
	  END DO
C
C Compute the first derivatives at node I. 
C
          D(1)=S(1) +(S(1)-S(2))*H(1)/(H(1)+H(2))
	  DO I=2,MIN(ND-1,ND_SM)
            D(I)=(S(I-1)*H(I)+S(I)*H(I-1))/(H(I-1)+H(I))
	  END DO
C             
C Adjust first derivatives so that function is monotonic  in each interval.
C
	  D(1)=( SIGN(ONE,S(1))+SIGN(ONE,D(1)) )*MIN(ABS(S(1)),0.5*ABS(D(1)))
	  DO I=2,MIN(ND-1,ND_SM)
	    D(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5*ABS(D(I)))
	  END DO
C
C Special treatment for computing slope at end point of R.
C
	  IF(ND .EQ. ND_SM)THEN
            D(ND)=S(ND-1)+(S(ND-1)-S(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
	    D(ND)=( SIGN(ONE,S(ND-1))+SIGN(ONE,D(ND)) )*
	1      MIN(ABS(S(ND-1)),0.5*ABS(D(ND)))
	  END IF
C
C Determine the ciefficients of the monotonic cubic polynomial.
C
C If T1=X-R(I) then
C             Y=A(I)*T1^3 + B(I)*T1^3 + C(I)*T1 + E(I)
C
	  DO I=1,ND_SM-1
            A(I)=(D(I)+D(I+1)-2.0*S(I))/H(I)/H(I)
	    B(I)=(3.0*S(I)-2.0*D(I)-D(I+1))/H(I)
	    C(I)=D(I)
	    E(I)=VARRAY(I,ML)
	  END DO
C
C Perform the interpolations.
C
	  DO J=1,NX
	    I=IVEC(J)
	    T1=(QZR(J)-R(I))
            QZ(J,ML)=((A(I)*T1+B(I))*T1+C(I))*T1+E(I)
	  END DO
C
	END DO
C
	RETURN
	END
