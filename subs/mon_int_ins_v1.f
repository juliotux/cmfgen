C
C Subroutine to obtain interpolated values at intermediate points as a function
C a gerneral monotonically increasing or decreasing absisca R. A modified cubic
C polynomial is used to do the interpolation. Instead of using
C the excact cubic estimates for the first derivative at the two nodes,
C we use revised estimates which insure that the interpolating function
C is mononotonic in the interpolating interval.
C
C A special format is used for the points at which the interpolations are
C made. Specifcally design for FG_J_CMF_V7. Use MON_INTERP for other
C formats. The first derivatives can also be computed and returned, at 
C both the nodes, ans the inserted points.
C
C Both the interpolated function and its first derivatives are continuous.
C
C The techniques is somewhat similar to that suggested by Nordulund.
C
C
C Ref: Steffen. M, 1990, A/&A, 239, 443-450
C
	SUBROUTINE MON_INT_INS_V1(CHI_INS,R_INS,NINS,CHI,R,ND,
	1                   LOGX,LOGY,dCHIdR,dCHIdR_INS,DERIV)
	IMPLICIT NONE
C
C Created 26-Sep-1997 - Based on MON_INTERP
C
	INTEGER ND,NINS
C
C NB: These arrays are dimensioned ND-1 since there are ONLY ND-1
C intervals. By doing this we can use MATRIX operations, and operate
C on the whole array.
C
	REAL*8 CHI_INS(ND-1,NINS)
	REAL*8 R_INS(ND-1,NINS)
	REAL*8 dCHIdR_INS(ND-1,NINS)

	REAL*8 R(ND)
	REAL*8 CHI(ND)
	REAL*8 dCHIdR(ND)
C
	LOGICAL LOGX			!Indicates interpolation in LOG(R)
	LOGICAL LOGY			!Indicate interpolation in LOG(Y)
	LOGICAL DERIV			!Indicates to compute derivatives.
C
	REAL*8 X(ND)			!Revised R array
	REAL*8 Y(ND)			!Revised CHI array
C
	REAL*8 H(ND)			!Delta R [ R(I+1)-R(I) ]
	REAL*8 S(ND)			!Slope in interval (I to I+1)
	REAL*8 D(ND)			!First derivative at node I
	REAL*8 COEF(ND,4)
C
	REAL*8 ONE
	REAL*8 DELR
	PARAMETER (ONE=1.0D0)
	INTEGER I,K
C
C Choose the correct plane for the interpolations.
C
	IF(LOGX)THEN
	  X=LOG(R)
	ELSE
	  X=R
	END IF
	IF(LOGY)THEN
	  Y=LOG(CHI)
	ELSE
	  Y=CHI
	END IF
C
C The array X may be either monotonically increasing, or decreasing.
C
	DO I=1,ND-1
	  H(I)=X(I+1)-X(I)
	  S(I)=(Y(I+1)-Y(I))/H(I)
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
	D(1)=( SIGN(ONE,S(1))+SIGN(ONE,D(1)) )*MIN(ABS(S(1)),0.5*ABS(D(1)))
	DO I=2,ND-1
	  D(I)=( SIGN(ONE,S(I-1))+SIGN(ONE,S(I)) )*
	1          MIN(ABS(S(I-1)),ABS(S(I)),0.5*ABS(D(I)))
	END DO
	D(ND)=( SIGN(ONE,S(ND-1))+SIGN(ONE,D(ND)) )*
	1        MIN(ABS(S(ND-1)),0.5*ABS(D(ND)))
C
C Determine the coeffients of the monotonic cubic polynmial.
C
C If T1=X-X(I) then
C             Y=COEF(I,1)*T1^4 + COEF(I,2)*T1^3 +COEF(I,3)*T1^2 +COEF(I,4)
C
	DO I=1,ND-1
          COEF(I,1)=(D(I)+D(I+1)-2.0*S(I))/H(I)/H(I)
	  COEF(I,2)=(3.0*S(I)-2.0*D(I)-D(I+1))/H(I)
	  COEF(I,3)=D(I)
	  COEF(I,4)=Y(I)
	END DO
C
C Now do the interpolations.
C
	DO K=1,NINS
	  DO I=1,ND-1
	    IF(LOGX)THEN
	      DELR=LOG(R_INS(I,K))-X(I)
	    ELSE
	      DELR=R_INS(I,K)-X(I)
	    END IF
	    CHI_INS(I,K)=((COEF(I,1)*DELR+COEF(I,2))*DELR+COEF(I,3))*DELR +
	1                   COEF(I,4)
	  END DO
	END DO
C
	IF(LOGY)CHI_INS=EXP(CHI_INS)
C
	IF(DERIV)THEN
	  DO K=1,NINS
	    DO I=1,ND-1
	      IF(LOGX)THEN
	        DELR=LOG(R_INS(I,K))-X(I)
	      ELSE
	        DELR=R_INS(I,K)-X(I)
	      END IF
	      dCHIdR_INS(I,K)=((3.0D0*COEF(I,1)*DELR+2.0D0*COEF(I,2))*DELR+
	1                   COEF(I,3))
	    END DO
	  END DO
C
C Now need to interpolated quantities if we did not perform the interpolations
C in a LINEAR-LINEAR plane. The instructions invloving INS variables are 2D 
C ARRRAY operations. The other instructions are 1D operations.
C
	  IF(LOGX .AND. LOGY)THEN
	    dCHIdR_INS=CHI_INS*dCHIdR_INS/R_INS
	    dCHIdR=CHI*D/R
	  ELSE IF(LOGX)THEN
	    dCHIdR_INS=dCHIdR_INS/R_INS
	    dCHIdR=D/R
	  ELSE IF(LOGY)THEN
	    dCHIdR_INS=CHI_INS*dCHIdR_INS
	    dCHIdR=CHI*D
	  ELSE
	    dCHIdR=D
	  END IF
	END IF
C
	RETURN
	END
