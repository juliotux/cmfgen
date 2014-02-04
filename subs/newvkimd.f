C
C Subroutine to compute the coefficient matrix (dimension NI*NX) for
C the variation in opacity. If the trapazoidal rule is being used
C compute the opacity, the matrix is tridiagonal, and is of dimension
C NI*NI. If we are correcting the trapazoidal rule using the first
C derivatives as indicated by Nordulund, the matrix is "pentadiagonal".
C Note that although dCHIdR at NI (LS > NC) depends on CHI(NI+1), we
C do not need to include the varaiation since it is multiplied by
C Z(NI) which is identically zero.
C
C Note that
C	            A(I)=d(dChIdr)/d[CHI(I-1)]  at I
C	            B(I)=d(dChIdr)/d[CHI(I)]    "  "
C	            C(I)=d(dChIdr)/d[CHI(I+1)]  "  "
C are contained in MOD_TRAP_DERIVATIVES.
C
	SUBROUTINE NEWVKIMD(W,DTAU,RKI,S,U,R,Z
	1                          ,DIFF,DBC,LS,NC,ND,NI)
	USE MOD_TRAP_DERIVATIVES
	IMPLICIT NONE
C
C Altered 02-Mar-1999 - Module MOD_TRAP_DERIVATIVES replaces COMMON block
C                          TRAPDERIVATIVES. Variable names remain same.
C Altered 24-May-1996 - DOUBLE PRECISION declaration replaced.
C                       CALL to DP_ZERO removed.
C Altered 14-Apr-1988 - Bug fix. W(I,K) had incorrect trapazoidal correction
C                          term (ALPHA term). Dimension of trapderivatives
C                          arrays extended to 200.
C Created 21-Jan-1988 - Based on VKIMD. Program does full linearization
C                           allowing for the corrections to the trapazoidal rule.
C
	INTEGER LS,NC,ND,NI
	REAL*8 W(NI,NI),DTAU(NI),S(NI),U(NI),Z(NI)
	REAL*8 R(ND),RKI(ND),DBC
	LOGICAL DIFF
C
	INTEGER I,J,K
	REAL*8 ALPHA,BETA
C
	W(:,:)=0.0D0
C
	ALPHA=0.5D0*(Z(1)-Z(2)) * (U(2)-U(1)) /DTAU(1)/DTAU(1)
	W(1,1)=ALPHA*(  1.0D0+(Z(1)-Z(2))/6.0D0*( A(2)*Z(2)/R(2)
	1                     -B(1)*Z(1)/R(1) )  )
	W(1,2)=ALPHA*(  1.0D0+(Z(1)-Z(2))/6.0D0*( B(2)*Z(2)/R(2)
	1                     -C(1)*Z(1)/R(1) )  )
	W(1,3)=ALPHA*(Z(1)-Z(2))*Z(2)*C(2)/R(2)/6.0D0
C
	DO 20 I=2,NI-1
	  K=I+1
	  J=I-1
	  ALPHA=( 0.5D0*(U(I)-S(I)) + (U(J)-U(I))/DTAU(J)/DTAU(J) )
	1           *(Z(J)-Z(I))*0.5D0
	  BETA=( 0.5D0*(U(I)-S(I)) + (U(K)-U(I))/DTAU(I)/DTAU(I) )
	1           *(Z(I)-Z(K))*0.5D0
	  IF(J .NE. 1)W(I,J-1)=-ALPHA*A(J)*Z(J)/R(J)*
	1                   (Z(J)-Z(I))/6.0D0
	  W(I,J)=ALPHA*(  1.0D0+(Z(J)-Z(I))/6.0D0*
	1            ( A(I)*Z(I)/R(I) - B(J)*Z(J)/R(J) )  )
	1            -BETA/6.0D0* (Z(I)-Z(K)) *A(I)*Z(I)/R(I)
	  W(I,K)=BETA*(  1.0D0+(Z(I)-Z(K))/6.0D0*
	1            ( B(K)*Z(K)/R(K) - C(I)*Z(I)/R(I) )  )
	1            +ALPHA/6.0D0* (Z(J)-Z(I)) *C(I)*Z(I)/R(I)
	  W(I,I)=ALPHA*(  1.0D0+(Z(J)-Z(I))/6.0D0*
	1            ( B(I)*Z(I)/R(I) - C(J)*Z(J)/R(J) )  )
	1       +BETA*(  1.0D0+(Z(I)-Z(K))/6.0D0*
	1            ( A(K)*Z(K)/R(K) - B(I)*Z(I)/R(I) )  )
	1       +0.5D0*S(I)*(DTAU(J)+DTAU(I))/RKI(I)
	  IF(K .NE. NI)W(I,K+1)=BETA*C(K)*Z(K)/R(K)*
	1                   (Z(I)-Z(K))/6.0D0
20	CONTINUE
C
	ALPHA=(U(NI)-U(NI-1))/DTAU(NI-1)/DTAU(NI-1)
	IF(LS .GT. NC)ALPHA=ALPHA+0.5D0*( S(NI)-U(NI) )
	ALPHA=0.5D0*(Z(NI-1)-Z(NI))*ALPHA
 
	W(NI,NI-2)=-ALPHA*(Z(NI-1)-Z(NI))*A(NI-1)*Z(NI-1)/R(NI-1)/6.0D0
	W(NI,NI-1)=ALPHA*(  1.0D0 + (Z(NI-1)-Z(NI))*( A(NI)*Z(NI)/R(NI)
	1                 - B(NI-1)*Z(NI-1)/R(NI-1) )/6.0D0  )
	W(NI,NI)=ALPHA*(  1.0D0 + (Z(NI-1)-Z(NI))*( B(NI)*Z(NI)/R(NI)
	1                 - C(NI-1)*Z(NI-1)/R(NI-1) )/6.0D0  )
C
	IF (LS .GT. NC)THEN
	  W(NI,NI)=W(NI,NI)-DTAU(NI-1)*S(NI)/RKI(NI)*0.5D0
	ELSE IF(DIFF)THEN
	  W(NI,NI)=W(NI,NI)-DBC/RKI(NI)
	END IF
C
	RETURN
	END
