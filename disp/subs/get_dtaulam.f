C
C Subroutine to compute the derivative of DTAU(I) with respect to 
C CHI(I) AND CHI(K). Routine is to be used only when computing the
C diagonal LAMBDA operator, and is not a general linearization
C routine. Only the required subset of the derivatives are computed.
C Handles the general integration rule.
C
C Note that
C	            A(I)=d(dChIdr)/d[CHI(I-1)]  at I
C	            B(I)=d(dChIdr)/d[CHI(I)]    "  "
C	            C(I)=d(dChIdr)/d[CHI(I+1)]  "  "
C are contained in module MOD_TRAP_DERIVATIVES.
C
	SUBROUTINE GET_DTAULAM(WA,WB,R,Z,NI)
	USE MOD_TRAP_DERIVATIVES
	IMPLICIT NONE
C
C Altered 02-Mar-1999. Module MOD_TRAP_DERIVARTIVES replace COMMON block
C                        TRAPDERIVATIVES. Same variable names (A, B and C).
C Created 04-May-1989.
C
	INTEGER NI
	REAL*8 WA(NI),WB(NI),Z(NI),R(NI)
C
	INTEGER I,K
	REAL*8 ALPHA
C
	DO I=1,NI
	  WA(I)=0.0D0
	  WB(I)=0.0D0
	END DO
C
C NB - WB(I) = dDTAU(I)/dCHI(I)
C    - WA(I) = dDTAU(I-1)/dCHI(I)
C
	DO I=1,NI-1
	  K=I+1
	  ALPHA=0.5D0*(Z(I)-Z(K))
	  WB(I)=ALPHA*(  1.0D0+(Z(I)-Z(K))*
	1            ( A(K)*Z(K)/R(K) - B(I)*Z(I)/R(I) )/6.0D0  )
	  WA(K)=ALPHA*(  1.0D0+(Z(I)-Z(K))*
	1            ( B(K)*Z(K)/R(K) - C(I)*Z(I)/R(I) )/6.0D0  )
	END DO
C
	RETURN
	END
