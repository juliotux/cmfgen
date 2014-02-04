C
C Subroutine to compute the coefficient matrix (dimension ND*NX) for
C the variation in DTAU with opacity. If the trapazoidal rule is being used
C compute the opacity, the matrix is upper-bi-diagonal, and is of dimension
C ND*ND. If we are correcting the trapazoidal rule using the first
C derivatives as indicated by Nordulund, the matrix has four non zero
C "diagonals" - 2 above the center diagonal, and one below.
C
C Note that
C	            A(I)=d(dChIdr)/d[CHI(I-1)]  at I
C	            B(I)=d(dChIdr)/d[CHI(I)]    "  "
C	            C(I)=d(dChIdr)/d[CHI(I+1)]  "  "
C
C and Q is the spherical optical depth fcator.
C
C
	SUBROUTINE dSPHEREdCHI(W,DTAU,R,Q,ND)
	USE MOD_TRAP_DERIVATIVES
	IMPLICIT NONE
	INTEGER ND
	REAL*8 W(ND,ND),DTAU(ND),Q(ND),R(ND)
C
C Altered 02-Mar-1999 - Module MOD_TRAP_DERIVATIVES replaces COMMON block
C                          TRAPDERIVATIVES. Variable names remain same.
C Altered 24-May-1996 - Call to DP_ZERO removed.
C Created 27-Apr-1989 - Based on NEWVKIMD
C
	INTEGER I,K
	REAL*8 ALPHA,BETA
C
	W(:,:)=0.0D0
C
	DO I=1,ND-1
	  K=I+1
	  ALPHA=0.5D0*(R(I)-R(K))
	  BETA=ALPHA*(R(I)-R(K))/6.0D0
	  IF(I .NE. 1)W(I,I-1)=-Q(I-1)*BETA*A(I)
	  W(I,I)=Q(I)*( ALPHA + BETA*(A(K)-B(I)) )
	  W(I,K)=Q(K)*( ALPHA + BETA*(B(K)-C(I)) )
	  IF(K .NE. ND)W(I,K+1)=Q(K+1)*BETA*C(K)
	END DO
C
	RETURN
	END
