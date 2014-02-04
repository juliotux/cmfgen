C
C Subroutine to compute the coefficient matrix (dimension ND*ND) for
C the variation in opacity. The matrix is tridiagonal. (This routine
C needs to be checked to see that the access is not too inefficient)
C Uses Schuster of diffusion approximation at lower boundary.
C
	SUBROUTINE VKIFEAU_IBC(W,DTAU,CHI,RJ,TA,TB,TC,R,Q,F,
	1                      ZETA,THETA,HBC_S,DIFF,DBB,ND)
	USE MOD_TRAP_DERIVATIVES
	IMPLICIT NONE
C
C Altered 02-Mar-1999 - Module MOD_TRAP_DERIVATIVES replaces COMMON block
C                          TRAPDERIVATIVES. Variable names (A, B and C) remain 
C                          same.
C Altered 28-MAy-1996 - Call to DP_ZERO removed.
C Altered 12-Jun-1991 - HBC_S boundary condition inserted. Name changed
C                       from VKIFEAUNEW.
C
C Created 15-Apr-1988 - Major rewrite. Based on VKIFEAU and NEWVKIMD.
C                       Corrections trapazoidal rule based on the
C                       first derivative has been installed. Arrays
C                       in trapderivitives extended to 200.
C Altered 26-Feb-1986 (Bug fixed)
C Created 18-Feb-1986
C
	INTEGER ND
	REAL*8 W(ND,ND),DTAU(ND),RJ(ND),R(ND),F(ND),Q(ND),CHI(ND)
	REAL*8 ZETA(ND),THETA(ND),TA(ND),TB(ND),TC(ND)
	REAL*8 HBC_S,DBB
	LOGICAL DIFF
C
C Local varaibles.
C
	INTEGER I,J,K
	REAL*8 ALPHA,BETA,T1,UIJ,UIK
C
	W(:,:)=0.0D0
C
C NB. It is plus HBC_S*... as XVEC (i.e. RHS) is -HBC_S*S*R(1)*R(1)
C
	UIK=(RJ(2)*F(2)*Q(2)*R(2)*R(2)-RJ(1)*F(1)*Q(1)*R(1)*R(1))/DTAU(1)
	ALPHA=0.5D0*(R(1)-R(2))*UIK/DTAU(1)
	W(1,1)=ALPHA*Q(1)*(  1.0D0+(R(1)-R(2))/6.0D0*(A(2)-B(1))  )
	1       +HBC_S*R(1)*R(1)*(ZETA(1)+THETA(1)*RJ(1))/CHI(1)
	W(1,2)=ALPHA*Q(2)*(  1.0D0+(R(1)-R(2))/6.0D0*(B(2)-C(1))  )
	W(1,3)=ALPHA*Q(3)*(R(1)-R(2))*C(2)/6.0D0
C
	DO I=2,ND-1
	  K=I+1
	  J=I-1
	  T1=DTAU(J)+DTAU(I)
C
C Use TA(I) and TC(I) which have been previously computed by TFEAU.
C
C 	  TA(I)=-2.0D0*R(J)*R(J)*F(J)*Q(J)/DTAU(J)/T1
C	  TC(I)=-2.0D0*R(K)*R(K)*F(K)*Q(K)/DTAU(I)/T1
C
C We use UIK for the sum of RHS terms containing 1/DTAU(I)
C We use UIJ for the sum of RHS terms containing 1/DTAU(I)
C Note that UIK+UIJ is the sum of the RHS terms containin 1/T1.
C
	  UIJ=2.0D0*RJ(I)*R(I)*R(I)*F(I)*Q(I)/T1
	  UIK=UIJ/DTAU(I)+TC(I)*RJ(K)
	  UIJ=UIJ/DTAU(J)+TA(I)*RJ(J)
C
C Here ALPHA=d(-RHS)/dDTAU(J) and ALPHA=d(-RHS)/dDTAU(I)
C
	  ALPHA=( UIJ/DTAU(J)+(UIJ+UIK)/T1 )*(R(J)-R(I))*0.5D0
	  BETA= ( UIK/DTAU(I)+(UIJ+UIK)/T1 )*(R(I)-R(K))*0.5D0
C
	  IF(J .NE. 1)W(I,J-1)=-ALPHA*A(J)*(R(J)-R(I))*Q(J-1)/6.0D0
	  W(I,J)=(  ALPHA*( 1.0D0+(R(J)-R(I))/6.0D0*(A(I)-B(J)) )
	1            -BETA/6.0D0* (R(I)-R(K))*A(I)  )*Q(J)
	  W(I,K)=(  BETA*( 1.0D0+(R(I)-R(K))/6.0D0*(B(K)-C(I)) )
	1            +ALPHA/6.0D0*(R(J)-R(I))*C(I)  )*Q(K)
	  W(I,I)=(  ALPHA*( 1.0D0+(R(J)-R(I))/6.0D0*(B(I)-C(J)) )
	1       +BETA*( 1.0D0+(R(I)-R(K))/6.0D0*(A(K)-B(I)) )  )*Q(I)
	1       -R(I)*R(I)*(ZETA(I)+THETA(I)*RJ(I))/CHI(I)/Q(I)
	  IF(K .NE. ND)W(I,K+1)=BETA*C(K)*(R(I)-R(K))*Q(K+1)/6.0D0
	END DO
C
	UIJ=F(ND)*R(ND)*R(ND)*RJ(ND)/DTAU(ND-1)+TA(ND)*RJ(ND-1)
	ALPHA=0.5D0*(R(ND-1)-R(ND))*UIJ/DTAU(ND-1)
	W(ND,ND-2)=-ALPHA*Q(ND-2)*(R(ND-1)-R(ND))*A(ND-1)/6.0D0
	W(ND,ND-1)=ALPHA*Q(ND-1)*
	1            ( 1.0D0 + (R(ND-1)-R(ND))*(A(ND)-B(ND-1))/6.0D0 )
	W(ND,ND)=ALPHA*Q(ND)*
	1            ( 1.0D0 + (R(ND-1)-R(ND))*(B(ND)-C(ND-1))/6.0D0 )
C
	IF(DIFF)THEN
	  W(ND,ND)=W(ND,ND)-R(ND)*R(ND)*DBB/CHI(ND)/CHI(ND)/3.0D0
	END IF
C
	RETURN
	END
