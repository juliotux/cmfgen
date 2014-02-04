C
C Subroutine to compute T matrix for continuum solution of
C the transfer equation in an extended atmosphere using the FEAUTRIER
C technique. For the continuum the effect of the velocity gradient can be
C ignored. Uses the diffusion approximation for the lower boundary condition
C or a Schuster bounary condition.
C
C Created 17-Feb-1986
C Altered 31-Oct-1986 - Schuster B.C. installed. (Calling changed)
C
	SUBROUTINE TFEAU(TA,TB,TC,R,Q,F,THETA,DTAU,HBC,INBC,DIFF,ND)
	IMPLICIT NONE
	INTEGER ND,I
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 R(ND),Q(ND),F(ND),THETA(ND),DTAU(ND),HBC,T1,INBC
	LOGICAL DIFF
C
	TA(1)=0.0D0
	TC(1)=F(2)*Q(2)*R(2)*R(2)/DTAU(1)
	TB(1)=-R(1)*R(1)*(F(1)*Q(1)/DTAU(1)+HBC)
C
	DO I=2,ND-1
	  T1=0.5D0*(DTAU(I-1)+DTAU(I))
	  TA(I)=-F(I-1)*Q(I-1)*R(I-1)*R(I-1)/DTAU(I-1)/T1
	  TC(I)=-F(I+1)*Q(I+1)*R(I+1)*R(I+1)/DTAU(I)/T1
	  TB(I)=R(I)*R(I)*((1.0D0-THETA(I))/Q(I)
	1        +F(I)*Q(I)*(1.0D0/DTAU(I)+1.0D0/DTAU(I-1))/T1)
	END DO
C
C Note that Q(ND)=1.0d0 by definition.
C
	IF(DIFF)THEN
	  TA(ND)=-F(ND-1)*Q(ND-1)*R(ND-1)*R(ND-1)/DTAU(ND-1)
	  TB(ND)=F(ND)*R(ND)*R(ND)/DTAU(ND-1)
	  TC(ND)=0.0D0
	ELSE
	  TA(ND)=-F(ND-1)*Q(ND-1)*R(ND-1)*R(ND-1)/DTAU(ND-1)
	  TB(ND)=R(ND)*R(ND)*(F(ND)/DTAU(ND-1)+INBC)
	  TC(ND)=0.0D0
	END IF
C
	RETURN
	END
