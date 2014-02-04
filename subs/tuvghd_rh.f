C
C Routine to compute the arrays which characterize the
C Comoving-frame transfer equation.
C
C NB: Equations are
C
C       TA(i).Y(i-1) - [ TA(i)+TB(i)+TC(i) ].Y(i) + TC(i).Y(i) = RHS(I) + ...
C
	SUBROUTINE TUVGHD_RH(TA,TB,TC,U,VB,VC,GB,H,RHS,
	1                   Q,QH,DTAU,SOURCE,DIFF,DBC,IC,LS,NC,NI)
	IMPLICIT NONE
C
C Created 14-Feb-1994 - Based on TUVGHD and XVECD
C
	INTEGER LS,NC,NI
	REAL*8 TA(NI),TB(NI),TC(NI),VB(NI),VC(NI),GB(NI),RHS(NI)
	REAL*8 H(NI),U(NI),Q(NI),QH(NI),DTAU(NI),SOURCE(NI)
	REAL*8 DBC,IC
	LOGICAL DIFF
C
C Local variables.
C
	INTEGER I
C
	TA(1)=0.0D0
	VB(1)=0.0D0
	VC(1)=0.0D0
	GB(1)=-1.0D0/((1.0D0+QH(1))*DTAU(1))
	H(1)=QH(1)/(1.0D0+QH(1))
	TC(1)=1.0D0/DTAU(1)
	TB(1)=1.0D0+Q(1)
	U(1)=-Q(1)
	RHS(1)=0.0D0
C
	DO I=2,NI-1
	  H(I)=QH(I)/(1.0D0+QH(I))
	  GB(I)=-1.0D0/((1.0D0+QH(I))*DTAU(I))
	END DO
C
	DO I=2,NI-1
	  TA(I)=-GB(I-1)
	  TC(I)=-GB(I)
	  VB(I)=-H(I-1)
	  VC(I)=H(I)
	  U(I)=-0.5D0*Q(I)*(DTAU(I)+DTAU(I-1))
	  TB(I)=0.5D0*(1.0D0+Q(I))*(DTAU(I)+DTAU(I-1))
	  RHS(I)=-0.5D0*SOURCE(I)*(DTAU(I-1)+DTAU(I))
	END DO
C
	IF(LS .GT. NC)THEN
	  TA(NI)=GB(NI-1)
	  TB(NI)=-0.5D0*(1.0D0+Q(NI))*DTAU(NI-1)
	  U(NI)=0.5D0*Q(NI)*DTAU(NI-1)
	  VB(NI)=H(NI-1)
	  RHS(NI)=0.5D0*DTAU(NI-1)*SOURCE(NI)
	ELSE IF(DIFF)THEN
	  TB(NI)=0.0D0
	  TA(NI)=-1.0D0/DTAU(NI-1)
	  U(NI)=0.0D0
	  VB(NI)=0.0D0
	  RHS(NI)=DBC
	ELSE
	  TA(NI)=-1.0D0/DTAU(NI-1)
	  TB(NI)=-1.0D0-Q(NI)
	  U(NI)=Q(NI)
	  VB(NI)=0.0D0
	  RHS(NI)=IC
	END IF
	TC(NI)=0.0D0
	VC(NI)=0.0D0
	GB(NI)=0.0D0
C
	RETURN
	END
