C
C Routine to compute the arrays which characterize the Comoving-frame
C transfer equation. This is based on the second order differencing
C scheme of Hamman.
C
C Altered 3-May-1989 - Cleaned.
C Altered 22-May-1989 - Diffusion bug fixed. Now set RHS for ML>1
C                       equal to the flux. Bug was for ML=1.
C                       Schuster boundary condition fixed.
C
	SUBROUTINE HAMTUVGH(TA,TB,TC,UA,UB,UC,VB,VC,GB,H,XM,Q,QH,
	1                     SRCE,DTAU,DBC,IC,DIFF,LS,NC,NI,ML)
	IMPLICIT NONE
C
	INTEGER LS,NC,NI,ML
	REAL*8 TA(NI),TB(NI),TC(NI),VB(NI),VC(NI),GB(NI),H(NI)
	REAL*8 UA(NI),UB(NI),UC(NI),XM(NI)
	REAL*8 Q(NI),QH(NI),SRCE(NI),DTAU(NI),DBC,IC
	LOGICAL DIFF
C
C local Variables.
C
	INTEGER I
	REAL*8 T1
C
C 
C
	IF(ML .EQ. 1)THEN
C
	  DO I=1,NI
	    H(I)=0.0D0
	    UA(I)=0.0D0
	    UB(I)=0.0D0
	    UC(I)=0.0D0
	    VB(I)=0.0D0
	    VC(I)=0.0D0
	  END DO
C
	  GB(1)=1.0D0/DTAU(1)
	  TC(1)=GB(1)
	  TB(1)=-TC(1)-1.0D0
	  XM(1)=0.0D0
C
	  DO I=2,NI-1
	    T1=0.5D0*(DTAU(I)+DTAU(I-1))
	    GB(I)=1.0D0/DTAU(I)
	    TA(I)=TC(I-1)
	    TC(I)=GB(I)
	    TB(I)=-TA(I)-TC(I)-T1
	    XM(I)=-SRCE(I)*T1
	  END DO
C
	  IF(LS .GT. NC)THEN
	    T1=DTAU(NI-1)
	    TA(NI)=2.0D0/DTAU(NI-1)
	    TB(NI)=-TA(NI)-T1
	    XM(NI)=-SRCE(NI)*T1
	  ELSE IF(DIFF)THEN
	    TB(NI)=1.0D0/DTAU(NI-1)
	    TA(NI)=-TB(NI)
	    XM(NI)=DBC
	  ELSE
	    TA(NI)=-1.0D0/DTAU(NI-1)
	    TB(NI)=-TA(NI)+1.0D0
	    XM(NI)=IC
	  END IF
	  TC(NI)=0.0D0
	  GB(NI)=0.0D0
	  RETURN
	END IF
C
C 
C
C****************************************************************************
C
	DO I=1,NI-1
	  GB(I)=0.5D0/((0.5D0+QH(I))*DTAU(I))
	  H(I)=(QH(I)-0.5D0)/(QH(I)+0.5D0)
	END DO
	GB(NI)=0.0D0
	H(NI)=0.0D0
C
	TA(1)=0.0D0
	TC(1)=0.5D0/DTAU(1)
	TB(1)=-TC(1)-0.5D0-Q(1)
	UA(1)=0.0D0
	UC(1)=-TC(1)
	UB(1)=TC(1)+0.5-Q(1)
	VB(1)=0.0D0
	VC(1)=0.0D0
	XM(1)=0.0D0
C
C For clarity with differencing derivation, T1 could be defined as
C 0.5D0*(DTAU(I)+DTAU(I-1)). To speed execution, we have included the
C factor of 0.5 directly into the defintions. Saves several multiplies.
C
	DO I=2,NI-1
	  T1=DTAU(I)+DTAU(I-1)
	  TA(I)=GB(I-1)
	  TC(I)=GB(I)
	  TB(I)=-GB(I-1)-GB(I)-(0.5D0+Q(I))*T1
	  VB(I)=H(I-1)+1.0D0
	  VC(I)=-H(I)-1.0D0
	  UA(I)=-GB(I-1)
	  UC(I)=-GB(I)
	  UB(I)=GB(I)+GB(I-1)+(0.5D0-Q(I))*T1
	  XM(I)=-SRCE(I)*T1
	END DO
C
	IF(LS .GT. NC)THEN
	  T1=2.0D0*DTAU(NI-1)
	  TA(NI)=2.0D0*GB(NI-1)
	  TB(NI)=-TA(NI)-(0.5D0+Q(NI))*T1
	  VB(NI)=2.0D0*(H(NI-1)+1.0D0)
	  UA(NI)=-TA(NI)
	  UB(NI)=(0.5D0-Q(NI))*T1+TA(NI)
	  XM(NI)=-SRCE(NI)*T1
	ELSE IF(DIFF)THEN
	  T1=0.5D0/DTAU(NI-1)
	  TB(NI)=T1
	  TA(NI)=-T1
	  UA(NI)=T1
	  UB(NI)=-T1
	  VB(NI)=0.0D0
	  VC(NI)=0.0D0
	  XM(NI)=DBC
	ELSE
	  TA(NI)=-0.5D0/DTAU(NI-1)
	  TB(NI)=-TA(NI)+0.5D0+Q(NI)
	  UA(NI)=-TA(NI)
	  UB(NI)=TA(NI)-0.5D0+Q(NI)
	  VB(NI)=0.0D0
	  XM(NI)=IC
	END IF
	TC(NI)=0.0D0
	UC(NI)=0.0D0
	VC(NI)=0.0D0
C
	RETURN
	END
