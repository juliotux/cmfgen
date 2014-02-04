C
C Routine to compute the opacity & emissivity variation matrices for
C the case with lines. Also computes the upper bidiagonal matrix
C (RKB,RKC) which multiply's %KI in the V equation. X is the line
C profile. It is assumed that :-
C				VK( , ,1)=%CHIL
C				VK( , ,2)=%ETAL
C				VK( , ,3)=%CHI (if needed)
C				VK( , ,4)=%ETA (if needed)
C
	SUBROUTINE KIVARNM(VK,RKB,RKC,X,DTAU,Z,Q,QH,UK,UKM1,VKM1,
	1             TCHI,SOURCE,ETAL,CHIL,CHI,TERF,LS,NC,NI,NM)
	IMPLICIT NONE
C
C Altered 24-May-1996 - IMPLCIT NONE installed.
C                         Initialization of VK now done using F90.
C Altered 7-SEP-1982 (variable NM)
C
	INTEGER LS,NC,NI,NM
C
	REAL*8 VK(NI,3,NM),RKB(NI),RKC(NI),X(NI)
	REAL*8 DTAU(NI),Q(NI),QH(NI),UK(NI),UKM1(NI),Z(NI)
	REAL*8 VKM1(NI),TCHI(NI),SOURCE(NI)
	REAL*8 CHI(NI),CHIL(NI),ETAL(NI)
	REAL*8 TERF
C
C Local variables.
C
	REAL*8 T1,T2,T3,T4,T5,T6,T7
	INTEGER I,J,K
C
	VK(:,:,:)=0.0D0
C
C Note that the d=1 equation contains no emission term.
C
	T4=1.0/(DTAU(1)*(1.0+QH(1)))
	T5=T4/DTAU(1)*0.5
	T6=QH(1)/((1.0+QH(1))*(TCHI(1)+TCHI(2)))
	VK(1,3,1)=(Z(1)-Z(2))*(UK(2)-UK(1))*0.5/(DTAU(1)*DTAU(1))
	VK(1,2,1)=VK(1,3,1)+(UKM1(1)-UK(1))*Q(1)/TCHI(1)
C
	RKB(1)=-VK(1,3,1)/(1.0+QH(1))-T6*T4*(UK(1)-UK(2)
	1+VKM1(1)*DTAU(1))
	RKC(1)=RKB(1)
C
C
	DO I=2,NI-1
	J=I-1
	K=I+1
	T1=T4
	T2=T5
	T3=T6
	T4=1.0/(DTAU(I)*(QH(I)+1.0))
	T5=T4/DTAU(I)*0.5
	T6=QH(I)/((1+QH(I))*(TCHI(I)+TCHI(I+1)))
	T7=0.25*(UK(I)*(1.0+Q(I))-SOURCE(I)-Q(I)*UKM1(I))
C
C OPACITY TERM (ie. ,1) )
C
	VK(I,1,1)=(Z(J)-Z(I))*((UK(J)-UK(I))*T2+T7)
	1-T3*((UK(J)-UK(I))*T1+VKM1(J)/(1.0+QH(J)))
C
	VK(I,3,1)=(Z(I)-Z(K))*((UK(K)-UK(I))*T5+T7)
	1-T6*((UK(K)-UK(I))*T4-VKM1(I)/(1.0+QH(I)))
C
	VK(I,2,1)=VK(I,1,1)+VK(I,3,1)+(DTAU(I)+DTAU(J))
	1/TCHI(I)*(SOURCE(I)+(UKM1(I)-UK(I))*Q(I))*0.5
C
C
	RKB(I)=(Z(I)-Z(K))*(UK(I)-UK(K))*T5-T6*T4
	1*((UK(I)-UK(K))+VKM1(I)*DTAU(I))
	RKC(I)=RKB(I)
C
C EMISSION TERM (ie. ,2)  )
C
	VK(I,2,2)=-0.5*(DTAU(J)+DTAU(I))/TCHI(I)
C
	END DO
C
C
	RKB(NI)=0.0
	RKC(NI)=0.0
			IF(LS .LE. NC)THEN
	VK(NI,1,1)=0.5*(Z(NI-1)-Z(NI))*(UK(NI)-UK(NI-1))/(DTAU
	1(NI-1)*DTAU(NI-1))
	VK(NI,2,1)=VK(NI,1,1)+Q(NI)/TCHI(NI)*(UK(NI)-UKM1(NI))
			ELSE
	VK(NI,1,1)=(UK(NI)-UK(NI-1))*T4*((Z(NI-1)-Z(NI))*0.5
	1/DTAU(NI-1)-T6)+VKM1(NI-1)*T6/(1.0+QH(NI-1))+0.25*(Z(NI-1)-
	1Z(NI))*(SOURCE(NI)-(1.0+Q(NI))*UK(NI)+Q(NI)*UKM1(NI))
C
	VK(NI,2,1)=VK(NI,1,1)+0.5*((UK(NI)-UKM1(NI))*Q(NI)
	1-SOURCE(NI))*DTAU(NI-1)/TCHI(I)
C
	VK(NI,2,2)=0.5*DTAU(NI-1)/TCHI(NI)
C
			END IF
C
	 IF(NM .EQ. 4)THEN
	  DO J=1,3
	   DO I=1,NI
	    VK(I,J,3)=VK(I,J,1)
	    VK(I,J,4)=VK(I,J,2)
	  END DO			!I
	 END DO			!J
	 END IF
C
C Multiply line opacity and emissivity variation by line profile.
C
	VK(1,2,1)=VK(1,2,1)*X(1)
	VK(1,3,1)=VK(1,3,1)*X(2)
	VK(1,2,2)=VK(1,2,2)*X(1)
C
	DO I=2,NI-1
	 VK(I,2,1)=VK(I,2,1)*X(I)
	 VK(I,3,1)=VK(I,3,1)*X(I+1)
	 VK(I,1,1)=VK(I,1,1)*X(I-1)
	 VK(I,2,2)=VK(I,2,2)*X(I)
	END DO
C
	VK(NI,2,1)=VK(NI,2,1)*X(NI)
	VK(NI,1,1)=VK(NI,1,1)*X(NI-1)
	VK(NI,2,2)=VK(NI,2,2)*X(NI)
C
C
C Boundary condition at outer boundary using the SOLOBOV
C aprroximation.
C
	T1=EXP(TERF*CHIL(1))
	T2=(CHI(1)/TCHI(1)*T1-1.0D0)
	VK(1,2,2)=VK(1,2,2)+T2/CHIL(1)
	  IF(NM .EQ. 4)THEN
	    VK(1,2,3)=VK(1,2,3)+ETAL(1)*T1*X(1)
	1   /TCHI(1)**2.0
	  END IF
	T1=T1*(TERF*CHI(1)-CHI(1)*X(1)/TCHI(1))/TCHI(1)
	VK(1,2,1)=VK(1,2,1)+(T1-T2/CHIL(1))*ETAL(1)/CHIL(1)
C
	RETURN
	END
