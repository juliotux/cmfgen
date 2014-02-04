C
C Routine to compute the line opacity variation vector, and the
C line source variation vector in the spirit of Schonberg and
C Hempe (1986, A&A, 163, 151).
C
C RKB multiply's the %KI in the V equation.
C LINE_PRO is the line profile.
C It is assumed that :-
C				KI( )=dCHIL(I)
C				VSRCE( )=dETAL(I)
C
	SUBROUTINE VARLAMKI(KI,VSRCE,RKB,
	1   LINE_PRO,DTAU,Z,Q,QH,UK,UKM1,VKM1,
	1   TCHI,SOURCE,ETAL,CHIL,CHI,TERF,THK_LINE,LS,NC,NI)
	IMPLICIT NONE
C
C Altered 28-May-1996 - Call to DP_ZERO removed.
C                       Generical calls for EXP
C                       ONE inserted.
C Created 02-Feb-1989 - Based on KIVARNM. In KIVARNM a vector RKC
C                       is also retrned. Since this is identically equal
C                       to RKB, we hav negelected it.
C
	LOGICAL THK_LINE
	INTEGER LS,NC,NI
	REAL*8 TERF
	REAL*8 KI(NI),VSRCE(NI),RKB(NI),LINE_PRO(NI)
	REAL*8 DTAU(NI),Q(NI),QH(NI),UK(NI),UKM1(NI),Z(NI)
	REAL*8 VKM1(NI),TCHI(NI),SOURCE(NI),CHI(NI),CHIL(NI),ETAL(NI)
C
C Local variables
C
	REAL*8, PARAMETER :: ONE=1.0D0
	INTEGER I,J,K
	REAL*8 T1,T2,T3,T4,T5,T6,T7,KIP,KIM
C
C 
C
C*********************************************************************
C                  Opacity Section
C*********************************************************************
C
	KI(:)=0.0D0
C
C Note that the D=1 equation contains no emission term.
C
	T4=1.0D0/(DTAU(1)*(1.0+QH(1)))
	T5=T4/DTAU(1)*0.5D0
	T6=QH(1)/((1.0+QH(1))*(TCHI(1)+TCHI(2)))
C
	KIP=(Z(1)-Z(2))*(UK(2)-UK(1))*0.5D0/(DTAU(1)*DTAU(1))
	KI(1)=KIP+(UKM1(1)-UK(1))*Q(1)/TCHI(1)
	RKB(1)=-KIP/(1.0+QH(1)) - T6*T4*( UK(1)-UK(2)
	1              +VKM1(1)*DTAU(1) )
C
	DO I=2,NI-1
	  J=I-1
	  K=I+1
	  T1=T4
	  T2=T5
	  T3=T6
	  T4=1.0D0/(DTAU(I)*(QH(I)+1.0D0))
	  T5=T4/DTAU(I)*0.5D0
	  T6=QH(I)/((1+QH(I))*(TCHI(I)+TCHI(I+1)))
	  T7=0.25D0*(UK(I)*(1.0D0+Q(I))-SOURCE(I)-Q(I)*UKM1(I))
C
	  KIM=(Z(J)-Z(I))*( (UK(J)-UK(I))*T2+T7 )
	1       -T3*((UK(J)-UK(I))*T1+VKM1(J)/( 1.0D0+QH(J)) )
C
	  KIP=(Z(I)-Z(K))*( (UK(K)-UK(I))*T5+T7 )
	1      -T6*((UK(K)-UK(I))*T4-VKM1(I)/(1.0+QH(I)))
C
C The SOURCE term has been drooped from:
C	1       /TCHI(I)*( SOURCE(I)+(UKM1(I)-UK(I))*Q(I) )*0.5D0
C as we are treating (CHIL,SOURCEL) as the independent variables.
C KIVARNM treats (CHIL,ETAL) as independent variables.
C
	  KI(I)=KIM+KIP+( DTAU(I)+DTAU(J) )
	1       *(UKM1(I)-UK(I))*Q(I)*0.5D0/TCHI(I)
C
	  RKB(I)=(Z(I)-Z(K))*(UK(I)-UK(K))*T5-T6*T4
	1      *( (UK(I)-UK(K))+VKM1(I)*DTAU(I) )
C
	END DO
C
	RKB(NI)=0.0
	IF(LS .LE. NC)THEN
	  KIM=0.5D0*(Z(NI-1)-Z(NI))*(UK(NI)-UK(NI-1))/
	1             (DTAU(NI-1)*DTAU(NI-1))
	  KI(NI)=KIM+Q(NI)/TCHI(NI)*(UK(NI)-UKM1(NI))
	ELSE
	  KIM=(UK(NI)-UK(NI-1))*T4*
	1       ( (Z(NI-1)-Z(NI))*0.5D0/DTAU(NI-1)-T6 )+
	1        VKM1(NI-1)*T6/(1.0+QH(NI-1))+
	1        0.25D0*(Z(NI-1)-Z(NI))*( SOURCE(NI)-(1.0D0+Q(NI))*
	1        UK(NI)+Q(NI)*UKM1(NI) )
C
C The SOURCE term has been drooped from:
C	  KI(NI)=KIM+0.5D0*( (UK(NI)-UKM1(NI))*Q(NI)-SOURCE(NI) )*
C as we are treating (CHIL,SOURCEL) as the independent variables.
C KIVARNM treats (CHIL,ETAL) as independent variables.
C
	  KI(NI)=KIM+0.5D0*(UK(NI)-UKM1(NI))*Q(NI)*
	1           DTAU(NI-1)/TCHI(I)
	END IF
C
C
C Multiply line opacity  variation by line profile. Unlike KIVARNM
C we also multiply RKB since it does not need to be carried for
C the continuum as well.
C
	DO I=1,NI
	 KI(I)=KI(I)*LINE_PRO(I)
	 RKB(I)=RKB(I)*LINE_PRO(I)
	END DO
C
C Boundary condition at outer boundary using the SOLOBOV
C aprroximation.
C
	T1=EXP(TERF*CHIL(1))
	T2=(CHI(1)/TCHI(1)*T1-ONE)
	KI(1)=KI(1)+T2/CHIL(1)
	T1=T1*(TERF*CHI(1)-CHI(1)*LINE_PRO(1)/TCHI(1))/TCHI(1)
	KI(1)=KI(1)+(T1-T2/CHIL(1))*ETAL(1)/CHIL(1)
C
C 
C
C This section of the routine computes the approximate diagonal
C operator. Alternatively, it can be considred to give the variation
C of U with S assuming dUdS is non zero only when U and S are at the
C same grid point.
C
	CALL DP_ZERO(VSRCE,NI)
	DO I=2,NI-1
	  VSRCE(I)=-0.5D0*(DTAU(I-1)+DTAU(I))
	END DO
C
C Both the diffusion and Schuster boundary conditions are independent
C of the SOURCE function. Only for rays not striking the core do we
C need to worry about the source function dependance.
C
	 IF(LS .GT. NC)THEN
	   VSRCE(NI)=0.5D0*DTAU(NI-1)
	 END IF
C
C Correct SOURCE function variation for continuum terms in
C soucre function at this frequency. Note that S= ( BETA*Sc+PRO*SL )
C /(BETA + PRO) . Thus dS/dSL= PRO/( BETA+PRO) = PRO*CHIL/TCHI
C
	DO I=1,NI
	  VSRCE(I)=VSRCE(I)*LINE_PRO(I)*CHIL(I)/TCHI(I)
	END DO
C
C Boundary condition at outer boundary using the SOLOBOV
C aprroximation. Note that the outer boundary condition is
C directly proportional to the line SOURCE function.
C
	IF(THK_LINE)THEN
	  VSRCE(1)=VSRCE(1)+( CHI(1)/TCHI(1)*EXP(TERF*CHIL(1))-ONE )
	END IF
C
	RETURN
	END
