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
C Created 04-May-1989 - Based on VARLAMKI
C Altered 22-May-1989 - Diffusion terms corrected. Schuster boundary corrected.
C
	SUBROUTINE HAM_VARLAM(KI,VSRCE,RKB,
	1   LINE_PRO,DTAU,R,Z,
	1   GB,Q,QH,UK,UKM1,VKM1,
	1   TCHI,SOURCE,ETAL,CHIL,CHI,
	1   IBOUND,TERF,THK_LINE,
	1   DBC,DIF,LS,NC,NI)
	IMPLICIT NONE
C
	LOGICAL THK_LINE,DIF
	INTEGER LS,NC,NI
	REAL*8 TERF,IBOUND,DBC
	REAL*8 KI(NI),VSRCE(NI),RKB(NI),LINE_PRO(NI)
	REAL*8 DTAU(NI),R(NI),Z(NI)
	REAL*8 GB(NI),Q(NI),QH(NI),UK(NI),UKM1(NI)
	REAL*8 VKM1(NI),TCHI(NI),SOURCE(NI),CHI(NI),CHIL(NI),ETAL(NI)
C
C Local variables
C
	INTEGER NV,I,J,K
	PARAMETER (NV=100)
	REAL*8 T1,dTAdCHI,dTBdCHI,dTCdCHI,dVBdCHI,dVCdCHI
	REAL*8 WA(NV),WB(NV),dGBdTAU(NV),dGBdCHI(NV),dHdCHI(NV)
C
	IF(NI .GT. NV)THEN
	  WRITE(6,*)'Error in HAM_VARLAM - NV too small'
	  STOP
	END IF
C
	CALL DP_ZERO(KI,NI)
	CALL DP_ZERO(VSRCE,NI)
	CALL DP_ZERO(RKB,NI)
C
C 
C
C*********************************************************************
C                  Opacity Section
C*********************************************************************
C
	CALL GET_DTAULAM(WA,WB,R,Z,NI)
C
C
C NB dGBdCHI is a partial derivative - it ignores the DTAU term.
C    All other d..dCHI variables are total derivatives.
C
	DO I=1,NI-1
	 dGBdTAU(I)=-GB(I)/DTAU(I)
	 T1=(0.5D0+QH(I))*(TCHI(I)+TCHI(I+1))
	 dGBdCHI(I)=GB(I)*QH(I)/T1
	 dHdCHI(I)= -QH(I)/(QH(I)+0.5D0)/T1
	END DO
C
C Note that the d=1 equation contains no emission term.
C
	
	T1=0.5D0/DTAU(1)/DTAU(1)
        KI(1)= Q(1)*(UKM1(1)-UK(1))/TCHI(1)
	1           + T1*( UK(2)+UKM1(2)-UK(1)-UKM1(1) ) * WB(1)
C
C NB dTdCHI=dGBdCHI(1).
C
	T1=dGBdCHI(1)+WB(1)*dGBdTAU(1)
	RKB(1)=T1*( UK(2)-UK(1)+UKM1(1)-UKM1(2) ) + dHdCHI(1)*VKM1(1)
C
	DO I=2,NI-1
	  J=I-1
	  K=I+1
	  dTAdCHI=dGBdCHI(J)+WA(I)*dGBdTAU(J)
	  dTCdCHI=dGBdCHI(I)+WB(I)*dGBdTAU(I)
	  dTBdCHI=Q(I)*(DTAU(J)+DTAU(I))/TCHI(I)
	1              - (0.5D0+Q(I))*(WA(I)+WB(I))
	1              -  dTAdCHI - dTBdCHI
	  dVBdCHI=dHdCHI(J)
	  dVCdCHI=-dHdCHI(I)
C
	  KI(I)= 2.0D0*UKM1(I)*Q(I)*(DTAU(I)+DTAU(J))/TCHI(I)
	1        - (SOURCE(I)+2.0D0*Q(I)*UKM1(I))*(WA(I)+WB(I))
	1        - dTAdCHI*(UK(J)+UKM1(J))
	1        - dTBdCHI*(UK(I)+UKM1(I))
	1        - dTCdCHI*(UK(K)+UKM1(K))
	1        + dVBdCHI*VKM1(J) + dVCdCHI*VKM1(I)
C
C NB TC(I)==GB(I) and VC(I)==(-1.0-H(I)) FOR I=2,...,NI-1
C
	  RKB(I)= dTCdCHI*(UK(K)-UK(I)+UKM1(K)-UKM1(I)) - dVCdCHI*VKM1(I)
C
	END DO
C
	J=NI-1
	IF(LS .GT. NC)THEN
	  dTAdCHI=2.0D0*( dGBdCHI(NI-1)+dGBdTAU(NI-1)*WA(NI) )
	  dTBdCHI=2.0D0*Q(NI)*( DTAU(NI-1)/TCHI(NI) - WA(NI) )-dTAdCHI
	  dVBdCHI=2.0D0*dHdCHI(NI-1)
	  KI(NI)= 4.0D0*UKM1(NI)*Q(NI)*DTAU(NI)/TCHI(NI)
	1        - (SOURCE(NI)+2.0D0*Q(NI)*UKM1(NI))*WA(NI)*2.0D0
	1        - dTAdCHI*(UK(J)+UKM1(J))
	1        - dTBdCHI*(UK(NI)+UKM1(NI))
	1        + dVBdCHI*VKM1(J)
	ELSE IF(DIF)THEN
	  dTAdCHI=0.5D0/DTAU(NI-1)/DTAU(NI-1)
	  dTBdCHI=-dTAdCHI
	  KI(NI)= -DBC/TCHI(NI)
	1            - dTAdCHI*(UK(J)+UKM1(J))
	1            - dTBdCHI*(UK(NI)+UKM1(NI))
	ELSE
	  dTAdCHI=0.5D0/DTAU(NI-1)/DTAU(NI-1)
	  dTBdCHI=-Q(NI)/TCHI(NI)-dTAdCHI
	  KI(NI)= -2.0D0*Q(NI)*UKM1(NI)/TCHI(NI)
	1            - dTAdCHI*(UK(J)+UKM1(J))
	1            - dTBdCHI*(UK(NI)+UKM1(NI))
	END IF
	RKB(NI)=0.0D0
C
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
C aprroximation. This depends both on TCHI and CHIL, hence we
C do it here. NB. The derivative of the ETAL(1)/CHIL(1) term
C is in the SOURCE section.
C
	IF(THK_LINE)THEN
	  T1=DEXP(TERF*CHIL(1))*CHI(1)/TCHI(1)*(ETAL(1)/CHIL(1)-IBOUND)
	  KI(1)=KI(1)+T1*( TERF-LINE_PRO(1)/TCHI(1) )
	END IF
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
	  VSRCE(I)=-(DTAU(I-1)+DTAU(I))
	END DO
C
C Both the diffusion and Schuster boundary conditions are independent
C of the SOURCE function. Only for rays not striking the core do we
C need to worry about the source function dependance.
C
	 IF(LS .GT. NC)THEN
	   VSRCE(NI)=-2.0D0*DTAU(NI-1)
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
	  VSRCE(1)=VSRCE(1)+( CHI(1)/TCHI(1)*DEXP(TERF*CHIL(1))-1.0D0 )
	END IF
C
	RETURN
	END
