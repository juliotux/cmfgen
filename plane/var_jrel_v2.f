!
! This subroutine computes the variation of J as a function of the emissivity
! ETA and the opacity CHI.
!
! TX has dimensions (ND,ND,NM) and TVX has dimension (ND-1,ND,NM).
!
! NM reflects the total number of matrices required to compute the variation
!    of J. All NM matrices are modified by a call to THOMAS in UP_TX_TVX. The
!    first 2 should represent dJ/dCHI and the second dJ/dETA --- the others
!    can be in any order.
!
! Typically they will be CHI_C, and ETA_C, CHIL_, AND ETAL of othe lines.
!
! NB TX(i,m, ) =dJ(i)/d(CHI(m),ETA(m),...)
! NB TVX(i,m, ) =dRSQH(i)/d(CHI(m),ETA(m),...)
!
! NM_KI reflects the 3rd dimension of KI. For this routine only the first 2 
!   are important and are used to compute the variation of J with respect to 
!   the CURRENT opacity and the CURRENT emissivity. 
!
! TX_DIF_d_T and TX_DIFF_d_dTdR are used to describe the variations in J
! caused by the DIFFUSION approximation at the inner boundary.
!
! The particular choice of the outer boundary condition adopted is irrelevant
! for this routine. Such information is incorporated by the outer boundary
! Eddington factors HBC, and NBC.
!
! Note TA, TBC, TC, HU etc have been defined so that the program computes
! JNU and r^2 HNU.
!
	SUBROUTINE VAR_JREL_V2(ETA,CHI,ESEC,THETA,V,SIGMA,R,
	1                  TX,TVX,TX_DIF_d_T,TX_DIF_d_dTdR,
	1                  TVX_DIF_d_T,TVX_DIF_d_dTdR,KI,WORKMAT,RHS_dHdCHI,
	1                  INIT,FREQ,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1	           INCL_ADVEC_TERMS,INCL_REL_TERMS,
	1                  DO_THIS_TX_MATRIX,METHOD,ND,NM,NM_KI)
	USE MOD_VAR_JREL_V2
	USE MOD_RAY_MOM_STORE
	USE MOD_TRAP_DERIVATIVES
	IMPLICIT NONE
!
! Altered 27-Dec-2008: 
! Created:  
!
	INTEGER ND
	INTEGER NM
	INTEGER NM_KI
!
	REAL*8 ETA(ND)
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 THETA(ND)
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
!
! Variation arrays and vectors.
!
	REAL*8 TX(ND,ND,NM)
	REAL*8 TVX(ND-1,ND,NM)
	REAL*8 KI(ND,ND,NM_KI)
	REAL*8 WORKMAT(ND,ND)
	REAL*8 RHS_dHdCHI(ND-1,ND)
	REAL*8 TX_DIF_d_T(ND),TX_DIF_d_dTdR(ND)
	REAL*8 TVX_DIF_d_T(ND),TVX_DIF_d_dTdR(ND)
!
	LOGICAL DO_THIS_TX_MATRIX(NM)
!
	REAL*8 dLOG_NU,dTdR,DBB,dDBBdT,IC
	CHARACTER*6 METHOD
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
	LOGICAL DIF,INIT
	LOGICAL INCL_ADVEC_TERMS
	LOGICAL INCL_REL_TERMS
	LOGICAL J_AT_INB_EQ_B
!
! Local variables.
!
	REAL*8 T1
	REAL*8 FREQ
!
	INTEGER LUER,ERROR_LU
        EXTERNAL ERROR_LU
	INTEGER I,J
	INTEGER IFAIL
! 
!
	IF(INIT)CALL ALLOC_MOD_VAR_JREL_V2(ND)
!
! Zero relevant vectors and matrices.
!
	JNU=0.0D0			!1:ND
	GAM_RSQHNU=0.0D0		!1:ND
!
	IF(INIT)THEN
	  TX(:,:,:)=0.0D0		!ND*ND*NM
	  TVX(:,:,:)=0.0D0		!(ND-1)*ND*NM )
	  JNU_PREV=0.0D0
	  GAM_RSQHNU_PREV=0.0D0
	  W=0.0D0
	  WPREV=0.0D0
	  PSI=0.0D0
	  PSIPREV=0.0D0
	  TX_DIF_d_T=0.0D0
	  TX_DIF_d_dTdR=0.0D0
	  EPS_A=0.0D0
	  EPS_B=0.0D0
	  EPS_PREV_A=0.0D0
	  EPS_PREV_B=0.0D0
	  DELTAH(:)=0.0D0
	  DELTA(:)=0.0D0
	  PSIPREV(:)=0.0D0
!
	  H_ON_J_PREV(:)=0.0D0
	  K_ON_J_PREV(:)=0.0D0
	  N_ON_J_PREV(:)=0.0D0
	  NMID_ON_HMID_PREV(:)=0.0D0
	  NMID_ON_J_PREV(:)=0.0D0
	  KMID_ON_J_PREV(:)=0.0D0
	END IF
!
	IF(INIT)THEN
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	  
	  BETA_FREQ(1:ND)=V(1:ND)*3.33564D-06       !/2.99794D+05
	  IF(INCL_REL_TERMS)THEN
	    BETA(1:ND)=BETA_FREQ(1:ND)
	  ELSE
	    BETA(1:ND)=0.0D0
	  END IF
!
	  IF(INCL_REL_TERMS)THEN
	    GAM_REL_SQ(1:ND)=1.0D0/(1.0D0-BETA(1:ND)*BETA(1:ND))
	    GAM_REL(1:ND)=SQRT(GAM_REL_SQ(1:ND))
	  ELSE
	    GAM_REL_SQ(1:ND)=1.0D0
	    GAM_REL(1:ND)=1.0D0
	  END IF
	  GAM_RSQ(1:ND)=GAM_REL(1:ND)*R(1:ND)*R(1:ND)
	  CON_DELTA(1:ND)=BETA_FREQ(1:ND)/R(1:ND)
	  CON_dKdNU(1:ND)=GAM_REL_SQ(1:ND)*(SIGMA(1:ND)+1.0D0)-1.0D0
	  CON_dHdNU(1:ND)=BETA(1:ND)*GAM_REL_SQ(1:ND)*(SIGMA(1:ND)+1.0D0)
!
	  DO I=1,ND-1
	    T1=0.5D0*(BETA(I)+BETA(I+1))
	    AV_SIGMA(I)=0.5D0*(SIGMA(I)+SIGMA(I+1))
	    CON_DELTAH(I)=2.0D0*(BETA_FREQ(I)+BETA_FREQ(I+1))/(R(I)+R(I+1))
	    CON_dKdNUH(I)=T1*(AV_SIGMA(I)+1.0D0)/(1.0D0-T1*T1)
	    CON_dNdNUH(I)=(AV_SIGMA(I)+1.0D0)/(1.0D0-T1*T1) -1.0D0
	  END DO
!
! These are set to zero to insure all velocity terms are neglected.
!
	  WRITE(171,*)' FREQ',FREQ
	  WRITE(171,*)'  DIF',DIF
	  WRITE(171,*)'  DBB',DBB
	  WRITE(171,*)'  ETA',ETA(1),ETA(ND)
	  WRITE(171,*)'  CHI',CHI(1),CHI(ND)
	  WRITE(171,*)' ESEC',ESEC(1),ESEC(ND)
	  WRITE(171,*)'THETA',THETA(1),THETA(ND)
	  WRITE(171,*)' SRCE',ETA(1)/(CHI(1)-ESEC(1)),ETA(ND)/(CHI(ND)-ESEC(ND))
	  WRITE(171,*)'    R',R(1),R(ND)
	  WRITE(171,*)'    V',V(1),V(ND)
	  WRITE(171,*)'SIGMA',SIGMA(1),SIGMA(ND)
	  WRITE(171,*)' K_ON_J',K_ON_J(1),K_ON_J(ND)
	  WRITE(171,*)'  H/J',H_ON_J(1),H_ON_J(ND)
	  WRITE(171,*)'  N/J',N_ON_J(1),N_ON_J(ND)
	  WRITE(171,*)'RN/RJ',NMID_ON_J(1),NMID_ON_J(ND)
	  CLOSE(UNIT=171)
!
	  H_ON_J(1:ND)=0.0D0
	  N_ON_J(1:ND)=0.0D0
	  KMID_ON_J(1:ND)=0.0D0
          NMID_ON_J(1:ND)=0.0D0
!
	END IF
!
! Zero relevant vectors and matrices.
!
	JNU(:)=0.0D0
	GAM_RSQHNU(:)=0.0D0
!
	IF(INCL_ADVEC_TERMS)THEN
	  VdHdR_TERM(1:ND)=H_ON_J(1:ND)*BETA(1:ND)
	ELSE
	  VdHdR_TERM(1:ND)=0.0D0
	END IF
!
!*****************************************************************************
!
! CHI_H refers to the modified CHI term that multiplies H in the 1ST
! moment equation. We evaluate it at the grid points, and then use the
! averaging procedure used in the non-relativistic case.
!
! CHI_J refers the the modified CHI term that multiplies J in the 0th
! moment equation.
!
	IF(INCL_REL_TERMS)THEN
	  DO I=1,ND
	    CHI_H(I)=CHI(I)/GAM_REL(I)+CON_DELTA(I)*(
	1       1.0D0+2.0D0*GAM_REL_SQ(I)*(SIGMA(I)+1.0D0) )
	    P_H(I)=1.0D0
	  END DO
	  IF(INCL_ADVEC_TERMS)THEN
	    DO I=1,ND
	      CHI_J(I)=CHI(I)/GAM_REL(I)+CON_DELTA(I)*(1.0D0+SIGMA(I))
	    END DO
	  ELSE
	    DO I=1,ND
	      CHI_J(I)=CHI(I)/GAM_REL(I)+CON_DELTA(I)*(
	1                   2.0D0+GAM_REL_SQ(I)*(1.0D0+SIGMA(I)) )
	    END DO
	  END IF
	ELSE
	  DO I=1,ND
	    CHI_H(I)=CHI(I)
	    P_H(I)=1.0D0
	    CHI_J(I)=CHI(I)
	  END DO
	END IF
!
	SOURCE(1:ND)=ETA(1:ND)/CHI_J(1:ND)
	COH_VEC(1:ND)=THETA(1:ND)*CHI(1:ND)/CHI_J(1:ND)/GAM_REL(1:ND)
!
! Compute the Q factors from F. The Q is used to make the J and K terms
! appearing in the first moment equation a perfect differential.
! Note that we integrate from the core outwards, and normalize Q to the
! value of unity at the core.
!
	DO I=1,ND
	  TA(ND-I+1)=( 3.0D0*K_ON_J(I)-1.0D0+
	1           BETA(I)*N_ON_J(I)-(SIGMA(I)+1.0D0)*VdHdR_TERM(I)+
	1      GAM_REL_SQ(I)*BETA(I)*(SIGMA(I)+1.0D0)*
	1       (BETA(I)*(1.0D0-K_ON_J(I)-VdHdR_TERM(I))-N_ON_J(I))
	1            )/(K_ON_J(I)+VdHdR_TERM(I))/R(I)
	  TB(I)=R(ND-I+1)
	END DO
	CALL INTEGRATE(TB,TA,Q,IFAIL,ND)
!
	DO I=1,ND-1		!Q(ND) undefined after exiting INTEGRATE
	  TB(I)=EXP(Q(ND-I))
	END DO
!
! Scale by r^2
!
	DO I=1,ND-1
	  Q(I)=TB(I)/(R(I)/R(ND))**2
	END DO
	Q(ND)=1.0D0
!
! Compute optical depth scales. We also compute the dTAUdCHI matrices
! for use in EDD_JREL_VAR.
!
	TA(1:ND)=CHI_H(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL d_DERIVCHI_dCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU_H,TA,R,R,TB,ND)
	CALL dSPHEREdCHI(dTAUdCHI_H,DTAU_H,R,Q,ND)
!
	TA(1:ND)=CHI_J(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL d_DERIVCHI_dCHI(TB,TA,R,ND,METHOD)
        CALL NORDTAU(DTAU_J,TA,R,R,TB,ND)
	CALL dSPHEREdCHI(dTAUdCHI_J,DTAU_J,R,Q,ND)
!
	IF(.NOT. INIT)THEN
!
! We are integrating from blue to red. dLOG_NU is define as vd / dv which is 
! the same as d / d ln v.
!
! EPS is used if we define N in terms of J rather than H, This is sometimes
! useful as H can approach zero, and hence N/H is undefined.
!
	  DO I=1,ND-1
	    DELTAH(I)=CON_DELTAH(I)/dLOG_NU/(CHI_H(I)+CHI_H(I+1))
	    W(I)=DELTAH(I)*(1.0D0+CON_dNdNUH(I)*NMID_ON_HMID(I))
	    WPREV(I)=DELTAH(I)*(1.0D0+CON_dNdNUH(I)*NMID_ON_HMID_PREV(I))
	    EPS_A(I)=DELTAH(I)*(CON_dNdNUH(I)*NMID_ON_J(I)+
	1            CON_dKdNUH(I)*KMID_ON_J(I))/(P_H(I)+W(I))
	    EPS_B(I)=EPS_A(I)*GAM_RSQ(I+1)
	    EPS_A(I)=EPS_A(I)*GAM_RSQ(I)
	    EPS_PREV_A(I)=DELTAH(I)*(CON_dNdNUH(I)*NMID_ON_J_PREV(I)+
	1            CON_dKdNUH(I)*KMID_ON_J_PREV(I))/(P_H(I)+W(I))
	    EPS_PREV_B(I)=EPS_PREV_A(I)*GAM_RSQ(I+1)
	    EPS_PREV_A(I)=EPS_PREV_A(I)*GAM_RSQ(I)
	  END DO
!
	  DO I=2,ND
	    DELTA(I)=CON_DELTA(I)/CHI_J(I)/dLOG_NU
	  END DO
	  DELTA(1)=CON_DELTA(1)/CHI_H(1)/dLOG_NU
	END IF
!
	DO I=2,ND-1
	  GAM_RSQ_DTAUONQ(I)=0.5D0*GAM_RSQ(I)*(DTAU_J(I)+DTAU_J(I-1))/Q(I)
	  PSI(I)=GAM_RSQ_DTAUONQ(I)*DELTA(I)*
	1            (1.0D0+CON_dKdNU(I)*K_ON_J(I)+CON_dHdNU(I)*H_ON_J(I) )
	  PSIPREV(I)=GAM_RSQ_DTAUONQ(I)*DELTA(I)*
	1            (1.0D0+CON_dKdNU(I)*K_ON_J_PREV(I)+
	1                       CON_dHdNU(I)*H_ON_J_PREV(I) )
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=GAM_RSQ(I+1)*(K_ON_J(I+1)+VdHdR_TERM(I+1))*Q(I+1)/(P_H(I)+W(I))/DTAU_H(I)
	  HL(I)=GAM_RSQ(I)*(K_ON_J(I)+VdHdR_TERM(I))*Q(I)/(P_H(I)+W(I))/DTAU_H(I)
	  HS(I)=WPREV(I)/(P_H(I)+W(I))
	END DO
!
! 
!
! We know dlnJdlnR, so we don't need to iterate.
!
	IF(INCL_ADVEC_TERMS)THEN
!	  VdJdR_TERM(1:ND)=GAM_RSQ(1:ND)*CON_DELTA(1:ND)*dlnGRSQJdlnR(1:ND)/CHI_J(1:ND)
	  VdJdR_TERM(1:ND)=CON_DELTA(1:ND)*dlnGRSQJdlnR(1:ND)/CHI_J(1:ND)
	  P_J(1:ND)=1.0D0+VdJdR_TERM(1:ND)
	ELSE
	  VdJdR_TERM(1:ND)=0.0D0
	  P_J(1:ND)=1.0D0
	END IF
!
! Compute the TRIDIAGONAL operators, and the RHS source vector. These
! vectors either depend  on dlnJdlnR etc, or are corrupted in the
! solution.
!
	DO I=2,ND-1
	  TA(I)=-HL(I-1)-EPS_A(I-1)
	  TC(I)=-HU(I)+EPS_B(I)
	  TB(I)=GAM_RSQ_DTAUONQ(I)*(P_J(I)-COH_VEC(I)) + PSI(I) + HL(I) + 
	1             HU(I-1)-EPS_B(I-1)+EPS_A(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=GAM_RSQ_DTAUONQ(I)*SOURCE(I)/GAM_REL(I)
	END DO
	XM(1)=0.0D0; XM(ND)=0.0D0
!
! Evaluate TA,TB,TC for boundary conditions
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	PSI(1)=GAM_RSQ(1)*DELTA(1)*( HBC-NBC+(NBC+BETA(1)*K_ON_J(1))*
	1             GAM_REL_SQ(1)*(SIGMA(1)+1.0D0) )
	PSIPREV(1)=GAM_RSQ(1)*DELTA(1)*( HBC_PREV-NBC_PREV+(NBC_PREV+
	1             BETA(1)*K_ON_J_PREV(1))*GAM_REL_SQ(1)*(SIGMA(1)+1.0D0) )
!
	TA(1)=0.0D0
	TC(1)= -GAM_RSQ(2)*( K_ON_J(2)+VdHdR_TERM(2) )*Q(2)/DTAU_H(1)
	TB(1)=  GAM_RSQ(1)*( K_ON_J(1)+VdHdR_TERM(1) )*Q(1)/DTAU_H(1) + 
	1                    PSI(1) + HBC*GAM_RSQ(1)*P_H(1)
	XM(1)=0.0D0
	VB(1)=0.0D0
	VC(1)=0.0D0
!
! Need to include relativistic terms.
!
! CHECK INHBC TERM.
!
	J_AT_INB_EQ_B=.FALSE.     !TRUE.
	IF(J_AT_INB_EQ_B)THEN
	  TA(ND)=0.0D0; TC(ND)=0.0D0;
	  TB(ND)=1.0D0
	  XM(ND)=IC
	ELSE IF(DIF)THEN
	  TA(ND)=-GAM_RSQ(ND-1)*Q(ND-1)*(K_ON_J(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	  TB(ND)=GAM_RSQ(ND)*(K_ON_J(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)
	  XM(ND)=GAM_RSQ(ND)*DBB/3.0D0/CHI(ND)
	ELSE
	  TA(ND)=-GAM_RSQ(ND-1)*Q(ND-1)*(K_ON_J(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	  TB(ND)=GAM_RSQ(ND)*(K_ON_J(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)+IN_HBC*GAM_RSQ(ND)
	  XM(ND)=GAM_RSQ(ND)*IC*(0.25D0+0.5D0*IN_HBC)
	END IF
	TC(ND)=0.0D0
	VB(ND)=0.0D0
	VC(ND)=0.0D0
	PSIPREV(ND)=0.0D0
!
! We create PSIPREV_MOD to save multiplications in the UP_TX_TVX routine/
! It is only different from PSIPREV when N_ON_J is non zero.
!
	PSIPREV_MOD(1)=PSIPREV(1)
	PSIPREV_MOD(ND)=PSIPREV(ND)
	DO I=2,ND-1
	  PSIPREV_MOD(I)=(EPS_PREV_A(I)-EPS_PREV_B(I-1)) + PSIPREV(I)
	END DO
!
	XM(1)=XM(1) + PSIPREV(1)*JNU_PREV(1)
	DO I=2,ND-1
	  XM(I)=XM(I) + VB(I)*GAM_RSQHNU_PREV(I-1) + VC(I)*GAM_RSQHNU_PREV(I)
	1          + PSIPREV_MOD(I)*JNU_PREV(I)
	1          + ( EPS_PREV_B(I)*JNU_PREV(I+1) - EPS_PREV_A(I-1)*JNU_PREV(I-1) )
	END DO
	XM(ND)=XM(ND)
!
!	IF(INIT)THEN
	IF(ABS(FREQ-49.8654D0) .LT. 0.001D0)THEN
	  WRITE(172,*)'TA'
	  WRITE(172,*)TA
	  WRITE(172,*)'TB'
	  WRITE(172,*)TB
	  WRITE(172,*)'TC'
	  WRITE(172,*)TC
	  WRITE(172,*)'XM'
	  WRITE(172,*)XM
	  WRITE(172,*)'HL'
	  WRITE(172,*)HL
	  WRITE(172,*)'HU'
	  WRITE(172,*)HU
	  WRITE(172,*)'HS'
	  WRITE(172,*)HS
	  CLOSE(UNIT=172)
	 END IF
!
! Solve for the radiation field along ray for this frequency.
!
	TA_SAV=TA; TB_SAV=TB; TC_SAV=TC; XM_SAV=XM
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
	DO I=1,ND
	  IF(XM(I) .GT. 1.0D+20)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in VAR_JREL_V2: RJ blowing up'
	    WRITE(LUER,'(6ES12.4)')(XM(J),J=1,ND)
	    STOP
	  END IF
	END DO
!
!	WRITE(211,'(8ES14.4)')FREQ,GAM_RSQ(1:7)*XM(1:7)
!
! Check that no negative mean intensities have been computed.
!
	DO I=1,ND
	  IF(XM(I) .LT. 0.0D0)THEN
	    XM(I)=ABS(XM(I))/10.0D0
	  END IF
	END DO
!
! Store J and evaluate GAM.R^2.H
!
	JNU(1:ND)=XM(1:ND)
	DO I=1,ND-1
	  GAM_RSQHNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*GAM_RSQHNU_PREV(I) +
	1              (EPS_PREV_A(I)*JNU_PREV(I)  -EPS_A(I)*XM(I)) +
	1              (EPS_PREV_B(I)*JNU_PREV(I+1)-EPS_B(I)*XM(I+1))
	END DO
!	WRITE(157,'(7ES16.6)')FREQ,XM(1),XM(1)*GAM_RSQ(1),GAM_RSQHNU(1),
!	1                          XM(ND),XM(ND)*GAM_RSQ(ND),GAM_RSQHNU(ND-1)
!
! NB: A & C are d[dx/dr]/dx where x=gam.r^2.J
!
	IF(INCL_ADVEC_TERMS)THEN
	  TC(1:ND)=GAM_RSQ(1:ND)*XM(1:ND)
	  CALL DERIVCHI(TB,TC,R,ND,'LINMON')
	  CALL d_DERIVCHI_dCHI(TB,TC,R,ND,'LINMON')
!	  DO I=1,ND
!	   WRITE(216,'(I3,5ES18.8)')I,FREQ,GAM_RSQ(I)*XM(I),dlnGRSQJdlnR(I),R(I)*TB(I)/GAM_RSQ(I)/XM(I)
!	  END DO
	  TA=TA_SAV; TB=TB_SAV; TC=TC_SAV
!	  WRITE(215,*)'New freq'
!	  DO I=1,ND
!	    WRITE(215,'(I3,5ES18.8)')I,FREQ,TA(I),TB(I),TC(I),(TA(I)+TB(I)+TC(I))/TB(I)
!	  END DO
	  DO I=2,ND-1
	    TA(I)=TA(I)+GAM_RSQ_DTAUONQ(I)*BETA(I)*A(I)*GAM_RSQ(I-1)/GAM_RSQ(I)/CHI_J(I)
	    TB(I)=TB(I)+GAM_RSQ_DTAUONQ(I)*BETA(I)*(B(I)-dlnGRSQJdlnR(I)/R(I))/CHI_J(I)
	    TC(I)=TC(I)+GAM_RSQ_DTAUONQ(I)*BETA(I)*C(I)*GAM_RSQ(I+1)/GAM_RSQ(I)/CHI_J(I)
	  END DO
!	  DO I=1,ND
!	    T1=R(I)*( A(I)*GAM_RSQ(I-1)*XM(I-1) + B(I)*GAM_RSQ(I)*XM(I) + C(I)*GAM_RSQ(I+1)*XM(I+1))/XM(I)/GAM_RSQ(I)
!	    WRITE(215,'(I3,7ES18.8)')I,FREQ,TA(I),TB(I),TC(I),(TA(I)+TB(I)+TC(I))/TB(I),dlnGRSQJdlnR(I),T1
!	  END DO
	  CALL THOMAS(TA,TB,TC,XM_SAV,ND,1)
	END IF
!	WRITE(211,'(8ES14.4)')FREQ,GAM_RSQ(1:7)*XM_SAV(1:7)
!
! 
!
! The J and H components of the radiation field have been found. We can thus
! begin the variation computation.
!
! Compute d{non-radiation field}/dchi matrix.
!
	CALL TUNE(1,'EDD_JREL_V1')
	CALL EDD_JREL_VAR_V2(KI,RHS_dHdCHI,
	1          R,SIGMA,CHI,ESEC,K_ON_J,
	1          DBB,DIF,J_AT_INB_EQ_B,METHOD,ND,NM_KI)
	CALL TUNE(2,'EDD_JREL_V1')
!
!	IF( ABS(FREQ-9.487047E-01) .LT. 1.0D-05)THEN
!	  OPEN(UNIT=7,FILE='KI_TEST',STATUS='UNKNOWN')
!	    CALL WR2D_V2(KI(1,1,1),ND,ND,'dCHI','*',.TRUE.,7) 
!	    CALL WR2D_V2(KI(1,1,2),ND,ND,'dETA','#',.TRUE.,7) 
!	    CALL WR2D_V2(RHS_dHdCHI,ND,ND,'dH','%',.TRUE.,7) 
!	  CLOSE(UNIT=7)
!	END IF
!
! Evaluate the intensity variations.
!                                  TX=dJ/d(chi,eta,....)
!                          and     TVX=dRSQH/d(chi,eta,...)
!
! WORKMAT is dimension (ND,ND) is is used to temporarily save TX( , ,K) for
! each K.
!
	CALL TUNE(1,'UP_TX')
	CALL UP_TX_TVX(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                       WORKMAT,ND,NM,NM_KI,
	1                       INIT,DO_THIS_TX_MATRIX)
	CALL TUNE(2,'UP_TX')
!
! 
!
! Evaluate diffusion variation. TXD is initially the value from the
! previous frequency.
!
	DO I=1,ND
	  TX_OLD_d_T(I)=TX_DIF_d_T(I)
	  TX_OLD_d_dTdR(I)=TX_DIF_d_dTDR(I)
	END DO
	IF(DIF)THEN
	  TX_DIF_d_T(1)=PSIPREV(1)*TX_DIF_d_T(1)
	  TX_DIF_d_dTdR(1)=PSIPREV(1)*TX_DIF_d_dTdR(1)
	  DO I=2,ND-1
	    TX_DIF_d_T(I)= PSIPREV_MOD(I)*TX_DIF_d_T(I)
	1             + VB(I)*TVX_DIF_d_T(I-1) + VC(I)*TVX_DIF_d_T(I)
	1          + ( EPS_PREV_B(I)*TX_OLD_d_T(I+1)
	1               - EPS_PREV_A(I-1)*TX_OLD_d_T(I-1) )
	    TX_DIF_d_dTdR(I)= PSIPREV_MOD(I)*TX_DIF_d_dTdR(I)
	1             + VB(I)*TVX_DIF_d_dTdR(I-1) + VC(I)*TVX_DIF_d_dTdR(I)
	1          + ( EPS_PREV_B(I)*TX_OLD_d_dTdR(I+1)
	1               - EPS_PREV_A(I-1)*TX_OLD_d_dTdR(I-1) )
	  END DO
	  TX_DIF_d_T(ND)=GAM_RSQ(ND)*dDBBdT/3.0D0/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_T(ND)
	  TX_DIF_d_dTdR(ND)=GAM_RSQ(ND)*DBB/dTdR/3.0D0/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_dTdR(ND)
!
! Solve for the radiation field along ray for this frequency.
!
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_T,ND,1)
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_dTdR,ND,1)
!
	  DO I=1,ND-1
	    TVX_DIF_d_T(I)=HU(I)*TX_DIF_d_T(I+1) - HL(I)*TX_DIF_d_T(I) +
	1        HS(I)*TVX_DIF_d_T(I) +
	1        (EPS_PREV_A(I)*TX_OLD_d_T(I)-EPS_A(I)*TX_DIF_d_T(I)) +
	1        (EPS_PREV_B(I)*TX_OLD_d_T(I+1)-EPS_B(I)*TX_DIF_d_T(I+1))
	    TVX_DIF_d_dTdR(I)=HU(I)*TX_DIF_d_dTdR(I+1) -
	1     HL(I)*TX_DIF_d_dTdR(I) +
	1     HS(I)*TVX_DIF_d_dTdR(I) +
	1     (EPS_PREV_A(I)*TX_OLD_d_dTdR(I)-EPS_A(I)*TX_DIF_d_dTdR(I)) +
	1     (EPS_PREV_B(I)*TX_OLD_d_dTdR(I+1)-EPS_B(I)*TX_DIF_d_dTdR(I+1))
	  END DO
!
	END IF	    	    !DIF END
!
! Save JNU and GAM_RSQHNU for next frequency integration.
!
	JNU_PREV(:)=JNU(:)
	GAM_RSQHNU_PREV(:)=GAM_RSQHNU(:)

!
	RETURN
	END
