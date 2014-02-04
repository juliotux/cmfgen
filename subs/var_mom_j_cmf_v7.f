C
C This subroutine computes the variation of J as a function of the emissivity
C ETA and the opacity CHI.
C
C TX has dimensions (ND,ND,NM) and TVX has dimension (ND-1,ND,NM).
C
C NM reflects the total number of matrices required to compute the variation
C    of J. All NM matrices are modified by a call to THOMAS in UP_TX_TVX. The
C    first 2 should represent dJ/dCHI and the second dJ/dETA --- the others
C    can be in any order.
C
C Typically they will be CHI_C, and ETA_C, CHIL_, AND ETAL of othe lines.
C
C NB TX(i,m, ) =dJ(i)/d(CHI(m),ETA(m),...)
C NB TVX(i,m, ) =dRSQH(i)/d(CHI(m),ETA(m),...)
C
C NM_KI reflects the 3rd dimension of KI. For this routine only the first 2 
C   are important and are used to compute the variation of J with respect to 
C   the CURRENT opacity and the CURRENT emissivity. 
C
C TX_DIF_d_T and TX_DIFF_d_dTdR are used to describe the variations in J
C caused by the DIFFUSION approximation at the inner boundary.
C
C The particular choice of the outer boundary condition adopted is irrelevant
C for this routine. Such information is incorporated by the outer boundary
C Eddington factors HBC, and NBC.
C
C Note TA, TBC, TC, HU etc have been defined so that the program computes
C JNU and r^2 HNU.
C
	SUBROUTINE VAR_MOM_J_CMF_V7(ETA,CHI,ESEC,
	1                  THETA,V,SIGMA,R,
	1                  TX,TVX,TX_DIF_d_T,TX_DIF_d_dTdR,
	1                  TVX_DIF_d_T,TVX_DIF_d_dTdR,
	1                  KI,WORKMAT,RHS_dHdCHI,
	1                  F,G,RSQN_ON_RSQJ,HBC,IN_HBC,NBC,
	1                  F_PREV,G_PREV,RSQN_ON_RSQJ_PREV,
	1                  HBC_PREV,IN_HBC_PREV,NBC_PREV,
	1                  INIT,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1                  DO_THIS_TX_MATRIX,METHOD,ND,NM,NM_KI)
	IMPLICIT NONE
C
C Altered:   04-Feb-1997 : THETA now passed. Passed to allow coherent
C                            scattering at some depths, and non-coherent
C                            scattering at other depths. Only introduced
C                            for the variation calculation. Logical variable
C                            COHERENT no longer required.
C Altered:   12-Dec-1996 : NM_KI installed. Changed to version V6.
C Altered:   05-Dec-1996 : PROGDESC set to REAL*8 value, PROG_ID installed
C                             ERROR_LU installed. TUNE installed.
C Altered:   25-Jan-1996 : HBC, NBC (and HBC_PREV, NBC_PREV) are now
C                            scalers. Other quantities were not needed
C                            with THK option processed by EXTENSION.
C                            Several lines deleted with HBC(2) etc.
C                          Changed to V5
C
C Altered:   11-Jan-1996 : Bug in calculation of DT_DIF_d_dTDR when using
C                             RSQN_ON_RSQJ. Initiliazation improved.
 
C Altered:   11-Mar-1995 : RSQN_ON_RSQJ installed so that N may be written in
C                          terms of H (using G) or J.
C                          Call modified, as were subroutines UP_TX_TVX and
C                          EDD_J_VAR.
C                          _V4 append to name.
C
C Finalized: 04-Nov-1994
C Created:   27-Sep-1995 : Diffusion approximation not tested yet.
C
	INTEGER ND
	INTEGER NM
	INTEGER NM_KI
	REAL*8 ETA(ND),CHI(ND),ESEC(ND),THETA(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND)
C
C Variation arrays and vectors.
C
	REAL*8 TX(ND,ND,NM),TVX(ND-1,ND,NM)
	REAL*8 KI(ND,ND,NM_KI)
	REAL*8 WORKMAT(ND,ND),RHS_dHdCHI(ND-1,ND)
	REAL*8 TX_DIF_d_T(ND),TX_DIF_d_dTdR(ND)
	REAL*8 TVX_DIF_d_T(ND),TVX_DIF_d_dTdR(ND)
C
	LOGICAL DO_THIS_TX_MATRIX(NM)
C
C "Eddington factors"
C
	REAL*8 F(ND),G(ND),RSQN_ON_RSQJ(ND)
	REAL*8 HBC,NBC,IN_HBC
	REAL*8 F_PREV(ND),G_PREV(ND),RSQN_ON_RSQJ_PREV(ND)
	REAL*8 HBC_PREV,NBC_PREV,IN_HBC_PREV
C
	REAL*8 dLOG_NU,dTdR,DBB,dDBBdT,IC
	CHARACTER*6 METHOD
C
C INIT is used to indicate that there is no coupling to the previous frequency.
C We are thus solving the normal continuum transfer equation (i.e. the absence
C of velocity fields)
C
	LOGICAL DIF,INIT
C
C Vectors required by future calls to VAR_MOM_J_CMF.
C
	INTEGER NV
	PARAMETER (NV=200)
	REAL*8 JNUM1(NV),RSQ_HNUM1(NV)
	SAVE JNUM1,RSQ_HNUM1
C
C Work vectors.
C
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,DTAU,RSQ_DTAUONQ,
	1                   XM,SOURCE,Q,JNU,RSQ_HNU,
	1                   VB,VC,HU,HL,HS,
	1                   EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                   TX_OLD_d_T,TX_OLD_d_dTdR,
	1                   GAM,GAMH,W,WPREV,PSI,PSIPREV_MOD,PSIPREV
C
	REAL*8 TA(NV),TB(NV),TC(NV),DTAU(NV),RSQ_DTAUONQ(NV)
	REAL*8 XM(NV),SOURCE(NV),Q(NV),JNU(NV),RSQ_HNU(NV)
	REAL*8 VB(NV),VC(NV),HU(NV),HL(NV),HS(NV)
	REAL*8 GAM(NV),GAMH(NV),W(NV),WPREV(NV)
	REAL*8 PSI(NV),PSIPREV_MOD(NV),PSIPREV(NV)
	REAL*8 EPS_A(NV),EPS_B(NV)
	REAL*8 EPS_PREV_A(NV),EPS_PREV_B(NV)
	REAL*8 TX_OLD_d_T(NV),TX_OLD_d_dTdR(NV)
C
	REAL*8 PROGDESC	
	REAL*8, PARAMETER :: PROG_ID=2.2281463D+08  !Must be unique (VAR_MOM_)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER I
	REAL*8 AV_SIGMA
C
C PROGDESC is a variable use to confirm that the scratch block is not
C being used by some other routine.
C
	PROGDESC=PROG_ID
	IF(ND .GT. NV)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in VAR_MOM_J_CMF - NV smaller than ND'
	  WRITE(I,*)'ND=',ND,'NV',NV
	  STOP
	END IF
C
C Zero common block. There are currently 29 vectors in the common block.
C TA must be the first vector, and PSIPREV the last.
C
	PSIPREV(NV-1)=1.0D0
	PSIPREV(NV)=1.0D0
	I=(NV*28)-1
	CALL DP_ZERO(TA,I)
	IF(PSIPREV(NV-1) .NE. 0 .OR. PSIPREV(NV) .NE. 1)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in zeroing SCRATCH block in VAR_MOM_J_CMF'
	  STOP
	ELSE
	  PSIPREV(NV)=0.0D0
	END IF
C
C 
C
C Zero relevant vectors and matrices.
C
	DO I=1,ND
	  JNU(I)=0.0D0
	  RSQ_HNU(I)=0.0D0
	END DO
C
	IF(INIT)THEN
	  CALL DP_ZERO(TX, ND*ND*NM )
	  CALL DP_ZERO(TVX, (ND-1)*ND*NM )
	  DO I=1,ND
	    JNUM1(I)=0.0D0
	    RSQ_HNUM1(I)=0.0D0
	    GAMH(I)=0.0D0
	    GAM(I)=0.0D0
	    W(I)=0.0D0
	    WPREV(I)=0.0D0
	    PSI(I)=0.0D0
	    PSIPREV(I)=0.0D0
	    TX_DIF_d_T(I)=0.0D0
	    TX_DIF_d_dTdR(I)=0.0D0
	    EPS_A(I)=0.0D0
	    EPS_B(I)=0.0D0
	    EPS_PREV_A(I)=0.0D0
	    EPS_PREV_B(I)=0.0D0
 	  END DO
	END IF
C
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
C
C Compute the Q factors from F. Then compute optical depth scale.
C
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA work vector
	DO I=1,ND
	  TA(I)=CHI(I)*Q(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
C
C We need to call d_DERIVCHI_dCHI to set the TRAP derivatives.
C
	CALL d_DERIVCHI_dCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
	DO I=2,ND
	  RSQ_DTAUONQ(I)=0.5D0*R(I)*R(I)*(DTAU(I)+DTAU(I-1))/Q(I)
	END DO
C
C 
C
C Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
C 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
C Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
C
C NB - By definition, G is defined at the mesh midpoints.
C
C We evaluate and store the constant terms in the computation of GAMH and
C GAM, since number of operations only proportional to ND. Later on scaling is
C proportional to NM*ND*ND.
C
	IF(.NOT. INIT)THEN
	  DO I=1,ND-1
	    AV_SIGMA=0.5D0*(SIGMA(I)+SIGMA(I+1))
	    GAMH(I)=2.0D0*3.33564D-06*(V(I)+V(I+1))/(R(I)+R(I+1))
	1         /dLOG_NU/( CHI(I)+CHI(I+1) )
	    W(I)=GAMH(I)*( 1.0D0+AV_SIGMA*G(I) )
	    WPREV(I)=GAMH(I)*( 1.0D0+AV_SIGMA*G_PREV(I) )
	    EPS_A(I)=GAMH(I)*AV_SIGMA*RSQN_ON_RSQJ(I)/(1.0D0+W(I))
	    EPS_B(I)=EPS_A(I)*R(I+1)*R(I+1)
	    EPS_A(I)=EPS_A(I)*R(I)*R(I)
	    EPS_PREV_A(I)=GAMH(I)*AV_SIGMA*RSQN_ON_RSQJ_PREV(I)/(1.0D0+W(I))
	    EPS_PREV_B(I)=EPS_PREV_A(I)*R(I+1)*R(I+1)
	    EPS_PREV_A(I)=EPS_PREV_A(I)*R(I)*R(I)
	  END DO
	  DO I=1,ND
	    GAM(I)=3.33564D-06*V(I)/R(I)/CHI(I)/dLOG_NU
	  END DO
C
C PSIPREV is equivalent to the U vector of FORMSOL.
C
	  PSI(1)=R(1)*R(1)*GAM(1)*( HBC+NBC*SIGMA(1) )
	  PSIPREV(1)=R(1)*R(1)*GAM(1)*( HBC_PREV+NBC_PREV*SIGMA(1) )
	  DO I=2,ND
	    PSI(I)=RSQ_DTAUONQ(I)*GAM(I)*( 1.0D0+SIGMA(I)*F(I) )
	    PSIPREV(I)=RSQ_DTAUONQ(I)*GAM(I)*
	1                   ( 1.0D0+SIGMA(I)*F_PREV(I) )
	  END DO
	END IF
C
C 
C
C Compute vectors used to compute the flux vector H.
C
	DO I=1,ND-1
	  HU(I)=R(I+1)*R(I+1)*F(I+1)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=R(I)*R(I)*F(I)*Q(I)/(1.0D0+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0D0+W(I))
	END DO
C
C Compute the TRIDIAGONAL operators, and the RHS source vector.
C
	DO I=2,ND-1
	  TA(I)=-HL(I-1)-EPS_A(I-1)
	  TC(I)=-HU(I)+EPS_B(I)
	  TB(I)=RSQ_DTAUONQ(I)*(1.0D0-THETA(I)) + PSI(I) +HU(I-1) +HL(I)
	1             -EPS_B(I-1)+EPS_A(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=RSQ_DTAUONQ(I)*SOURCE(I)
	END DO
C
C Evaluate TA,TB,TC for boundary conditions
C
	TC(1)=-R(2)*R(2)*F(2)*Q(2)/DTAU(1)
	TB(1)=R(1)*R(1)*( F(1)*Q(1)/DTAU(1) + HBC ) + PSI(1)
	XM(1)=0.0D0
	TA(1)=0.0D0
	VB(1)=0.0D0
	VC(1)=0.0D0
C
	TA(ND)=-R(ND-1)*R(ND-1)*F(ND-1)*Q(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  TB(ND)=R(ND)*R(ND)*F(ND)/DTAU(ND-1)
	  XM(ND)=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	ELSE
	  TB(ND)=R(ND)*R(ND)*F(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)
	END IF
	TC(ND)=0.0D0
	VB(ND)=0.0D0
	VC(ND)=0.0D0
	PSIPREV(ND)=0.0D0
C
C We create PSIPREV_MOD to save multiplications in the UP_TX_TVX routine/
C It is only different from PSIPREV when N_ON_J is non zero.
C
	PSIPREV_MOD(1)=PSIPREV(1)
	PSIPREV_MOD(ND)=PSIPREV(ND)
	DO I=2,ND-1
	  PSIPREV_MOD(I)=(EPS_PREV_A(I)-EPS_PREV_B(I-1)) + PSIPREV(I)
	END DO
C
	XM(1)=XM(1) + PSIPREV_MOD(1)*JNUM1(1)
	DO I=2,ND-1
	  XM(I)=XM(I) + VB(I)*RSQ_HNUM1(I-1) + VC(I)*RSQ_HNUM1(I)
	1          + PSIPREV_MOD(I)*JNUM1(I)
	1          + ( EPS_PREV_B(I)*JNUM1(I+1)
	1               - EPS_PREV_A(I-1)*JNUM1(I-1) )
	END DO
	XM(ND)=XM(ND) + PSIPREV_MOD(ND)*JNUM1(ND)
C
C Solve for the radiation field along ray for this frequency.
C
	CALL THOMAS(TA,TB,TC,XM,ND,1)
C
	DO I=1,ND
	  JNU(I)=XM(I)
	END DO
C
	DO I=1,ND-1
	  RSQ_HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQ_HNUM1(I) +
	1              (EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*XM(I)) +
	1              (EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*XM(I+1))
	END DO
C 
C
C The J and H components of the radiation field have been found. We can thus
C begin the variation computation.
C
C Compute d{non-radiation field}/dchi matrix.
C
	CALL TUNE(1,'MOM_EDD')
	CALL EDD_J_VAR_V5(KI,RHS_dHdCHI,WORKMAT,
	1                SOURCE,CHI,ESEC,THETA,DTAU,R,SIGMA,
	1                F,Q,HU,HL,HS,RSQ_DTAUONQ,
	1                W,WPREV,PSI,PSIPREV,
	1                EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                JNU,JNUM1,RSQ_HNUM1,
	1                DBB,DIF,ND,NM_KI)
	CALL TUNE(2,'MOM_EDD')
C
C Evaluate the intensity variations.
C                                  TX=dJ/d(chi,eta,....)
C                          and     TVX=dRSQH/d(chi,eta,...)
C
C WORKMAT is dimension (ND,ND) is is used to temporarily save TX( , ,K) for
C each K.
C
	CALL TUNE(1,'UP_TX')
	CALL UP_TX_TVX(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                       WORKMAT,ND,NM,NM_KI,
	1                       INIT,DO_THIS_TX_MATRIX)
	CALL TUNE(2,'UP_TX')
C
C 
C
C Evaluate diffusion variation. TXD is initially the value from the
C previous frequency.
C
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
	  TX_DIF_d_T(ND)=dDBBdT*R(ND)*R(ND)/3.0D0/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_T(ND)
	  TX_DIF_d_dTdR(ND)=DBB/dTdR*R(ND)*R(ND)/3.0D0/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_dTdR(ND)
C
C Solve for the radiation field along ray for this frequency.
C
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_T,ND,1)
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_dTdR,ND,1)
C
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
C
	END IF	    	    !DIF END
C
C 
C
C Save JNU and RSQ_HNU for next frequency integration.
C
	DO I=1,ND
	  JNUM1(I)=JNU(I)
	  RSQ_HNUM1(I)=RSQ_HNU(I)
	END DO
C
C 
C
	IF(PROGDESC .NE. PROG_ID)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error - SCRATCH block corrupted in VAR_MOM_J_CMF'
	  STOP
	END IF
C
	RETURN
	END
