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
	SUBROUTINE VAR_MOM_J_CMF_V8(ETA,CHI,ESEC,
	1               THETA,V,SIGMA,R,
	1               TX,TVX,TX_DIF_d_T,TX_DIF_d_dTdR,
	1               TVX_DIF_d_T,TVX_DIF_d_dTdR,
	1               KI,WORKMAT,RHS_dHdCHI,
	1               F,G,RSQN_ON_RSQJ,HBC,IN_HBC,NBC,
	1               F_PREV,G_PREV,RSQN_ON_RSQJ_PREV,
	1               HBC_PREV,IN_HBC_PREV,NBC_PREV,
	1               INIT,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1               FREQ,OUT_BC_TYPE,
	1               DO_THIS_TX_MATRIX,METHOD,ND,NM,NM_KI)
	IMPLICIT NONE
!
! Altered:   06-Aug-2007 : Installed integer OUT_BC_TYPE
!                            so that could have choice of boundary conditions
!                            This allows other options to be included.
! Altered:   04-Feb-1997 : THETA now passed. Passed to allow coherent
!                            scattering at some depths, and non-coherent
!                            scattering at other depths. Only introduced
!                            for the variation calculation. Logical variable
!                            COHERENT no longer required.
! Altered:   12-Dec-1996 : NM_KI installed. Changed to version V6.
! Altered:   05-Dec-1996 : PROGDESC set to REAL*8 value, PROG_ID installed
!                             ERROR_LU installed. TUNE installed.
! Altered:   25-Jan-1996 : HBC, NBC (and HBC_PREV, NBC_PREV) are now
!                            scalers. Other quantities were not needed
!                            with THK option processed by EXTENSION.
!                            Several lines deleted with HBC(2) etc.
!                          Changed to V5
!
! Altered:   11-Jan-1996 : Bug in calculation of DT_DIF_d_dTDR when using
!                             RSQN_ON_RSQJ. Initiliazation improved.
 
! Altered:   11-Mar-1995 : RSQN_ON_RSQJ installed so that N may be written in
!                          terms of H (using G) or J.
!                          Call modified, as were subroutines UP_TX_TVX and
!                          EDD_J_VAR.
!                          _V4 append to name.
!
! Finalized: 04-Nov-1994
! Created:   27-Sep-1995 : Diffusion approximation not tested yet.
!
	INTEGER ND
	INTEGER NM
	INTEGER NM_KI
	REAL*8 ETA(ND),CHI(ND),ESEC(ND),THETA(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND)
!
! Variation arrays and vectors.
!
	REAL*8 TX(ND,ND,NM),TVX(ND-1,ND,NM)
	REAL*8 KI(ND,ND,NM_KI)
	REAL*8 WORKMAT(ND,ND),RHS_dHdCHI(ND-1,ND)
	REAL*8 TX_DIF_d_T(ND),TX_DIF_d_dTdR(ND)
	REAL*8 TVX_DIF_d_T(ND),TVX_DIF_d_dTdR(ND)
!
	LOGICAL DO_THIS_TX_MATRIX(NM)
!
! "Eddington factors"
!
	REAL*8 F(ND),G(ND),RSQN_ON_RSQJ(ND)
	REAL*8 HBC,NBC,IN_HBC
	REAL*8 F_PREV(ND),G_PREV(ND),RSQN_ON_RSQJ_PREV(ND)
	REAL*8 HBC_PREV,NBC_PREV,IN_HBC_PREV
!
	REAL*8 dLOG_NU,dTdR,DBB,dDBBdT,IC
	REAL*8 FREQ
	CHARACTER*6 METHOD
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
	LOGICAL DIF,INIT
	INTEGER OUT_BC_TYPE
!
! Vectors required by future calls to VAR_MOM_J_CMF.
!
	INTEGER NV
	PARAMETER (NV=1000)
	REAL*8 JNUM1(NV),RSQ_HNUM1(NV)
	SAVE JNUM1,RSQ_HNUM1
!
! Work vectors.
!
	COMMON /SCRATCH/  PROGDESC,TA,TB,TC,DTAU,RSQ_DTAUONQ,
	1                   XM,SOURCE,Q,JNU,RSQ_HNU,
	1                   VB,VC,HU,HL,HS,
	1                   EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                   TX_OLD_d_T,TX_OLD_d_dTdR,
	1                   GAM,GAMH,W,WPREV,PSI,PSIPREV_MOD,PSIPREV
!
	REAL*8 TA(NV),TB(NV),TC(NV),DTAU(NV),RSQ_DTAUONQ(NV)
	REAL*8 XM(NV),SOURCE(NV),Q(NV),JNU(NV),RSQ_HNU(NV)
	REAL*8 VB(NV),VC(NV),HU(NV),HL(NV),HS(NV)
	REAL*8 GAM(NV),GAMH(NV),W(NV),WPREV(NV)
	REAL*8 PSI(NV),PSIPREV_MOD(NV),PSIPREV(NV)
	REAL*8 EPS_A(NV),EPS_B(NV)
	REAL*8 EPS_PREV_A(NV),EPS_PREV_B(NV)
	REAL*8 TX_OLD_d_T(NV),TX_OLD_d_dTdR(NV)
!
	REAL*8 PROGDESC	
	REAL*8 T1
	REAL*8 DTAU_BND
	REAL*8, PARAMETER :: PROG_ID=2.2281463D+08  !Must be unique (VAR_MOM_)
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER I
	REAL*8 AV_SIGMA
!
! PROGDESC is a variable use to confirm that the scratch block is not
! being used by some other routine.
!
	PROGDESC=PROG_ID
	IF(ND .GT. NV)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in VAR_MOM_J_CMF - NV smaller than ND'
	  WRITE(I,*)'ND=',ND,'NV',NV
	  STOP
	END IF
!
! Zero common block. There are currently 29 vectors in the common block.
! TA must be the first vector, and PSIPREV the last.
!
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
!
! 
!
! Zero relevant vectors and matrices.
!
	DO I=1,ND
	  JNU(I)=0.0D0
	  RSQ_HNU(I)=0.0D0
	END DO
!
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
!
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
!
! Compute the Q factors from F. Then compute optical depth scale.
!
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA work vector
	DO I=1,ND
	  TA(I)=CHI(I)*Q(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
!
! We need to call d_DERIVCHI_dCHI to set the TRAP derivatives.
!
	CALL d_DERIVCHI_dCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
	DO I=2,ND
	  RSQ_DTAUONQ(I)=0.5D0*R(I)*R(I)*(DTAU(I)+DTAU(I-1))/Q(I)
	END DO
!
! 
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
! NB - By definition, G is defined at the mesh midpoints.
!
! We evaluate and store the constant terms in the computation of GAMH and
! GAM, since number of operations only proportional to ND. Later on scaling is
! proportional to NM*ND*ND.
!
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
!
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	  DO I=2,ND
	    PSI(I)=RSQ_DTAUONQ(I)*GAM(I)*( 1.0D0+SIGMA(I)*F(I) )
	    PSIPREV(I)=RSQ_DTAUONQ(I)*GAM(I)*
	1                   ( 1.0D0+SIGMA(I)*F_PREV(I) )
	  END DO
	END IF
!
! 
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=R(I+1)*R(I+1)*F(I+1)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=R(I)*R(I)*F(I)*Q(I)/(1.0D0+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0D0+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	DO I=2,ND-1
	  TA(I)=-HL(I-1)-EPS_A(I-1)
	  TC(I)=-HU(I)+EPS_B(I)
	  TB(I)=RSQ_DTAUONQ(I)*(1.0D0-THETA(I)) + PSI(I) +HU(I-1) +HL(I)
	1             -EPS_B(I-1)+EPS_A(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=RSQ_DTAUONQ(I)*SOURCE(I)
	END DO
!
! Evaluate TA,TB,TC for boundary conditions
!
	IF(OUT_BC_TYPE .LE. 1)THEN
	  PSI(1)=R(1)*R(1)*GAM(1)*( HBC+NBC*SIGMA(1) )
	  PSIPREV(1)=R(1)*R(1)*GAM(1)*( HBC_PREV+NBC_PREV*SIGMA(1) )
	  TC(1)=-R(2)*R(2)*F(2)*Q(2)/DTAU(1)
	  TB(1)=R(1)*R(1)*( F(1)*Q(1)/DTAU(1) + HBC ) + PSI(1)
	  XM(1)=XM(1) + PSIPREV(1)*JNUM1(1)
	  TA(1)=0.0D0
	  VB(1)=0.0D0
	  VC(1)=0.0D0
	ELSE
	  PSI(1)=R(1)*R(1)*GAM(1)*(1.0D0+SIGMA(1)*F(1))
	  PSIPREV(1)=R(1)*R(1)*GAM(1)*(1.0D0+SIGMA(1)*F_PREV(1))
	  T1=0.25D0*(CHI(2)+CHI(1))*(R(2)-R(1))
	  DTAU_BND=T1
	  TC(1)=(HU(1)-EPS_B(1))/T1
	  TB(1)=-(HL(1)+R(1)*R(1)*HBC+EPS_A(1))/T1-PSI(1)-R(1)*R(1)*(1.0D0-THETA(1))
	  XM(1)=-SOURCE(1)*R(1)*R(1)- HS(1)*RSQ_HNUM1(1)/T1- PSIPREV(1)*JNUM1(1)
	  XM(1)=XM(1)-(EPS_PREV_A(1)*JNUM1(1)+EPS_PREV_B(1)*JNUM1(2))/T1
	END IF
!
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
	DO I=2,ND-1
	  XM(I)=XM(I) + VB(I)*RSQ_HNUM1(I-1) + VC(I)*RSQ_HNUM1(I)
	1          + PSIPREV_MOD(I)*JNUM1(I)
	1          + ( EPS_PREV_B(I)*JNUM1(I+1)
	1               - EPS_PREV_A(I-1)*JNUM1(I-1) )
	END DO
	XM(ND)=XM(ND) + PSIPREV_MOD(ND)*JNUM1(ND)
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
	DO I=1,ND
	  JNU(I)=XM(I)
	END DO
!
	DO I=1,ND-1
	  RSQ_HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQ_HNUM1(I) +
	1              (EPS_PREV_A(I)*JNUM1(I)-EPS_A(I)*XM(I)) +
	1              (EPS_PREV_B(I)*JNUM1(I+1)-EPS_B(I)*XM(I+1))
	END DO
! 
!
! The J and H components of the radiation field have been found. We can thus
! begin the variation computation.
!
! Compute d{non-radiation field}/dchi matrix.
!
	CALL TUNE(1,'MOM_EDD')
	CALL EDD_J_VAR_V6(KI,RHS_dHdCHI,WORKMAT,
	1           SOURCE,CHI,ESEC,THETA,DTAU,R,SIGMA,
	1           F,Q,HU,HL,HS,RSQ_DTAUONQ,
	1           W,WPREV,PSI,PSIPREV,
	1           EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1           JNU,JNUM1,RSQ_HNUM1,
	1           DBB,DIF,HBC,OUT_BC_TYPE,
	1           ND,NM_KI)
	CALL TUNE(2,'MOM_EDD')
!
! Evaluate the intensity variations.
!                                  TX=dJ/d(chi,eta,....)
!                          and     TVX=dRSQH/d(chi,eta,...)
!
! WORKMAT is dimension (ND,ND) is is used to temporarily save TX( , ,K) for
! each K.
!
	CALL TUNE(1,'UP_TX')
	CALL UP_TX_TVX_V2(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                       WORKMAT,ND,NM,NM_KI,
	1                       DTAU_BND,OUT_BC_TYPE,
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
	  TX_DIF_d_T(ND)=dDBBdT*R(ND)*R(ND)/3.0D0/CHI(ND) +
	1                       PSIPREV(ND)*TX_DIF_d_T(ND)
	  TX_DIF_d_dTdR(ND)=DBB/dTdR*R(ND)*R(ND)/3.0D0/CHI(ND) +
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
! 
!
! Save JNU and RSQ_HNU for next frequency integration.
!
	DO I=1,ND
	  JNUM1(I)=JNU(I)
	  RSQ_HNUM1(I)=RSQ_HNU(I)
	END DO
!
! 
!
	IF(PROGDESC .NE. PROG_ID)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error - SCRATCH block corrupted in VAR_MOM_J_CMF'
	  STOP
	END IF
!
	RETURN
	END
