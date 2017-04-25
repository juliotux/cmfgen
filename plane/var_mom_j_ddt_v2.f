 	MODULE MOD_VAR_HUB_J_V2
!
	REAL*8, ALLOCATABLE :: JNU(:)
	REAL*8, ALLOCATABLE :: JNUM1(:)
	REAL*8, ALLOCATABLE :: JNU_OLDT(:)
	REAL*8, ALLOCATABLE :: JNU_MOD(:)
	REAL*8, ALLOCATABLE :: RSQ_HNU(:)
	REAL*8, ALLOCATABLE :: RSQ_HNUM1(:)
	REAL*8, ALLOCATABLE :: RSQ_HNU_OLDT(:)
!
	REAL*8, ALLOCATABLE :: dH(:)
	REAL*8, ALLOCATABLE :: dH_OLDT(:)
	REAL*8, ALLOCATABLE :: DJDt(:)
	REAL*8, ALLOCATABLE :: DJDt_OLDT(:)
	REAL*8, ALLOCATABLE :: DTAU(:)
	REAL*8, ALLOCATABLE :: DUMMY_VEC(:)
	REAL*8, ALLOCATABLE :: HU(:)
	REAL*8, ALLOCATABLE :: HL(:)
	REAL*8, ALLOCATABLE :: HS(:)
	REAL*8, ALLOCATABLE :: HT(:)
	REAL*8, ALLOCATABLE :: GAM(:)
	REAL*8, ALLOCATABLE :: GAMH(:)
	REAL*8, ALLOCATABLE :: PSI(:)
	REAL*8, ALLOCATABLE :: PSIPREV(:)
	REAL*8, ALLOCATABLE :: Q(:)
	REAL*8, ALLOCATABLE :: RSQ_DTAUONQ(:)
	REAL*8, ALLOCATABLE :: SOURCE(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: TX_OLD_d_T(:)
	REAL*8, ALLOCATABLE :: TX_OLD_d_dTdR(:)
	REAL*8, ALLOCATABLE :: VB(:)
	REAL*8, ALLOCATABLE :: VC(:)
	REAL*8, ALLOCATABLE :: W(:)
	REAL*8, ALLOCATABLE :: WPREV(:)
	REAL*8, ALLOCATABLE :: XM(:)
!
	REAL*8 CHI_AT_INB_PREV
	REAL*8 DBB_PREV
	REAL*8 dDBBDT_PREV
!
	REAL*8 HONJ_OUTBC_PREV
	REAL*8 HONJ_OUTBC_OLDT
	REAL*8 RSQH_AT_IB_OLDT
	REAL*8 RSQH_AT_IB_PREV
!
	REAL*8 ROLD_ON_R
	REAL*8 DELTA_TIME_SECS
	REAL*8 RECIP_CDELTAT
!
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
 	END MODULE MOD_VAR_HUB_J_V2
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
! Typically they will be CHI_C, and ETA_C, CHIL_, AND ETAL of other lines.
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
	SUBROUTINE VAR_MOM_J_DDT_V2(ETA,CHI,ESEC,THETA,V,R,
	1                  TX,TVX,TX_DIF_d_T,TX_DIF_d_dTdR,
	1                  TVX_DIF_d_T,TVX_DIF_d_dTdR,
	1                  KI,WORKMAT,RHS_dHdCHI,F,
	1                  DO_TIME_VAR,RELAX_PARAM,INIT,FREQ,dLOG_NU,
	1                  dTdR,DBB,dDBBdT,
	1                  DO_THIS_TX_MATRIX,METHOD,
	1                  INNER_BND_METH,OUTER_BND_METH,
	1                  ND,NM,NM_KI)
 	USE MOD_VAR_HUB_J_V2
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
! Altered: 10-Feb-2010 - Minor correction to dRHSdCHI_IB for DIFFUSION approximation.
! Altered:  2-Feb-2010 - Bug fix with DJDt_OLDT(1).
! Altered: 18-Jan-2010 - Changed to V2.
!                        Installed INNER_BND_METH,OUTER_BND_METH and now use MOD_RAY_MOM_STORE.
!                        Changes to allow "THIN" inner boundary condition, and more flexible
!                          boundary conditions.
! Created: 16-July-2006
!
	INTEGER ND
	INTEGER NM
	INTEGER NM_KI
!
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 ETA(ND)
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 THETA(ND)
!
! Variation arrays and vectors.
!
	REAL*8 TX(ND,ND,NM)
	REAL*8 TVX(ND-1,ND,NM)
	REAL*8 KI(ND,ND,NM_KI)
	REAL*8 WORKMAT(ND,ND)
	REAL*8 RHS_dHdCHI(ND-1,ND)
	REAL*8 TX_DIF_d_T(ND)
	REAL*8 TX_DIF_d_dTdR(ND)
	REAL*8 TVX_DIF_d_T(ND)
	REAL*8 TVX_DIF_d_dTdR(ND)
!
	LOGICAL DO_THIS_TX_MATRIX(NM)
!
! "Eddington factors"
!
	REAL*8 F(ND)
	REAL*8 HONJ_OUTBC
	REAL*8 RSQH_AT_IB
	REAL*8 RSQH_AT_OB
!
	REAL*8 dLOG_NU,dTdR,DBB,dDBBdT,IC
	REAL*8 dNU_TERM_DIF_BC
	REAL*8 dRHSdCHI_IB
	REAL*8 dRHSdCHI_OB
	REAL*8 FREQ
	REAL*8 C_KMS
	REAL*8 T1,T2,T3
	REAL*8 GAM_INB
	REAL*8 MOD_DTAU
	REAL*8 HPLUS,FPLUS
	REAL*8 HMIN,FMIN
	REAL*8 RSQ
	REAL*8 RELAX_PARAM
	CHARACTER(LEN=*) METHOD
	CHARACTER(LEN=*) INNER_BND_METH
	CHARACTER(LEN=*) OUTER_BND_METH
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
	LOGICAL DO_TIME_VAR
	LOGICAL INIT
!
	REAL*8 SPEED_OF_LIGHT
	INTEGER ERROR_LU
	EXTERNAL SPEED_OF_LIGHT,ERROR_LU
	CHARACTER(LEN=6), PARAMETER :: OPTION='NORMAL'
!
! Local variables.
!
	INTEGER I
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
!
	IF(FIRST_TIME)THEN
	  ALLOCATE ( JNU(ND) )
	  ALLOCATE ( JNUM1(ND) )
	  ALLOCATE ( JNU_OLDT(ND) )
	  ALLOCATE ( JNU_MOD(ND) )
	  ALLOCATE ( RSQ_HNU(ND) )
	  ALLOCATE ( RSQ_HNUM1(ND) )
	  ALLOCATE ( RSQ_HNU_OLDT(ND) )
!
	  ALLOCATE ( dH(ND) )
	  ALLOCATE ( dH_OLDT(ND) )
	  ALLOCATE ( DJDt(ND) )
	  ALLOCATE ( DJDt_OLDT(ND) )
	  ALLOCATE ( DTAU(ND) )
	  ALLOCATE ( DUMMY_VEC(ND) )
	  ALLOCATE ( HU(ND) )
	  ALLOCATE ( HL(ND) )
	  ALLOCATE ( HS(ND) )
	  ALLOCATE ( HT(ND) )
	  ALLOCATE ( GAM(ND) )
	  ALLOCATE ( GAMH(ND) )
	  ALLOCATE ( PSI(ND) )
	  ALLOCATE ( PSIPREV(ND) )
	  ALLOCATE ( Q(ND) )
	  ALLOCATE ( RSQ_DTAUONQ(ND) )
	  ALLOCATE ( SOURCE(ND) )
	  ALLOCATE ( TA(ND) )
	  ALLOCATE ( TB(ND) )
	  ALLOCATE ( TC(ND) )
	  ALLOCATE ( TX_OLD_d_T(ND) )
	  ALLOCATE ( TX_OLD_d_dTdR(ND) )
	  ALLOCATE ( VB(ND) )
	  ALLOCATE ( VC(ND) )
	  ALLOCATE ( W(ND) )
	  ALLOCATE ( WPREV(ND) )
	  ALLOCATE ( XM(ND) )
!
	  FIRST_TIME=.FALSE.
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
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
!
! We use TA as a temporary variable for ROLD.
!
	IF(DO_TIME_VAR)THEN
	  CALL GET_JH_AT_PREV_TIME_STEP(JNU_OLDt,RSQ_HNU_OLDt,
	1        RSQH_AT_IB_OLDT,HONJ_OUTBC_OLDT,
	1        DELTA_TIME_SECS,FREQ,R,V,ND,INIT,'NORMAL')
	  TA(1:ND)=R(1:ND)-1.0D-05*V(1:ND)*DELTA_TIME_SECS
	  JNU_OLDT(1:ND)=JNU_OLDT(1:ND)/TA(1:ND)/TA(1:ND)
	ELSE
	  JNU_OLDT=0.0D0; RSQ_HNU_OLDt=0.0D0
	  RSQH_AT_IB_OLDT=0.0D0; HONJ_OUTBC_OLDT=0.0D0
	END IF
!
! NB: The factor of 10^10 occurs because c. /\t is a length, and R in
!     cmfgen is in units of 10^10 cm. NB: In the differenced equations
!     we always have terms like 1/(c . /\t . chi)
!
!     The factor of 10^{-5} in ROLD_ON_R occurs because V is in km/s, and
!      R is in units of 10^10cm..
!
	IF(INIT)THEN
	  IF(DO_TIME_VAR)THEN
	    RECIP_CDELTAT=1.0D+10*RELAX_PARAM/SPEED_OF_LIGHT()/DELTA_TIME_SECS
	    ROLD_ON_R=1.0D0-1.0D-05*V(ND)*DELTA_TIME_SECS/R(ND)
	  ELSE
	    RECIP_CDELTAT=0.0D0
	    ROLD_ON_R=0.0D0
	  END IF
	END IF
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
	DO I=2,ND-1
	  RSQ_DTAUONQ(I)=0.5D0*R(I)*R(I)*(DTAU(I)+DTAU(I-1))/Q(I)
	END DO
!
! 
!
! Assume (1) Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s.
!
! We evaluate and store the constant terms in the computation of GAMH and
! GAM, since number of operations only proportional to ND. Later on scaling is
! proportional to NM*ND*ND.
!
	IF(INIT)THEN
	  TX=0.0D0		! ND:ND:NM
	  TVX=0.0D0		! ND-1:ND:NM
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
	    DJDT(I)=0.0D0
	    DJDT_OLDT(I)=0.0D0
 	  END DO
	  CHI_AT_INB_PREV=1.0D0
	ELSE
	  DO I=1,ND-1
	    GAMH(I)=2.0D0*(V(I)+V(I+1))/(R(I)+R(I+1))/dLOG_NU/( CHI(I)+CHI(I+1) )/C_KMS
	    dH(I)=2.0D0*RECIP_CDELTAT/( CHI(I)+CHI(I+1) )
	    dH_OLDT(I)=dH(I)*ROLD_ON_R
	    W(I)=GAMH(I)+dH(I)
	    WPREV(I)=GAMH(I)
	  END DO
	  DO I=1,ND
	    GAM(I)=V(I)/R(I)/CHI(I)/dLOG_NU/C_KMS
	  END DO
	  T1=ROLD_ON_R**3
	  DO I=2,ND-1
	    PSI(I)=RSQ_DTAUONQ(I)*GAM(I)
	    PSIPREV(I)=RSQ_DTAUONQ(I)*GAM(I)
	    DJDT(I)=RSQ_DTAUONQ(I)*RECIP_CDELTAT/CHI(I)
	    DJDT_OLDT(I)=T1*RSQ_DTAUONQ(I)*RECIP_CDELTAT/CHI(I)
	  END DO
	END IF
!
! If if it is the first frequency, we still need to allow for the time variability
! terms.
!
	IF(INIT .AND. DO_TIME_VAR)THEN
	  DO I=1,ND-1
	    dH(I)=2.0D0*RECIP_CDELTAT/( CHI(I)+CHI(I+1) )
	    dH_OLDT(I)=dH(I)*ROLD_ON_R
	    W(I)=dH(I)
	  END DO
	  T1=ROLD_ON_R**3
	  DO I=2,ND-1
	    DJDT(I)=RSQ_DTAUONQ(I)*RECIP_CDELTAT/CHI(I)
	    DJDT_OLDT(I)=T1*RSQ_DTAUONQ(I)*RECIP_CDELTAT/CHI(I)
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
	  HT(I)=dH_OLDT(I)/(1.0D0+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	DO I=2,ND-1
	  TA(I)=-HL(I-1)
	  TC(I)=-HU(I)
	  TB(I)=RSQ_DTAUONQ(I)*(1.0D0-THETA(I)) + PSI(I) + DJDT(I) +HU(I-1) +HL(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=RSQ_DTAUONQ(I)*SOURCE(I)
	END DO
!
	DO I=2,ND-1
	  XM(I)=XM(I) + 
	1        (VB(I)*RSQ_HNUM1(I-1) + VC(I)*RSQ_HNUM1(I)) +
	1        (HT(I)*RSQ_HNU_OLDT(I) - HT(I-1)*RSQ_HNU_OLDT(I-1)) +
	1         PSIPREV(I)*JNUM1(I) + DJDT_OLDT(I)*JNU_OLDt(I)
	END DO
!
! Evaluate TA,TB,TC for boundary conditions
! NB: GAM(1) is zero if inital frequency. Likewise, RECIP_CDELTAT is zero if not doing
!            time variation.
!
! NB: dRHSdCHI_OB is only used for contributions directly related to CHI, but not
! those related to PSI and DJDt.
!
	TA(1)=0.0D0
	VB(1)=0.0D0
	VC(1)=0.0D0
        HONJ_OUTBC=(HPLUS_OB-HMIN_OB)/(JPLUS_OB+JMIN_OB)
        IF(OUTER_BND_METH .EQ. 'HONJ')THEN
	  DJDt(1)=R(1)*R(1)*RECIP_CDELTAT*HONJ_OUTBC/CHI(1)
	  DJDt_OLDT(1)=(ROLD_ON_R**3)*R(1)*R(1)*RECIP_CDELTAT*HONJ_OUTBC_OLDT/CHI(1)
	  PSI(1)=R(1)*R(1)*GAM(1)*HONJ_OUTBC
	  PSIPREV(1)=R(1)*R(1)*GAM(1)*HONJ_OUTBC_PREV
	  TC(1)=-R(2)*R(2)*F(2)*Q(2)/DTAU(1)
	  TB(1)=R(1)*R(1)*( F(1)*Q(1)/DTAU(1) + HONJ_OUTBC ) + PSI(1) + DJDt(1)
	  XM(1)=PSIPREV(1)*JNUM1(1) + DJDT_OLDT(1)*JNU_OLDt(1)
	  dRHSdCHI_OB=0.0D0
!
	ELSE IF(OUTER_BND_METH .EQ. 'HALF_MOM')THEN
	  MOD_DTAU=0.5D0*(CHI(1)+CHI(2))*(R(1)-R(2))
	  RSQ=R(1)*R(1)
	  HPLUS=HPLUS_OB/JPLUS_OB
	  FPLUS=KPLUS_OB/JPLUS_OB
	  PSI(1)=RSQ*HPLUS*GAM(1)
	  PSIPREV(1)=RSQ*HONJ_OUTBC_PREV*GAM(1)
	  T1=R(1)*R(1)*RECIP_CDELTAT/CHI(1)
	  DJDt(1)=HPLUS*T1
	  DJDt_OLDt(1)=(ROLD_ON_R**3)*HONJ_OUTBC_OLDT*T1
	  TC(1)=-R(2)*R(2)*F(2)/MOD_DTAU
	  TB(1)=RSQ*( FPLUS/MOD_DTAU -  (1.0D0-FPLUS)/R(1)/CHI(1) + HPLUS ) + PSI(1) +DJDT(1)
	  XM(1)=RSQ*( HMIN_OB - KMIN_OB/MOD_DTAU+(JMIN_OB-KMIN_OB)/R(1)/CHI(1)  +
	1             GAM(1)*(HMIN_OB+HONJ_OUTBC_PREV*JNUM1(1)) ) +
	1            (T1*HMIN_OB+DJDt_OLDt(1)*JNU_OLDt(1))
	  dRHSdCHI_OB=RSQ*( (JMIN_OB-KMIN_OB)/R(1)/CHI(1)+GAM(1)*HMIN_OB ) + T1*HMIN_OB
	  dRHSdCHI_OB=-dRHSdCHI_OB/CHI(1)
	  XM(2)=XM(2)-JMIN_OB*TA(2)
	ELSE
	  I=ERROR_LU()
	  WRITE(I,*)'Only HONJ & HALF_MD boundary conditions implemented at outer boundary'
	  WRITE(I,*)'Routine is VAR_MOM_J_DDT_V2'
	  STOP
	END IF
!
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  PSI(ND)=0.0D0
	  PSIPREV(ND)=0.0D0
	  DJDT(ND)=0.0D0
	  DJDT_OLDt(ND)=0.0D0
	  RSQH_AT_IB=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	  GAM_INB=0.0D0                               !or set to GAM(ND)
	  dNU_TERM_DIF_BC=GAM_INB*(RSQH_AT_IB-RSQH_AT_IB_PREV)
	  T1=RECIP_CDELTAT*(RSQH_AT_IB-RSQH_AT_IB_OLDt*ROLD_ON_R)/CHI(ND)
	  TA(ND)=-R(ND-1)*R(ND-1)*F(ND-1)*Q(ND-1)/DTAU(ND-1)
	  TB(ND)=R(ND)*R(ND)*F(ND)/DTAU(ND-1)
	  XM(ND)=RSQH_AT_IB+T1+dNU_TERM_DIF_BC
	  dRHSdCHI_IB=RSQH_AT_IB + T1 + dNU_TERM_DIF_BC + (RECIP_CDELTAT/CHI(ND)+GAM_INB)*RSQH_AT_IB
	  dRHSdCHI_IB=-dRHSdCHI_IB/CHI(ND)
!
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
!
	  RSQH_AT_IB=0.0D0
	  PSI(ND)=0.0D0
	  PSIPREV(ND)=0.0D0
	  DJDT(ND)=0.0D0
	  DJDT_OLDt(ND)=0.0D0
	  TA(ND)=-R(ND-1)*R(ND-1)*F(ND-1)*Q(ND-1)/DTAU(ND-1)
	  TB(ND)=R(ND)*R(ND)*F(ND)/DTAU(ND-1)
	  T1=RECIP_CDELTAt*(RSQH_AT_IB-ROLD_ON_R*RSQH_AT_IB_OLDt)/CHI(ND)
	  XM(ND)=RSQH_AT_IB + T1
	  dRHSdCHI_IB=-T1/CHI(ND)
!
	  TB(ND)=TB(ND)+0.1D0*R(ND)*R(ND)*F(ND)/DTAU(ND-1)
	  XM(ND)=XM(ND)+0.1D0*R(ND)*R(ND)*F(ND)*(JPLUS_IB+JMIN_IB)/DTAU(ND-1)
!
! Since Q(ND)=1, we can still use DTAU(ND-1) at the inner boundary.
! [Terms contain a /Q(ND)].
!
	ELSE IF(INNER_BND_METH(1:6) .EQ. 'HOLLOW')THEN
	  RSQ=R(ND)*R(ND)
	  HMIN=HMIN_IB/JMIN_IB
	  FMIN=KMIN_IB/JMIN_IB
	  TA(ND)=-R(ND-1)*R(ND-1)*F(ND-1)/DTAU(ND-1)
	  TB(ND)=RSQ*( FMIN/DTAU(ND-1) + (1.0D0-FMIN)/R(ND)/CHI(ND) + HMIN*(1.0D0+GAM(ND)) )
	  XM(ND)=RSQ*( HPLUS_IB-KPLUS_IB/DTAU(ND-1)-(JPLUS_IB-KPLUS_IB)/R(ND)/CHI(ND)
	1              + GAM(ND)*(HPLUS_IB-RSQH_AT_IB_PREV/RSQ) )
	  XM(ND-1)=XM(ND-1)-JPLUS_IB*TC(ND-1)
	  dRHSdCHI_IB=RSQ*( (JPLUS_IB-KPLUS_IB)/R(ND)/CHI(ND)-GAM(ND)*HPLUS_IB )/CHI(ND)
!
	  TB(ND)=TB(ND)+0.1D0*RSQ*FMIN/DTAU(ND-1)
	  XM(ND)=XM(ND)+0.1D0*RSQ*FMIN*(JPLUS_IB+JMIN_IB)/DTAU(ND-1)
!
! For consistency PSI is the term on the LHS, and PSIRPEV is the term on the RHS.
!
	  IF(.NOT. INIT)THEN
	    PSI(ND)=RSQ*HMIN*GAM(ND) 
	    PSIPREV(ND)=-GAM(ND)*RSQH_AT_IB_PREV/JNUM1(ND)
	  END IF
	  DJDt(ND)=0.0D0
	  DJDt_OLDt(ND)=0.0D0
	ELSE
	  I=ERROR_LU()
	  WRITE(I,*)'Only DIF (diffusion), ZERO_FLUX & HOLLOW boundary conditions currently implemented'
	  WRITE(I,*)'Routine is VAR_MOM_J_DDT_V2'
	  STOP
	END IF
	TC(ND)=0.0D0
	VB(ND)=0.0D0
	VC(ND)=0.0D0
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
	JNU_MOD(1:ND)=XM(1:ND)
        IF(INNER_BND_METH(1:6) .EQ. 'HOLLOW')THEN
          RSQH_AT_IB=R(ND)*R(ND)*(HPLUS_IB-HMIN*XM(ND))
          XM(ND)=XM(ND)+JPLUS_IB
        END IF
        IF(OUTER_BND_METH .EQ. 'HALF_MOM')THEN
          RSQH_AT_OB=R(1)*R(1)*(HPLUS*XM(1)-HMIN_OB)
          XM(1)=XM(1)+JMIN_OB
        ELSE
          RSQH_AT_OB=R(1)*R(1)*HONJ_OUTBC*XM(1)
        END IF
!
	JNU(1:ND)=XM(1:ND)
	DO I=1,ND-1
	  RSQ_HNU(I)=(HU(I)*XM(I+1)-HL(I)*XM(I))+HS(I)*RSQ_HNUM1(I)+HT(I)*RSQ_HNU_OLDT(I)
	END DO
!
! 
!
! The J and H components of the radiation field have been found. We can thus
! begin the variation computation.
!
! Compute d{non-radiation field}/dchi matrix.
!
	CALL TUNE(1,'MOM_EDD')
	CALL EDD_J_HUB_VAR_V2(KI,RHS_dHdCHI,WORKMAT,
	1                SOURCE,CHI,ESEC,THETA,DTAU,R,
	1                F,Q,HU,HL,HS,HT,RSQ_DTAUONQ,
	1                W,WPREV,PSI,PSIPREV,DJDt,DJDT_OLDT,
	1                JNU,JNUM1,JNU_OLDT,JNU_MOD,
	1                RSQ_HNUM1,RSQ_HNU_OLDT,
	1                dRHSdCHI_IB,dRHSdCHI_OB,
	1                JMIN_IB,KMIN_IB,JPLUS_IB,KPLUS_IB,
	1                JMIN_OB,KMIN_OB,JPLUS_OB,KPLUS_OB,
	1                INNER_BND_METH,OUTER_BND_METH,
	1                ND,NM_KI)
	CALL TUNE(2,'MOM_EDD')
!
! Evaluate the intensity variations.
!                                  TX=dJ/d(chi,eta,....)
!                          and     TVX=dRSQH/d(chi,eta,...)
!
! WORKMAT is dimension (ND,ND) is is used to temporarily save TX( , ,K) for
! each K. We don't need to pass HT, since JNU_OLDT is fixed --- it does
! not vary in the linearization.
!
	CALL TUNE(1,'UP_TX')
	DUMMY_VEC=0.0D0
	CALL UP_TX_TVX(TX,TVX,KI,TA,TB,TC,PSIPREV,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       DUMMY_VEC,DUMMY_VEC,DUMMY_VEC,DUMMY_VEC,
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
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  TX_DIF_d_T(1)=PSIPREV(1)*TX_DIF_d_T(1)
	  TX_DIF_d_dTdR(1)=PSIPREV(1)*TX_DIF_d_dTdR(1)
	  DO I=2,ND-1
	    TX_DIF_d_T(I)= PSIPREV(I)*TX_DIF_d_T(I)
	1             + VB(I)*TVX_DIF_d_T(I-1) + VC(I)*TVX_DIF_d_T(I)
	    TX_DIF_d_dTdR(I)= PSIPREV(I)*TX_DIF_d_dTdR(I)
	1             + VB(I)*TVX_DIF_d_dTdR(I-1) + VC(I)*TVX_DIF_d_dTdR(I)
	  END DO
	  T1=R(ND)*R(ND)/3.0D0/CHI(ND)
	  TX_DIF_d_T(ND)=T1*dDBBdT*(1.0D0+RECIP_CDELTAT/CHI(ND)) +
	1          T1*GAM_INB*(dDBBdT-dDBBdT_PREV*CHI(ND)/CHI_AT_INB_PREV)
	  T1=T1/dTdR
	  TX_DIF_d_dTdR(ND)=T1*DBB*(1.0D0+RECIP_CDELTAT/CHI(ND)) +
	1          T1*GAM_INB*(DBB-DBB_PREV*CHI(ND)/CHI_AT_INB_PREV)
!
! Solve for the radiation field along ray for this frequency.
!
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_T,ND,1)
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_dTdR,ND,1)
!
	  DO I=1,ND-1
	    TVX_DIF_d_T(I)=HU(I)*TX_DIF_d_T(I+1) - HL(I)*TX_DIF_d_T(I) +
	1        HS(I)*TVX_DIF_d_T(I)
	    TVX_DIF_d_dTdR(I)=HU(I)*TX_DIF_d_dTdR(I+1) -
	1     HL(I)*TX_DIF_d_dTdR(I) +
	1     HS(I)*TVX_DIF_d_dTdR(I)
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
	DBB_PREV=DBB
	dDBBdT_PREV=dDBBdT
	CHI_AT_INB_PREV=CHI(ND)
	HONJ_OUTBC_PREV=HONJ_OUTBC
	RSQH_AT_IB_PREV=RSQH_AT_IB
!
	RETURN
	END
