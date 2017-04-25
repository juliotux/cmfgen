	MODULE PP_VAR_MOM_CMF_MOD_V1
!
	REAL*8, ALLOCATABLE :: F_PREV(:)
	REAL*8, ALLOCATABLE :: G_PREV(:)
	REAL*8, ALLOCATABLE :: N_ON_J_PREV(:)
	REAL*8, ALLOCATABLE :: JNU(:)
	REAL*8, ALLOCATABLE :: HNU(:)
	REAL*8, ALLOCATABLE :: JNUM1(:)
	REAL*8, ALLOCATABLE :: HNUM1(:)
!
	REAL*8, ALLOCATABLE :: KI(:,:,:)
	REAL*8, ALLOCATABLE :: WORKMAT(:,:)
	REAL*8, ALLOCATABLE :: RHS_dHdCHI(:,:)
!
	REAL*8, ALLOCATABLE :: TA(:),TB(:),TC(:),DTAU(:),MID_DTAU(:)
	REAL*8, ALLOCATABLE :: XM(:),SOURCE(:),THETA(:)
	REAL*8, ALLOCATABLE :: VB(:),VC(:),HU(:),HL(:),HS(:)
	REAL*8, ALLOCATABLE :: GAM(:),GAMH(:),W(:),WPREV(:)
	REAL*8, ALLOCATABLE :: PSI(:),PSIPREV_MOD(:),PSIPREV(:)
	REAL*8, ALLOCATABLE :: EPS(:),EPS_PREV(:)
	REAL*8, ALLOCATABLE :: TX_OLD_d_T(:),TX_OLD_d_dTdR(:)
!
	REAL*8 HBC_PREV
	REAL*8 IN_HBC_PREV
	REAL*8 NBC_PREV
	REAL*8 NBC_INCID_PREV
!
	INTEGER, PARAMETER :: NM_KI=2
!
	END MODULE PP_VAR_MOM_CMF_MOD_V1
!
! 
!
! This subroutine computes the variation of J as a function of the emissivity
! ETA and the opacity CHI for a plane-parralel atmoshere. To be used in conjunction
! with PP_MOM_CMF_V1.F
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
	SUBROUTINE PP_VAR_MOM_CMF_V1(ETA,CHI,ESEC,V,SIGMA,R,
	1              TX,TVX,TX_DIF_d_T,TX_DIF_d_dTdR,
	1              TVX_DIF_d_T,TVX_DIF_d_dTdR,F,G,N_ON_J,
	1              IN_HBC,HBC,HBC_INCID,NBC,NBC_INCID,
	1              INIT,REALLOCATE_ARRAYS,
	1              dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1              DO_THIS_TX_MATRIX,METHOD,COHERENT,ND,NM)
	USE PP_VAR_MOM_CMF_MOD_V1
	IMPLICIT NONE
!
! Created: 24-Mar-2006 : Based on VAR_MOM_J_CMF_V7.F
!
	INTEGER ND
	INTEGER NM
	REAL*8 ETA(ND)
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
	REAL*8 R(ND)
!
! Variation arrays and vectors.
!
	REAL*8 TX(ND,ND,NM)
	REAL*8 TVX(ND-1,ND,NM)
	REAL*8 TX_DIF_d_T(ND)
	REAL*8 TX_DIF_d_dTdR(ND)
	REAL*8 TVX_DIF_d_T(ND)
	REAL*8 TVX_DIF_d_dTdR(ND)
!
	LOGICAL DO_THIS_TX_MATRIX(NM)
!
! "Eddington factors"
!
	REAL*8 F(ND),G(ND),N_ON_J(ND)
	REAL*8 HBC,HBC_INCID
	REAL*8 NBC,NBC_INCID
	REAL*8 IN_HBC
!
	REAL*8 dLOG_NU
	REAL*8 dTdR
	REAL*8 DBB
	REAL*8 dDBBdT
	REAL*8 IC
	CHARACTER*6 METHOD
	LOGICAL COHERENT
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
	LOGICAL INIT
	LOGICAL REALLOCATE_ARRAYS
	LOGICAL DIF
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER I
	INTEGER IOS
	INTEGER LUER
	REAL*8 dN_INCID
	REAL*8 AV_SIGMA
!
! Deallocate arrays if requested. This might be done because we are testing
! with a new number of depth points, for example.
!
	IF(INIT .AND. ALLOCATED(JNU) .AND. REALLOCATE_ARRAYS)THEN
	  DEALLOCATE (JNU,HNU,JNUM1,HNUM1,F_PREV,G_PREV,N_ON_J_PREV)
	  DEALLOCATE (KI,WORKMAT,RHS_dHdCHI)
	  DEALLOCATE (TA,TB,TC,DTAU,MID_DTAU)
	  DEALLOCATE (XM,SOURCE,THETA)
	  DEALLOCATE (VB,VC,HU,HL,HS)
	  DEALLOCATE (GAM,GAMH,W,WPREV)
	  DEALLOCATE (PSI,PSIPREV_MOD,PSIPREV)
	  DEALLOCATE (EPS,EPS_PREV)
	  DEALLOCATE (TX_OLD_d_T,TX_OLD_d_dTdR)
	END IF
!
	IF(.NOT. ALLOCATED(JNU))THEN
	  ALLOCATE (JNU(ND),HNU(ND),JNUM1(ND),HNUM1(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (F_PREV(ND),G_PREV(ND),N_ON_J_PREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TA(ND),TB(ND),TC(ND),DTAU(ND),MID_DTAU(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (XM(ND),SOURCE(ND),THETA(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VB(ND),VC(ND),HU(ND),HL(ND),HS(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM(ND),GAMH(ND),W(ND),WPREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PSI(ND),PSIPREV_MOD(ND),PSIPREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EPS(ND),EPS_PREV(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TX_OLD_d_T(ND),TX_OLD_d_dTdR(ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (KI(ND,ND,2), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (WORKMAT(ND,ND), STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (RHS_dHdCHI(ND-1,ND), STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in PP_VAR_MOM_CMF_V1'
	    WRITE(LUER,*)'Unable to allocate dynamic memory: STAT=',IOS
	    STOP
	  END IF
	END IF
!
! 
!
! Zero relevant vectors and matrices.
!
	JNU=0.0D0; HNU=0.0D0
	IF(INIT)THEN
	  TX=0.0D0; TVX=0.0D0
	  JNUM1=0.0D0; HNUM1=0.0D0
	  GAM=0.0D0; GAMH=0.0D0
	  W=0.0D0; WPREV=0.0D0
	  PSI=0.0D0; PSIPREV=0.0D0
	  TX_DIF_d_T=0.0D0; TX_DIF_d_dTdR=0.0D0
	  EPS=0.0D0; EPS_PREV=0.0D0
	  dN_INCID=0.0D0
	END IF
!
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
	THETA=0.0D0
	IF(COHERENT)THETA=ESEC/CHI
!
! Compute dCHIdZ and then compute the optical depth scale.
! We also need to call d_DERIVCHI_dCHI to set the TRAP derivatives 
! which are use in PP_EDD_VAR_CMF_V1.
!
	CALL DERIVCHI(TB,CHI,R,ND,METHOD)
	CALL NORDTAU(DTAU,CHI,R,R,TB,ND)
	DO I=2,ND-1
	  MID_DTAU(I)=0.5D0*(DTAU(I-1)+DTAU(I))
	END DO
	CALL d_DERIVCHI_dCHI(TB,CHI,R,ND,METHOD)
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
	1         /dLOG_NU/(CHI(I)+CHI(I+1))
	    W(I)=GAMH(I)*(1.0D0+AV_SIGMA)*G(I)
	    WPREV(I)=GAMH(I)*(1.0D0+AV_SIGMA)*G_PREV(I)
	    EPS(I)=GAMH(I)*(1.0D0+AV_SIGMA)*N_ON_J(I)/(1.0D0+W(I))
	    EPS_PREV(I)=GAMH(I)*(1.0D0+AV_SIGMA)*N_ON_J_PREV(I)/(1.0D0+W(I))
	  END DO
	  DO I=1,ND
	    GAM(I)=3.33564D-06*V(I)/R(I)/CHI(I)/dLOG_NU
	  END DO
!
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	  PSI(1)=GAM(1)*NBC*(1.0D0+SIGMA(1))
	  PSIPREV(1)=GAM(1)*NBC_PREV*(1.0D0+SIGMA(1))
	  dN_INCID=GAM(1)*(1.0D0+SIGMA(1))*(NBC_INCID-NBC_INCID_PREV)
	  DO I=2,ND
	    PSI(I)=MID_DTAU(I)*GAM(I)*(1.0D0+SIGMA(I))*F(I)
	    PSIPREV(I)=MID_DTAU(I)*GAM(I)*(1.0D0+SIGMA(I))*F_PREV(I)
	  END DO
	END IF
!
! 
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=F(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=F(I)/(1.0D0+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0D0+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	DO I=2,ND-1
	  TA(I)=-HL(I-1)-EPS(I-1)
	  TC(I)=-HU(I)+EPS(I)
	  TB(I)=MID_DTAU(I)*(1.0D0-THETA(I)) + PSI(I) +HU(I-1) +HL(I)
	1             -EPS(I-1)+EPS(I)
	  VB(I)=-HS(I-1)
	  VC(I)=HS(I)
	  XM(I)=MID_DTAU(I)*SOURCE(I)
	END DO
!
! Evaluate TA,TB,TC for boundary conditions
!
	TC(1)=-F(2)/DTAU(1)
	TB(1)=F(1)/DTAU(1) + HBC + PSI(1)
	XM(1)=HBC_INCID + dN_INCID
	TA(1)=0.0D0
	VB(1)=0.0D0
	VC(1)=0.0D0
!
	TA(ND)=-F(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  TB(ND)=F(ND)/DTAU(ND-1)
	  XM(ND)=DBB/3.0D0/CHI(ND)
	ELSE
	  TB(ND)=F(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=IC*(0.25D0+0.5D0*IN_HBC)
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
	  PSIPREV_MOD(I)=(EPS_PREV(I)-EPS_PREV(I-1)) + PSIPREV(I)
	END DO
!
	XM(1)=XM(1) + PSIPREV_MOD(1)*JNUM1(1)
	DO I=2,ND-1
	  XM(I)=XM(I) + VB(I)*HNUM1(I-1) + VC(I)*HNUM1(I)
	1        + PSIPREV_MOD(I)*JNUM1(I)
	1        + (EPS_PREV(I)*JNUM1(I+1)-EPS_PREV(I-1)*JNUM1(I-1))
	END DO
	XM(ND)=XM(ND) + PSIPREV_MOD(ND)*JNUM1(ND)
!
! Solve for the radiation field along ray for this frequency, and store
! in JNU.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
	JNU=XM
!
	DO I=1,ND-1
	  HNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*HNUM1(I) +
	1              (EPS_PREV(I)*JNUM1(I)-EPS(I)*XM(I)) +
	1              (EPS_PREV(I)*JNUM1(I+1)-EPS(I)*XM(I+1))
	END DO
! 
!
! The J and H components of the radiation field have been found. We can thus
! begin the variation computation.
!
! Compute d{non-radiation field}/dchi matrix.
!
	CALL TUNE(1,'MOM_EDD')
	CALL PP_EDD_VAR_CMF_V1(KI,RHS_dHdCHI,WORKMAT,
	1                SOURCE,CHI,ESEC,THETA,DTAU,R,SIGMA,
	1                F,HU,HL,HS,MID_DTAU,
	1                W,WPREV,PSI,PSIPREV,
	1                EPS,EPS_PREV,
	1                JNU,JNUM1,HNUM1,
	1                DBB,DIF,ND,NM_KI)
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
	CALL UP_TX_TVX(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       EPS,EPS,EPS_PREV,EPS_PREV,
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
	1          + ( EPS_PREV(I)*TX_OLD_d_T(I+1)
	1               - EPS_PREV(I-1)*TX_OLD_d_T(I-1) )
	    TX_DIF_d_dTdR(I)= PSIPREV_MOD(I)*TX_DIF_d_dTdR(I)
	1             + VB(I)*TVX_DIF_d_dTdR(I-1) + VC(I)*TVX_DIF_d_dTdR(I)
	1          + ( EPS_PREV(I)*TX_OLD_d_dTdR(I+1)
	1               - EPS_PREV(I-1)*TX_OLD_d_dTdR(I-1) )
	  END DO
	  TX_DIF_d_T(ND)=dDBBdT/3.0D0/CHI(ND) +
	1                    PSIPREV(ND)*TX_DIF_d_T(ND)
	  TX_DIF_d_dTdR(ND)=DBB/dTdR/3.0D0/CHI(ND) +
	1                    PSIPREV(ND)*TX_DIF_d_dTdR(ND)
!
! Solve for the radiation field along ray for this frequency.
!
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_T,ND,1)
	  CALL SIMPTH(TA,TB,TC,TX_DIF_d_dTdR,ND,1)
!
	  DO I=1,ND-1
	    TVX_DIF_d_T(I)=HU(I)*TX_DIF_d_T(I+1) - HL(I)*TX_DIF_d_T(I) +
	1        HS(I)*TVX_DIF_d_T(I) +
	1        (EPS_PREV(I)*TX_OLD_d_T(I)-EPS(I)*TX_DIF_d_T(I)) +
	1        (EPS_PREV(I)*TX_OLD_d_T(I+1)-EPS(I)*TX_DIF_d_T(I+1))
	    TVX_DIF_d_dTdR(I)=HU(I)*TX_DIF_d_dTdR(I+1) -
	1     HL(I)*TX_DIF_d_dTdR(I) +
	1     HS(I)*TVX_DIF_d_dTdR(I) +
	1     (EPS_PREV(I)*TX_OLD_d_dTdR(I)-EPS(I)*TX_DIF_d_dTdR(I)) +
	1     (EPS_PREV(I)*TX_OLD_d_dTdR(I+1)-EPS(I)*TX_DIF_d_dTdR(I+1))
	  END DO
!
	END IF	    	    !DIF END
!
! 
!
! Save JNU,  HNU  and the Eddington factors for the next frequency. As we only
! call the variation subroutine once for each frequency, these are saved
! on each call.
!
	JNUM1=JNU; HNUM1=HNU
	F_PREV=F; G_PREV=G; N_ON_J_PREV=N_ON_J
	HBC_PREV=HBC; IN_HBC_PREV=IN_HBC; NBC_PREV=NBC; NBC_INCID_PREV=NBC_INCID
!
	RETURN
	END
