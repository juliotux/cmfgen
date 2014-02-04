	MODULE MOD_J_REL
	IMPLICIT NONE
!
	INTEGER MOM_ERR_CNT
	INTEGER, PARAMETER :: N_ERR_MAX=1000
	REAL*8 MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
! These vectors must be saved as they will be used on subsequent iterations.
!
	REAL*8, ALLOCATABLE :: AV_SIGMA(:)
	REAL*8, ALLOCATABLE :: BETA(:)
	REAL*8, ALLOCATABLE :: BETA_FREQ(:)
	REAL*8, ALLOCATABLE :: GAM_REL(:)
	REAL*8, ALLOCATABLE :: GAM_REL_SQ(:)
	REAL*8, ALLOCATABLE :: CON_DELTA(:)
	REAL*8, ALLOCATABLE :: CON_DELTAH(:)
	REAL*8, ALLOCATABLE :: CON_dKdNUH(:)
	REAL*8, ALLOCATABLE :: CON_dNdNUH(:)
	REAL*8, ALLOCATABLE :: CON_dKdNU(:)
	REAL*8, ALLOCATABLE :: CON_dHdNU(:)
	REAL*8, ALLOCATABLE :: GAM_RSQHNU(:)
!
        REAL*8, ALLOCATABLE :: FEDD_PREV(:)
        REAL*8, ALLOCATABLE :: GEDD_PREV(:)
        REAL*8, ALLOCATABLE :: H_ON_J_PREV(:)
        REAL*8, ALLOCATABLE :: N_ON_J_PREV(:)
        REAL*8, ALLOCATABLE :: RSQN_ON_RSQJ_PREV(:)
        REAL*8, ALLOCATABLE :: KMID_ON_J_PREV(:)
        REAL*8, ALLOCATABLE :: JNU_PREV(:)
        REAL*8, ALLOCATABLE :: GAM_RSQHNU_PREV(:)
!
        REAL*8, ALLOCATABLE :: FEDD_SAVE(:)
        REAL*8, ALLOCATABLE :: GEDD_SAVE(:)
        REAL*8, ALLOCATABLE :: H_ON_J_SAVE(:)
        REAL*8, ALLOCATABLE :: N_ON_J_SAVE(:)
        REAL*8, ALLOCATABLE :: RSQN_ON_RSQJ_SAVE(:)
        REAL*8, ALLOCATABLE :: KMID_ON_J_SAVE(:)
        REAL*8, ALLOCATABLE :: JNU_SAVE(:)
        REAL*8, ALLOCATABLE :: GAM_RSQHNU_SAVE(:)
!
	REAL*8, ALLOCATABLE :: TA(:),TB(:),TC(:)
	REAL*8, ALLOCATABLE :: CHI_H(:),CHI_J(:)
	REAL*8, ALLOCATABLE :: DTAU_H(:),DTAU_J(:),DTAUONQ(:)
	REAL*8, ALLOCATABLE :: Q(:),XM(:),SOURCE(:)
	REAL*8, ALLOCATABLE :: VB(:),VC(:),COH_VEC(:)
	REAL*8, ALLOCATABLE :: HU(:),HL(:),HS(:),HD(:)
	REAL*8, ALLOCATABLE :: P_H(:),P_J(:),JOLD(:)
	REAL*8, ALLOCATABLE :: VdJdR_TERM(:),VdHdR_TERM(:)
	REAL*8, ALLOCATABLE :: DELTA(:),DELTAH(:),W(:),WPREV(:)
	REAL*8, ALLOCATABLE :: PSI(:),PSIPREV(:)
	REAL*8, ALLOCATABLE :: EPS(:),EPS_PREV(:)
	REAL*8, ALLOCATABLE :: GAM_RSQJNU_PREV(:)
!
	REAL*8 HBC_PREV,HBC_SAVE
	REAL*8 NBC_PREV,NBC_SAVE
	REAL*8 IN_HBC_PREV,IN_HBC_SAVE
!
	REAL*8 FREQ_SAVE
	INTEGER ND_SAV
!
	END MODULE MOD_J_REL
!
!
! Subroutine to allocate the required vectors.
!
	SUBROUTINE ALLOC_MOD_J_REL(ND)
	USE MOD_J_REL
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER IOS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
!
	IF(.NOT. ALLOCATED(AV_SIGMA))THEN
!
	  FREQ_SAVE=0.0D0
	  ND_SAV=ND
!
	  ALLOCATE (AV_SIGMA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (BETA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (BETA_FREQ(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_REL(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_REL_SQ(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_DELTA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_DELTAH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dKdNUH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dNdNUH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dKdNU(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CON_dHdNU(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(1) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
          IF(IOS .EQ. 0)ALLOCATE (FEDD_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GEDD_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (H_ON_J_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (N_ON_J_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RSQN_ON_RSQJ_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (KMID_ON_J_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (JNU_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU_PREV(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(2) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
          IF(IOS .EQ. 0)ALLOCATE (FEDD_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GEDD_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (H_ON_J_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (N_ON_J_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (RSQN_ON_RSQJ_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (KMID_ON_J_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (JNU_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU_SAVE(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(3) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
	  ALLOCATE (TA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TB(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (TC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_H(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI_J(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAU_H(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAU_J(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DTAUONQ(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (Q(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (XM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (SOURCE(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VB(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (COH_VEC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HU(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HL(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HS(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HD(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (P_H(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (P_J(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (JOLD(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VdJdR_TERM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (VdHdR_TERM(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DELTA(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (DELTAH(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (W(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (WPREV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PSI(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PSIPREV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EPS(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (EPS_PREV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_RSQJNU_PREV(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(3) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
	END IF
!
        FEDD_SAVE(:)=0.0D0
        GEDD_SAVE(:)=0.0D0
        H_ON_J_SAVE(:)=0.0D0
        N_ON_J_SAVE(:)=0.0D0
        RSQN_ON_RSQJ_SAVE(:)=0.0D0
        KMID_ON_J_SAVE(:)=0.0D0
        JNU_SAVE(:)=0.0D0
        GAM_RSQHNU_SAVE(:)=0.0D0
	HBC_SAVE=0.0D0
	NBC_SAVE=0.0D0
	IN_HBC_SAVE=0.0D0
!
	IF(ND_SAV .EQ. ND)THEN
	  RETURN
	ELSE
	  WRITE(LUER,*)'Error in ALLOC_MOD_J_REL'
	  WRITE(LUER,*)'At present allocation size cannot be changed'
	  STOP
	END IF
!
	RETURN
	END
!
! Routine to compute the mean intensity J at a single frequency in the
! Comoving-Frame for an EXPANDING atmosphere with a monotonic velocity law.
! The computed intensity thus depends on the intensity computed for the 
! previous (bluer) frequency.
!
! FULLY RELATIVISTIC SOLUTION:
!          Ref: Mihalas, ApJ, 237, 574
!               Solution of the Comoving-Frame Equation of Transfer in
!               Spherically Symmetric Flows. VI Relativistic flows.
!               See notes for implementation.
!
! The solution involves the use of several moment ratios. These moment
! ratios should be computed using a formal ray by ray solution. Routine
! should be called in a loop to converge the "MOMENT" ratios (or
! EDDINGTON) factors.
!
! Required MOMENT ratios:
!
!                	     F = K / J    (d=1,2,..., N)
!	                H_ON_J = H/J      (d=1,2,...,N)
!	                N_ON_J = N/J      (d=1,2,...,N)
!
!	                     G =  N / H (d=1.5, 2.5, ..., N-0.5)
!	       RSQN_ON_RSQJ(I) = GAM RSQ_N(I)/( GAM RSQ_J(I)+ GAM RSQ_J(I+1))
!	          KMID_ON_J(I) = GAM RSQ_K(I)/( GAM RSQ_J(I)+ GAM RSQ_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!	RSQ_K=0.25*K(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
!                           RSQN_ON_RSQJ=0 at all depths.
!     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) RSQN_ON_RSQJ is defined at all
!                           depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or RSQN_ON_RSQJ is non-zero,
!                           and is the value to be used in MOM_J_CMF
!
	SUBROUTINE MOM_JREL_V3(ETA,CHI,ESEC,V,SIGMA,R,
	1               JNU,RSQHNU,dlnJdlnR,
	1		FEDD,GEDD,H_ON_J,N_ON_J,RSQN_ON_RSQJ,KMID_ON_J,
	1               HBC,IN_HBC,NBC,FREQ,dLOG_NU,
	1               DIF,DBB,IC,METHOD,COHERENT,
	1               INCL_ADVEC_TERMS,INCL_REL_TERMS,INIT,ND)
	USE MOD_J_REL
	IMPLICIT NONE
!
! Created: 31-Dec-2004 : Based on MOM_J_REL_V1
!                        Removed all *PREV quantities from call and placed 
!                          in module.
!
	INTEGER NC
	INTEGER NP
	INTEGER ND
!
	REAL*8 ETA(ND)
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 V(ND)			!in km/s
	REAL*8 SIGMA(ND)		!dlnV/dlnR
	REAL*8 R(ND)			!in units of 10^10 cm
!
! Moment ratio variables. All must be supplied.
!
	REAL*8 FEDD(ND)			!J/K at nodes
	REAL*8 GEDD(ND)			!N/G at midpoints
	REAL*8 H_ON_J(ND)		!H/J at nodes
	REAL*8 N_ON_J(ND)		!N/J at nodes
	REAL*8 RSQN_ON_RSQJ(ND)		!N/J at nodes
	REAL*8 KMID_ON_J(ND)		!
!
! These values are computed, and returned.
!
	REAL*8 JNU(ND)
	REAL*8 RSQHNU(ND)
	REAL*8 dlnJdlnR(ND)
!
! Boundary conditions: Must be supplied.
!
	REAL*8 HBC,NBC,IN_HBC
!
	REAL*8 DBB,IC
	REAL*8 FREQ,dLOG_NU
	CHARACTER*6 METHOD
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields). When INIT is true, we also calcualte quantities
! that will be used for other frequencies.
!
	LOGICAL INIT
!
! COHERENT indicates whether the scattering is coherent. If it is, we
! explicitly take it into account. If COHERENT is FALSE, any electron
! scattering term should be included directly in the ETA that is passed
! to the routine.
!
	LOGICAL COHERENT
	LOGICAL DIF			!Use diffusion approximation?
	LOGICAL J_AT_INB_EQ_B
!
! If INCL_REL_TERMS is true, terms of order v/c, which are usually
!         are included. These terms are usually excluded from stellar
!         wind models, but included in SN models.
!
	LOGICAL INCL_REL_TERMS
!
! If INCL_ADVEC_TERMS is true, the advection terms (dJ/dr and dH/dr)
!   are included. For WIND models these terms are usually neglected.
!   For supernova models, the ADVECTION terms are sometimes neglected.
!
	LOGICAL INCL_ADVEC_TERMS
!
! Local variables.
!
	REAL*8 DAMP_FAC 
	REAL*8 T1
	REAL*8 MAX_ER
	INTEGER COUNT
	INTEGER IFAIL
	INTEGER I,J
	INTEGER LUER,ERROR_LU
	LOGICAL ACCURATE
!
! 
!
	LUER=ERROR_LU()
	IF(INIT)THEN
!
! Allocate storage
!
	  CALL ALLOC_MOD_J_REL(ND)
!
! Zero all vectors.
!
	  DELTAH(:)=0.0D0
	  DELTA(:)=0.0D0
	  W(:)=0.0D0
	  WPREV(:)=0.0D0
	  PSI(:)=0.0D0
	  PSIPREV(:)=0.0D0
	  JNU_PREV(:)=0.0D0
	  GAM_RSQHNU_PREV(:)=0.0D0
	  EPS(:)=0.0D0
	  EPS_PREV(:)=0.0D0
!
	  H_ON_J_PREV(:)=0.0D0
	  FEDD_PREV(:)=0.0D0
	  N_ON_J_PREV(:)=0.0D0
	  GEDD_PREV(:)=0.0D0
	  RSQN_ON_RSQJ_PREV(:)=0.0D0
	  KMID_ON_J_PREV(:)=0.0D0
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	  
	  BETA_FREQ(1:ND)=V(1:ND)*3.33564D-06                  !/2.99794D+05
	  IF(INCL_ADVEC_TERMS .OR. INCL_REL_TERMS)THEN
	    BETA(1:ND)=BETA_FREQ(1:ND)
	  ELSE
	    BETA(1:ND)=0.0D0
	  END IF
!
	  GAM_REL_SQ(1:ND)=1.0D0/(1.0D0-BETA(1:ND)*BETA(1:ND))
	  GAM_REL(1:ND)=SQRT(GAM_REL_SQ(1:ND))
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
	  H_ON_J(1:ND)=0.0D0
	  N_ON_J(1:ND)=0.0D0
	  KMID_ON_J(1:ND)=0.0D0
          RSQN_ON_RSQJ(1:ND)=0.0D0
!
	  WRITE(LUER,*)'Using MOM_JREL_V3 for the solution o the radiative transfer equation'
	  WRITE(LUER,*)'INCL_ADVEC_TERMS in MOM_JREL_V3 is',INCL_ADVEC_TERMS
	  WRITE(LUER,*)'  INCL_REL_TERMS in MOM_JREL_V3 is',INCL_REL_TERMS
!
	END IF
!
! If new frequency, we need to update PREV vectors which refer to the previous
! frequency. We can't update these on exit, since we iterate for each freqyency
! in this routine.
!
	IF(FREQ .NE. FREQ_SAVE)THEN
          JNU_PREV(:)=JNU_SAVE(:)
          GAM_RSQHNU_PREV(:)=GAM_RSQHNU_SAVE(:)
          FEDD_PREV(:)=FEDD_SAVE(:)
          GEDD_PREV(:)=GEDD_SAVE(:)
          H_ON_J_PREV(:)=H_ON_J_SAVE(:)
          N_ON_J_PREV(:)=N_ON_J_SAVE(:)
          RSQN_ON_RSQJ_PREV(:)=RSQN_ON_RSQJ_SAVE(:)
          KMID_ON_J_PREV(:)=KMID_ON_J_SAVE(:)
	  HBC_PREV=HBC_SAVE
	  NBC_PREV=NBC_SAVE
	  IN_HBC_PREV=IN_HBC_SAVE
	END IF
!
! 
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
	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0D0
	  END DO
	  MOM_ERR_CNT=0
	END IF
!
	SOURCE(1:ND)=ETA(1:ND)/CHI_J(1:ND)
	IF(COHERENT)THEN
	  COH_VEC(1:ND)=ESEC(1:ND)/CHI_J(1:ND)/GAM_REL(1:ND)
	ELSE
	  COH_VEC(1:ND)=0.0D0
	END IF
!
! NB: We actually solve for r^2 J, not J.
!
! Compute the Q factors from F. The Q is used to make the J and K terms
! appearing in the first moment equation a perfect differential.
! Note that we integrate from the core outwards, and normalize Q to the
! value of unity at the core.
!
	DO I=1,ND
	  TA(ND-I+1)=( 3.0D0*FEDD(I)-1.0D0+
	1           BETA(I)*N_ON_J(I)-(SIGMA(I)+1.0D0)*VdHdR_TERM(I)+
	1      GAM_REL_SQ(I)*BETA(I)*(SIGMA(I)+1.0D0)*
	1       (BETA(I)*(1.0D0-FEDD(I)-VdHdR_TERM(I))-N_ON_J(I))
	1            )/(FEDD(I)+VdHdR_TERM(I))/R(I)
	  TB(I)=R(ND-I+1)
	END DO
	CALL INTEGRATE(TB,TA,Q,IFAIL,ND)
!
	OPEN(UNIT=147,FILE='Q_CHECK')
	  DO I=ND,1,-1
	    WRITE(147,'(I5,10ES14.4)')I,TB(ND-I+1),FEDD(I),H_ON_J(I),N_ON_J(I),TA(ND-I+1),Q(ND-I+1),
	1                 BETA(I),SIGMA(I),GAM_REL_SQ(I),VdHdR_TERM(I)
	  END DO
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
!	OPEN(UNIT=147,FILE='Q_CHECK')
	  DO I=1,ND
	    WRITE(147,'(I5,10ES14.4)')I,FEDD(I),H_ON_J(I),N_ON_J(I),TA(I),Q(I),BETA(I),SIGMA(I),VdHdR_TERM(I)
	  END DO
	CLOSE(UNIT=147)
!
! Compute optical depth scales.
!
	TA(1:ND)=CHI_H(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU_H,TA,R,R,TB,ND)
!
	TA(1:ND)=CHI_J(1:ND)*Q(1:ND)
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU_J,TA,R,R,TB,ND)
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
	    W(I)=DELTAH(I)*(1.0D0+CON_dNdNUH(I)*GEDD(I))
	    WPREV(I)=DELTAH(I)*(1.0D0+CON_dNdNUH(I)*GEDD_PREV(I))
	    EPS(I)=DELTAH(I)*(CON_dNdNUH(I)*RSQN_ON_RSQJ(I)+
	1            CON_dKdNUH(I)*KMID_ON_J(I))/(P_H(I)+W(I))
	    EPS_PREV(I)=DELTAH(I)*(CON_dNdNUH(I)*RSQN_ON_RSQJ_PREV(I)+
	1            CON_dKdNUH(I)*KMID_ON_J_PREV(I))/(P_H(I)+W(I))
	  END DO
!
	  DO I=2,ND
	    DELTA(I)=CON_DELTA(I)/CHI_J(I)/dLOG_NU
	  END DO
	  DELTA(1)=CON_DELTA(1)/CHI_H(1)/dLOG_NU
!
! PSIPREV is equivalent to the U vector of FORMSOL.
!
	  PSI(1)=DELTA(1)*( HBC-NBC+(NBC+BETA(1)*FEDD(1))*
	1             GAM_REL_SQ(1)*(SIGMA(1)+1.0D0) )
	  PSIPREV(1)=DELTA(1)*( HBC_PREV-NBC_PREV+(NBC_PREV+
	1             BETA(1)*FEDD_PREV(1))*
	1             GAM_REL_SQ(1)*(SIGMA(1)+1.0D0) )
	END IF
!
	DO I=2,ND
	  DTAUONQ(I)=0.5D0*(DTAU_J(I)+DTAU_J(I-1))/Q(I)
	  PSI(I)=DTAUONQ(I)*DELTA(I)*
	1            (1.0D0+CON_dKdNU(I)*FEDD(I)+CON_dHdNU(I)*H_ON_J(I) )
	  PSIPREV(I)=DTAUONQ(I)*DELTA(I)*
	1            (1.0D0+CON_dKdNU(I)*FEDD_PREV(I)+
	1                       CON_dHdNU(I)*H_ON_J_PREV(I) )
	END DO
!
! NB: We are initially computing GAM_REL R^2 J. We need to multiply the
!     original JNU_PREV by GAM_REL R^2, since is was divided by 
!     GAM_REL R^2 before it was stored.
!
	DO I=1,ND
	  GAM_RSQJNU_PREV(I)=GAM_REL(I)*R(I)*R(I)*JNU_PREV(I)
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=(FEDD(I+1)+VdHdR_TERM(I+1))*Q(I+1)/
	1                (P_H(I)+W(I))/DTAU_H(I)
	  HL(I)=(FEDD(I)+VdHdR_TERM(I))*Q(I)/
	1                (P_H(I)+W(I))/DTAU_H(I)
	  HS(I)=WPREV(I)/(P_H(I)+W(I))
	END DO
!
! 
!
! As we don't know dlnJdlnR and dHdlnR we need to iterate.
!
	ACCURATE=.FALSE.
	COUNT=0
	IF(INIT)dlnJdlnR(1:ND)=0.0D0
	JOLD(1:ND)=0.0D0
!
	DO WHILE(.NOT. ACCURATE)
!
	  COUNT=COUNT+1
!
	  IF(INCL_ADVEC_TERMS)THEN
	    VdJdR_TERM(1:ND)=CON_DELTA(1:ND)*dlnJdlnR(1:ND)/CHI_J(1:ND)
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
	    TA(I)=-HL(I-1)-EPS(I-1)
	    TC(I)=-HU(I)+EPS(I)
	    TB(I)=DTAUONQ(I)*(P_J(I)-COH_VEC(I)) + PSI(I) + HL(I) + 
	1             HU(I-1)-EPS(I-1)+EPS(I)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	  END DO
!
! Evaluate TA,TB,TC for boundary conditions
!
	  TC(1)=-( FEDD(2)+VdHdR_TERM(2) )*Q(2)/DTAU_H(1)
	  TB(1)= ( FEDD(1)+VdHdR_TERM(1) )*Q(1)/DTAU_H(1) + PSI(1) + HBC*P_H(1)
	  XM(1)=0.0D0
	  TA(1)=0.0D0
	  VB(1)=0.0D0
	  VC(1)=0.0D0
!
! Need to include relativistic terms.
!
	  J_AT_INB_EQ_B=.TRUE.
	  IF(J_AT_INB_EQ_B)THEN
	    TA(ND)=0.0D0; TC(ND)=0.0D0
	    TB(ND)=1.0D0
	    XM(ND)=GAM_REL(ND)*R(ND)*R(ND)*IC
	  ELSE IF(DIF)THEN
	    TA(ND)=-Q(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	    TB(ND)=(FEDD(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)
	    XM(ND)=GAM_REL(ND)*DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	  ELSE
	    TA(ND)=-Q(ND-1)*(FEDD(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	    TB(ND)=(FEDD(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)+IN_HBC*GAM_REL(ND)
	    XM(ND)=GAM_REL(ND)*R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)
	  END IF
	  TC(ND)=0.0D0
	  VB(ND)=0.0D0
	  VC(ND)=0.0D0
	  PSIPREV(ND)=0.0D0
!
	  XM(1)=XM(1) + PSIPREV(1)*GAM_RSQJNU_PREV(1)
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*GAM_RSQHNU_PREV(I-1) + VC(I)*GAM_RSQHNU_PREV(I)
	1          + PSIPREV(I)*GAM_RSQJNU_PREV(I)
	1          - EPS_PREV(I-1)*(GAM_RSQJNU_PREV(I-1)+GAM_RSQJNU_PREV(I))
	1          + EPS_PREV(I)*(GAM_RSQJNU_PREV(I)+GAM_RSQJNU_PREV(I+1))
	  END DO
	  XM(ND)=XM(ND)
!
!	IF(INIT)THEN
        IF(ABS(FREQ-49.8654D0) .LT. 0.001)THEN
	  OPEN(UNIT=173,STATUS='UNKNOWN')
          WRITE(173,*)'TA'
          WRITE(173,*)0.0D0,(TA(I)*R(I-1)*R(I-1)*GAM_REL(I-1),I=2,ND)
          WRITE(173,*)'TB'
          WRITE(173,*)(TB(I)*R(I)*R(I)*GAM_REL(I),I=1,ND)
          WRITE(173,*)'TC'
          WRITE(173,*)(TC(I)*R(I+1)*R(I+1)*GAM_REL(I+1),I=1,ND-1),0.0D0
          WRITE(173,*)'XM'
          WRITE(173,*)XM
          WRITE(173,*)'HL'
          WRITE(173,*)(HL(I)*R(I)*R(I)*GAM_REL(I),I=1,ND-1),0.0D0
          WRITE(173,*)'HU'
          WRITE(173,*)(HU(I)*R(I+1)*R(I+1)*GAM_REL(I+1),I=1,ND-1),0.0D0
          WRITE(173,*)'HS'
          WRITE(173,*)HS
          CLOSE(UNIT=173)
         END IF
!
! Solve for the radiation field along ray for this frequency.
!
	  CALL THOMAS(TA,TB,TC,XM,ND,1)
!
! Check that no negative mean intensities have been computed.
!
	  DO I=1,ND
	    IF(XM(I) .LT. 0)THEN
	      XM(I)=ABS(XM(I))/10.0D0
	      RECORDED_ERROR=.FALSE.
	      J=1
	      DO WHILE (J .LE. MOM_ERR_CNT .AND. .NOT. RECORDED_ERROR)
	        IF(MOM_ERR_ON_FREQ(J) .EQ. FREQ)RECORDED_ERROR=.TRUE.
	        J=J+1
	      END DO
	      IF(.NOT. RECORDED_ERROR .AND. 
	1                     MOM_ERR_CNT .LT. N_ERR_MAX)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	        MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	      END IF	
	    END IF
	  END DO
!
! Store J, correcting for the fact that we actually compute gamma r^2 J
!
	  DO I=1,ND
	    JNU(I)=XM(I)/R(I)/R(I)/GAM_REL(I)
	  END DO
!
	  DO I=1,ND-1
	    GAM_RSQHNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*GAM_RSQHNU_PREV(I) +
	1        ( EPS_PREV(I)*(GAM_RSQJNU_PREV(I)+GAM_RSQJNU_PREV(I+1)) -
	1          EPS(I)*(XM(I)+XM(I+1)) )
	  END DO
!
	  IF(.NOT. INCL_ADVEC_TERMS)THEN
	     ACCURATE=.TRUE.
	  ELSE
!
! Compute new derivatives, and decide whether to iterate further.
!
	    MAX_ER=0.0D0
	    DO I=1,ND
	      MAX_ER=MAX(MAX_ER, ABS((JOLD(I)-JNU(I))/JNU(I)))
	    END DO
	    IF(MAX_ER .LT. 1.0D-08)ACCURATE=.TRUE.
	    JOLD(1:ND)=JNU(1:ND)
!
	    DAMP_FAC=0.8D0
	    IF(.NOT. ACCURATE .AND. INCL_ADVEC_TERMS)THEN
!	      CALL DERIVCHI(TB,JNU,R,ND,'LINMON')
!	      TB(1:ND)=R(1:ND)*TB(1:ND)/JNU(1:ND)
	      CALL DERIVCHI(TB,XM,R,ND,'LINMON')
	      TB(1:ND)=R(1:ND)*TB(1:ND)/XM(1:ND)
	      IF(MAX_ER .LT. 0.01D0)THEN
	        dlnJdlnR(1:ND)=DAMP_FAC*TB(1:ND)+(1.0D0-DAMP_FAC)*dlnJdlnR(1:ND)
	      ELSE
	        dlnJdlnR(1:ND)=0.1D0*TB(1:ND)+0.9D0*dlnJdlnR(1:ND)
	      END IF
	    END IF
	    IF(COUNT .EQ.  100)THEN
	      WRITE(LUER,*)'WARNINNG: Error in MOM_J_REL_V2: excessive iteration count (> 100).'
	      WRITE(LUER,*)'FREQ=',FREQ
	      WRITE(LUER,*)'Maximum error =',MAX_ER
	      ACCURATE=.TRUE.
	    END IF
	  END IF
!
	END DO
!
! Save variables for next frequency
!
	FREQ_SAVE=FREQ
        FEDD_SAVE(:)=FEDD(:)
        GEDD_SAVE(:)=GEDD(:)
        H_ON_J_SAVE(:)=H_ON_J(:)
        N_ON_J_SAVE(:)=N_ON_J(:)
        RSQN_ON_RSQJ_SAVE(:)=RSQN_ON_RSQJ(:)
        KMID_ON_J_SAVE(:)=KMID_ON_J(:)
        JNU_SAVE(:)=JNU(:)
        GAM_RSQHNU_SAVE(:)=GAM_RSQHNU(:)
	HBC_SAVE=HBC
	NBC_SAVE=NBC
	IN_HBC_SAVE=IN_HBC
!
! Return RSQHNU for consistency with other programs.
!
	DO I=1,ND-1
	  T1=0.5D0*(BETA(I)+BETA(I+1))
	  T1=SQRT(1.0D0-T1*T1)
	  RSQHNU(I)=T1*GAM_RSQHNU(I)
	END DO
!
	RETURN
	END
