	MODULE MOD_JREL_V7
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER ND_SAV
	INTEGER MOM_ERR_CNT
	INTEGER, PARAMETER :: N_ERR_MAX=1000
	REAL*8 MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
! These vectors must be supplied by the calling routine. We use these
! vectors to contain the data on a finer grid.
!
	REAL*8, ALLOCATABLE ::  ETA(:)
	REAL*8, ALLOCATABLE ::  CHI(:)
	REAL*8, ALLOCATABLE ::  ESEC(:)
	REAL*8, ALLOCATABLE ::  V(:)			!in km/s
	REAL*8, ALLOCATABLE ::  SIGMA(:)		!dlnV/dlnR
	REAL*8, ALLOCATABLE ::  R(:)			!in units of 10^10 cm
!
! These values are computed, and returned.
!
! R_PNT(K) defines the interpolation for the variable at depth K.
! ?_INDX are used to indicate the location of J and H on the small grid
! in the larger array.
!
        INTEGER, ALLOCATABLE :: R_PNT(:)
        INTEGER, ALLOCATABLE :: J_INDX(:)
        INTEGER, ALLOCATABLE :: H_INDX(:)
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
	REAL*8, ALLOCATABLE :: GAM_RSQJNU(:)
	REAL*8, ALLOCATABLE :: GAM_RSQHNU(:)
!
        REAL*8, ALLOCATABLE :: GAM_RSQJNU_PREV(:)
        REAL*8, ALLOCATABLE :: GAM_RSQHNU_PREV(:)
!
        REAL*8, ALLOCATABLE :: GAM_RSQJNU_SAVE(:)
        REAL*8, ALLOCATABLE :: GAM_RSQHNU_SAVE(:)
	REAL*8, ALLOCATABLE :: GAM_RSQJOLD(:)
!
	REAL*8, ALLOCATABLE :: TA(:),TB(:),TC(:)
	REAL*8, ALLOCATABLE :: CHI_H(:),CHI_J(:)
	REAL*8, ALLOCATABLE :: DTAU_H(:),DTAU_J(:),DTAUONQ(:)
	REAL*8, ALLOCATABLE :: Q(:),XM(:),SOURCE(:)
	REAL*8, ALLOCATABLE :: VB(:),VC(:),COH_VEC(:)
	REAL*8, ALLOCATABLE :: HU(:),HL(:),HS(:),HD(:)
	REAL*8, ALLOCATABLE :: P_H(:),P_J(:)
	REAL*8, ALLOCATABLE :: VdJdR_TERM(:),VdHdR_TERM(:)
	REAL*8, ALLOCATABLE :: DELTA(:),DELTAH(:),W(:),WPREV(:)
	REAL*8, ALLOCATABLE :: PSI(:),PSIPREV(:)
	REAL*8, ALLOCATABLE :: EPS(:),EPS_PREV(:)
!
	REAL*8 FREQ_SAVE
	REAL*8 IN_NBC_SAVE
	REAL*8 IN_NBC_PREV
!
	END MODULE MOD_JREL_V7
!
!
! Subroutine to allocate the required vectors.
!
	SUBROUTINE ALLOC_MOD_JREL_V7()
	USE MOD_JREL_V7
	IMPLICIT NONE
!
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
	  ALLOCATE (ETA(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHI(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ESEC(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (V(ND),STAT=IOS)			!in km/s
	  IF(IOS .EQ. 0)ALLOCATE (SIGMA(ND),STAT=IOS)			!dlnV/dlnR
	  IF(IOS .EQ. 0)ALLOCATE (R(ND),STAT=IOS)			!in units of 10^10 cm
!
	  IF(IOS .EQ. 0)ALLOCATE (R_PNT(ND),STAT=IOS)			!
	  IF(IOS .EQ. 0)ALLOCATE (J_INDX(ND),STAT=IOS)			!
	  IF(IOS .EQ. 0)ALLOCATE (H_INDX(ND),STAT=IOS)			!
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(0) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
	  IF(IOS .EQ. 0)ALLOCATE (AV_SIGMA(ND),STAT=IOS)
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
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(1) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
!
	  IF(IOS .EQ. 0)ALLOCATE (GAM_RSQJNU(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQJNU_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU_PREV(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQJNU_SAVE(ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (GAM_RSQHNU_SAVE(ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(2) in ALLOC_MOD_J_REL'
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
	  IF(IOS .EQ. 0)ALLOCATE (GAM_RSQJOLD(ND),STAT=IOS)
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
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Unable to allocate memory(3) in ALLOC_MOD_J_REL'
	    WRITE(LUER,*)'STAT=',IOS
	    STOP
	  END IF
	END IF
!
        GAM_RSQJNU_SAVE(:)=0.0D0
        GAM_RSQHNU_SAVE(:)=0.0D0
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
!	       NMID_ON_J(I) = GAM RSQ_N(I)/( GAM RSQ_J(I)+ GAM RSQ_J(I+1))
!	          KMID_ON_J(I) = GAM RSQ_K(I)/( GAM RSQ_J(I)+ GAM RSQ_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!	RSQ_K=0.25*K(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and NMID_ON_J is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
!                           NMID_ON_J=0 at all depths.
!     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) NMID_ON_J is defined at all
!                           depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or NMID_ON_J is non-zero,
!                           and is the value to be used in MOM_J_CMF
!
	SUBROUTINE MOM_JREL_V7(ETA_SM,CHI_SM,ESEC_SM,V_SM,SIGMA_SM,R_SM,
	1               JNU_SM,RSQHNU_SM,HFLUX_AT_IB,HFLUX_AT_OB,
	1               VDOP_VEC,VDOP_FRAC,
	1               FREQ,dLOG_NU,DBB,IB_STAB_FACTOR,
	1               INNER_BND_METH,OUTER_BND_METH,METHOD,COHERENT,N_TYPE,
	1               INCL_ADVEC_TERMS,INCL_REL_TERMS,INIT,ND_SM)
	USE MOD_JREL_V7
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
! Alteerd: 15-Feb-2014 : Changed to V7. Added IB_STAB_FACTOR to call.
! Altered: 28-Jan-2012 : Minor bug fix. dLOG_NU was being used with HOLLOW option when INIT was true.
! Altered: 25-Aug-2010 : Bug fix with IN_NBC_SAVE/ IN_HBC_SAVE. Both values
!                          incorrectly set because of type'o.
! Altered: 17-Dec-2009 : Call changed: changed to V6
!                        INNER_BND_METH and OUTER_BND_METH inserted.
! Altered: 15-Nov-2009 : Bug fix: LOG(CHI_SM(1:ND))--> LOG(CHI_SM(1:ND_SM))
!                          Also output warning when ND is different.
! Created: 31-Dec-2004 : Based on MOM_J_REL_V1
!                        Removed all *PREV quantities from call and placed 
!                          in module.
!
	INTEGER ND_SM
!
	REAL*8 ETA_SM(ND_SM)
	REAL*8 CHI_SM(ND_SM)
	REAL*8 ESEC_SM(ND_SM)
	REAL*8 V_SM(ND_SM)			!in km/s
	REAL*8 SIGMA_SM(ND_SM)			!dlnV/dlnR
	REAL*8 R_SM(ND_SM)			!in units of 10^10 cm
!
! These values are computed, and returned.
!
	REAL*8 JNU_SM(ND_SM)
	REAL*8 RSQHNU_SM(ND_SM)
	REAL*8 VDOP_VEC(ND_SM)
	REAL*8 VDOP_FRAC
!
! These are returned
!
	REAL*8 HFLUX_AT_OB
	REAL*8 HFLUX_AT_IB
!
! Boundary conditions: Must be supplied.
!
	REAL*8 DBB
	REAL*8 IB_STAB_FACTOR
	CHARACTER(LEN=*) INNER_BND_METH
	CHARACTER(LEN=*) OUTER_BND_METH
!
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
	CHARACTER(LEN=*) N_TYPE
!
! Local variables.
!
	REAL*8 DAMP_FAC 
	REAL*8 T1,T2
	REAL*8 DTAU
	REAL*8 RSQ_JP,RSQ_HP,RSQ_NP
	REAL*8 FMIN,FPLUS,HMIN,NMIN
	REAL*8 MAX_ER
	REAL*8 DELTA_R
	REAL*8 VDOP_MIN
	REAL*8 TA_SAV,TB_SAV,XM_SAV
!
	INTEGER COUNT
	INTEGER IFAIL
	INTEGER I,J,K
	INTEGER IT1
	INTEGER LUER,ERROR_LU
	INTEGER, PARAMETER :: IONE=1
!
	LOGICAL ACCURATE
!
! 
!
	LUER=ERROR_LU()
	IF(INIT)THEN
!
! Determin # of points in the new grid. This will allow us to allocate
! the required memory. We always insert an EVEN number of points. This guarentees that
! H_SM (defined at the midpoints of the pass grid) has an exact correspondence
! with H defined on the extended gid.
!
	  K=1
	  VDOP_MIN=VDOP_FRAC*MINVAL(VDOP_VEC(1:ND_SM))
	  DO I=1,ND_SM-1
	    IT1=INT( (V_SM(I)-V_SM(I+1))/VDOP_MIN )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
!	    IF(I .EQ. ND_SM-1)IT1=20
	    IF(IT1 .GT. 0)K=K+IT1
	    K=K+1
	  END DO
	  ND=K
	  IF(ND .NE. ND_SM)THEN
	    WRITE(LUER,*)'Warning in MOM_JREL_V7: we have adjusted ND for more accuracy'
	    WRITE(LUER,*)'ND(adjusted)=',ND,'ND(original)=',ND_SM
	    WRITE(LUER,*)'This could effect convergence with CMFGEN models'
	  END IF
!
	  IF(N_TYPE .NE. 'G_ONLY' .AND. N_TYPE .NE. 'N_ON_J' .AND. ND .NE. ND_SM)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in MOM_J_CMF_V7'
	    WRITE(LUER,*)'Cannot use N_TYPE=''MIXED'' when inserting extra points'
	    STOP
	  END IF
!
! Allocate storage
!
	  CALL ALLOC_MOD_JREL_V7( ) 
!
! Define the revise R grid.
!
	  K=1
	  R(1)=R_SM(1)
	  R_PNT(1)=1
	  DO I=1,ND_SM-1
            IT1=INT( (V_SM(I)-V_SM(I+1))/VDOP_MIN )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
!	    IF(I .EQ. ND_SM-1)IT1=20
            IF(IT1 .GT. 0)THEN
              DELTA_R=(R_SM(I+1)-R_SM(I))/(IT1+1)
                DO J=1,IT1
                  K=K+1
                  R(K)=R(K-1)+DELTA_R
                  R_PNT(K)=I
                END DO
            END IF
            K=K+1
            R(K)=R_SM(I+1)
            R_PNT(K)=I
	  END DO
	  IF(ND .NE. K)THEN
	    WRITE(LUER,*)'Error in MOM_JREL_V7'
	    WRITE(LUER,*)'Inconsistent grid sizes: ND & K',ND,K
	    STOP
	  END IF
!
	  J_INDX(1:ND_SM)=0; H_INDX(1:ND_SM)=0
	  K=1
	  DO I=1,ND_SM
	    DO WHILE(J_INDX(I) .EQ. 0)
	      IF(R_SM(I) .LE. R(K) .AND. R_SM(I) .GE. R(K+1))THEN
	        IF( (R(K)-R_SM(I)) .LT. (R_SM(I)-R(K+1)) )THEN
	          J_INDX(I)=K
	        ELSE
	          J_INDX(I)=K+1
	        END IF
	      ELSE
	        K=K+1
	      END IF
	    END DO
	  END DO
!
	  K=1
	  DO I=1,ND_SM-1
	    T1=0.5D0*(R_SM(I)+R_SM(I+1))
	    DO WHILE(H_INDX(I) .EQ. 0)
	      IF(T1 .LT. R(K) .AND. T1 .GT. R(K+1))THEN
	        H_INDX(I)=K
	      ELSE
	        K=K+1
	      END IF
	    END DO
	  END DO
!
! Interpolate V & SIGMA.
!
	  CALL MON_INTERP(V,ND,IONE,R,ND,V_SM,ND_SM,R_SM,ND_SM)
	  CALL MON_INTERP(SIGMA,ND,IONE,R,ND,SIGMA_SM,ND_SM,R_SM,ND_SM)
!
! Zero all vectors.
!
	  DELTAH(:)=0.0D0
	  DELTA(:)=0.0D0
	  W(:)=0.0D0
	  WPREV(:)=0.0D0
	  PSI(:)=0.0D0
	  PSIPREV(:)=0.0D0
	  GAM_RSQJNU_PREV(:)=0.0D0
	  GAM_RSQHNU_PREV(:)=0.0D0
	  EPS(:)=0.0D0
	  EPS_PREV(:)=0.0D0
!
	  H_ON_J_PREV(:)=0.0D0
	  K_ON_J_PREV(:)=0.0D0
	  N_ON_J_PREV(:)=0.0D0
	  NMID_ON_HMID_PREV(:)=0.0D0
	  NMID_ON_J_PREV(:)=0.0D0
	  KMID_ON_J_PREV(:)=0.0D0
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	  
	  BETA_FREQ(1:ND)=V(1:ND)*3.33564D-06                  !/2.99794D+05
	  IF(INCL_ADVEC_TERMS .OR. INCL_REL_TERMS)THEN
	    BETA(1:ND)=BETA_FREQ(1:ND)
	    GAM_REL_SQ(1:ND)=1.0D0/(1.0D0-BETA(1:ND)*BETA(1:ND))
	    GAM_REL(1:ND)=SQRT(GAM_REL_SQ(1:ND))
	  ELSE
	    GAM_REL(1:ND)=1.0D0
	    GAM_REL_SQ(1:ND)=1.0D0
	    GAM_REL_STORE(1:ND_STORE)=1.0D0
	    BETA(1:ND)=0.0D0
	  END IF
!
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
	END IF
!
	IF(ND .GT. ND_SM)THEN
!
	  TA(1:ND_SM)=LOG(ETA_SM(1:ND_SM))
	  CALL MON_INTERP(ETA,ND,IONE,R,ND,TA,ND_SM,R_SM,ND_SM)
	  ETA=EXP(ETA)
!
	  TA(1:ND_SM)=LOG(CHI_SM(1:ND_SM))
	  CALL MON_INTERP(CHI,ND,IONE,R,ND,TA,ND_SM,R_SM,ND_SM)
	  CHI=EXP(CHI)
!
	  TA(1:ND_SM)=LOG(ESEC_SM(1:ND_SM))
	  CALL MON_INTERP(ESEC,ND,IONE,R,ND,TA,ND_SM,R_SM,ND_SM)
	  ESEC=EXP(ESEC)
!
	ELSE
	  ESEC(1:ND)=ESEC_SM(1:ND)
	  CHI(1:ND)=CHI_SM(1:ND)
	  ETA(1:ND)=ETA_SM(1:ND)
	END IF
!
! Get Eddington values on the revised grid.
!
	CALL GET_MOMS_REL(R, V, FREQ, N_TYPE, ND)
!
! These are set to zero to insure all velocity terms are neglected.
!
	IF(INIT)THEN
	  H_ON_J(1:ND)=0.0D0
	  N_ON_J(1:ND)=0.0D0
	  KMID_ON_J(1:ND)=0.0D0
          NMID_ON_J(1:ND)=0.0D0
	  dlnGRSQJdlnR(1:ND)=0.0D0
	END IF
	DO I=1,ND
	  IF(NMID_ON_HMID(I) .GT. 1.0D0)NMID_ON_HMID(I)=1.0D0
	END DO
!
! If new frequency, we need to update PREV vectors which refer to the previous
! frequency. We can't update these on exit, since we iterate for each freqyency
! in this routine.
!
	IF(FREQ .NE. FREQ_SAVE .AND. .NOT. INIT)THEN
          GAM_RSQJNU_PREV(:)=GAM_RSQJNU_SAVE(:)
          GAM_RSQHNU_PREV(:)=GAM_RSQHNU_SAVE(:)
          K_ON_J_PREV(:)=K_ON_J_SAVE(:)
          NMID_ON_HMID_PREV(:)=NMID_ON_HMID_SAVE(:)
          H_ON_J_PREV(:)=H_ON_J_SAVE(:)
          N_ON_J_PREV(:)=N_ON_J_SAVE(:)
          NMID_ON_J_PREV(:)=NMID_ON_J_SAVE(:)
          KMID_ON_J_PREV(:)=KMID_ON_J_SAVE(:)
	  HBC_PREV=HBC_SAVE
	  NBC_PREV=NBC_SAVE
	  IN_HBC_PREV=IN_HBC_SAVE
	  IN_NBC_PREV=IN_NBC_SAVE
	END IF
!
! 
!
! Zero relevant vectors and matrices.
!
	GAM_RSQJNU(:)=0.0D0
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
! Compute optical depth scales.
!
	TA(1:ND)=CHI_H(1:ND)*Q(1:ND)
	DO I=1,ND
	  IF(TA(I) .LT. 0.0D0)TA(I)=0.1D0*ABS(TA(I))
	  IF(TA(I) .EQ. 0.0D0)TA(I)=0.1D0*ABS(TA(I-1))
	END DO
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
	    W(I)=DELTAH(I)*(1.0D0+CON_dNdNUH(I)*NMID_ON_HMID(I))
	    WPREV(I)=DELTAH(I)*(1.0D0+CON_dNdNUH(I)*NMID_ON_HMID_PREV(I))
	    EPS(I)=DELTAH(I)*(CON_dNdNUH(I)*NMID_ON_J(I)+
	1            CON_dKdNUH(I)*KMID_ON_J(I))/(P_H(I)+W(I))
	    EPS_PREV(I)=DELTAH(I)*(CON_dNdNUH(I)*NMID_ON_J_PREV(I)+
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
	  PSI(1)=DELTA(1)*( HBC-NBC+(NBC+BETA(1)*K_ON_J(1))*
	1             GAM_REL_SQ(1)*(SIGMA(1)+1.0D0) )
	  PSIPREV(1)=DELTA(1)*( HBC_PREV-NBC_PREV+(NBC_PREV+
	1             BETA(1)*K_ON_J_PREV(1))*
	1             GAM_REL_SQ(1)*(SIGMA(1)+1.0D0) )
	END IF
!
	DO I=2,ND
	  DTAUONQ(I)=0.5D0*(DTAU_J(I)+DTAU_J(I-1))/Q(I)
	  PSI(I)=DTAUONQ(I)*DELTA(I)*
	1            (1.0D0+CON_dKdNU(I)*K_ON_J(I)+CON_dHdNU(I)*H_ON_J(I) )
	  PSIPREV(I)=DTAUONQ(I)*DELTA(I)*
	1            (1.0D0+CON_dKdNU(I)*K_ON_J_PREV(I)+
	1                       CON_dHdNU(I)*H_ON_J_PREV(I) )
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=(K_ON_J(I+1)+VdHdR_TERM(I+1))*Q(I+1)/(P_H(I)+W(I))/DTAU_H(I)
	  HL(I)=(K_ON_J(I)+VdHdR_TERM(I))*Q(I)/(P_H(I)+W(I))/DTAU_H(I)
	  HS(I)=WPREV(I)/(P_H(I)+W(I))
	END DO
!
! 
!
! As we don't know dlnGRSQJdlnR and dHdlnR we need to iterate.
!
	ACCURATE=.FALSE.
	COUNT=0
	GAM_RSQJOLD(1:ND)=0.0D0
!
	DO WHILE(.NOT. ACCURATE)
!
	  COUNT=COUNT+1
!
	  IF(INCL_ADVEC_TERMS)THEN
	    VdJdR_TERM(1:ND)=CON_DELTA(1:ND)*dlnGRSQJdlnR(1:ND)/CHI_J(1:ND)
	    P_J(1:ND)=1.0D0+VdJdR_TERM(1:ND)
	  ELSE
	    VdJdR_TERM(1:ND)=0.0D0
	    P_J(1:ND)=1.0D0
	  END IF
!
! Compute the TRIDIAGONAL operators, and the RHS source vector. These
! vectors either depend  on dlnGRSQJdlnR etc, or are corrupted in the
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
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*GAM_RSQHNU_PREV(I-1) + VC(I)*GAM_RSQHNU_PREV(I)
	1          + PSIPREV(I)*GAM_RSQJNU_PREV(I)
	1          - EPS_PREV(I-1)*(GAM_RSQJNU_PREV(I-1)+GAM_RSQJNU_PREV(I))
	1          + EPS_PREV(I)*(GAM_RSQJNU_PREV(I)+GAM_RSQJNU_PREV(I+1))
	  END DO
!
! Evaluate TA,TB,TC for boundary conditions
!
	  IF(OUTER_BND_METH .EQ. 'HONJ')THEN
	    TC(1)=-( K_ON_J(2)+VdHdR_TERM(2) )*Q(2)/DTAU_H(1)
	    TB(1)= ( K_ON_J(1)+VdHdR_TERM(1) )*Q(1)/DTAU_H(1) + PSI(1) + HBC*P_H(1)
	    XM(1)= PSIPREV(1)*GAM_RSQJNU_PREV(1)
	    TA(1)=0.0D0
	    VB(1)=0.0D0
	    VC(1)=0.0D0
	  ELSE
	    WRITE(LUER,*)'Error in MOM_JREL_V7: Only outer boundary method implemented is HONJ'
	    WRITE(LUER,*)'Passed outer boundary method is',TRIM(OUTER_BND_METH)
	  END IF
!
	  IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	    TA(ND)=-Q(ND-1)*(K_ON_J(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	    TB(ND)=(K_ON_J(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)
	    XM(ND)=GAM_REL(ND)*DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
!
	  ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	    TA(ND)=-Q(ND-1)*(K_ON_J(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	    TB(ND)=(K_ON_J(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)
	    XM(ND)=0.0D0
!
! Done to stablize solution: 8-Feb-2014
!
	    XM(ND)=XM(ND)+IB_STAB_FACTOR*TB(ND)*GAM_REL(ND)*R(ND)*R(ND)*(JPLUS_IB+JMIN_IB)
	    TB(ND)=(1.0D0+IB_STAB_FACTOR)*TB(ND)
!
! This first order inner boundary conditions assumes a hollow core, and that JPLUS_IB, HPLUS_IB,
! KPLUS_IB, and HMIN_IB/JMIN_IB & KMIN_IB/JMIN_IB are known. These later quantities must be computed 
! by the ray-ray solution, and are supplied to this routine by the module MOD_RAY_MOM_STORE
!
! This boundary condtion is only stable when JPLUS_IB etc are detemined by the radiation field at an
! earlier freqiency, as is the case for an expanding hollow core.

	  ELSE IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	    RSQ_JP=GAM_REL(ND)*R(ND)*R(ND)*JPLUS_IB
	    RSQ_HP=GAM_REL(ND)*R(ND)*R(ND)*HPLUS_IB
	    RSQ_NP=GAM_REL(ND)*R(ND)*R(ND)*NPLUS_IB
	    HMIN=HMIN_IB/JMIN_IB
	    NMIN=NMIN_IB/JMIN_IB
	    FPLUS=KPLUS_IB/JPLUS_IB
	    FMIN=KMIN_IB/JMIN_IB
	    DTAU=0.5D0*(R(ND-1)-R(ND))*(CHI(ND)+CHI(ND-1))
	    IF(INIT)THEN
	      TA(ND)=-K_ON_J(ND-1)/DTAU
	      TB(ND)=FMIN/DTAU + (1.0D0-FMIN)/R(ND)/CHI(ND) + HMIN
	      XM(ND)=RSQ_HP-FPLUS*RSQ_JP/DTAU-(1.0D0-FPLUS)*RSQ_JP/R(ND)/CHI(ND)
	      XM(ND-1)=XM(ND-1)-RSQ_JP*TC(ND-1)
	    ELSE
	      T1=CON_DELTA(ND)/dLOG_NU/CHI(ND)
	      TA(ND)=-K_ON_J(ND-1)/DTAU
	      TB(ND)=FMIN/DTAU + (1.0D0-FMIN)/R(ND)/CHI(ND) + HMIN*(1.0D0+T1)
	      XM(ND)=RSQ_HP-FPLUS*RSQ_JP/DTAU- (1.0D0-FPLUS)*RSQ_JP/R(ND)/CHI(ND)
	1                    +T1*(RSQ_HP-IN_HBC_PREV)
	      XM(ND-1)=XM(ND-1)-RSQ_JP*TC(ND-1)
!
! Correction to improve stability.
!
	      XM(ND)=XM(ND)+IB_STAB_FACTOR*TB(ND)*GAM_REL(ND)*R(ND)*R(ND)*JMIN_IB
	      TB(ND)=(1.0D0+IB_STAB_FACTOR)*TB(ND)
	    END IF
!
	  ELSE IF(INNER_BND_METH(1:3) .EQ. 'OLD')THEN
	    TA(ND)=-Q(ND-1)*(K_ON_J(ND-1)+VdHdR_TERM(ND-1))/DTAU_H(ND-1)
	    TB(ND)=(K_ON_J(ND)+VdHdR_TERM(ND))/DTAU_H(ND-1)
	    XM(ND)=GAM_REL(ND)*R(ND)*R(ND)*HNU_AT_IB
	    IF(.NOT. INIT)THEN
	      T1=BETA(ND)*BETA(ND)*GAM_REL_SQ(ND)*(SIGMA(ND)+1.0D0)/CHI_H(ND)/R(ND)/DLOG_NU
	      XM(ND)=XM(ND)-T1*K_ON_J_PREV(ND)*GAM_RSQJNU_PREV(ND)
	      TB(ND)=TB(ND)-T1*K_ON_J(ND)
	      T1=GAM_REL_SQ(ND)*(SIGMA(ND)+1.0D0)-1.0D0
	      XM(ND)=XM(ND)+GAM_REL(ND)*R(ND)*BETA(ND)*( (HNU_AT_IB-HNU_AT_IB_PREV) +
	1                   T1*(NNU_AT_IB-NNU_AT_IB_PREV) )/CHI_H(ND)/DLOG_NU
	    END IF
	  END IF
	  TC(ND)=0.0D0
	  VB(ND)=0.0D0
	  VC(ND)=0.0D0
	  PSIPREV(ND)=0.0D0
!
! Solve for the radiation field along ray for this frequency.
!
	  TA_SAV=TA(ND);TB_SAV=TB(ND); XM_SAV=XM(ND)
	  CALL THOMAS(TA,TB,TC,XM,ND,1)
!
!	  IF(FREQ .EQ. FREQ_SAVE)BACKSPACE(UNIT=178)
!	  WRITE(178,'(ES15.7,9ES14.6)')FREQ,RSQ_JP,RSQ_HP,XM(ND)+RSQ_JP,XM(ND-1),XM(ND-2),RSQ_HP-HMIN*XM(ND)
!
! We use IN_HBC_SAVE to save H at the previous freuqency.
!
	IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	  IN_HBC_SAVE=RSQ_HP-HMIN*XM(ND)
	  IN_NBC_SAVE=RSQ_NP-NMIN*XM(ND)
	  XM(ND)=XM(ND)+RSQ_JP
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	  IN_HBC_SAVE=0.0D0
	  IN_NBC_SAVE=0.0D0
	ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  IN_HBC_SAVE=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	  IN_NBC_SAVE=DBB*R(ND)*R(ND)/5.0D0/CHI(ND)
	END IF
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
	      IF(.NOT. RECORDED_ERROR .AND. MOM_ERR_CNT .LT. N_ERR_MAX)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	        MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	      END IF	
	    END IF
	  END DO
	  GAM_RSQJNU(1:ND)=XM(1:ND)
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
	      MAX_ER=MAX(MAX_ER, ABS((GAM_RSQJOLD(I)-GAM_RSQJNU(I))/GAM_RSQJNU(I)))
	    END DO
	    IF(MAX_ER .LT. 1.0D-08)ACCURATE=.TRUE.
	    GAM_RSQJOLD(1:ND)=GAM_RSQJNU(1:ND)
!
! NB: DERIVCHI computes d(GAM.R^2J)/dR. We then mulitply that derivative by
! R/[GAM.R^2.J] to get dln(GAM.R^2J)/dlnR).
!
	    DAMP_FAC=0.8D0
	    IF(.NOT. ACCURATE .AND. INCL_ADVEC_TERMS)THEN
	      CALL DERIVCHI(TB,XM,R,ND,'LINMON')
	      TB(1:ND)=R(1:ND)*TB(1:ND)/XM(1:ND)
	      IF(MAX_ER .LT. 0.01D0)THEN
	        dlnGRSQJdlnR(1:ND)=DAMP_FAC*TB(1:ND)+(1.0D0-DAMP_FAC)*dlnGRSQJdlnR(1:ND)
	      ELSE
	        dlnGRSQJdlnR(1:ND)=0.1D0*TB(1:ND)+0.9D0*dlnGRSQJdlnR(1:ND)
	      END IF
	    END IF
	    IF(COUNT .EQ.  100)THEN
	      WRITE(LUER,*)'Error in MOM_J_REL_V7: excessive iteration count.'
	      WRITE(LUER,'(A,ES15.8,4X,A,ES9.2)')' FREQ= ',FREQ,'Error =',MAX_ER
	      ACCURATE=.TRUE.
	    END IF
	  END IF
!
	END DO
!
! Save variables for next frequency
!
	FREQ_SAVE=FREQ
        K_ON_J_SAVE(:)=K_ON_J(:)
        NMID_ON_HMID_SAVE(:)=NMID_ON_HMID(:)
        H_ON_J_SAVE(:)=H_ON_J(:)
        N_ON_J_SAVE(:)=N_ON_J(:)
        NMID_ON_J_SAVE(:)=NMID_ON_J(:)
        KMID_ON_J_SAVE(:)=KMID_ON_J(:)
        GAM_RSQJNU_SAVE(:)=GAM_RSQJNU(:)
        GAM_RSQHNU_SAVE(:)=GAM_RSQHNU(:)
	HBC_SAVE=HBC
	NBC_SAVE=NBC
!
! Regrid derived J and RSQH values onto small grid. We devide RSQJ by R^2 so that
! we return J.
!
        DO I=1,ND_SM
          K=J_INDX(I)
          JNU_SM(I)=GAM_RSQJNU(K)/GAM_REL(K)/R_SM(I)/R_SM(I)
        END DO
!
! Return RSQHNU for consistency with other programs.
!
        DO I=1,ND_SM-1
          K=H_INDX(I)
	  T1=0.5D0*(BETA(K)+BETA(K+1))
	  T1=SQRT(1.0D0-T1*T1)
          RSQHNU_SM(I)=T1*GAM_RSQHNU(K)
        END DO
!
	HFLUX_AT_OB=IN_HBC_SAVE/R_SM(ND)/R_SM(ND)
	IF(OUTER_BND_METH .EQ. 'HONJ')THEN
	  HFLUX_AT_OB=HBC*JNU_SM(1)
	END IF
!
	RETURN
	END
