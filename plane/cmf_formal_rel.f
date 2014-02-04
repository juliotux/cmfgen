	MODULE EXT_REL_GRID
	  REAL*8, ALLOCATABLE :: R_EXT(:)
	  REAL*8, ALLOCATABLE :: LOG_R_EXT(:)
	  REAL*8, ALLOCATABLE :: Z_EXT(:)
	  REAL*8, ALLOCATABLE :: V_EXT(:)
	  REAL*8, ALLOCATABLE :: VDOP_VEC_EXT(:)
	  REAL*8, ALLOCATABLE :: SIGMA_EXT(:)
	  REAL*8, ALLOCATABLE :: ETA_EXT(:)
	  REAL*8, ALLOCATABLE :: CHI_EXT(:)
	  REAL*8, ALLOCATABLE :: LOG_ETA_EXT(:)
	  REAL*8, ALLOCATABLE :: LOG_CHI_EXT(:)
	  REAL*8, ALLOCATABLE :: TMP_VEC(:)
!
	  REAL*8, ALLOCATABLE :: JNU_MID(:)
	  REAL*8, ALLOCATABLE :: HNU_MID(:)
	  REAL*8, ALLOCATABLE :: KNU_MID(:)
	  REAL*8, ALLOCATABLE :: NNU_MID(:)
!
	  REAL*8, ALLOCATABLE :: MIDR(:)
	  REAL*8, ALLOCATABLE :: GAM_REL(:)
!
	  REAL*8, ALLOCATABLE :: CHI_COEF(:,:)
	  REAL*8, ALLOCATABLE :: ETA_COEF(:,:)
!
	  REAL*8, ALLOCATABLE :: CHI_RAY(:)
	  REAL*8, ALLOCATABLE :: ETA_RAY(:)
!
	  INTEGER ND_EXT
	  INTEGER ND_ADD
	  INTEGER NP_MAX
	  INTEGER IDMIN,IDMAX
	  INTEGER, PARAMETER :: ND_ADD_MAX=24
	  LOGICAL, SAVE ::  FIRST_TIME=.TRUE.
!
	END MODULE EXT_REL_GRID
!
	SUBROUTINE CMF_FORMAL_REL
	1                 (ETA,CHI,ESEC,V,SIGMA,R,P,
	1                  JNU,HNU,KNU,NNU,H_ON_J,K_ON_J,N_ON_J,
	1                  K_ON_J_MID,RSQN_MID_ON_RSQJ,N_ON_H_MID,
	1                  IN_HBC,HBC,NBC,
	1                  IPLUS,MU_AT_RMAX,HQW_AT_RMAX,
	1                  FREQ,dLOG_NU,DIF,B_PLANCK,DBB,IC,THICK_OB,
	1                  N_TYPE,VDOP_VEC,VDOP_FRAC,
	1                  METHOD,INITIALIZE,NEW_FREQ,NC,NP,ND)
        USE EXT_REL_GRID
	USE MOD_SPACE_GRID_V2
	IMPLICIT NONE
!
        INTEGER NC
        INTEGER ND
        INTEGER NP
!
	REAL*8 R(ND)
        REAL*8 V(ND)
	REAL*8 SIGMA(ND)
	REAL*8 P(NP)
	REAL*8 ETA(ND)	
	REAL*8 CHI(ND)
	REAL*8 ESEC(ND)
	REAL*8 VDOP_VEC(ND)
	REAL*8 VDOP_FRAC
	REAL*8 DELV_FRAC_FG
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!
        REAL*8 JNU(ND)
	REAL*8 HNU(ND)
	REAL*8 KNU(ND)
	REAL*8 NNU(ND)
	REAL*8 H_ON_J(ND)
	REAL*8 K_ON_J(ND)
	REAL*8 N_ON_J(ND)
!
	REAL*8 K_ON_J_MID(ND)
	REAL*8 N_ON_H_MID(ND)
	REAL*8 RSQN_MID_ON_RSQJ(ND)
!
        REAL*8 IPLUS(NP)
        REAL*8 MU_AT_RMAX(NP)
        REAL*8 HQW_AT_RMAX(NP)
        REAL*8 IN_HBC
	REAL*8 HBC
	REAL*8 NBC
!
	REAL*8 B_PLANCK
	REAL*8 DBB
	REAL*8 IC
	REAL*8 FREQ
	REAL*8 dLOG_NU
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE
	LOGICAL DIF		!Use diffusion approximation
!
! Use "Thick" boundary condition. at outer boundary. Only noted when INITIALIZE 
! is true. All subsequent frequencies will use the same boundary condition
! independent of the passed value (Until INITIALIZE is set to TRUE again).
!
	LOGICAL THICK_OB
!
! First frequency -- no frequency coupling.

	LOGICAL INITIALIZE
!
! Upon leaving this routine the radiation field along each ray is stored. This
! will provide the blue wing information necessary for the next frequency.
! This routine may, however, be used in an iterative loop. In this case the
! "blue wing" information should remain unaltered between calls.
! NEW_FREQ indicates that a new_frequency is being passed, and hence the "blue
! wing" information should be updated.
!
	LOGICAL NEW_FREQ	
!
! Local variables.
!
	INTEGER, PARAMETER :: IONE=1
	LOGICAL, PARAMETER :: LFALSE=.FALSE.
	LOGICAL, PARAMETER :: LTRUE=.TRUE.
!                
	REAL*8 DBC
	REAL*8 I_CORE
	REAL*8 T1,T2
	REAL*8 ALPHA
	REAL*8 ESEC_POW
	REAL*8 BETA
	REAL*8 VINF
	REAL*8 RMAX
	REAL*8 DEL_R_FAC
	REAL*8 NU_ON_dNU
	CHARACTER(LEN=20) BOUNDARY
!
	INTEGER NDM1
	INTEGER I,J,K,IP,ID
	INTEGER NP_LIMIT
	INTEGER NRAY
!
	REAL*8 C_KMS
	REAL*8 SPEED_OF_LIGHT
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU, SPEED_OF_LIGHT
	LOGICAL NEW_R_GRID
!
!
!
	IF(INITIALIZE)THEN
	  NU_ON_dNU=0.0D0
	ELSE
	  NU_ON_dNU=1.0D0/dLOG_NU
	END IF
!	IF(dLOG_NU .NE. 0.0D0)NU_ON_dNU=1.0D0/dLOG_NU
!
! Check to see whether we have a new R grid, or the solution options have
! change. This can only happen when INIT is TRUE.
!
	LUER=ERROR_LU()
	NEW_R_GRID=.FALSE.
	IF(INITIALIZE .AND. .NOT. FIRST_TIME)THEN
	  DO I=1,ND
	    IF(R(I) .NE. R_EXT(ND_ADD+I))THEN
	      NEW_R_GRID=.TRUE.
	      WRITE(LUER,*)'Updating RGRID in FG_J_CMF_V10'
	      EXIT
	    END IF
	  END DO
	ELSE IF(INITIALIZE)THEN
	  NEW_R_GRID=.TRUE.
	END IF
	
! Deallocate all allocated rays if we are using a diferent solution technique.
! This option will only be used when testing, since in CMFGEN we will always use
! the same atmospheric structure.
!
	IF(ALLOCATED(R_EXT) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R_EXT )
	  DEALLOCATE ( LOG_R_EXT )
	  DEALLOCATE ( V_EXT )
	  DEALLOCATE ( Z_EXT )
	  DEALLOCATE ( SIGMA_EXT )
	  DEALLOCATE ( ETA_EXT )
	  DEALLOCATE ( CHI_EXT )
	  DEALLOCATE ( LOG_ETA_EXT )
	  DEALLOCATE ( LOG_CHI_EXT )
	  DEALLOCATE ( TMP_VEC )
	  DEALLOCATE ( GAM_REL )
	  DEALLOCATE ( JNU_MID, HNU_MID, KNU_MID, NNU_MID)
	  DEALLOCATE ( MIDR )
	END IF
!
! Set up the revised grid to improve computational accuracy. Unles we are
! carrying out tests, these need only be constructed once.
!
	IF(FIRST_TIME .OR. .NOT. ALLOCATED(R_EXT) )THEN
!
          ND_ADD=0
          IF(THICK_OB)ND_ADD=ND_ADD_MAX
          ND_EXT=ND+ND_ADD
!
	  ALLOCATE ( R_EXT(ND_EXT) )
	  ALLOCATE ( LOG_R_EXT(ND_EXT) )
	  ALLOCATE ( V_EXT(ND_EXT) )
	  ALLOCATE ( VDOP_VEC_EXT(ND_EXT) )
	  ALLOCATE ( Z_EXT(ND_EXT) )
	  ALLOCATE ( SIGMA_EXT(ND_EXT) )
	  ALLOCATE ( ETA_EXT(ND_EXT) )
	  ALLOCATE ( CHI_EXT(ND_EXT) )
	  ALLOCATE ( LOG_ETA_EXT(ND_EXT) )
	  ALLOCATE ( LOG_CHI_EXT(ND_EXT) )
	  ALLOCATE ( TMP_VEC(ND_EXT) )
!
	  ALLOCATE ( JNU_MID(ND), HNU_MID(ND), KNU_MID(ND), NNU_MID(ND))
	  ALLOCATE ( GAM_REL(ND) )
	  ALLOCATE ( MIDR(ND) )
!
	  ALLOCATE ( CHI_COEF(ND_EXT,4) )
	  ALLOCATE ( ETA_COEF(ND_EXT,4) )
!
! Compute the extended R grid, excluding inserted points.
!
	  DO I=1,ND
	    R_EXT(ND_ADD+I)=R(I)
	  END DO
	  IF(THICK_OB)THEN
	    IF(V(ND) .LT. 10.0D0 .AND. R(1)/R(ND) .GE. 9.99D0)THEN
	      RMAX=10.0D0*R(1)		!Stellar wind
	    ELSE IF(V(ND) .GT. 10.0D0 .OR. V(1) .GT. 2.0D+04)THEN
	      RMAX=1.5D0*R(1)		!SN model
	    ELSE
	      RMAX=MIN(10.0D0,SQRT(R(1)/R(ND)))*R(1)
	    END IF
	    ALPHA=R(1)+(R(1)-R(2))                               
	    DEL_R_FAC=EXP( LOG(RMAX/ALPHA)/(ND_ADD-4) )
	    R_EXT(1)=RMAX
	    R_EXT(5)=RMAX/DEL_R_FAC
	    R_EXT(2)=R_EXT(1)-0.001*(R_EXT(1)-R_EXT(5))
	    R_EXT(3)=R_EXT(1)-0.1*(R_EXT(1)-R_EXT(5))
	    R_EXT(4)=R_EXT(1)-0.4*(R_EXT(1)-R_EXT(5))
	    DO I=5,ND_ADD-1
	      R_EXT(I)=R_EXT(I-1)/DEL_R_FAC
	    END DO
	    R_EXT(ND_ADD)=ALPHA
!
	  END IF
!
	  WRITE(6,'(5ES14.4)')(R_EXT(I),I=1,ND_EXT)
!
	  C_KMS=1.0D-05*SPEED_OF_LIGHT()
	  GAM_REL(1:ND)=1.0D0/SQRT(1.0D0-(V(1:ND)/C_KMS)**2)
	  DO I=1,ND-1
	    MIDR(I)=0.5D0*(R(I)+R(I+1))
	  END DO
!
!
! Compute VEXT and R_EXT. We assume a BETA velocity law at large R.
!
	  V_EXT(ND_ADD+1:ND_EXT)=V(1:ND)
	  VDOP_VEC_EXT(ND_ADD+1:ND_EXT)=VDOP_VEC(1:ND)
	  SIGMA_EXT(ND_ADD+1:ND_EXT)=SIGMA(1:ND)
	  IF(THICK_OB)THEN
	    BETA=(SIGMA(1)+1.0D0)*(R(1)/R(ND)-1.0D0)
            VINF=V(1)/(1-R(ND)/R(1))**BETA
	    DO I=1,ND_ADD
	      V_EXT(I)=VINF*(1.0D0-R_EXT(ND_EXT)/R_EXT(I))**BETA
	      SIGMA_EXT(I)=BETA/(R_EXT(I)/R_EXT(ND_EXT)-1.0D0)-1.0D0
	    END DO
	    VDOP_VEC_EXT(1:ND_ADD)=VDOP_VEC(1)
	    WRITE(LUER,'(A)')' '
	    WRITE(LUER,*)'Using thick boundary condition in CMF_FORM_REL'
	    WRITE(LUER,'(2(A,ES16.8,3X))')' R(1)=',R(1),'RMAX=',RMAX
	    WRITE(LUER,'(2(A,ES16.8,3X))')' V(1)=',V(1),'VMAX=',V_EXT(1)
	  END IF
	  LOG_R_EXT(1:ND_EXT)=LOG(R_EXT(1:ND_EXT))
!
	END IF
!
! Compute CHI_EXT, and ETA_EXT. CHI_EXT could be saved, as it doesn't change
! during the iteration procedure. ETA does however change (since it depends
! of J) and thus ETA_EXT must be re-computed on each entry.
!
! The first checks whether we may have negative line opacities due
! to stimulated emission. In such a case we simply assume an 1/r^2
! extrapolation.
!
! We also interpolate in ESEC, since ESEC (in the absence of negative 
! absorption) provides a lower bound to the opacity. NB: When CHI is much
! larger then ESEC its variation with r dominates, and it is possible to
! extrapolate CHI below ESEC.
!
	IF(ND_ADD .NE. 0)THEN
	  IDMIN=1; IDMAX=4
	  IF(CHI(IDMIN) .LE. ESEC(IDMIN) .OR. CHI(IDMAX) .LE. ESEC(IDMAX))THEN
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.)ESEC_POW=2
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(IDMAX)-ESEC(IDMAX)) / (CHI(IDMIN)-ESEC(IDMIN)) )
	1          /LOG(R(IDMIN)/R(IDMAX))
	    IF(ALPHA .LT. 2.)ALPHA=2.0
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.)ESEC_POW=2
   	    DO I=1,ND_ADD
	      T1=(CHI(IDMIN)-ESEC(IDMIN))*(R(IDMIN)/R_EXT(I))**ALPHA
	      T2=ESEC(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	      CHI_EXT(I)=T1+T2
	    END DO
	  END IF
	  DO I=ND_ADD+1,ND_EXT
	    CHI_EXT(I)=CHI(I-ND_ADD)
	  END DO
!
! We limit alpha to 3.5 to avoid excess envelope emission. If alpha were
! 3 we would get a logarithmic flux divergence as we increase the volume.
!
	  ALPHA=LOG(ETA(4)/ETA(1))/LOG(R(1)/R(4))
	  IF(ALPHA .LT. 3.5)ALPHA=3.5
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(1)*(R(1)/R_EXT(I))**ALPHA
	    IF(ETA_EXT(I) .LE. 1.0D-280)ETA_EXT(I)=1.0D-280
	  END DO
	  DO I=ND_ADD+1,ND_EXT
	    ETA_EXT(I)=ETA(I-ND_ADD)
	  END DO
	ELSE
	  CHI_EXT(1:ND)=CHI(1:ND)		!NB: In this can ND=ND_EXT
	  ETA_EXT(1:ND)=ETA(1:ND)
	END IF
!
	IF(FIRST_TIME .OR. NEW_R_GRID)THEN
	  CALL DEFINE_GRID_V2(R_EXT,V_EXT,VDOP_VEC_EXT,VDOP_FRAC,ND_EXT,R,P,ND,NC,NP)
	  J=0
!	  OPEN(UNIT=7,FILE='MU_VALUES',STATUS='UNKNOWN')
	  T1=0.0D0
	  DO IP=1,NP
	    J=MAX(J,RAY(IP)%NZ)
!	    MU_AT_RMAX(IP)=RAY(IP)%MU_P(RAY(IP)%LNK(1))
!	    HQW_AT_RMAX(IP)=2.0D0*HQW_P(1,IP)
!	    T1=T1+MU_AT_RMAX(IP)*HQW_AT_RMAX(IP)
!	    WRITE(7,'(I4,4ES16.6)')IP,MU_AT_RMAX(IP),HQW_AT_RMAX(IP),
!	1                 (MU_AT_RMAX(IP)+V(1)/C_KMS)/(1.0D0+MU_AT_RMAX(IP)*V(1)/C_KMS),
!	1                 SQRT( (R(1)-P(IP))*(R(1)+P(IP)) )/R(1)
	  END DO
!	  WRITE(7,*)'Sum of HQW is ',SUM(HQW_AT_RMAX(1:NP)),T1
!	  CLOSE(UNIT=7)
!	  IF(ALLOCATED(CHI_RAY))DEALLOCATE(CHI_RAY,ETA_RAY)
	  ALLOCATE (CHI_RAY(J),ETA_RAY(J))
	END IF
!
! This must be after the call to DEFINE_GRID so that RAY_POINTS_INSERTED is defined.
!
	IF(RAY_POINTS_INSERTED)THEN
	  LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	  LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	END IF
!
!
	JNU(:)=0.0D0
	HNU(:)=0.0D0
	KNU(:)=0.0D0
	NNU(:)=0.0D0
!
	IF(INITIALIZE)THEN
	  DO IP=1,NP
	    RAY(IP)%I_P=0.0D0; RAY(IP)%I_M=0.0D0
	    RAY(IP)%I_P_PREV=0.0D0; RAY(IP)%I_M_PREV=0.0D0
	    RAY(IP)%I_P_SAVE=0.0D0; RAY(IP)%I_M_SAVE=0.0D0
	  END DO	
	ELSE IF(NEW_FREQ)THEN
	  DO IP=1,NP
	    RAY(IP)%I_P=0.0D0; RAY(IP)%I_M=0.0D0
	    RAY(IP)%I_P_PREV=RAY(IP)%I_P_SAVE
	    RAY(IP)%I_M_PREV=RAY(IP)%I_M_SAVE
	  END DO	
	END IF
!
! If no points have been inserted on the rays, CHI and ETA are the same for 
! all rays.
!
	IF(.NOT. RAY_POINTS_INSERTED)THEN
	  CHI_RAY(1:ND_EXT)=CHI_EXT(1:ND_EXT)
	  ETA_RAY(1:ND_EXT)=ETA_EXT(1:ND_EXT)
	END IF
!
! Determine radiative transfer along each p-ray
!
	NP_LIMIT=NP-1
	IF(THICK_OB)NP_LIMIT=NP
	T1=DBB/CHI(ND)			!dB/dTAU
	BOUNDARY='diffusion'
	IF(.not. DIF)T1=0.0D0
	DO IP=1,NP_LIMIT
!
	  NRAY=RAY(IP)%NZ
	  IF(RAY_POINTS_INSERTED)THEN
            K=1
            DO I=1,RAY(IP)%NZ
100	      CONTINUE
              IF( RAY(IP)%R_RAY(I) .EQ. R_EXT(K))THEN
                CHI_RAY(I)=CHI_EXT(K)
                ETA_RAY(I)=ETA_EXT(K)
	      ELSE IF(RAY(IP)%R_RAY(I) .GT. R_EXT(K+1))THEN
                T1=LOG(RAY(IP)%R_RAY(I)/R_EXT(K))
                T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
                CHI_RAY(I)=EXP(T2)
                T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
                ETA_RAY(I)=EXP(T2)
	      ELSE 
	        K=K+1
                GOTO 100
	      END IF
            END DO
	  END IF
!
! Solve using Relativistic Formal Integral
!
	  T1=DBB/CHI(ND)
          CALL SOLVE_CMF_FORMAL_V2(CHI_RAY,ETA_RAY,IP,NU_ON_dNU,BOUNDARY,b_planck,t1,NRAY,NP,NC)
!
	  DO ID=1,MIN(ND,NP-IP+1)
	    I_P_GRID(ID)=RAY(IP)%I_P(RAY(IP)%LNK(ID))
	    I_M_GRID(ID)=RAY(IP)%I_M(RAY(IP)%LNK(ID))
	    IF( I_P_GRID(ID) .LT. 0 .OR. I_M_GRID(ID) .LT. 0.0D0)THEN
	      WRITE(145,*)'Error: invalid intensities is SOVE_CMF_REL'
	      WRITE(145,*)IP,NRAY
	      WRITE(145,'(4ES14.4)')FREQ,dLOG_NU,B_PLANCK,DBB
	      WRITE(145,'(A)')' '
	      DO I=1,NRAY
	        WRITE(145,'(4ES14.4)')CHI_RAY(I),ETA_RAY(I),RAY(IP)%I_P(I),RAY(IP)%I_M(I)
	      END DO
	      STOP
	    END IF
	  END DO
!
! Integrate over p to get J and K. NB: nz(ip) is the number of
! data points along each segment assocaited with the 2 directions.
!
           DO ID=1,MIN(ND,NP-IP+1)
             Jnu(ID)=Jnu(ID)+I_p_grid(ID)*Jqw_p(ID,ip)+I_m_grid(ID)*Jqw_m(ID,ip)
             Hnu(ID)=Hnu(ID)+I_p_grid(ID)*Hqw_p(ID,ip)+I_m_grid(ID)*Hqw_m(ID,ip)
             Knu(ID)=Knu(ID)+I_p_grid(ID)*Kqw_p(ID,ip)+I_m_grid(ID)*Kqw_m(ID,ip)
             Nnu(ID)=Nnu(ID)+I_p_grid(ID)*Nqw_p(ID,ip)+I_m_grid(ID)*Nqw_m(ID,ip)
	  END DO
!
! Store intensity at the outer boundary.
!
	  IPLUS(IP)=I_P_GRID(1)-I_M_GRID(1)
!
	END DO
!
	DO IP=1,NC+1
	END DO
!
! Save intnesity for integration at next frequency.
!
	DO IP=1,NP
	  RAY(IP)%I_P_SAVE=RAY(IP)%I_P
	  RAY(IP)%I_M_SAVE=RAY(IP)%I_M
	END DO
!
! Compute {H,K,N}/J at the nodes.
!
	H_ON_J(1:ND)=HNU(1:ND)/JNU(1:ND)
	K_ON_J(1:ND)=KNU(1:ND)/JNU(1:ND)
	N_ON_J(1:ND)=NNU(1:ND)/JNU(1:ND)
	HBC=H_ON_J(1)
	NBC=N_ON_J(1)
!
! Compute r^2. Gam. {J,H,K,N} are mid points. These are then used to compute
! Eddington-like factor at the grid mid points.
!
	NDM1=ND-1
	TMP_VEC(1:ND)=GAM_REL(1:ND)*JNU(1:ND)*R(1:ND)*R(1:ND)
	CALL MON_INTERP(JNU_MID,NDM1,IONE,MIDR,NDM1,TMP_VEC,ND,R,ND)
	TMP_VEC(1:ND)=GAM_REL(1:ND)*HNU(1:ND)*R(1:ND)*R(1:ND)
	CALL MON_INTERP(HNU_MID,NDM1,IONE,MIDR,NDM1,TMP_VEC,ND,R,ND)
	TMP_VEC(1:ND)=GAM_REL(1:ND)*KNU(1:ND)*R(1:ND)*R(1:ND)
	CALL MON_INTERP(KNU_MID,NDM1,IONE,MIDR,NDM1,TMP_VEC,ND,R,ND)
	TMP_VEC(1:ND)=GAM_REL(1:ND)*NNU(1:ND)*R(1:ND)*R(1:ND)
	CALL MON_INTERP(NNU_MID,NDM1,IONE,MIDR,NDM1,TMP_VEC,ND,R,ND)
!
! Now compute the "G" Eddington factors. Because H may be zero we have to
! be careful. Depending on N_TYPE we can specify N in terms of J, N in
! terms of H, or N in terms of J and H.
!
	TMP_VEC(1:ND)=GAM_REL(1:ND)*JNU(1:ND)*R(1:ND)*R(1:ND)
	DO I=1,NDM1
	    K_ON_J_MID(I)=KNU_MID(I)/(TMP_VEC(I)+TMP_VEC(I+1))
	END DO
	IF(N_TYPE .EQ. 'N_ON_J')THEN
	  DO I=1,NDM1
	    RSQN_MID_ON_RSQJ(I)=NNU_MID(I)/(TMP_VEC(I)+TMP_VEC(I+1))
	  END DO
	ELSE IF(N_TYPE .EQ. 'MIXED')THEN
	  RSQN_MID_ON_RSQJ(1:NDM1)=0.0D0
!	  N_ON_H_MID(1)=NNU_MID(1)/HNU_MID(1)
	  DO I=1,NDM1
	    IF(HNU_MID(I) .NE. 0)THEN
	      T1=NNU_MID(I)/HNU_MID(I)
	    ELSE
	      T1=100.0D0
	    END IF
	    IF(T1 .GT. 1.1 .OR. T1 .LT. 0.05)THEN
	      RSQN_MID_ON_RSQJ(I)=NNU_MID(I)/(TMP_VEC(I)+TMP_VEC(I+1))
	      N_ON_H_MID(I)=0.0D0
	    ELSE
	      N_ON_H_MID(I)=T1
	    END IF
	  END DO
	ELSE IF(N_TYPE .EQ. 'G_ONLY')THEN
!
! Compute G Eddington factor storing in N. We also check the validity of
! the Eddington factor in case strange values are occurring because H
! is near zero (switching sign?).
!
	  RSQN_MID_ON_RSQJ(1:NDM1)=0.0D0
	  DO I=1,NDM1
	    IF(HNU_MID(I) .NE. 0)THEN
	       N_ON_H_MID(I)=NNU_MID(I)/HNU_MID(I)
	       IF(N_ON_H_MID(I) .GT. 1.0D0)THEN
	         N_ON_H_MID(I)=1.0D0
	       ELSE IF( NNU_MID(I) .LT. 0.01) THEN
	         N_ON_H_MID(I)=0.01D0
	       END IF
	    ELSE
	      N_ON_H_MID(I)=0.0D0              !Later replaced by average.
	      WRITE(LUER,'(1X,A,1PE16.8)')'HNU zero for frequency:',FREQ
	      WRITE(LUER,'(1X,A,I4)')'Error occurred at depth:',I
	     END IF
	  END DO
	END IF
!
	FIRST_TIME=.FALSE.
	RETURN
	END
