	MODULE EXT_REL_GRID_V3
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
	END MODULE EXT_REL_GRID_V3
!
	SUBROUTINE CMF_FORMAL_REL_V3
	1           (ETA,CHI,ESEC,V,SIGMA,R,P,
	1            JNU,FEDD,RET_HNU_AT_IB,RET_HNU_AT_OB,IPLUS,
	1            FREQ,dLOG_NU,B_PLANCK,DBB,
	1            INNER_BND_METH,THICK_OB,
	1            VDOP_VEC,VDOP_FRAC,
	1            METHOD,INITIALIZE,NEW_FREQ,NC,NP,ND)
        USE EXT_REL_GRID_V3
	USE MOD_SPACE_GRID_V2
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
! Altered 14-May-2009: Altered value of IDMAX (for ETA and CHI interpolaton).
!                        IDMAX & IDMIN now stored in FG_J_CMF_MOD_V11.
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
	REAL*8 RET_HNU_AT_IB
	REAL*8 RET_HNU_AT_OB
!
! NB: J,H,K,N refer to the first 4 moments of the radiation field.
!
        REAL*8 JNU(ND)
	REAL*8 FEDD(ND)
        REAL*8 IPLUS(NP)
!
	REAL*8 B_PLANCK
	REAL*8 DBB
	REAL*8 FREQ
	REAL*8 dLOG_NU
!
	CHARACTER*6 METHOD
!
	CHARACTER(LEN=*) INNER_BND_METH
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
	REAL*8 dBdTAU
	REAL*8 ALPHA
	REAL*8 ESEC_POW
	REAL*8 BETA
	REAL*8 VINF
	REAL*8 RMAX
	REAL*8 DEL_R_FAC
	REAL*8 NU_ON_dNU
	REAL*8 MU_AT_RMAX
	REAL*8 HQW_AT_RMAX
	CHARACTER(LEN=20) BOUNDARY
!
	INTEGER NDM1
	INTEGER I,J,K,IP,ID,IOS
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
!
! Allocate data for moments which will be used to construct the Eddington 
! factors.
!
	IF(.NOT. ALLOCATED(JNU_STORE))THEN
	  ND_STORE=ND
	  ALLOCATE (JNU_STORE(ND))
	  ALLOCATE (HNU_STORE(ND))
	  ALLOCATE (KNU_STORE(ND))
	  ALLOCATE (NNU_STORE(ND))
	  ALLOCATE (R_STORE(ND))
	  ALLOCATE (GAM_REL_STORE(ND))
	  ALLOCATE (RMID_STORE(ND))
	  ALLOCATE (EXT_RMID_STORE(ND+1))
	END IF
!
	IF(INITIALIZE)THEN
	  R_STORE(1:ND)=R(1:ND)
	  GAM_REL_STORE(1:ND)=1.0D0/SQRT(1.0D0-(V(1:ND)/2.99792458D+05)**2)
	  DO I=1,ND-1
	    RMID_STORE(I)=0.5D0*(R(I)+R(I+1))
	    EXT_RMID_STORE(I+1)=0.5D0*(R(I)+R(I+1))
	  END DO
	  EXT_RMID_STORE(1)=R(1); EXT_RMID_STORE(ND+1)=R(ND)
	END IF
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
!
! Deallocate all allocated rays if we are using a diferent solution technique.
! This option will only be used when testing, since in CMFGEN we will always use
! the same atmospheric structure.
!
	IF(ALLOCATED(R_EXT) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( LOG_R_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( V_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( Z_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( SIGMA_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( ETA_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( CHI_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( LOG_ETA_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( LOG_CHI_EXT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( TMP_VEC, STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error deallocating R_EXT etc in CMF_FORMA_REL_V3'
	    STOP
	  END IF
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
	    R_EXT(2)=R_EXT(1)-0.001D0*(R_EXT(1)-R_EXT(5))
	    R_EXT(3)=R_EXT(1)-0.1D0*(R_EXT(1)-R_EXT(5))
	    R_EXT(4)=R_EXT(1)-0.4D0*(R_EXT(1)-R_EXT(5))
	    DO I=5,ND_ADD-1
	      R_EXT(I)=R_EXT(I-1)/DEL_R_FAC
	    END DO
	    R_EXT(ND_ADD)=ALPHA
!
	  END IF
	  C_KMS=1.0D-05*SPEED_OF_LIGHT()
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
! Define zone used to extrapolate opacities and emissivities.
!
	  IF(ND_ADD .NE. 0)THEN
	    IDMIN=1
	    T1=R_EXT(1)/R(1)
	    IF(T1 .LT. R(1)/R(ND/3))THEN
	      T1=MIN(3.0D0,T1)
	      DO I=1,ND
	        IF(R(1)/R(I) .GT. T1)THEN
	          IDMAX=MAX(4,I)
	          EXIT
	        END IF
	      END DO
	      IDMAX=MIN(IDMAX,ND/6)
	    ELSE
	      IDMAX=ND/6
	    END IF
	    WRITE(6,'(A,I4,A)')' Using depth 1 and',IDMAX,' to extrapolate opacities'
	  END IF
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
	  IF(CHI(IDMIN) .LE. ESEC(IDMIN) .OR. CHI(IDMAX) .LE. ESEC(IDMAX))THEN
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0D0)ESEC_POW=2.0D0
	    DO I=1,ND_ADD
	      CHI_EXT(I)=CHI(IDMIN)*(R(IDMIN)/R_EXT(I))**ESEC_POW
	    END DO
	  ELSE
	    ALPHA=LOG( (CHI(IDMAX)-ESEC(IDMAX)) / (CHI(IDMIN)-ESEC(IDMIN)) )
	1          /LOG(R(IDMIN)/R(IDMAX))
	    IF(ALPHA .LT. 2.0D0)ALPHA=2.0D0
	    ESEC_POW=LOG(ESEC(IDMAX)/ESEC(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	    IF(ESEC_POW .LT. 2.0D0)ESEC_POW=2.0D0
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
	  ALPHA=LOG(ETA(IDMAX)/ETA(IDMIN))/LOG(R(IDMIN)/R(IDMAX))
	  IF(ALPHA .LT. 3.5)ALPHA=3.5
	  DO I=1,ND_ADD
	    ETA_EXT(I)=ETA(IDMIN)*(R(IDMIN)/R_EXT(I))**ALPHA
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
! This must be after the call to DEFINE_GRID so that RAY_POINTS_INSERTED is defined.
!
	IF(FIRST_TIME .OR. NEW_R_GRID)THEN
	  CALL DEFINE_GRID_V2(R_EXT,V_EXT,VDOP_VEC_EXT,VDOP_FRAC,ND_EXT,R,P,ND,NC,NP)
	  J=0
	  OPEN(UNIT=7,FILE='MU_VALUE_CHK',STATUS='UNKNOWN')
	  WRITE(7,'(A)')' '
	  WRITE(7,'(A)')' Comparison of MU(cmf) and MU(obs) at outer boundary (CMF_FORMAL_REL_V3)'
	  WRITE(7,'(A)')' The  first MU(obs) is the transformed value of MU(cmf)'
	  WRITE(7,'(A)')' The second MU(obs) is simply computed from P(ip) and RMAX'
	  WRITE(7,'(A)')' '
	  WRITE(7,'(2X,A,4(6X,A,2X))')'IP',' MU(cmf)',' HQW(cmf)',' MU(obs)',' MU(obs)'
	  T1=0.0D0
	  DO IP=1,NP
	    J=MAX(J,RAY(IP)%NZ)
	    MU_AT_RMAX=RAY(IP)%MU_P(RAY(IP)%LNK(1))
	    HQW_AT_RMAX=2.0D0*HQW_P(1,IP)
	    WRITE(7,'(I4,4ES16.6)')IP,MU_AT_RMAX,HQW_AT_RMAX,
	1                 (MU_AT_RMAX+V(1)/C_KMS)/(1.0D0+MU_AT_RMAX*V(1)/C_KMS),
	1                 SQRT( (R(1)-P(IP))*(R(1)+P(IP)) )/R(1)
	  END DO
	  CLOSE(UNIT=7)
	  ALLOCATE (CHI_RAY(J),ETA_RAY(J))
	END IF
!
! For SN we can have a hollow core. In this case we need to store the
! inward radiation field for use in calculating the outward radiation field.
! Because of the expansion, the comoving frequencies do not match, hence
! the need for storage of results at earlier frequencies.
!
	IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	  IF(.NOT. ALLOCATED(FREQ_STORE))THEN
	    N_STORE=2.0D0*V(ND)/5.0D0
	    WRITE(6,*)'N_STORE=',N_STORE
	    ALLOCATE (FREQ_STORE(0:N_STORE-1))
	    DO IP=1,NC
	      ALLOCATE (RAY(IP)%I_IN_BND_STORE(0:N_STORE-1))
	    END DO
	    IF(NEW_R_GRID)THEN
	      BETA=V(ND)/2.99792458D+05
	      DO IP=1,NC
	        T1=SQRT( (R(ND)-P(IP))*(R(ND)+P(IP)) )/R(ND)
	        T2=1.0D0-BETA*(T1+BETA)/(1.0D0+BETA*T1)
	        RAY(IP)%FREQ_CONV_FAC=1.0D0/GAM_REL_STORE(ND)/GAM_REL_STORE(ND)/T2/(1.0D0-BETA*T1)
	        WRITE(6,'(5ES14.4)')RAY(IP)%FREQ_CONV_FAC,GAM_REL_STORE(ND),T1,T2,BETA
	      END DO
	    END IF
	  END IF
!
	  IF(INITIALIZE)THEN
	    CUR_LOC=-1
	    FREQ_STORE=0.0D0
	    DO IP=1,NP
	      DO I=1,NC
	        RAY(IP)%I_IN_BND_STORE=0.0D0
	      END DO
	    END DO
	  END IF
	END IF 
!
	IF(RAY_POINTS_INSERTED)THEN
	  LOG_CHI_EXT(1:ND_EXT)=LOG(CHI_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(CHI_COEF,LOG_CHI_EXT,LOG_R_EXT,ND_EXT)
	  LOG_ETA_EXT(1:ND_EXT)=LOG(ETA_EXT(1:ND_EXT))
	  CALL MON_INT_FUNS_V2(ETA_COEF,LOG_ETA_EXT,LOG_R_EXT,ND_EXT)
	END IF
!
	IF(INITIALIZE)THEN
	  DO IP=1,NP
	    RAY(IP)%I_P=0.0D0; RAY(IP)%I_M=0.0D0
	    RAY(IP)%I_P_PREV=0.0D0; RAY(IP)%I_M_PREV=0.0D0
	    RAY(IP)%I_P_SAVE=0.0D0; RAY(IP)%I_M_SAVE=0.0D0
	    HNU_AT_OB_PREV=0.0D0; NNU_AT_OB_PREV=0.0D0
	    HNU_AT_IB_PREV=0.0D0; NNU_AT_IB_PREV=0.0D0
	  END DO	
	ELSE IF(NEW_FREQ)THEN
	  DO IP=1,NP
	    RAY(IP)%I_P=0.0D0; RAY(IP)%I_M=0.0D0
	    RAY(IP)%I_P_PREV=RAY(IP)%I_P_SAVE
	    RAY(IP)%I_M_PREV=RAY(IP)%I_M_SAVE
	  END DO	
	  HNU_AT_OB_PREV=HNU_AT_OB; NNU_AT_OB_PREV=NNU_AT_OB
	  HNU_AT_IB_PREV=HNU_AT_IB; NNU_AT_IB_PREV=NNU_AT_IB
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
	Jnu_store=0.0D0; Hnu_store=0.0D0
	Knu_store=0.0D0; Nnu_store=0.0D0
	JPLUS_IB=0.0D0; HPLUS_IB=0.0D0; KPLUS_IB=0.0D0
	JMIN_IB=0.0D0;  HMIN_IB=0.0D0; KMIN_IB=0.0D0
	JPLUS_OB=0.0D0;  HPLUS_OB=0.0D0; KPLUS_OB=0.0D0
	JMIN_OB=0.0D0;   HMIN_OB=0.0D0;  KMIN_OB=0.0D0
!
! If using the HOLLOW core option, we need to determine location to
! store inner boundary intensity. We do it here, since the storage 
! location hase the same pointer for all rays.
!
	IF(INNER_BND_METH .eq. 'HOLLOW')THEN
          IF(CUR_LOC .EQ. -1)THEN
            CUR_LOC=0
            FREQ_STORE(CUR_LOC)=FREQ
          ELSE IF(FREQ_STORE(CUR_LOC) .NE. FREQ)THEN
            CUR_LOC=MOD(CUR_LOC+1,N_STORE)
            FREQ_STORE(CUR_LOC)=FREQ
	  END IF
	END IF
!
! Determine radiative transfer along each p-ray
!
	NP_LIMIT=NP-1
	IF(THICK_OB)NP_LIMIT=NP
	dBdTAU=DBB/CHI(ND)			!dB/dTAU
	IPLUS=0.0D0
!
	IF(RAY_POINTS_INSERTED)THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(NRAY,I,K,T1,T2,CHI_RAY,ETA_RAY)
	  DO IP=1,NP_LIMIT
!
	    NRAY=RAY(IP)%NZ
	    IF(RAY_POINTS_INSERTED)THEN
              K=1
              DO I=1,RAY(IP)%NZ
100	        CONTINUE
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
            CALL SOLVE_CMF_FORMAL_V2(CHI_RAY,ETA_RAY,IP,FREQ,NU_ON_dNU,INNER_BND_METH,b_planck,dBdTAU,NRAY,NP,NC)
	  END DO
!$OMP END PARALLEL DO
!
	ELSE 
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(NRAY)
	  DO IP=1,NP_LIMIT
	    NRAY=RAY(IP)%NZ
            CALL SOLVE_CMF_FORMAL_V2(CHI_RAY,ETA_RAY,IP,FREQ,NU_ON_dNU,INNER_BND_METH,b_planck,dBdTAU,NRAY,NP,NC)
	  END DO
!$OMP END PARALLEL DO
	END IF
!
	DO IP=1,NP_LIMIT
	  NRAY=RAY(IP)%NZ
	  DO ID=1,MIN(ND,NP-IP+1)
	    I_P_GRID(ID)=RAY(IP)%I_P(RAY(IP)%LNK(ID))
	    I_M_GRID(ID)=RAY(IP)%I_M(RAY(IP)%LNK(ID))
	    IF( I_P_GRID(ID) .LT. 0 .OR. I_M_GRID(ID) .LT. 0.0D0)THEN
	      WRITE(LUER,*)'Error: invalid intensities in CMF_FORMAL_REL_V3'
	      WRITE(LUER,*)'Check file CMF_FORMAL_REL_ERRORS'
	      OPEN(UNIT=7,FILE='CMF_FORMAL_REL_ERRORS',STATUS='UNKNOWN')
	        WRITE(7,*)IP,NRAY
	        WRITE(7,'(6ES14.4)')FREQ,dLOG_NU,NU_ON_dNU,B_PLANCK,DBB
	        WRITE(7,'(A)')' '
	        DO I=1,NRAY
	          WRITE(7,'(8ES14.4)')CHI_RAY(I),ETA_RAY(I),
	1                     RAY(IP)%S_P(I),RAY(IP)%B_P(I),RAY(IP)%I_P(I),
	1                     RAY(IP)%S_M(I),RAY(IP)%B_M(I),RAY(IP)%I_M(I)
	        END DO
	      CLOSE(UNIT=7)
	      STOP
	    END IF
	  END DO
!
! Integrate over p to get J and K. NB: nz(ip) is the number of
! data points along each segment associated with the 2 directions.
!
           DO ID=1,MIN(ND,NP-IP+1)
             Jnu_store(ID)=Jnu_store(ID)+I_p_grid(ID)*Jqw_p(ID,ip)+I_m_grid(ID)*Jqw_m(ID,ip)
             Hnu_store(ID)=Hnu_store(ID)+I_p_grid(ID)*Hqw_p(ID,ip)+I_m_grid(ID)*Hqw_m(ID,ip)
             Knu_store(ID)=Knu_store(ID)+I_p_grid(ID)*Kqw_p(ID,ip)+I_m_grid(ID)*Kqw_m(ID,ip)
             Nnu_store(ID)=Nnu_store(ID)+I_p_grid(ID)*Nqw_p(ID,ip)+I_m_grid(ID)*Nqw_m(ID,ip)
	  END DO
!
          IF(IP .LE. NC+1)THEN
	    JPLUS_IB=JPLUS_IB+I_p_grid(ND)*Jqw_p(ND,ip)
	    HPLUS_IB=HPLUS_IB+I_p_grid(ND)*Hqw_p(ND,ip)
	    KPLUS_IB=KPLUS_IB+I_p_grid(ND)*Kqw_p(ND,ip)
	    NPLUS_IB=NPLUS_IB+I_p_grid(ND)*Nqw_p(ND,ip)
	    JMIN_IB =JMIN_IB +I_m_grid(ND)*Jqw_m(ND,ip)
	    HMIN_IB =HMIN_IB -I_m_grid(ND)*Hqw_m(ND,ip)    !- to make +ve
	    KMIN_IB =KMIN_IB +I_m_grid(ND)*Kqw_m(ND,ip)
	    NMIN_IB =NMIN_IB -I_m_grid(ND)*Nqw_m(ND,ip)
	  END IF
	  JPLUS_OB=JPLUS_OB+I_p_grid(1)*Jqw_p(1,ip)
	  HPLUS_OB=HPLUS_OB+I_p_grid(1)*Hqw_p(1,ip)
	  KPLUS_OB=KPLUS_OB+I_p_grid(1)*Kqw_p(1,ip)
	  NPLUS_OB=NPLUS_OB+I_p_grid(1)*Nqw_p(1,ip)
	  JMIN_OB =JMIN_OB +I_m_grid(1)*Jqw_m(1,ip)
	  HMIN_OB =HMIN_OB -I_m_grid(1)*Hqw_m(1,ip)    !- to make +ve
	  KMIN_OB =KMIN_OB +I_m_grid(1)*Kqw_m(1,ip)
	  NMIN_OB =NMIN_OB -I_m_grid(1)*Nqw_m(1,ip)
!
! Store intensity at the outer boundary.
!
	  IPLUS(IP)=I_P_GRID(1)-I_M_GRID(1)
!
	END DO
!
! Save intensity for integration at next frequency.
!
	DO IP=1,NP
	  RAY(IP)%I_P_SAVE=RAY(IP)%I_P
	  RAY(IP)%I_M_SAVE=RAY(IP)%I_M
	END DO
!
	JNU=JNU_STORE
	FEDD=KNU_STORE/JNU_STORE
	HN_DEF_ON_NODES=.TRUE.
	HNU_AT_OB=Hnu_store(1);  NNU_AT_OB=Nnu_store(1)
	HNU_AT_IB=Hnu_store(ND); NNU_AT_IB=Nnu_store(ND)
!
	HBC=HNU_AT_OB/JNU_STORE(1)
	NBC=NNU_AT_OB/JNU_STORE(1)
	IF(HBC .LT. 0.0D0)HBC=0.0D0
	IF(NBC .LT. 0.0D0)NBC=0.0D0
!
	RET_HNU_AT_OB=HNU_AT_OB
	RET_HNU_AT_IB=HNU_AT_IB
!
!	T1=1.0D-200
!	WRITE(180,'(ES16.6,4(ES16.6,2F8.4))')FREQ,
!	1          JPLUS_IB, HPLUS_IB/MAX(JPLUS_IB,T1), KPLUS_IB/MAX(JPLUS_IB,T1),
!	1          JMIN_IB,  HMIN_IB/MAX(JMIN_IB,T1),   KMIN_IB/MAX(JMIN_IB,T1),
!	1          JPLUS_OB, HPLUS_OB/MAX(JPLUS_OB,T1), KPLUS_OB/MAX(JPLUS_OB,T1),
!	1          JMIN_OB,  HMIN_OB/MAX(JMIN_OB,T1),   KMIN_OB/MAX(JMIN_OB,T1)
!	FLUSH(UNIT=180)
!
	FIRST_TIME=.FALSE.
	RETURN
	END
