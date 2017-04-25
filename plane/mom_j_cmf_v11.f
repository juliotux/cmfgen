!
! Data module for MOM_J_CMF_V8. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE MOD_MOM_J_V11
	IMPLICIT NONE
!
! To be dimenensioned ND_SM where ND_SM is the size of the R grid
! as passed to MOM_J_CMF.
!
	REAL*8, ALLOCATABLE :: LOG_R_SM(:)
!
! Dimensioned ND_SM,4
!
	REAL*8, ALLOCATABLE :: V_COEF(:,:)
	REAL*8, ALLOCATABLE :: SIGMA_COEF(:,:)
	REAL*8, ALLOCATABLE :: ESEC_COEF(:,:)
	REAL*8, ALLOCATABLE :: CHI_COEF(:,:)
	REAL*8, ALLOCATABLE :: ETA_COEF(:,:)
!
! All the following rays have dimension ND, where ND >= ND_SM.
! Some of the data in the arrays is need in subsequent calls.
! 
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: ESEC(:)
	REAL*8, ALLOCATABLE :: CHI(:)
	REAL*8, ALLOCATABLE :: ETA(:)
!
	REAL*8, ALLOCATABLE :: RSQJNU(:)
	REAL*8, ALLOCATABLE :: RSQHNU(:)
!
	REAL*8, ALLOCATABLE :: RSQJNU_PREV(:)
	REAL*8, ALLOCATABLE :: RSQHNU_PREV(:)
!
	REAL*8, ALLOCATABLE :: CON_GAM(:)
	REAL*8, ALLOCATABLE :: CON_GAMH(:)
	REAL*8, ALLOCATABLE :: AV_SIGMA(:)
!
	REAL*8, ALLOCATABLE :: DTAU(:)
	REAL*8, ALLOCATABLE :: DTAUONQ(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: Q(:)
	REAL*8, ALLOCATABLE :: XM(:)
	REAL*8, ALLOCATABLE :: SOURCE(:)
	REAL*8, ALLOCATABLE :: VB(:)
	REAL*8, ALLOCATABLE :: VC(:)
	REAL*8, ALLOCATABLE :: HU(:)
	REAL*8, ALLOCATABLE :: HL(:)
	REAL*8, ALLOCATABLE :: HS(:)
	REAL*8, ALLOCATABLE :: COH_VEC(:)
	REAL*8, ALLOCATABLE :: GAM(:)
	REAL*8, ALLOCATABLE :: GAMH(:)
	REAL*8, ALLOCATABLE :: W(:)
	REAL*8, ALLOCATABLE :: WPREV(:)
	REAL*8, ALLOCATABLE :: PSI(:)
	REAL*8, ALLOCATABLE :: PSIPREV(:)
	REAL*8, ALLOCATABLE :: EPS(:)
	REAL*8, ALLOCATABLE :: EPS_PREV(:)
!
	REAL*8 VDOP_FRAC_SAVE
!
	INTEGER ND
!
! R_PNT(K) defines the interpolation for the variable at depth K.
!
	INTEGER, ALLOCATABLE :: R_PNT(:)
!
! ?_INDX are used to indicate the location of J and H on the small grid
! in the larger array.
!
	INTEGER, ALLOCATABLE :: J_INDX(:)
	INTEGER, ALLOCATABLE :: H_INDX(:)
!
	LOGICAL VERBOSE
	LOGICAL FIRST_TIME
	DATA FIRST_TIME/.TRUE./
	DATA VDOP_FRAC_SAVE/-10001.1D0/    !Absurd value
!
	END MODULE MOD_MOM_J_V11
!
!
!
!
! Routine to compute the mean intensity J at a single frequency in the
! Comoving-Frame. The computed intensity thus depends on the intensity
! computed for the previous (bluer) frequency.
!
! The F, G, and NMID_ON_J Eddingto factors must be supplied.
!
! NB:
!	F = K / J
!	G=  N / H
!	NMID_ON_J(I) = RSQ_N(I)/( RSQ_J(I)+ RQS_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and NMID_ON_J is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
!       NMID_ON_J=0 at all depths.
!     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) NMID_ON_J is defined at all
!       depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or NMID_ON_J is
!       non-zero, and is the value to be used in MOM_J_CMF
!
	SUBROUTINE MOM_J_CMF_V11(ETA_SM,CHI_SM,ESEC_SM,
	1                  V_SM,SIGMA_SM,R_SM,
	1                  JNU,RSQHNU_SM,
	1                  VDOP_VEC,VDOP_FRAC,
	1                  FREQ,dLOG_NU,
	1                  INNER_BND_METH,DBB,IC,IB_STAB_FACTOR,
	1                  N_TYPE,H_CHK_OPTION,METHOD,COHERENT,OUT_BC_TYPE,
	1                  INIT,NEW_FREQ,NC,NP,ND_SM)
	USE MOD_MOM_J_V11
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
! Altered: 17-Oct-2016 : Changed to V11: H_CHK_OPTION inserted into call.
!                          Replaces H_ON_J_CHECK option.
!                          Will allow greater flexibility in testing etc.
! Altered: 17-Feb-2016: Added ZERO_FLUX option. Now pass INNER_BND_METH (instead of DIF).
!                         IB_STAB_FACTOR added for ZERO_FLUX option.
! Altered: 24-Feb-2015: Added CHECK_H_ON_J to call.
! Created: 17-Dec-2008: Based on MOM_J_CMF_V7
!
!
	INTEGER NC,NP,ND_SM
	REAL*8 ETA_SM(ND_SM),CHI_SM(ND_SM),ESEC_SM(ND_SM)
	REAL*8 V_SM(ND_SM),SIGMA_SM(ND_SM),R_SM(ND_SM)
!
! Radiation field variables. F, G, JNU_PREV, and RSQHNU_PREV must be supplied.
! JNU and RSQHNU recomputed.
!
	REAL*8 JNU(ND_SM),RSQHNU_SM(ND_SM)
	REAL*8 VDOP_VEC(ND_SM),VDOP_FRAC
!
	INTEGER N_ERR_MAX,MOM_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL*8 MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
!
! Boundary conditions.
!
	INTEGER OUT_BC_TYPE
	REAL*8 DBB,IC,FREQ,dLOG_NU
	REAL*8 IB_STAB_FACTOR
	CHARACTER(LEN=*) INNER_BND_METH
	CHARACTER(LEN=*) H_CHK_OPTION
	CHARACTER*6 METHOD
	CHARACTER*6 N_TYPE 
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
! NEW_FREQ is used to indicatae that we are computing J for a new frequency.
! If we were iterating between computing J and the Eddington factors, NEW_FREQ
! would be set to false.
!
! COHERENT indicates whether the scattering is coherent. If it is, we
! explicitly take it into account. If COHERENT is FALSE, any electron
! scattering term should be included directly in the ETA that is passed
! to the routine.
!
	LOGICAL INIT,COHERENT,NEW_FREQ
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL*8 T1,T2
	REAL*8  DELTA_R
	INTEGER IT1,I,J,K
	INTEGER IOS
	LOGICAL NEW_R_GRID
!
!
!
	NEW_R_GRID=.FALSE.
	IF(INIT .AND. ALLOCATED(R))THEN
	  DO I=1,ND_SM
	    IF(LOG(R_SM(I)) .NE. LOG_R_SM(I))THEN
	      NEW_R_GRID=.TRUE.
	      WRITE(171,*)'Updating RGRID in MOM_J_CMF_V8'
	      EXIT
	    END IF
	  END DO
          IF(VDOP_FRAC .NE. VDOP_FRAC_SAVE)NEW_R_GRID=.TRUE.
	END IF
	IF(FIRST_TIME)THEN
          OPEN(UNIT=47,FILE='MOM_J_ERRORS',STATUS='UNKNOWN')
          CALL GET_VERBOSE_INFO(VERBOSE)
	ELSE IF(INIT)THEN
	  REWIND(47)
	END IF
!
! Deallocate all arrayes if we have changed VDOP_FRAC. This will only
! be done in testing this routine (e.g., using DISPGEN).
!
	IF(ALLOCATED(R) .AND. NEW_R_GRID)THEN
	  DEALLOCATE ( R )
	  DEALLOCATE ( R_PNT )
	  DEALLOCATE ( LOG_R_SM )
	  DEALLOCATE ( V_COEF )
	  DEALLOCATE ( SIGMA_COEF )
	  DEALLOCATE ( ETA_COEF )
	  DEALLOCATE ( ESEC_COEF )
	  DEALLOCATE ( CHI_COEF )
!
	  DEALLOCATE ( V )
	  DEALLOCATE ( SIGMA )
	  DEALLOCATE ( ETA )
	  DEALLOCATE ( ESEC )
	  DEALLOCATE ( CHI )
!
	  DEALLOCATE ( RSQJNU )
	  DEALLOCATE ( RSQHNU )
	  DEALLOCATE ( RSQJNU_PREV )
	  DEALLOCATE ( RSQHNU_PREV )
!
	  DEALLOCATE ( CON_GAM )
	  DEALLOCATE ( CON_GAMH )
	  DEALLOCATE ( AV_SIGMA )
!
	  DEALLOCATE ( DTAU )
	  DEALLOCATE ( DTAUONQ )
	  DEALLOCATE ( TA )
	  DEALLOCATE ( TB )
	  DEALLOCATE ( TC )
	  DEALLOCATE ( Q )
	  DEALLOCATE ( XM )
	  DEALLOCATE ( SOURCE )
	  DEALLOCATE ( VB )
	  DEALLOCATE ( VC )
	  DEALLOCATE ( HU )
	  DEALLOCATE ( HL )
	  DEALLOCATE ( HS )
	  DEALLOCATE ( COH_VEC )
	  DEALLOCATE ( GAM )
	  DEALLOCATE ( GAMH )
	  DEALLOCATE ( W )
	  DEALLOCATE ( WPREV )
	  DEALLOCATE ( PSI )
	  DEALLOCATE ( PSIPREV )
	  DEALLOCATE ( EPS )
	  DEALLOCATE ( EPS_PREV )
	  DEALLOCATE ( J_INDX )
	  DEALLOCATE ( H_INDX )
	END IF
	VDOP_FRAC_SAVE=VDOP_FRAC
!
! On the very first entry, we define the improved R grid, and allocate all
! data arrays.
!
	IF( FIRST_TIME .OR. .NOT. ALLOCATED(R) )THEN
!
! Determine the number of points for the expanded R grid.
! We always insert an EVEN number of points. This guarentees that
! H_SM (defined at the midpoints of the pass grid) has an exact correspondence
! with H defined on the extended gid.
!
          K=1
          T2=VDOP_FRAC*MINVAL(VDOP_VEC(1:ND_SM))
          DO I=1,ND_SM-1
            IT1=INT( (V_SM(I)-V_SM(I+1))/T2 )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
            IF(IT1 .GT. 0)K=K+IT1
            K=K+1
          END DO
          ND=K
!
	  ALLOCATE ( R(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( R_PNT(ND),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate R and R_PNT in MOM_J_CMF_V8'
	     STOP
	  END IF
          K=1
	  R(1)=R_SM(1)
          R_PNT(1)=1
	  DO I=1,ND_SM-1
            IT1=INT( (V_SM(I)-V_SM(I+1))/T2 )
	    IF( MOD(IT1,2) .NE. 0)IT1=IT1+1
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
!
	  ALLOCATE ( LOG_R_SM(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( V_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ESEC_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate COEF memory in MOM_J_CMF_V8'
	     STOP
	  END IF
!
	  IF(IOS .EQ. 0)ALLOCATE ( V(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SIGMA(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ETA(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( ESEC(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CHI(ND),STAT=IOS )
!
	  IF(IOS .EQ. 0)ALLOCATE ( RSQJNU(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( RSQHNU(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( RSQJNU_PREV(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( RSQHNU_PREV(ND),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate V_BLOCK in MOM_J_CMF_V8'
	     STOP
	  END IF
	  RSQJNU(1:ND)=0.0D0; RSQJNU_PREV(1:ND)=0.0D0
	  RSQHNU(1:ND)=0.0D0; RSQHNU_PREV(1:ND)=0.0D0
!
	  IF(IOS .EQ. 0)ALLOCATE ( CON_GAM(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( CON_GAMH(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( AV_SIGMA(ND),STAT=IOS )
!
	  IF(IOS .EQ. 0)ALLOCATE ( DTAU(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( DTAUONQ(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( TA(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( TB(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( TC(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( Q(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( XM(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( SOURCE(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( VB(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( VC(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( HU(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( HL(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( HS(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( COH_VEC(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( GAM(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( GAMH(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( W(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( WPREV(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( PSI(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( PSIPREV(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( EPS(ND),STAT=IOS )
	  IF(IOS .EQ. 0)ALLOCATE ( EPS_PREV(ND),STAT=IOS )
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to allocate CON_GAM block in MOM_J_CMF_V8'
	     STOP
	  END IF
!
!
	  LOG_R_SM(1:ND_SM)=LOG(R_SM(1:ND_SM))
          CALL MON_INT_FUNS_V2(V_COEF,V_SM,LOG_R_SM,ND_SM)
          CALL MON_INT_FUNS_V2(SIGMA_COEF,SIGMA_SM,LOG_R_SM,ND_SM)
          DO I=1,ND
            K=R_PNT(I)
            T1=LOG(R(I)/R_SM(K))
            V(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
            SIGMA(I)=((SIGMA_COEF(K,1)*T1+SIGMA_COEF(K,2))*T1+SIGMA_COEF(K,3))*T1+SIGMA_COEF(K,4)
          END DO
!
	  ALLOCATE ( J_INDX(ND_SM) );       J_INDX(1:ND_SM)=0
	  ALLOCATE ( H_INDX(ND_SM) );       H_INDX(1:ND_SM)=0
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
	  IF(N_TYPE .NE. 'G_ONLY' .AND. N_TYPE .NE. 'N_ON_J')THEN
	    IF(ND .NE. ND_SM)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in MOM_J_CMF_V8'
	      WRITE(LUER,*)'Cannot use N_TYPE=''MIXED'' when inserting extra points'
	    END IF
	  END IF
!
	  FIRST_TIME=.FALSE.
	END IF
!
!
!
! Get Eddingto values on the revised grid.
!
	CALL GET_MOMS_NON_REL(R, V, FREQ, N_TYPE, ND)
	DO I=1,ND
	  IF(NMID_ON_HMID(I) .GT. 1.0D0)NMID_ON_HMID(I)=1.0D0
	END DO
!
	IF(ND .GT. ND_SM)THEN
!
! Interpolate quantities onto revised grid.
!
	  TA(1:ND_SM)=LOG(CHI_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(CHI_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ESEC_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ESEC_COEF,TA,LOG_R_SM,ND_SM)
	  TA(1:ND_SM)=LOG(ETA_SM(1:ND_SM))
	  CALL MON_INT_FUNS_V2(ETA_COEF,TA,LOG_R_SM,ND_SM)

	  DO I=1,ND
	    K=R_PNT(I)
	    T1=LOG(R(I)/R_SM(K))
	    T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	    CHI(I)=EXP(T2)
	    T2=((ESEC_COEF(K,1)*T1+ESEC_COEF(K,2))*T1+ESEC_COEF(K,3))*T1+ESEC_COEF(K,4)
	    ESEC(I)=EXP(T2)
	    T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	    ETA(I)=EXP(T2)
	  END DO
!
	ELSE
	  ESEC(1:ND)=ESEC_SM(1:ND)
	  CHI(1:ND)=CHI_SM(1:ND)
	  ETA(1:ND)=ETA_SM(1:ND)
	END IF
!
! 
!

	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0D0
	  END DO
	  MOM_ERR_CNT=0
	END IF
!
!*****************************************************************************
!
	DO I=1,ND
	  SOURCE(I)=ETA(I)/CHI(I)
	END DO
	IF(COHERENT)THEN
	  DO I=1,ND
	    COH_VEC(I)=ESEC(I)/CHI(I)
	  END DO
	ELSE
	  DO I=1,ND
	    COH_VEC(I)=0.0D0
	  END DO
	END IF
!
! NB: We actually solve for r^2 J, not J.
!
! Compute the Q factors from F. Then compute optical depth scale.
!
	CALL QFROMF(K_ON_J,Q,R,TA,TB,ND)	!TA work vector
	DO I=1,ND
	  TA(I)=CHI(I)*Q(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
!
	IF(NEW_FREQ)THEN
	  K_ON_J_PREV(1:ND)=K_ON_J_SAVE(1:ND)
	  NMID_ON_HMID_PREV(1:ND)=NMID_ON_HMID_SAVE(1:ND)
	  NMID_ON_J_PREV(1:ND)=NMID_ON_J_SAVE(1:ND)
	  RSQJNU_PREV(1:ND)=RSQJNU(1:ND)
	  RSQHNU_PREV(1:ND)=RSQHNU(1:ND)
	  HBC_PREV=HBC_SAVE;  IN_HBC_PREV=IN_HBC_SAVE
	  NBC_PREV=NBC_SAVE
	END IF
!
!
	IF(INIT)THEN
	  DO I=1,ND
	    GAMH(I)=0.0D0
	    GAM(I)=0.0D0
	    W(I)=0.0D0
	    WPREV(I)=0.0D0
	    PSI(I)=0.0D0
	    PSIPREV(I)=0.0D0
	    RSQJNU_PREV(I)=0.0D0
	    RSQHNU_PREV(I)=0.0D0
	    EPS(I)=0.0D0
	    EPS_PREV(I)=0.0D0
	  END DO
	  HBC_PREV=0.0D0;  IN_HBC_PREV=0.0D0; NBC_PREV=0.0D0
	ELSE
!
! Assume (1)	SIGMAd+1/2 = 0.5*( SIGMAd+1+SIGMAd )
! 	 (2)	Vd+1/2=0.5*( Vd + Vd+1 )
! Note that V is in km/s and SIGMA=(dlnV/dlnR-1.0)
!
	  DO I=1,ND-1
	    CON_GAMH(I)=2.0D0*3.33564D-06*(V(I)+V(I+1))/(R(I)+R(I+1))
	    AV_SIGMA(I)=0.5D0*(SIGMA(I)+SIGMA(I+1))
	    CON_GAM(I)=3.33564D-06*V(I)/R(I)
	  END DO
	  CON_GAM(ND)=3.33564D-06*V(ND)/R(ND)
!
! Since we are intgerating from blue to red, FL_PREV is always larger than
! FL. dLOG_NU is define as vd / dv which is the same as d / d ln v.
!
! EPS is used if we define N in terms of J rather than H, This is sometimes
! useful as H can approach zero, and hence N/H is undefined.
!
	  IF(N_TYPE .EQ. 'G_ONLY')THEN
	    DO I=1,ND-1
	      GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	      W(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*NMID_ON_HMID(I) )
	      WPREV(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*NMID_ON_HMID_PREV(I) )
	    END DO
	  ELSE
	    DO I=1,ND-1
	      GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	      W(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*NMID_ON_HMID(I) )
	      WPREV(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*NMID_ON_HMID_PREV(I) )
	      EPS(I)=GAMH(I)*AV_SIGMA(I)*NMID_ON_J(I)/(1.0D0+W(I))
	      EPS_PREV(I)=GAMH(I)*AV_SIGMA(I)*NMID_ON_J_PREV(I)/(1.0D0+W(I))
	    END DO
	  END IF
!
	  DO I=1,ND
	    GAM(I)=CON_GAM(I)/CHI(I)/dLOG_NU
	  END DO
!
	END IF
!
! 
!
	DO I=2,ND-1
	  DTAUONQ(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
	  PSI(I)=DTAUONQ(I)*GAM(I)*( 1.0D0+SIGMA(I)*K_ON_J(I) )
	  PSIPREV(I)=DTAUONQ(I)*GAM(I)*(  1.0D0+SIGMA(I)*K_ON_J_PREV(I) )
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=K_ON_J(I+1)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=K_ON_J(I)*Q(I)/(1.0D0+W(I))/DTAU(I)
	  HS(I)=WPREV(I)/(1.0D0+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)
	    TC(I)=-HU(I)
	    TB(I)=DTAUONQ(I)*(1.0D0-COH_VEC(I)) + PSI(I) + HL(I) + HU(I-1)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	  END DO
	ELSE
	  DO I=2,ND-1
	    TA(I)=-HL(I-1)-EPS(I-1)
	    TC(I)=-HU(I)+EPS(I)
	    TB(I)=DTAUONQ(I)*(1.0D0-COH_VEC(I)) + PSI(I) + HL(I) + HU(I-1)
	1               -EPS(I-1)+EPS(I)
	    VB(I)=-HS(I-1)
	    VC(I)=HS(I)
	    XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	  END DO
	END IF
!
	TA(ND)=-K_ON_J(ND-1)*Q(ND-1)/DTAU(ND-1)
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  TB(ND)=K_ON_J(ND)/DTAU(ND-1)
	  XM(ND)=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	  TB(ND)=K_ON_J(ND)/DTAU(ND-1)
	  XM(ND)=IB_STAB_FACTOR*TB(ND)*R(ND)*R(ND)*(JPLUS_IB+JMIN_IB)
	  TB(ND)=(1.0D0+IB_STAB_FACTOR)*TB(ND)
	ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	  TB(ND)=K_ON_J(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid INNER_BND_METH:', INNER_BND_METH
	  STOP
	END IF
	TC(ND)=0.0D0
	VB(ND)=0.0D0
	VC(ND)=0.0D0
	PSIPREV(ND)=0.0D0
!
! Note that EPS and EPS_PREV will be identically zero hence when N_TYPE is
! G_ONLY.
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*RSQHNU_PREV(I-1) + VC(I)*RSQHNU_PREV(I)
	1          + PSIPREV(I)*RSQJNU_PREV(I)
	  END DO
	  XM(ND)=XM(ND)
	ELSE
	  DO I=2,ND-1
	    XM(I)=XM(I) + VB(I)*RSQHNU_PREV(I-1) + VC(I)*RSQHNU_PREV(I)
	1          + PSIPREV(I)*RSQJNU_PREV(I)
	1          - EPS_PREV(I-1)*(RSQJNU_PREV(I-1)+RSQJNU_PREV(I))
	1          + EPS_PREV(I)*(RSQJNU_PREV(I)+RSQJNU_PREV(I+1))
	  END DO
	  XM(ND)=XM(ND)
	END IF
!
!
! Evaluate TA,TB,TC for boudary conditions
!
	IF(OUT_BC_TYPE .LE. 1)THEN 		!Old (def) BC.
	  PSI(1)=GAM(1)*(HBC+NBC*SIGMA(1))
	  PSIPREV(1)=GAM(1)*( HBC_PREV+NBC_PREV*SIGMA(1) )
	  TC(1)=-K_ON_J(2)*Q(2)/DTAU(1)
	  TB(1)=K_ON_J(1)*Q(1)/DTAU(1) + PSI(1) + HBC
	  XM(1)=PSIPREV(1)*RSQJNU_PREV(1)
	  TA(1)=0.0D0
	  VB(1)=0.0D0
	  VC(1)=0.0D0
	ELSE
	  PSI(1)=GAM(1)*(1.0D0+SIGMA(1)*K_ON_J(1))
	  PSIPREV(1)=GAM(1)*(1.0D0+SIGMA(1)*K_ON_J_PREV(1))
	  T1=0.25D0*(CHI(2)+CHI(1))*(R(2)-R(1))
	  TC(1)=(HU(1)-EPS(1))/T1
	  TB(1)= (COH_VEC(1)-1.0D0) -(HL(1)+HBC+EPS(1))/T1 -PSI(1)
	  XM(1)=-SOURCE(1)*R(1)*R(1) -HS(1)*RSQHNU_PREV(1)/T1 -PSIPREV(1)*RSQJNU_PREV(1)
	  IF(N_TYPE .NE. 'G_ONLY')THEN
	   XM(1)=XM(1)-EPS_PREV(1)*(RSQJNU_PREV(1)+RSQJNU_PREV(2))/T1
	  END IF
	END IF
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
! Check that no negative mean intensities have been computed.
!
	IF(MINVAL(XM(1:ND)) .LE. 0.0D0 .AND. VERBOSE)THEN
	  WRITE(47,'(/,A,E16.8)')'Freq=',FREQ
	  TA(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
	  CALL WRITE_VEC(TA,ND,'XM Vec',47)
	  CALL WRITE_VEC(K_ON_J,ND,'K_ON_J (FEDD) Vec',47)
	  IF(N_TYPE .NE. 'G_ONLY')CALL WRITE_VEC(NMID_ON_J,ND,'NMID_ON_J Vec',47)
	  IF(N_TYPE .NE. 'N_ON_J')CALL WRITE_VEC(NMID_ON_HMID,ND,'NMID_ON_HMID (G) Vec',47)
	  CALL WRITE_VEC(ETA,ND,'ETA Vec',47)
	  CALL WRITE_VEC(ESEC,ND,'ESEC Vec',47)
	  CALL WRITE_VEC(CHI,ND,'CHI Vec',47)
	END IF
!
	RECORDED_ERROR=.FALSE.
	DO I=1,ND
	  IF(XM(I) .LT. 0.0D0)THEN
	    IF(.NOT. VERBOSE)THEN
	      WRITE(47,'(I5,ES16.8,12ES13.4)')I,FREQ,XM(I),ETA(I),CHI(I),ESEC(I),K_ON_J(I),
	1                         NMID_ON_J(I),NMID_ON_HMID(I),XM(MAX(1,I-2):MIN(I+2,ND))
	    END IF
	    XM(I)=ABS(XM(I))/10.0D0
	    IF(.NOT. RECORDED_ERROR)THEN
	      IF(MOM_ERR_CNT .GT. N_ERR_MAX)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	      ELSE IF(MOM_ERR_ON_FREQ(MOM_ERR_CNT) .NE. FREQ)THEN
	        MOM_ERR_CNT=MOM_ERR_CNT+1
	        IF(MOM_ERR_CNT .LT. N_ERR_MAX)MOM_ERR_ON_FREQ(MOM_ERR_CNT)=FREQ
	      END IF
	      RECORDED_ERROR=.TRUE.
	    END IF
	  END IF
	END DO
!
! Save R^2 J for next iteration.
!
	RSQJNU(1:ND)=XM(1:ND)
!
	IF(N_TYPE .EQ. 'G_ONLY')THEN
	  DO I=1,ND-1
	    RSQHNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQHNU_PREV(I)
	  END DO
	ELSE
	  DO I=1,ND-1
	    RSQHNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQHNU_PREV(I)+
	1              ( EPS_PREV(I)*(RSQJNU_PREV(I)+RSQJNU_PREV(I+1)) -
	1                  EPS(I)*(XM(I)+XM(I+1)) )
	  END DO
	END IF
!
! Make sure H satisfies the basic requirement that it is less than J.
!
	IF(H_CHK_OPTION .EQ. 'AV_VAL')THEN
	  DO I=1,ND-1
	    T1=(RSQJNU(I)+RSQJNU(I+1))/2.0D0
	    IF(RSQHNU(I) .GT. T1)THEN
	      RSQHNU(I)=0.99999D0*T1
	    ELSE IF(RSQHNU(I) .LT. -T1)THEN
	      RSQHNU(I)=-0.99999D0*T1
	    END IF
	  END DO
	ELSE IF(H_CHK_OPTION .EQ. 'MAX_VAL')THEN
	  DO I=1,ND-1
	    T1=MAX(RSQJNU(I),RSQJNU(I+1))
	    IF(RSQHNU(I) .GT. T1)THEN
	      RSQHNU(I)=0.99999D0*T1
	    ELSE IF(RSQHNU(I) .LT. -T1)THEN
	      RSQHNU(I)=-0.99999D0*T1
	    END IF
	  END DO
	END IF
!
! Regrid derived J and RSQH values onto small grid. We devide RSQJ by R^2 so that
! we return J.
!
	DO I=1,ND_SM
	  K=J_INDX(I)
	  JNU(I)=RSQJNU(K)/R_SM(I)/R_SM(I)
	END DO
!
	DO I=1,ND_SM-1
	  K=H_INDX(I)
	  RSQHNU_SM(I)=RSQHNU(K)
	END DO
!
	K_ON_J_SAVE(1:ND)=K_ON_J(1:ND)
	NMID_ON_HMID_SAVE(1:ND)=NMID_ON_HMID(1:ND)
	NMID_ON_J_SAVE(1:ND)=NMID_ON_J(1:ND)
	HBC_SAVE=HBC
	IN_HBC_SAVE=IN_HBC
	NBC_SAVE=NBC
!
	IF(NEW_R_GRID)WRITE(171,*)'Exiting MOM_J_CMF_V8'
!
	RETURN
	END
