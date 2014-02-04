
! Data module for MOM_J_CMF_V7. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE MOD_MOM_J_V7
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
	REAL*8, ALLOCATABLE :: F_COEF(:,:)
	REAL*8, ALLOCATABLE :: G_COEF(:,:)
	REAL*8, ALLOCATABLE :: RSQN_ON_RSQJ_COEF(:,:)
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
	REAL*8, ALLOCATABLE :: F(:)
	REAL*8, ALLOCATABLE :: G(:)
	REAL*8, ALLOCATABLE :: RSQN_ON_RSQJ(:)
	REAL*8, ALLOCATABLE :: F_SAV(:)
	REAL*8, ALLOCATABLE :: G_SAV(:)
	REAL*8, ALLOCATABLE :: RSQN_ON_RSQJ_SAV(:)
!
	REAL*8, ALLOCATABLE :: RSQJNU_PREV(:)
	REAL*8, ALLOCATABLE :: RSQHNU_PREV(:)
	REAL*8, ALLOCATABLE :: F_PREV(:)
	REAL*8, ALLOCATABLE :: G_PREV(:)
	REAL*8, ALLOCATABLE :: RSQN_ON_RSQJ_PREV(:)
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
! Boundary conditions
!
        REAL*8 HBC_PREV
	REAL*8 IN_HBC_PREV
	REAL*8 NBC_PREV
!
	REAL*8 HBC_SAV
	REAL*8 IN_HBC_SAV
	REAL*8 NBC_SAV
	REAL*8 VDOP_FRAC_SAV
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
	LOGICAL FIRST_TIME
	DATA FIRST_TIME/.TRUE./
	DATA VDOP_FRAC_SAV/-10001.1D0/    !Absurd value
!
	END MODULE MOD_MOM_J_V7
!
!
!
!
! Routine to compute the mean intensity J at a single frequency in the
! Comoving-Frame. The computed intensity thus depends on the intensity
! computed for the previous (bluer) frequency.
!
! The F, G, and RSQN_ON_RSQJ Eddingto factors must be supplied.
!
! NB:
!	F = K / J
!	G=  N / H
!	RSQN_ON_RSQJ(I) = RSQ_N(I)/( RSQ_J(I)+ RQS_J(I+1))
!
! where
!	RSQ_N=0.25*N(I)*(R(I)+R(I+1))**2
!
! NB: Only one of G and RSQN_ON_RSQJ is defined at a given depth. This
!     avoids having to test what mode I am using for the Eddington factors.
!
!     IF N_TYPE='G_ONLY' (in FG_J_CMF_V4) G is defined at all depths, and
!       RSQN_ON_RSQJ=0 at all depths.
!     IF N_TYPE='N_ON_J' (in FG_J_CMF_V4) RSQN_ON_RSQJ is defined at all
!       depths, and G=0 at all depths.
!     IF N_TYPE='MIXED' (in FG_J_CMF_V4) one of G or RSQN_ON_RSQJ is
!       non-zero, and is the value to be used in MOM_J_CMF
!
	SUBROUTINE MOM_J_CMF_V7(ETA_SM,CHI_SM,ESEC_SM,
	1                  V_SM,SIGMA_SM,R_SM,
	1		   F_SM,G_SM,RSQN_ON_RSQJ_SM,
	1                  JNU_SM,RSQHNU_SM,
	1                  VDOP_VEC,VDOP_FRAC,
	1                  HBC,IN_HBC,NBC,
	1                  FREQ,dLOG_NU,DIF,DBB,IC,
	1                  N_TYPE,METHOD,COHERENT,OUT_BC_TYPE,
	1                  INIT,NEW_FREQ,NC,NP,ND_SM)
	USE MOD_MOM_J_V7
	IMPLICIT NONE
!
! Altered :  27-Feb-2004: Bug fix with RSQHNU_PREV (not sure effect) & access to DTAU(ND)
!                           which doesn't effect results. 
! Altered :  26-Jan-2004: Changed to allow it to work with changes in the R grid.
! Altered :  01-Jan-2004: Important update: SP constants to DP constants;
!                          Effected error messgaes to fort.47 
! Altered:   28-Feb-2002: Extensive modifications and cleaning.
!                         Routine now allow variable VDOP_FRAC for testing 
!                           purposes. This option shoul NOT be used with CMFGEN.
!                              
! Created:   14-Nov-2001: Based on MOM_J_CMF_V5
!
	INTEGER NC,NP,ND_SM
	REAL*8 ETA_SM(ND_SM),CHI_SM(ND_SM),ESEC_SM(ND_SM)
	REAL*8 V_SM(ND_SM),SIGMA_SM(ND_SM),R_SM(ND_SM)
!
! Radiation field variables. F, G, JNU_PREV, and RSQHNU_PREV must be supplied.
! JNU and RSQHNU recomputed.
!
	REAL*8 F_SM(ND_SM),G_SM(ND_SM),RSQN_ON_RSQJ_SM(ND_SM)
	REAL*8 JNU_SM(ND_SM),RSQHNU_SM(ND_SM)
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
	REAL*8 HBC,NBC,IN_HBC
	REAL*8 DBB,IC,FREQ,dLOG_NU
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
	LOGICAL DIF,INIT,COHERENT,NEW_FREQ
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
	      WRITE(171,*)'Updating RGRID in MOM_J_CMF_V7'
	      EXIT
	    END IF
	  END DO
          IF(VDOP_FRAC .NE. VDOP_FRAC_SAV)NEW_R_GRID=.TRUE.
	END IF
	IF(FIRST_TIME)THEN
          OPEN(UNIT=47,FILE='MOM_J_ERRORS',STATUS='UNKNOWN')
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
	  DEALLOCATE ( F_COEF )
	  DEALLOCATE ( G_COEF )
	  DEALLOCATE ( RSQN_ON_RSQJ_COEF )
!
	  DEALLOCATE ( V )
	  DEALLOCATE ( SIGMA )
	  DEALLOCATE ( ETA )
	  DEALLOCATE ( ESEC )
	  DEALLOCATE ( CHI )
!
	  DEALLOCATE ( RSQJNU )
	  DEALLOCATE ( RSQHNU )
	  DEALLOCATE ( F )
	  DEALLOCATE ( G )
	  DEALLOCATE ( RSQN_ON_RSQJ )
	  DEALLOCATE ( F_SAV )
	  DEALLOCATE ( G_SAV )
	  DEALLOCATE ( RSQN_ON_RSQJ_SAV )
!
	  DEALLOCATE ( RSQJNU_PREV )
	  DEALLOCATE ( RSQHNU_PREV )
	  DEALLOCATE ( F_PREV )
	  DEALLOCATE ( G_PREV )
	  DEALLOCATE ( RSQN_ON_RSQJ_PREV )
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
	VDOP_FRAC_SAV=VDOP_FRAC
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
	  ALLOCATE ( R(ND) )
	  ALLOCATE ( R_PNT(ND) )
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
	  IF(IOS .EQ. 0)ALLOCATE ( F_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( G_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( RSQN_ON_RSQJ_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(2,*)'Unable to allocate COEF memory in MOM_J_CMF_V7'
	     STOP
	  END IF
!
	  ALLOCATE ( V(ND) )
	  ALLOCATE ( SIGMA(ND) )
	  ALLOCATE ( ETA(ND) )
	  ALLOCATE ( ESEC(ND) )
	  ALLOCATE ( CHI(ND) )
!
	  ALLOCATE ( RSQJNU(ND) )             ; RSQJNU(1:ND)=0.0D0
	  ALLOCATE ( RSQHNU(ND) )             ; RSQHNU(1:ND)=0.0D0
	  ALLOCATE ( F(ND) )
	  ALLOCATE ( G(ND) )
	  ALLOCATE ( RSQN_ON_RSQJ(ND) )
	  ALLOCATE ( F_SAV(ND) )              ; F_SAV(1:ND)=0.0D0
	  ALLOCATE ( G_SAV(ND) )              ; G_SAV(1:ND)=0.0D0
	  ALLOCATE ( RSQN_ON_RSQJ_SAV(ND) )   ; RSQN_ON_RSQJ_SAV(1:ND)=0.0D0
!
	  ALLOCATE ( RSQJNU_PREV(ND) )        ; RSQJNU_PREV(1:ND)=0.0D0
	  ALLOCATE ( RSQHNU_PREV(ND) )        ; RSQHNU_PREV(1:ND)=0.0D0
	  ALLOCATE ( F_PREV(ND) )
	  ALLOCATE ( G_PREV(ND) )
	  ALLOCATE ( RSQN_ON_RSQJ_PREV(ND) )
!
	  ALLOCATE ( CON_GAM(ND) )
	  ALLOCATE ( CON_GAMH(ND) )
	  ALLOCATE ( AV_SIGMA(ND) )
!
	  ALLOCATE ( DTAU(ND) )
	  ALLOCATE ( DTAUONQ(ND) )
	  ALLOCATE ( TA(ND) )
	  ALLOCATE ( TB(ND) )
	  ALLOCATE ( TC(ND) )
	  ALLOCATE ( Q(ND) )
	  ALLOCATE ( XM(ND) )
	  ALLOCATE ( SOURCE(ND) )
	  ALLOCATE ( VB(ND) )
	  ALLOCATE ( VC(ND) )
	  ALLOCATE ( HU(ND) )
	  ALLOCATE ( HL(ND) )
	  ALLOCATE ( HS(ND) )
	  ALLOCATE ( COH_VEC(ND) )
	  ALLOCATE ( GAM(ND) )
	  ALLOCATE ( GAMH(ND) )
	  ALLOCATE ( W(ND) )
	  ALLOCATE ( WPREV(ND) )
	  ALLOCATE ( PSI(ND) )
	  ALLOCATE ( PSIPREV(ND) )
	  ALLOCATE ( EPS(ND) )
	  ALLOCATE ( EPS_PREV(ND) )
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
	      WRITE(LUER,*)'Error in MOM_J_CMF_V7'
	      WRITE(LUER,*)'Cannot use N_TYPE=''MIXED'' when inserting extra points'
	    END IF
	  END IF
!
	  FIRST_TIME=.FALSE.
	END IF
!
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
	  CALL MON_INT_FUNS_V2(F_COEF,F_SM,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(G_COEF,G_SM,LOG_R_SM,ND_SM)
	  CALL MON_INT_FUNS_V2(RSQN_ON_RSQJ_COEF,RSQN_ON_RSQJ_SM,LOG_R_SM,ND_SM)

	  DO I=1,ND
	    K=R_PNT(I)
	    T1=LOG(R(I)/R_SM(K))
	    T2=((CHI_COEF(K,1)*T1+CHI_COEF(K,2))*T1+CHI_COEF(K,3))*T1+CHI_COEF(K,4)
	    CHI(I)=EXP(T2)
	    T2=((ESEC_COEF(K,1)*T1+ESEC_COEF(K,2))*T1+ESEC_COEF(K,3))*T1+ESEC_COEF(K,4)
	    ESEC(I)=EXP(T2)
	    T2=((ETA_COEF(K,1)*T1+ETA_COEF(K,2))*T1+ETA_COEF(K,3))*T1+ETA_COEF(K,4)
	    ETA(I)=EXP(T2)
	    F(I)=((F_COEF(K,1)*T1+F_COEF(K,2))*T1+F_COEF(K,3))*T1+F_COEF(K,4)
	    G(I)=((G_COEF(K,1)*T1+G_COEF(K,2))*T1+G_COEF(K,3))*T1+G_COEF(K,4)
	    RSQN_ON_RSQJ(I)=((RSQN_ON_RSQJ_COEF(K,1)*T1+RSQN_ON_RSQJ_COEF(K,2))*T1+
	1                               RSQN_ON_RSQJ_COEF(K,3))*T1+RSQN_ON_RSQJ_COEF(K,4)
	  END DO
!
	ELSE
	  ESEC(1:ND)=ESEC_SM(1:ND)
	  CHI(1:ND)=CHI_SM(1:ND)
	  ETA(1:ND)=ETA_SM(1:ND)
	  F(1:ND)=F_SM(1:ND)
	  G(1:ND)=G_SM(1:ND)
	  RSQN_ON_RSQJ(1:ND)=RSQN_ON_RSQJ_SM(1:ND)
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
	DO I=1,ND
	  IF(G(I) .GT. 1.0D0)G(I)=1.0D0
	END DO
!
! NB: We actually solve for r^2 J, not J.
!
! Compute the Q factors from F. Then compute optical depth scale.
!
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA work vector
	DO I=1,ND
	  TA(I)=CHI(I)*Q(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
!
	IF(NEW_FREQ)THEN
	  F_PREV(1:ND)=F_SAV(1:ND)
	  G_PREV(1:ND)=G_SAV(1:ND)
	  RSQN_ON_RSQJ_PREV(1:ND)=RSQN_ON_RSQJ_SAV(1:ND)
	  RSQJNU_PREV(1:ND)=RSQJNU(1:ND)
	  RSQHNU_PREV(1:ND)=RSQHNU(1:ND)
	  HBC_PREV=HBC_SAV;  IN_HBC_PREV=IN_HBC_SAV
	  NBC_PREV=NBC_SAV
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
	      W(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*G(I) )
	      WPREV(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*G_PREV(I) )
	    END DO
	  ELSE
	    DO I=1,ND-1
	      GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	      W(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*G(I) )
	      WPREV(I)=GAMH(I)*( 1.0D0+AV_SIGMA(I)*G_PREV(I) )
	      EPS(I)=GAMH(I)*AV_SIGMA(I)*RSQN_ON_RSQJ(I)/(1.0D0+W(I))
	      EPS_PREV(I)=GAMH(I)*AV_SIGMA(I)*RSQN_ON_RSQJ_PREV(I)/(1.0D0+W(I))
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
	  PSI(I)=DTAUONQ(I)*GAM(I)*( 1.0D0+SIGMA(I)*F(I) )
	  PSIPREV(I)=DTAUONQ(I)*GAM(I)*(  1.0D0+SIGMA(I)*F_PREV(I) )
	END DO
!
! Compute vectors used to compute the flux vector H.
!
	DO I=1,ND-1
	  HU(I)=F(I+1)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=F(I)*Q(I)/(1.0D0+W(I))/DTAU(I)
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
	TA(ND)=-F(ND-1)*Q(ND-1)/DTAU(ND-1)
	IF(DIF)THEN
	  TB(ND)=F(ND)/DTAU(ND-1)
	  XM(ND)=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	ELSE
	  TB(ND)=F(ND)/DTAU(ND-1)+IN_HBC
	  XM(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*IN_HBC)
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
	  TC(1)=-F(2)*Q(2)/DTAU(1)
	  TB(1)=F(1)*Q(1)/DTAU(1) + PSI(1) + HBC
	  XM(1)=PSIPREV(1)*RSQJNU_PREV(1)
	  TA(1)=0.0D0
	  VB(1)=0.0D0
	  VC(1)=0.0D0
	ELSE
	  PSI(1)=GAM(1)*(1.0D0+SIGMA(1)*F(1))
	  PSIPREV(1)=GAM(1)*(1.0D0+SIGMA(1)*F_PREV(1))
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
	IF(MINVAL(XM(1:ND)) .LE. 0.0D0)THEN
	   WRITE(47,*)'Freq=',FREQ
	   TA(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
	   CALL WRITV(TA,ND,'XM Vec',47)
	   CALL WRITV(F,ND,'F Vec',47)
	   CALL WRITV(G,ND,'G Vec',47)
	   CALL WRITV(ETA,ND,'ETA Vec',47)
	   CALL WRITV(ESEC,ND,'ESEC Vec',47)
	   CALL WRITV(CHI,ND,'CHI Vec',47)
	END IF
!
	DO I=1,ND
	  IF(XM(I) .LT. 0.0D0)THEN
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
! Regrid derived J and RSQH values onto small grid. We devide RSQJ by R^2 so that
! we return J.
!
	DO I=1,ND_SM
	  K=J_INDX(I)
	  JNU_SM(I)=RSQJNU(K)/R_SM(I)/R_SM(I)
	END DO
!
	DO I=1,ND_SM-1
	  K=H_INDX(I)
	  RSQHNU_SM(I)=RSQHNU(K)
	END DO
!
	F_SAV(1:ND)=F(1:ND)
	G_SAV(1:ND)=G(1:ND)
	RSQN_ON_RSQJ_SAV(1:ND)=RSQN_ON_RSQJ(1:ND)
	HBC_SAV=HBC
	IN_HBC_SAV=IN_HBC
	NBC_SAV=NBC
!
	IF(NEW_R_GRID)WRITE(171,*)'Exiting MOM_J_CMF_V7'
!
	RETURN
	END
