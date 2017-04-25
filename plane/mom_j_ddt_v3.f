!
! Data module for MOM_J_DDT_V3. Data placed in this module is automatically
! saved between subroutine calls..
!
	MODULE MOD_MOM_J_DDT_V3
	IMPLICIT NONE
!
! Dimensioned ND_SM,4
! ND_SM is the depth dimension passed to the routine, and hence is the same
! as ND in CMFGEN.
!
	REAL*8, ALLOCATABLE :: V_COEF(:,:)
	REAL*8, ALLOCATABLE :: ESEC_COEF(:,:)
	REAL*8, ALLOCATABLE :: CHI_COEF(:,:)
	REAL*8, ALLOCATABLE :: ETA_COEF(:,:)
	REAL*8, ALLOCATABLE :: F_COEF(:,:)
!
! All the following rays have dimension ND, where ND >= ND_SM.
! Some of the data in the arrays is need in subsequent calls.
! 
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: ESEC(:)
	REAL*8, ALLOCATABLE :: CHI(:)
	REAL*8, ALLOCATABLE :: ETA(:)
!
	REAL*8, ALLOCATABLE :: RSQJNU(:)
	REAL*8, ALLOCATABLE :: RSQHNU(:)
	REAL*8, ALLOCATABLE :: RSQJNU_PREV(:)
	REAL*8, ALLOCATABLE :: RSQHNU_PREV(:)
	REAL*8, ALLOCATABLE :: RSQJNU_OLDT(:)
	REAL*8, ALLOCATABLE :: RSQHNU_OLDT(:)
!
	REAL*8, ALLOCATABLE :: F(:)
	REAL*8, ALLOCATABLE :: F_SAV(:)
	REAL*8, ALLOCATABLE :: F_PREV(:)
!
	REAL*8, ALLOCATABLE :: CON_GAM(:)
	REAL*8, ALLOCATABLE :: CON_GAMH(:)
!
	REAL*8, ALLOCATABLE :: COH_VEC(:)
	REAL*8, ALLOCATABLE :: DJDT(:)
	REAL*8, ALLOCATABLE :: DTAU(:)
	REAL*8, ALLOCATABLE :: DTAUONQ(:)
	REAL*8, ALLOCATABLE :: GAM(:)
	REAL*8, ALLOCATABLE :: GAMH(:)
	REAL*8, ALLOCATABLE :: dH(:)
	REAL*8, ALLOCATABLE :: dH_OLDT(:)
	REAL*8, ALLOCATABLE :: HL(:)
	REAL*8, ALLOCATABLE :: HU(:)
	REAL*8, ALLOCATABLE :: HS(:)
	REAL*8, ALLOCATABLE :: HT(:)
	REAL*8, ALLOCATABLE :: PSI(:)
	REAL*8, ALLOCATABLE :: PSIPREV(:)
	REAL*8, ALLOCATABLE :: Q(:)
	REAL*8, ALLOCATABLE :: SOURCE(:)
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: W(:)
	REAL*8, ALLOCATABLE :: WPREV(:)
	REAL*8, ALLOCATABLE :: XM(:)
!
! Boundary conditions: INBC denotes the inner boundary.
! OUTBC denotes the outer boundary. _PREV refers to the
! boundary condition at the previous frequency, _SAV
! refers the boundary condition just computed (which will
! become _PREV when NEW_FREW is TRUE), and _OLD refers to
! the boundary condition at the previous time step.
!
        REAL*8 HONJ_OUTBC_PREV
	REAL*8 RSQH_AT_OB_PREV
	REAL*8 RSQH_AT_OB_SAV
	REAL*8 HONJ_OUTBC_SAV
!
	REAL*8 RSQH_AT_IB_PREV
	REAL*8 RSQH_AT_IB_SAV
	REAL*8 RSQH_AT_IB_OLDT
	REAL*8 HONJ_OUTBC_OLDT
!
	REAL*8 DELTA_TIME_SECS
	REAL*8 RECIP_CDELTAT
	REAL*8 ROLD_ON_R
	REAL*8 R_RAT_FOR_J
	REAL*8 C_KMS
	REAL*8 VDOP_FRAC_SAV
!
! To be dimensioned ND_SM where ND_SM is the size of the R grid
! as passed to MOM_J_CMF.
!
	REAL*8, ALLOCATABLE :: R_SM_SAV(:)
	REAL*8, ALLOCATABLE :: LOG_R_SM(:)
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
	INTEGER ND
	INTEGER, SAVE :: ND_SAV=0
	INTEGER, SAVE :: ND_SM_SAV=0
!
	LOGICAL, SAVE::  FIRST_TIME=.TRUE.
	DATA VDOP_FRAC_SAV/-10001.1D0/    !Absurd value
!
	END MODULE MOD_MOM_J_DDT_V3
!
!
!
! Routine to compute the mean intensity J at a single frequency in the
! comoving-frame. Time dependence is taken into account. The time
! derivative is in the Lagrangian frame. This routine ASSUMES a
! Hubble flow (i.e., v proportional to r).
!
! The computed mean intensity depends on:
!
!       J & H computed for the previous (bluer) frequency.
!       J & H computed at the previous time step.
!
! The F Eddington factor (=K/J)  must be supplied. J & H  for the
! previous frequency are stored internally, while J and H for
! the old time step are obtained via a subroutine call.
!
	SUBROUTINE MOM_J_DDT_V3(ETA_SM,CHI_SM,ESEC_SM,
	1            V_SM,R_SM,F_SM,
	1            JNU_SM,RSQHNU_SM,DJDt_TERM,
	1            HFLUX_AT_IB,HFLUX_AT_OB,
	1            VDOP_VEC,VDOP_FRAC,FREQ,dLOG_NU,DBB,
	1            INNER_BND_METH,OUTER_BND_METH,
	1            METHOD,COHERENT,INIT,NEW_FREQ,
	1            DO_TIME_VAR,USE_DR4JDT,RELAX_PARAM,NC,NP,ND_SM,NCF)
	USE MOD_MOM_J_DDT_V3
	USE MOD_RAY_MOM_STORE
	IMPLICIT NONE
!
! Altered: 21-Apr-2016 : Changed to V3: USE_DR4JDT option installed.
!                          Old method can still be used.
! Altered 17-Feb-2015 : Check that |r^2.H| < r^2.J
!                         Modified error output to MOM_J_ERRORS
!                         [OSPREY/cur_cmf_gam: 24-Jan-2015]
! Altered 31-Jan-2010 : Cleaned (removed write statements).
!                         For ZERO_FLUX option we add an additional term containing
!                         J computed by the ray-ray solution to improve stability.
! Altered 18-Jan-2010 : More flexible options installed in CALL to allow improved
!                         boundary conditions. Routine now handles "HOLLOW" and 
!                         "ZERO_FLUX" option at inner boundary, and "HALF_MOMJ" at outer
!                         boundary. Other options could also easily be installed.
!                         Based on MOM_J_DDT_V1 (and a modified test routine).
! Altered 31-July-2006: Routine now only deallocates/allocates arrays when 
!                         absolutely necessary. 
!
	INTEGER NC
	INTEGER NP
	INTEGER ND_SM
	INTEGER NCF
!
! _SM denotes the original grid as supplied, for example, by CMFGEN.
! This routine has the option of using a finer grid for calculating
! the radiation field (primarily for use with CMF_FLUX).
!
	REAL*8 ETA_SM(ND_SM)
	REAL*8 CHI_SM(ND_SM)
	REAL*8 ESEC_SM(ND_SM)
	REAL*8 V_SM(ND_SM)
	REAL*8 R_SM(ND_SM)
!
! Radiation field variables. F supplied. JNU_SM and RSQHNU_SM recomputed.
!
	REAL*8 F_SM(ND_SM)
	REAL*8 JNU_SM(ND_SM)
	REAL*8 RSQHNU_SM(ND_SM)
	REAL*8 DJDt_TERM(ND_SM)
	REAL*8 VDOP_VEC(ND_SM)
	REAL*8 VDOP_FRAC
!
	INTEGER N_ERR_MAX,MOM_ERR_CNT
	PARAMETER (N_ERR_MAX=1000)
	REAL*8 MOM_ERR_ON_FREQ
	COMMON /MOM_J_CMF_ERR/MOM_ERR_ON_FREQ(N_ERR_MAX),MOM_ERR_CNT
	LOGICAL RECORDED_ERROR
	LOGICAL, PARAMETER :: VERBOSE=.FALSE.
!
! Boundary conditions.
!
	REAL*8 HFLUX_AT_OB
	REAL*8 HFLUX_AT_IB
!
	REAL*8 RELAX_PARAM
	REAL*8 DBB
	REAL*8 FREQ
	REAL*8 dLOG_NU
	CHARACTER*6 METHOD
!
	CHARACTER(LEN=*) OUTER_BND_METH
	CHARACTER(LEN=*) INNER_BND_METH
!
! INIT is used to indicate that there is no coupling to the previous frequency.
! We are thus solving the normal continuum transfer equation (i.e. the absence
! of velocity fields)
!
! NEW_FREQ is used to indicate that we are computing J for a new frequency.
! If we were iterating between computing J and the Eddington factors, NEW_FREQ
! would be set to false.
!
! COHERENT indicates whether the scattering is coherent. If it is, we
! explicitly take it into account. If COHERENT is FALSE, any electron
! scattering term should be included directly in the ETA that is passed
! to the routine.
!
! When TRUE, DO_TIME_VAR indicates that the Lagrangian time derivative should
! be included into the Radiation Transfer Equation.
!
	LOGICAL DIF
	LOGICAL INIT
	LOGICAL COHERENT
	LOGICAL NEW_FREQ
	LOGICAL DO_TIME_VAR
	LOGICAL USE_DR4JDT
!
	REAL*8 SPEED_OF_LIGHT
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	EXTERNAL SPEED_OF_LIGHT
!
! Local variables.
!
	REAL*8 RSQH_AT_OB
	REAL*8 RSQH_AT_IB
	REAL*8 HONJ_OUTBC
	REAL*8 HPLUS,HMIN
	REAL*8 FMIN,FPLUS
	REAL*8 T1,T2,T3,T4
	REAL*8 MOD_DTAU
	REAL*8 DELTA_R
	INTEGER IT1,I,J,K
	INTEGER IOS
	LOGICAL NEW_R_GRID
	CHARACTER(LEN=6), PARAMETER :: OPTION='NORMAL'
!
!
	LUER=ERROR_LU()
	NEW_R_GRID=.FALSE.
	IF(ND_SM .NE. ND_SM_SAV)THEN
	  NEW_R_GRID=.TRUE.
	  IF(.NOT. FIRST_TIME)THEN
	    WRITE(LUER,*)'Updating R grid in MOM_J_DDT_V2 as ND_SM has changed'
	  END IF
	  ELSE IF(ALLOCATED(R))THEN
!
! We only use R_SM_SAV here. We do the comparison with R_SM_SAV, rather than
! LOG_R_SM, because it avoids precision problems.
!
	  DO I=1,ND_SM
	    IF(R_SM(I) .NE. R_SM_SAV(I))THEN
	      NEW_R_GRID=.TRUE.
	      WRITE(LUER,*)'Updating R grid in MOM_J_DDT_V2 as R grid has changed'
	      EXIT
	    END IF
	  END DO
          IF(VDOP_FRAC .NE. VDOP_FRAC_SAV)THEN
	    NEW_R_GRID=.TRUE.
	    WRITE(LUER,*)'Updating R grid in MOM_J_DDT_V2 as VDOP_FRAC has changed'
	  END IF
	END IF
	IF(NEW_R_GRID .AND. .NOT. INIT)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in MOM_J_DDT_V2'
	  WRITE(LUER,*)'NEW_R_GRID is TRUE but INIT was set to false'
	  WRITE(LUER,*)'Major logic problem'
	  STOP
	END IF
!
	IF(FIRST_TIME)THEN
          OPEN(UNIT=47,FILE='MOM_J_ERRORS',STATUS='UNKNOWN')
	  WRITE(LUER,*)'The value of USE_DR4JDT is',USE_DR4JDT
	ELSE IF(INIT)THEN
	  REWIND(47)
	END IF
!
! On the very first entry, we define the improved R grid, and allocate all
! data arrays.
!
	IF(NEW_R_GRID)THEN
!
! Determine the number of points for the expanded R grid.
! We always insert an EVEN number of points. This guarantees that
! H_SM (defined at the midpoints of the pass grid) has an exact correspondence
! with H defined on the extended grid.
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
	  IF(ND .NE. ND_SAV)THEN
	    IOS=0; IF(ALLOCATED(R))DEALLOCATE(R,R_PNT,STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error deallocating R & R_PNT in MOM_J_DDT_V2'
	      STOP
	    END IF
	    ALLOCATE ( R(ND) )
	    ALLOCATE ( R_PNT(ND) )
	  END IF
	  VDOP_FRAC_SAV=VDOP_FRAC
!
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
	END IF
!
! Deallocate all arrays if we have changed VDOP_FRAC. This will only
! be done in testing this routine (e.g., using DISPGEN).
!
	IF(ND_SM .NE. ND_SM_SAV .AND. .NOT. FIRST_TIME)THEN
	  DEALLOCATE ( LOG_R_SM, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( R_SM_SAV, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( V_COEF, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( ETA_COEF, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( ESEC_COEF, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( CHI_COEF, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( F_COEF, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( J_INDX, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( H_INDX, STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error deallocating R_SM_AV, V_COEF, etc in MOM_J_DDT_V2'
	    STOP
	  END IF
	END IF
!
	IF(ND .NE. ND_SAV .AND. .NOT. FIRST_TIME)THEN
	  DEALLOCATE ( V, STAT=IOS )
	  IF(IOS .EQ. 0)DEALLOCATE ( ETA, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( ESEC, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( CHI, STAT=IOS)
!
	  IF(IOS .EQ. 0)DEALLOCATE ( RSQJNU, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( RSQHNU, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( RSQJNU_PREV, STAT=IOS)		!J at previous frequency
	  IF(IOS .EQ. 0)DEALLOCATE ( RSQHNU_PREV, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( RSQJNU_OLDt, STAT=IOS)		!J at previous time step
	  IF(IOS .EQ. 0)DEALLOCATE ( RSQHNU_OLDt, STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error deallocating ETA, ESEC, etc in MOM_J_DDT_V2'
	    STOP
	  END IF
!
	  DEALLOCATE ( F, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( F_SAV, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( F_PREV, STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error deallocating F, F_SAV, etc in MOM_J_DDT_V2'
	    STOP
	  END IF
!
	  DEALLOCATE ( CON_GAM, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( CON_GAMH, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( COH_VEC, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( DJDt, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( dH, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( dH_OLDT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( DTAU, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( DTAUONQ, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( GAM, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( GAMH, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( HU, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( HL, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( HS, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( HT, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( Q, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( SOURCE, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( TA, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( TB, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( TC, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( PSI, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( PSIPREV, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( W, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( WPREV, STAT=IOS)
	  IF(IOS .EQ. 0)DEALLOCATE ( XM, STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error deallocating CON_GAM etc in MOM_J_DDT_V2'
	    STOP
	  END IF
!
	END IF
!
!
! Allocate all arrays if not previously allocated (ND_SM_SAV=0, ND_SAV=0),
! or if grid size has changed.
!
	IF(ND_SM .NE. ND_SM_SAV)THEN
	  ALLOCATE ( R_SM_SAV(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( LOG_R_SM(ND_SM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( V_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ETA_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( ESEC_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( CHI_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE ( F_COEF(ND_SM,4),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(LUER,*)'Unable to allocate COEF memory in MOM_J_DDT_V2'
	     STOP
	  END IF
	  ALLOCATE ( J_INDX(ND_SM) );       J_INDX(1:ND_SM)=0
	  ALLOCATE ( H_INDX(ND_SM) );       H_INDX(1:ND_SM)=0
	  ND_SM_SAV=ND_SM
	END IF
!
! Allocate arrays on the large ND grid.
!
	IF(ND .NE. ND_SAV)THEN
	  ALLOCATE ( V(ND) )
	  ALLOCATE ( ETA(ND) )
	  ALLOCATE ( ESEC(ND) )
	  ALLOCATE ( CHI(ND) )
!
	  ALLOCATE ( RSQJNU(ND) )             ; RSQJNU(1:ND)=0.0D0
	  ALLOCATE ( RSQHNU(ND) )             ; RSQHNU(1:ND)=0.0D0
	  ALLOCATE ( RSQJNU_PREV(ND) )        ; RSQJNU_PREV(1:ND)=0.0D0
	  ALLOCATE ( RSQHNU_PREV(ND) )        ; RSQHNU_PREV(1:ND)=0.0D0
	  ALLOCATE ( RSQJNU_OLDt(ND) )        ; RSQJNU_OLDt(1:ND)=0.0D0
	  ALLOCATE ( RSQHNU_OLDt(ND) )        ; RSQHNU_OLDt(1:ND)=0.0D0
!
	  ALLOCATE ( F(ND) )
	  ALLOCATE ( F_SAV(ND) )              ; F_SAV(1:ND)=0.0D0
	  ALLOCATE ( F_PREV(ND) )
!
	  ALLOCATE ( CON_GAM(ND) )
	  ALLOCATE ( CON_GAMH(ND) )
	  ALLOCATE ( DJDt(ND) )
	  ALLOCATE ( COH_VEC(ND) )
	  ALLOCATE ( dH(ND) )
	  ALLOCATE ( dH_OLDT(ND) )
	  ALLOCATE ( DTAU(ND) )
	  ALLOCATE ( DTAUONQ(ND) )
	  ALLOCATE ( GAM(ND) )
	  ALLOCATE ( GAMH(ND) )
	  ALLOCATE ( HU(ND) )
	  ALLOCATE ( HL(ND) )
	  ALLOCATE ( HS(ND) )
	  ALLOCATE ( HT(ND) )
	  ALLOCATE ( Q(ND) )
	  ALLOCATE ( SOURCE(ND) )
	  ALLOCATE ( TA(ND) )
	  ALLOCATE ( TB(ND) )
	  ALLOCATE ( TC(ND) )
	  ALLOCATE ( PSI(ND) )
	  ALLOCATE ( PSIPREV(ND) )
	  ALLOCATE ( W(ND) )
	  ALLOCATE ( WPREV(ND) )
	  ALLOCATE ( XM(ND) )
	  ND_SAV=ND
!
	END IF
!
!
	IF(NEW_R_GRID)THEN
	  R_SM_SAV(1:ND_SM)=R_SM(1:ND_SM)
	  LOG_R_SM(1:ND_SM)=LOG(R_SM(1:ND_SM))
          CALL MON_INT_FUNS_V2(V_COEF,V_SM,LOG_R_SM,ND_SM)
          DO I=1,ND
            K=R_PNT(I)
            T1=LOG(R(I)/R_SM(K))
            V(I)=((V_COEF(K,1)*T1+V_COEF(K,2))*T1+V_COEF(K,3))*T1+V_COEF(K,4)
          END DO
!
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
! Note that V is in km/s. The factor of 2 in CON_GAMH is an allowance
! for a division by 0.5(CHI(I)+CHI(J)).
!
	  C_KMS=1.0D-05*SPEED_OF_LIGHT()
	  DO I=1,ND-1
	    CON_GAMH(I)=2.0D0*(V(I)+V(I+1))/(R(I)+R(I+1))/C_KMS
	    CON_GAM(I)=V(I)/R(I)/C_KMS
	  END DO
	  CON_GAM(ND)=V(ND)/R(ND)/C_KMS
!
	END IF
!
	IF(INIT)THEN
	  DO I=1,N_ERR_MAX
	    MOM_ERR_ON_FREQ(I)=0.0D0
	  END DO
	  MOM_ERR_CNT=0
	END IF
!
!
!
! We have no finished all initializations that need to be done on
! the first entry, or on the initial entry of a new frequency sequence.
! We can store the new opacities/emissivities for the upcoming calculation.
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
	  END DO
!
	ELSE
	  ESEC(1:ND)=ESEC_SM(1:ND)
	  CHI(1:ND)=CHI_SM(1:ND)
	  ETA(1:ND)=ETA_SM(1:ND)
	  F(1:ND)=F_SM(1:ND)
	END IF
!
!*****************************************************************************
!
	SOURCE(:)=ETA(:)/CHI(:)
	IF(COHERENT .AND. USE_DR4JDT)THEN
	  COH_VEC(:)=(ESEC(:)+V(:)/R(:)/C_KMS)/CHI(:)
	ELSE IF(COHERENT)THEN
	  COH_VEC(:)=ESEC(:)/CHI(:)
	ELSE IF(USE_DR4JDT)THEN
	  COH_VEC(:)=V(:)/R(:)/C_KMS/CHI(:)
	ELSE
	  COH_VEC(:)=0.0D0
	END IF
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
!*************************************************
! ******** FUDGE ****
!*************************************************
!
	IF(DTAU(ND-1) .LT. 1.0D-05)DTAU(ND-1)=1.0D-05
!
	IF(NEW_FREQ)THEN
	  F_PREV(1:ND)=F_SAV(1:ND)
	  RSQJNU_PREV(1:ND)=RSQJNU(1:ND)
	  RSQHNU_PREV(1:ND)=RSQHNU(1:ND)
	  HONJ_OUTBC_PREV=HONJ_OUTBC_SAV
	  RSQH_AT_IB_PREV=RSQH_AT_IB_SAV
	  RSQH_AT_OB_PREV=RSQH_AT_OB_SAV
	  IF(DO_TIME_VAR)THEN
	    CALL GET_JH_AT_PREV_TIME_STEP(RSQJNU_OLDt,RSQHNU_OLDt,
	1        RSQH_AT_IB_OLDT,HONJ_OUTBC_OLDT,DELTA_TIME_SECS,
	1        FREQ,R,V,ND,INIT,OPTION)
	  END IF
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
	    IF(FIRST_TIME)THEN
	       WRITE(LUER,'(4X,A,ES12.4)')' RECIP_CDELTAT=',RECIP_CDELTAT
	       WRITE(LUER,'(4X,A,ES12.4)')'     ROLD_ON_R=',ROLD_ON_R
	    END IF
	  ELSE
	    RECIP_CDELTAT=0.0D0
	    ROLD_ON_R=0.0D0
	  END IF
	END IF
	R_RAT_FOR_J=ROLD_ON_R
	IF(USE_DR4JDT)R_RAT_FOR_J=ROLD_ON_R*ROLD_ON_R
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
	    dH(I)=0.0D0
	    dH_OLDT(I)=0.0D0
	    DJDt(I)=0.0D0
	  END DO
	  HONJ_OUTBC_PREV=0.0D0;  RSQH_AT_IB_PREV=0.0D0; RSQH_AT_OB_PREV=0.0D0
	ELSE
!
!
! Since we are integrating from blue to red, FL_PREV is always larger than
! FL. dLOG_NU is define as vd / dv which is the same as d / d ln v.
!
	  DO I=1,ND-1
	    GAMH(I)=CON_GAMH(I)/dLOG_NU/( CHI(I)+CHI(I+1) )
	    dH(I)=2.0D0*RECIP_CDELTAT/( CHI(I)+CHI(I+1) )
	    dH_OLDT(I)=dH(I)*ROLD_ON_R
	    W(I)=GAMH(I)+dH(I)
	  END DO
	  GAM(:)=CON_GAM(:)/CHI(:)/dLOG_NU
	END IF
!
! Even if it is the first frequency, we still need to allow for the time 
! variability terms.
!
	IF(INIT .AND. DO_TIME_VAR)THEN
	  DO I=1,ND-1
	    dH(I)=2.0D0*RECIP_CDELTAT/( CHI(I)+CHI(I+1) )
	    dH_OLDT(I)=dH(I)*ROLD_ON_R
	    W(I)=dH(I)
	  END DO
	END IF
!
! 
!
	DO I=2,ND-1
	  DTAUONQ(I)=0.5D0*(DTAU(I)+DTAU(I-1))/Q(I)
!	  DTAUONQ(I)=0.5D0*(R(I-1)-R(I+1))*CHI(I)
	  PSI(I)=DTAUONQ(I)*GAM(I)
	  PSIPREV(I)=DTAUONQ(I)*GAM(I)
	  DJDt(I)=DTAUONQ(I)*RECIP_CDELTAT/CHI(I)
	END DO
!
! Compute vectors used to compute the flux vector H. I could have
! redefined W as 1+W, but keep the current definition for consistency
! with older versions, and with the compatible VAR_MOM routine.
!
	DO I=1,ND-1
	  HU(I)=F(I+1)*Q(I+1)/(1.0D0+W(I))/DTAU(I)
	  HL(I)=F(I)*Q(I)/(1.0D0+W(I))/DTAU(I)
	  HS(I)=GAMH(I)/(1.0D0+W(I))
	  HT(I)=dH_OLDT(I)/(1.0D0+W(I))
	END DO
!
! Compute the TRIDIAGONAL operators, and the RHS source vector.
!
	DO I=2,ND-1
	  TA(I)=-HL(I-1)
	  TC(I)=-HU(I)
	  TB(I)=DTAUONQ(I)*(1.0D0-COH_VEC(I)) + PSI(I) + DJDt(I) + HL(I) + HU(I-1)
	  XM(I)=DTAUONQ(I)*SOURCE(I)*R(I)*R(I)
	END DO
!	I=2
!	WRITE(401,'(ES14.7,2ES22.14,3ES14.7)')FREQ,HL(I),HU(I),DTAUONQ(I)*(1.0D0-COH_VEC(I)),PSI(I),DJDt(I)
!	I=ND-1
!	WRITE(402,'(ES14.7,2ES22.14,3ES14.7)')FREQ,HL(I),HU(I),DTAUONQ(I)*(1.0D0-COH_VEC(I)),PSI(I),DJDt(I)
!
	DO I=2,ND-1
	  XM(I)=XM(I) + 
	1        (HS(I)*RSQHNU_PREV(I) - HS(I-1)*RSQHNU_PREV(I-1)) + 
	1        (HT(I)*RSQHNU_OLDT(I) - HT(I-1)*RSQHNU_OLDT(I-1)) + 
	1        PSIPREV(I)*RSQJNU_PREV(I) + DJDt(I)*RSQJNU_OLDt(I)*R_RAT_FOR_J
	END DO
!
! Evaluate TA,TB,TC for boundary conditions. These will automatically
! handle INIT=.TRUE. (since GAM=0) and DO_TIME_VR=.FALSE. (since
! RECIP_DELTA=0)
!
	TA(1)=0.0D0
	DJDt(1)=RECIP_CDELTAT/CHI(1)
	IF(INIT)DJDt(1)=0.0D0
	IF(OUTER_BND_METH .EQ. 'HONJ')THEN
	  HONJ_OUTBC=(HPLUS_OB-HMIN_OB)/(JPLUS_OB+JMIN_OB)
	  PSI(1)=GAM(1)*HONJ_OUTBC
	  PSIPREV(1)=GAM(1)*HONJ_OUTBC_PREV
	  TC(1)=-F(2)*Q(2)/DTAU(1)
	  TB(1)=F(1)*Q(1)/DTAU(1) + HONJ_OUTBC*(1.0D0+GAM(1)+DJDt(1))
	  XM(1)=PSIPREV(1)*RSQJNU_PREV(1) + DJDt(1)*RSQJNU_OLDt(1)*ROLD_ON_R*HONJ_OUTBC_OLDT
!
!	  WRITE(156,'(5ES14.5,2E30.16)')FREQ,HONJ_OUTBC,GAM(1),DJDt(1),(TB(1)+TC(1))/TB(1),TB(1),TC(1)
!
	ELSE IF(OUTER_BND_METH .EQ. 'HALF_MOM')THEN
	  MOD_DTAU=0.5D0*(CHI(1)+CHI(2))*(R(1)-R(2))
	  T1=R(1)*R(1)
	  HPLUS=HPLUS_OB/JPLUS_OB
	  FPLUS=KPLUS_OB/JPLUS_OB
	  TC(1)=-F(2)/MOD_DTAU
	  TB(1)=FPLUS/MOD_DTAU -  (1.0D0-FPLUS)/R(1)/CHI(1) + HPLUS*(1.0D0+GAM(1)+DJDT(1))
	  XM(1)=T1*( HMIN_OB - KMIN_OB/MOD_DTAU+(JMIN_OB-KMIN_OB)/R(1)/CHI(1) ) +
	1       GAM(1)*(T1*HMIN_OB+RSQH_AT_OB_PREV) +
	1       DJDT(1)*(T1*HMIN_OB+RSQJNU_OLDt(1)*ROLD_ON_R*HONJ_OUTBC_OLDT)
	  XM(2)=XM(2)-TA(2)*JMIN_OB*T1
	ELSE
	  I=ERROR_LU()
          WRITE(I,*)'Only HONJ & HALF_MD boundary conditions implemented at outer boundary'
          WRITE(I,*)'Routine is MOM_J_DDT_V2'
          STOP
	END IF
!
! ***  INNER BOUNDARY CONDITION ****
!
	PSI(ND)=0.0D0
	PSIPREV(ND)=0.0D0
	DJDt(ND)=0.0D0
!
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  RSQH_AT_IB=DBB*R(ND)*R(ND)/3.0D0/CHI(ND)
	  TB(ND)=F(ND)/DTAU(ND-1)
	  TA(ND)=-F(ND-1)*Q(ND-1)/DTAU(ND-1)
	  XM(ND)=RSQH_AT_IB+RECIP_CDELTAt*(RSQH_AT_IB-ROLD_ON_R*RSQH_AT_IB_OLDt)/CHI(ND)
!
!	1              +GAM(ND)*(RSQH_AT_IB-RSQH_AT_IB_PREV)
!
! Boundry condition assumes a zero-flux conditiont at the the inner boundary. For
! this routine (but NOT the ray soluton), it is eqivalent to the diffussion approximation
! with DBB=0.0D0
!
	ELSE IF(INNER_BND_METH .EQ. 'ZERO_FLUX')THEN
	  RSQH_AT_IB=0.0D0
	  TA(ND)=-F(ND-1)*Q(ND-1)/DTAU(ND-1)
	  TB(ND)=F(ND)/DTAU(ND-1)
	  XM(ND)=RSQH_AT_IB+RECIP_CDELTAt*(RSQH_AT_IB-ROLD_ON_R*RSQH_AT_IB_OLDt)/CHI(ND)
	  XM(ND)=XM(ND)+0.10D0*F(ND)*R(ND)*R(ND)*(JPLUS_IB+JMIN_IB)/DTAU(ND-1)
	  TB(ND)=TB(ND)+0.10D0*F(ND)/DTAU(ND-1)
!
! With this boundary condition we specify JPLUS, HPLUS, and KPLUS at the
! inner boundary. These, in general, will be dependent on the radiation
! field previously computed at higher freqencies (as the inner
! core is expanding, and it sees the other side of the hollow core moving
! away).
!
	ELSE IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	  T1=R(ND)*R(ND)
	  HMIN=HMIN_IB/JMIN_IB
	  FPLUS=KPLUS_IB/JPLUS_IB
	  FMIN=KMIN_IB/JMIN_IB
	  TA(ND)=-F(ND-1)/DTAU(ND-1)
	  TB(ND)=FMIN/DTAU(ND-1) + (1.0D0-FMIN)/R(ND)/CHI(ND) + HMIN*(1.0D0+DjdT(ND)+GAM(ND))
	  XM(ND)=T1*( HPLUS_IB-FPLUS*JPLUS_IB/DTAU(ND-1)-(1.0D0-FPLUS)/R(ND)/CHI(ND)*JPLUS_IB )
	1              +GAM(ND)*(T1*HPLUS_IB-RSQH_AT_IB_PREV)
	1              +RECIP_CDELTAt*(T1*HPLUS_IB-ROLD_ON_R*RSQH_AT_IB_OLDt)/CHI(ND)
	  XM(ND-1)=XM(ND-1)-T1*JPLUS_IB*TC(ND-1)
!
! Modification designed to improve stability.
!
          TB(ND)=TB(ND)+0.1D0*FMIN/DTAU(ND-1)
          XM(ND)=XM(ND)+0.1D0*FMIN*R(ND)*R(ND)*JMIN_IB/DTAU(ND-1)
!
	ELSE
          I=ERROR_LU()
          WRITE(I,*)'Only DIFFUSION, HOLLOW, and ZERO_FLUX boundary conditions are currently implemented'
          WRITE(I,*)'Routine is MOM_J_DDT_V2'
          STOP
	END IF
	TC(ND)=0.0D0
!
! Solve for the radiation field along ray for this frequency.
!
	CALL THOMAS(TA,TB,TC,XM,ND,1)
!
! Correct XM at boundary if we have only comuted a half-mment.
!
	IF(INNER_BND_METH .EQ. 'HOLLOW')THEN
	  RSQH_AT_IB=R(ND)*R(ND)*HPLUS_IB-HMIN*XM(ND)
	  XM(ND)=XM(ND)+R(ND)*R(ND)*JPLUS_IB
	END IF
	IF(OUTER_BND_METH .EQ. 'HALF_MOM')THEN
	  RSQH_AT_OB=HPLUS*XM(1)-HMIN_OB*R(1)*R(1)
	  XM(1)=XM(1)+JMIN_OB*R(1)*R(1)
	ELSE
	  RSQH_AT_OB=HONJ_OUTBC*XM(1)
	END IF
!
! Check that no negative mean intensities have been computed.
!
	IF(MINVAL(XM(1:ND)) .LE. 0.0D0)THEN
	  IF(VERBOSE)THEN
	    WRITE(47,*)'Freq=',FREQ
	    TA(1:ND)=XM(1:ND)/R(1:ND)/R(1:ND)
	    CALL WRITE_VEC(TA,ND,'XM Vec',47)
	    CALL WRITE_VEC(F,ND,'F Vec',47)
	    CALL WRITE_VEC(ETA,ND,'ETA Vec',47)
	    CALL WRITE_VEC(ESEC,ND,'ESEC Vec',47)
	    CALL WRITE_VEC(CHI,ND,'CHI Vec',47)
	  ELSE
	    DO I=1,ND
	     IF(XM(I) .LE. 0.0D0)THEN
	       WRITE(47,'(I5,ES16.8,10ES13.4)')I,FREQ,XM(I),ETA(I),CHI(I),ESEC(I),F(I),XM(MAX(1,I-2):MIN(I+2,ND))
	     END IF
	    END DO
	  END IF
	END IF
!
	RECORDED_ERROR=.FALSE.
	DO I=1,ND
	  IF(XM(I) .LT. 0.0D0)THEN
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
	DO I=1,ND-1
	  RSQHNU(I)=HU(I)*XM(I+1)-HL(I)*XM(I)+HS(I)*RSQHNU_PREV(I)+HT(I)*RSQHNU_OLDT(I)
	END DO
!
! Make sure H satisfies the basic requirement that it is less than J.
!
	DO I=1,ND-1
	  T1=(XM(I)+XM(I+1))/2.0D0
	  IF(RSQHNU(I) .GT. T1)THEN
	    RSQHNU(I)=0.99D0*T1
	  ELSE IF(RSQHNU(I) .LT. -T1)THEN
	    RSQHNU(I)=-0.99D0*T1
	  END IF
	END DO
!
	CALL OUT_JH(RSQJNU,RSQHNU,RSQH_AT_IB,HONJ_OUTBC,FREQ,NCF,R,V,ND,INIT,OPTION)
!
! Re-grid derived J and RSQH values onto small grid. We divide RSQJ by R^2 so that
! we return J. We also compute DJDt on the small grid. This is returned to CMFGEN_SUB
! so that we can allow for the DJDt term when doing the luminosity check for 
! OBSFLUX. Note that the term is Dr^3J/Dt . 1/rc .
!
	DO I=1,ND_SM
	  K=J_INDX(I)
	  JNU_SM(I)=RSQJNU(K)/R_SM(I)/R_SM(I)
	  DJDt_TERM(I)=RELAX_PARAM*RECIP_CDELTAT*(RSQJNU(K)-R_RAT_FOR_J*RSQJNU_OLDT(K))
	END DO
!
	DO I=1,ND_SM-1
	  K=H_INDX(I)
	  RSQHNU_SM(I)=RSQHNU(K)
	END DO
!
	F_SAV(1:ND)=F(1:ND)
	HONJ_OUTBC_SAV=HONJ_OUTBC
	RSQH_AT_IB_SAV=RSQH_AT_IB
	RSQH_AT_OB_SAV=RSQH_AT_OB
	FIRST_TIME=.FALSE.
!
	HFLUX_AT_IB=RSQH_AT_IB/R(ND)/R(ND)
	HFLUX_AT_OB=RSQH_AT_OB/R(1)/R(1)
!
	RETURN
	END
