!
! Module for handling lines. Used for blanketing mode, Sobolev mode,
! and CMF mode.
!
	MODULE LINE_MOD
!
! Incorporated 02-Jan-2013 : Altered to allow depth dependent profiles (from cur_cmf_25jun13).
! Arrays and variables for treating lines simultaneously.
!
	REAL*8,  ALLOCATABLE :: EINA(:)                !MAX_SIM
	REAL*8,  ALLOCATABLE :: OSCIL(:)               !MAX_SIM
	REAL*8,  ALLOCATABLE :: GLDGU(:)               !MAX_SIM
	REAL*8,  ALLOCATABLE :: AMASS_SIM(:)           !MAX_SIM
	REAL*8,  ALLOCATABLE :: FL_SIM(:)              !MAX_SIM
	INTEGER, ALLOCATABLE :: SIM_NL(:)              !MAX_SIM
	INTEGER, ALLOCATABLE :: SIM_NUP(:)             !MAX_SIM
	LOGICAL, ALLOCATABLE :: WEAK_LINE(:)           !MAX_SIM
!
	REAL*8,  ALLOCATABLE :: CHIL_MAT(:,:)          !ND,MAX_SIM
	REAL*8,  ALLOCATABLE :: ETAL_MAT(:,:)          !ND,MAX_SIM
	REAL*8,  ALLOCATABLE :: BB_COR(:,:)            !ND,MAX_SIM
!
	REAL*8,  ALLOCATABLE :: VB_SIM(:,:)            !ND,MAX_SIM
	REAL*8,  ALLOCATABLE :: VC_SIM(:,:)                            !ND,MAX_SIM)
	REAL*8,  ALLOCATABLE :: VB_2(:)                !ND
	REAL*8,  ALLOCATABLE :: VC_2(:)                !ND
!
	REAL*8,  ALLOCATABLE :: BETA(:)                !ND : Sobolev escape probability
	REAL*8,  ALLOCATABLE :: BETAC(:)               !ND : Dilution factor weighted "escape probaility"
!
	REAL*8,  ALLOCATABLE :: BETAC_SIM(:,:)         !ND,MAX_SIM
	REAL*8,  ALLOCATABLE :: ZNET_SIM(:,:)          !ND,MAX_SIM
	REAL*8,  ALLOCATABLE :: JBAR_SIM(:,:)          !ND,MAX_SIM
!
	CHARACTER(LEN=50), ALLOCATABLE :: TRANS_NAME_SIM(:)          !MAX_SIM
!
	INTEGER NUM_SIM_LINES
	INTEGER SIM_INDX
	INTEGER TMP_MAX_SIM
!
! Temporary variables used only local to compute quantities associated with
! the opacity and emissivity, and the rate equations.
!
	REAL*8 OPAC_FAC
	REAL*8 EMIS_FAC
	REAL*8 STIM_FAC
	REAL*8 MUL_FAC
	REAL*8 dRATE_dT
	REAL*8 dRATE_dLOW
	REAL*8 dRATE_dUP
	REAL*8 RATE_FAC
C
C These pointers are used to indicate which locations are being used in
C the SIM variation arrays. LOW refers to the lower level of the transition,
C UP the upper level.
C 
	INTEGER, ALLOCATABLE :: LOW_POINTER(:)             !MAX_SIM
	INTEGER, ALLOCATABLE :: UP_POINTER(:)              !MAX_SIM
C
C Used to count how many different lines are using an individual storage
C location.
C
	INTEGER, ALLOCATABLE :: VAR_IN_USE_CNT(:)          !NM
C
C Indicates the variable (e.g. 1 to NT-2) occupying a particular storage
C location in the the variation arrays (e.g. TX)
C
	INTEGER, ALLOCATABLE :: VAR_LEV_ID(:)              !NM
!
	LOGICAL, ALLOCATABLE :: USE_THIS_VAR_MAT(:)          !NM
	LOGICAL, ALLOCATABLE :: IMP_TRANS_VEC(:)             !NM
	LOGICAL THIS_TRANS_IMP
C
C The following matrices are required to handle the treatment of super levels.
C
C ?_STAR_RATIO is defined by:
C
C       [LTE Pop. of uncombined level / LTE Pop. of SUPER level]
C                        
C L refers to the lower level, U to the upper level.
C
	REAL*8, ALLOCATABLE :: L_STAR_RATIO(:,:)    !ND,MAX_SIM
	REAL*8, ALLOCATABLE :: U_STAR_RATIO(:,:)    !ND,MAX_SIM
	REAL*8, ALLOCATABLE :: dL_RAT_dT(:,:)       !ND,MAX_SIM
	REAL*8, ALLOCATABLE :: dU_RAT_dT(:,:)       !ND,MAX_SIM
	REAL*8, ALLOCATABLE :: LOW_OCC_PROB(:)      !ND
C
C Variables, vectors and arrays for treating lines simultaneously with the
C continuum.
C
	REAL*8, ALLOCATABLE :: LINE_PROF_SIM(:,:)             !ND,MAX_SIM
	REAL*8, ALLOCATABLE :: LINE_QW_SIM(:,:)               !ND,MAX_SIM
	REAL*8, ALLOCATABLE :: LINE_QW_SUM(:,:)               !ND,MAX_SIM
	REAL*8, ALLOCATABLE :: NEG_OPAC_FAC(:)                !ND

	REAL*8, ALLOCATABLE :: LINE_OPAC_CON(:)               !MAX_SIM
	REAL*8, ALLOCATABLE :: LINE_EMIS_CON(:)               !MAX_SIM
C
	LOGICAL, ALLOCATABLE :: RESONANCE_ZONE(:)             !MAX_SIM
	LOGICAL, ALLOCATABLE :: END_RES_ZONE(:)               !MAX_SIM
	LOGICAL, ALLOCATABLE :: DO_THIS_TX_MATRIX(:)          !NM
	LOGICAL, ALLOCATABLE :: LINE_STORAGE_USED(:)          !MAX_SIM
	LOGICAL, ALLOCATABLE :: NEW_LINE_STORAGE(:)           !MAX_SIM
C
	INTEGER, ALLOCATABLE ::  SIM_LINE_POINTER(:)          !MAX_SIM
!
	REAL*8, ALLOCATABLE ::  AVE_ENERGY(:)                 !NT
!
	END MODULE LINE_MOD
!
! 
	SUBROUTINE SET_LINE_MOD(ND,NT,MAX_SIM,NM)
	USE LINE_MOD
	IMPLICIT NONE
!
	INTEGER ND,NT,MAX_SIM,NM
	INTEGER IOS
	INTEGER LUER, ERROR_LU
	EXTERNAL ERROR_LU
!
	IOS=0
!
	IF(IOS .EQ. 0)ALLOCATE( EINA(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( OSCIL(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( GLDGU(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( AMASS_SIM(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( FL_SIM(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( SIM_NL(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( SIM_NUP(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( WEAK_LINE(MAX_SIM) ,STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE( CHIL_MAT(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( ETAL_MAT(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( BB_COR(ND,MAX_SIM) ,STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE( VB_SIM(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VC_SIM(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VB_2(ND) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VC_2(ND) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( BETA(ND) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( BETAC(ND) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( BETAC_SIM(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( ZNET_SIM(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( JBAR_SIM(ND,MAX_SIM) ,STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE( TRANS_NAME_SIM(MAX_SIM) ,STAT=IOS)
!
! These pointers are used to indicate which locations are being used in
! the SIM variation arrays. LOW refers to the lower level of the transition,
! UP the upper level.
! 
	IF(IOS .EQ. 0)ALLOCATE( LOW_POINTER(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( UP_POINTER(MAX_SIM) ,STAT=IOS)
!
! Used to count how many different lines are using an individual storage
! location.
!
	IF(IOS .EQ. 0)ALLOCATE( VAR_IN_USE_CNT(NM) ,STAT=IOS)
!
! Indicates the variable (e.g. 1 to NT-2) occupying a particular storage
! location in the the variation arrays (e.g. TX)
!
	IF(IOS .EQ. 0)ALLOCATE( VAR_LEV_ID(NM) ,STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE( USE_THIS_VAR_MAT(NM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( IMP_TRANS_VEC(NM) ,STAT=IOS)
!
! The following matrices are required to handle the treatment of super levels.
!
! ?_STAR_RATIO is defined by:
!
!       [LTE Pop. of uncombined level / LTE Pop. of SUPER level]
!                        
! L refers to the lower level, U to the upper level.
!
	IF(IOS .EQ. 0)ALLOCATE( L_STAR_RATIO(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( U_STAR_RATIO(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( dL_RAT_dT(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( dU_RAT_dT(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( LOW_OCC_PROB(ND) ,STAT=IOS)
!
! Variables, vectors and arrays for treating lines simultaneously with the
! continuum.
!
	IF(IOS .EQ. 0)ALLOCATE( LINE_PROF_SIM(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( LINE_QW_SIM(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( LINE_QW_SUM(ND,MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( NEG_OPAC_FAC(ND) ,STAT=IOS)

	IF(IOS .EQ. 0)ALLOCATE( LINE_OPAC_CON(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( LINE_EMIS_CON(MAX_SIM) ,STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE( RESONANCE_ZONE(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( END_RES_ZONE(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( DO_THIS_TX_MATRIX(NM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( LINE_STORAGE_USED(MAX_SIM) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( NEW_LINE_STORAGE(MAX_SIM) ,STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE( SIM_LINE_POINTER(MAX_SIM), STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE(AVE_ENERGY(NT), STAT=IOS)
!
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Unable to allocate memory in SET_LINE_MOD'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
	RETURN
	END
