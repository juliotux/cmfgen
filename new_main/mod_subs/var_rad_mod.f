!
! This moduls declares vectors and arrays required to compute dJ in the continuum section of
! the code. Storage is allocated by a call to SET_VAR_RAD_MOD_V2.
!
	MODULE VAR_RAD_MOD
!
	REAL*8, ALLOCATABLE :: DIFFW(:)         !NT - Variation of diffusion approx. at inner boundary.
!
	REAL*8, ALLOCATABLE :: VK(:,:)    	!ND,ND - Coef. matrix of %CHI vector
	REAL*8, ALLOCATABLE :: FC(:,:)    	!ND,ND - Coef. of %EMIS vector in angular equ.
	REAL*8, ALLOCATABLE :: F2DA(:,:)    	!ND,ND - Coef. of %CHi in angular equ.
	REAL*8, ALLOCATABLE :: FA(:)    		!ND 
!
	REAL*8, ALLOCATABLE :: F2DAEXT(:,:)    	!NDMAX,NDMAX 
	REAL*8, ALLOCATABLE :: FCEXT(:,:)    	!NDMAX,NDMAX 
	REAL*8, ALLOCATABLE :: FAEXT(:)    	!NDMAX- 
!
! Variation arrays
! Variable,depth of variable,depth of J. If NUM_BNDS .ne. ND the
! the variable depth is given by [ VJ(I1,I2,I3) ] I3+I2-NDIAG.
!
	REAL*8, ALLOCATABLE :: VJ(:,:,:)    	!NT,NUM_BNDS,ND - 
!
! Variation line arrays
!
	REAL*8, ALLOCATABLE :: TX(:,:,:)    	!ND,ND,NM - 
	REAL*8, ALLOCATABLE :: TVX(:,:,:)    	!ND-1,ND,NM - 
!
! We make TX_EXT and TVX_EXT allocatable as they are accessed directly
! in VARCONT and thus must have the correct dimensions.
!
	REAL*8, ALLOCATABLE :: TX_EXT(:,:,:)
	REAL*8, ALLOCATABLE :: TVX_EXT(:,:,:)
	REAL*8, ALLOCATABLE :: KI(:,:,:)    		  !NDMAX,ND,NM_KI - 
!
        REAL*8, ALLOCATABLE :: dJ_LOC(:,:,:)              !NM,NUM_BNDS,ND
        REAL*8, ALLOCATABLE :: dZ(:,:,:,:)                !NM,NUM_BNDS,ND,MAX_SIM
        REAL*8, ALLOCATABLE :: dZ_POPS(:,:,:)             !NT,NUM_BNDS,ND
!
	REAL*8, ALLOCATABLE :: dJ_DIF_d_T_EXT(:)          !NDMAX -
	REAL*8, ALLOCATABLE :: dJ_DIF_d_dTdR_EXT(:)       !NDMAX -
	REAL*8, ALLOCATABLE :: dJ_DIF_d_T(:)              !NDMAX -
	REAL*8, ALLOCATABLE :: dJ_DIF_d_dTdR(:)           !NDMAX -
	REAL*8, ALLOCATABLE :: RHS_dHdCHI(:,:)            !NDMAX,ND -
	REAL*8, ALLOCATABLE :: dRSQH_DIF_d_T(:)           !NDMAX -
	REAL*8, ALLOCATABLE :: dRSQH_DIF_d_dTdR(:)        !NDMAX -
!
	END MODULE VAR_RAD_MOD
!
! 
! Subroutine to allocate vectors and arrays.
!
	SUBROUTINE SET_VAR_RAD_MOD_V2(ND,NDEXT,NT,NUM_BNDS,NM,MAX_SIM,NM_KI,
	1                ACCURATE,ALLOCATE_TX)
	USE VAR_RAD_MOD
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NDEXT
	INTEGER NT
	INTEGER NUM_BNDS
	INTEGER NM
	INTEGER NM_KI
	INTEGER MAX_SIM
	LOGICAL ACCURATE
	LOGICAL ALLOCATE_TX
!
	INTEGER IOS
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	IOS=0
	IF(IOS .EQ. 0)ALLOCATE( DIFFW(NT),STAT=IOS)       !Diffusion variation
	IF(IOS .EQ. 0)ALLOCATE( VK(ND,ND),STAT=IOS)       !Coef. matrix of %CHI vector
	IF(IOS .EQ. 0)ALLOCATE( FC(ND,ND),STAT=IOS)       !Coef. of %EMIS vector in angular equ.
	IF(IOS .EQ. 0)ALLOCATE( F2DA(ND,ND),STAT=IOS)     !Coef. of %CHi in angular equ.
	IF(IOS .EQ. 0)ALLOCATE( FA(ND),STAT=IOS)          !
!
	IF(ACCURATE)THEN
	  IF(IOS .EQ. 0)ALLOCATE( F2DAEXT(NDEXT,NDEXT),STAT=IOS)  	!These arrays don't need to be
	  IF(IOS .EQ. 0)ALLOCATE( FCEXT(NDEXT,NDEXT),STAT=IOS)    	!NDEXT,NDEXT - contiguous as for PERTJD.
	  IF(IOS .EQ. 0)ALLOCATE( FAEXT(NDEXT), STAT=IOS)
	END IF
!
! Variation arrays
! Variable,depth of variable,depth of J. If NUM_BNDS .ne. ND the
! the variable depth is given by [ VJ(I1,I2,I3) ] I3+I2-NDIAG.
!
	IF(IOS .EQ. 0)ALLOCATE( VJ(NT,NUM_BNDS,ND), STAT=IOS)
!
! Variation line arrays
!
	IF(ALLOCATE_TX)THEN
	  IF(IOS .EQ. 0)ALLOCATE( TX(ND,ND,NM),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE( TVX(ND-1,ND,NM),STAT=IOS)
!
! We make TX_EXT and TVX_EXT allocatable as they are accessed directly
! in VARCONT and thus must have the correct dimensions.
!
	  IF(ACCURATE)THEN
	    IF(IOS .EQ. 0)ALLOCATE( TX_EXT(NDEXT,ND,NM), STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE( TVX_EXT(NDEXT-1,ND,NM), STAT=IOS)
	  END IF
	END IF
!
! KI is assumed to be of dimension:
!                         (ND,3,NM_KI) in VAR_FORMSOL (NM >= 4)
!                         (ND,ND,NM_KI) in VAR_MOMHAM (NM >= 4)
!                         (ND,ND,NM_KI) in VAR_MOM_J_CMF_V6 (NM >= 2)
!
	IF(IOS .EQ. 0)ALLOCATE(KI(NDEXT,NDEXT,NM_KI), STAT=IOS)
!
        IF(IOS .EQ. 0)ALLOCATE( dJ_LOC(NM,NUM_BNDS,ND), STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE( dZ(NM,NUM_BNDS,ND,MAX_SIM), STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE( dZ_POPS(NT,NUM_BNDS,ND), STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_T_EXT(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_dTdR_EXT(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_T(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dJ_DIF_d_dTdR(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( RHS_dHdCHI(NDEXT,ND),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dRSQH_DIF_d_T(NDEXT),STAT=IOS )
	IF(IOS .EQ. 0)ALLOCATE ( dRSQH_DIF_d_dTdR(NDEXT),STAT=IOS )
!
	IF(IOS .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in SET_VAR_RAD_MOD_V2'
	  WRITE(LUER,*)'Unable to allocate required memory'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
	END
