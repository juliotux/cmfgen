!
! Subroutine to assign the 2-photon data to the appropriate program
! species.
!
! This routine must be executed for each iteration as the arrays
! FS_RAT_LOW and FS_RAT_UP need to be updated.
!
! This routine should be identical to SET_TWO_PHOT_V3 except we do not
! access HNST_F_ON_S (we assume it unity) and EQSPC. We also set 
! UP_LEV_TWO and LOW_LEV_TWO to the level in the FULL atom -- not in POPS.
!
	SUBROUTINE SET_TWO_PHOT_DISP_V3(SPECIES,ID,
	1            HNST_S,N_S,
	1            HNST_F_ON_S,LEVEL_NAME,EDGE_F,G_F,F_TO_S,N_F,
	1            ND,ZION,EQSPEC,SPECIES_PRESENT)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Altered 19-Aug-2015: Added options to improve two photon absorption (cur_hmi,14-Jul-2015)
! Altered 05-Apr-2011 - Changed to V3.
!                       HNST_F_ON_S (rather than HNST_F) is passed in call.
!                       HNST_F/HNST_S replaced by HNST_F_ON_S - done to facilitate
!                         modifications allowing lower temperatures.
! Created 26-Jun-1998
!
	INTEGER N_S
	INTEGER N_F
	INTEGER ND
	INTEGER ID
	INTEGER EQSPEC
!
	REAL*8 HNST_S(N_S,ND)
	REAL*8 HNST_F_ON_S(N_F,ND)
	REAL*8 EDGE_F(N_F)
	REAL*8 G_F(N_F)
	INTEGER F_TO_S(N_F)
!
	LOGICAL SPECIES_PRESENT
!
	CHARACTER*(*) SPECIES
	CHARACTER*(*) LEVEL_NAME(N_F)
!
	REAL*8 GION
	REAL*8 ZION
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local variables
!
	INTEGER I,J
	INTEGER I_S,I_F
!
! Allocate required data vectors. These vectors are conatined and described
! in the data module TWO_PHOT_MOD
!
	IF( .NOT. ALLOCATED(FREQ_TWO))THEN
	  ALLOCATE (FREQ_TWO(N_TWO))
	  ALLOCATE (Z_TWO(N_TWO))
	  ALLOCATE (G_LOW_TWO(N_TWO))
	  ALLOCATE (G_UP_TWO(N_TWO))
	  ALLOCATE (LOW_LEV_TWO(N_TWO))
	  ALLOCATE (UP_LEV_TWO(N_TWO))
	  ALLOCATE (TWO_PHOT_AVAILABLE(N_TWO))
	  ALLOCATE (ION_LOW_LEV_TWO(N_TWO))
	  ALLOCATE (ION_UP_LEV_TWO(N_TWO))
	  ALLOCATE (ION_ID_TWO(N_TWO))
	  ALLOCATE (LST_FREQ_INDX_TWO(N_TWO))
C
	  ALLOCATE (FS_RAT_LOW(ND,N_TWO))
	  ALLOCATE (FS_RAT_UP(ND,N_TWO))
	  ALLOCATE (UP_RATE_TWO(ND,N_TWO))
	  ALLOCATE (DOWN_RATE_TWO(ND,N_TWO))
	  ALLOCATE (PHOT_OC_TWO(ND,N_TWO))
C
	  INITIALIZE_TWO=.TRUE.
	END IF
C
	IF(INITIALIZE_TWO)THEN
	  INITIALIZE_TWO=.FALSE.
	  FREQ_TWO(:)=0.0D0
	  G_LOW_TWO(:)=0.0D0
	  G_UP_TWO(:)=0.0D0
	  LOW_LEV_TWO(:)=0.0D0
	  UP_LEV_TWO(:)=0.0D0
	  ION_LOW_LEV_TWO(:)=0.0D0
	  ION_UP_LEV_TWO(:)=0.0D0
	  ION_ID_TWO(:)=0.0D0
	  Z_TWO(:)=0.0D0
	  FS_RAT_LOW(:,:)=0.0D0
	  FS_RAT_UP(:,:)=0.0D0
	  DOWN_RATE_TWO(:,:)=0.0D0
	  UP_RATE_TWO(:,:)=0.0D0
	  LST_FREQ_INDX_TWO=0
	  TWO_PHOT_AVAILABLE(:)=.FALSE.
	END IF
!
	IF(.NOT. SPECIES_PRESENT)RETURN
!
	DO J=1,N_TWO
!
! Identify lower level of 2-photon transition,
!
	  IF(SPEC_ID_TWO(J) .EQ. SPECIES)THEN
	    Z_TWO(J)=ZION
	    ION_ID_TWO(J)=ID
	    DO I=1,N_F
	       IF( LEVEL_NAME(I) .EQ. LOW_NAME_TWO(J) .OR. 
	1              LEVEL_NAME(I) .EQ. A_LOW_NAME_TWO(J) )THEN
	          I_F=I
	          I_S=I                          !F_TO_S(I)
	          ION_LOW_LEV_TWO(J)=I_S
	          LOW_LEV_TWO(J)=I_S             !EQSPEC+I_S-1
	          FREQ_TWO(J)=EDGE_F(I_F)
	          G_LOW_TWO(J)=G_F(I_F)
	          FS_RAT_LOW(1:ND,J)=1.0D0             !HNST_F_ON_S(I_F,1:ND)
	       END IF
	    END DO
!
! Identify upper level of 2-photon transition
!
	    DO I=1,N_F
	       IF( LEVEL_NAME(I) .EQ. UP_NAME_TWO(J) .OR. 
	1                 LEVEL_NAME(I) .EQ. A_UP_NAME_TWO(J))THEN
	          I_F=I
	          I_S=I                                 !F_TO_S(I)
	          ION_UP_LEV_TWO(J)=I_S
	          UP_LEV_TWO(J)=I_S                     ! EQSPEC+I_S-1
	          G_UP_TWO(J)=G_F(I_F)
	          FREQ_TWO(J)=FREQ_TWO(J)-EDGE_F(I_F)
	          FS_RAT_UP(1:ND,J)=1.0D0               !HNST_F_ON_S(I_F,1:ND)
!
! The following treats the case when the 2s and 2p state are treated as a single level.
!
	          IF(LEVEL_NAME(I) .EQ. '2___' .AND. TWO_PHOT_FORMAT_DATE .EQ. ' ')THEN
	            COEF_TWO(2,J)=2
	            LUER=ERROR_LU()
	            WRITE(LUER,*)'Warning in SET_TWO_PHOT -- adjusting A as working with n=2 level'
	            WRITE(LUER,'(1X,A,2X,A)')TRIM(SPECIES),TRIM(LEVEL_NAME(I))
	          END IF
	          COEF_TWO(1,J)=COEF_TWO(1,J)*COEF_TWO(2,J)/G_UP_TWO(J)
	       END IF
	    END DO
!
! Verify that ordering of level names in 2-photon data file is correct.
!
	    IF(UP_LEV_TWO(J) .LE. LOW_LEV_TWO(J))THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in SET_TWO_PHOT --- invalid level ordering'
	      WRITE(LUER,*)SPEC_ID_TWO(J)
	      WRITE(LUER,*)LOW_LEV_TWO(J),UP_LEV_TWO(J)
	      STOP
	    END IF
	    TWO_PHOT_AVAILABLE(J)=.TRUE.
	  END IF
!
	END DO		!J: Which transition
!
	RETURN
	END
