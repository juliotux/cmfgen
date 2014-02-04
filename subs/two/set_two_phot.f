!
! Subroutine to assign the 2-photon data to the appropriate program
! species.
!
! This routine must be executed for each iteration as the arrays
! FS_RAT_LOW and FS_RAT_UP need to be updated.
!
	SUBROUTINE SET_TWO_PHOT(SPECIES,
	1            HNST_S,N_S,
	1            HNST_F,LEVEL_NAME,EDGE_F,G_F,F_TO_S,N_F,
	1            ND,ZION,EQSPEC,SPECIES_PRESENT)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Created 26-Jun-1998
!
	INTEGER*4 N_S
	INTEGER*4 N_F
	INTEGER*4 ND
	INTEGER*4 EQSPEC
!
	REAL*8 HNST_S(N_S,ND)
	REAL*8 HNST_F(N_F,ND)
	REAL*8 EDGE_F(N_F)
	REAL*8 G_F(N_F)
	INTEGER*4 F_TO_S(N_F)
!
	LOGICAL SPECIES_PRESENT
!
	CHARACTER*(*) SPECIES
	CHARACTER*(*) LEVEL_NAME(N_F)
!
	REAL*8 T(ND)
	REAL*8 GION
	REAL*8 ZION
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER*4 ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local variables
!
	INTEGER*4 I,J
	INTEGER*4 I_S,I_F
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
C
	  ALLOCATE (FS_RAT_LOW(ND,N_TWO))
	  ALLOCATE (FS_RAT_UP(ND,N_TWO))
	  ALLOCATE (UP_RATE_TWO(ND,N_TWO))
	  ALLOCATE (DOWN_RATE_TWO(ND,N_TWO))
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
	  Z_TWO(:)=0.0D0
	  FS_RAT_LOW(:,:)=0.0D0
	  FS_RAT_UP(:,:)=0.0D0
	  DOWN_RATE_TWO(:,:)=0.0D0
	  UP_RATE_TWO(:,:)=0.0D0
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
	    DO I=1,N_F
	       IF( LEVEL_NAME(I) .EQ. LOW_NAME_TWO(J) .OR. 
	1              LEVEL_NAME(I) .EQ. A_LOW_NAME_TWO(J) )THEN
	          I_F=I
	          I_S=F_TO_S(I)
	          LOW_LEV_TWO(J)=EQSPEC+I_S-1
	          FREQ_TWO(J)=EDGE_F(I_F)
	          G_LOW_TWO(J)=G_F(I_F)
	          FS_RAT_LOW(1:ND,J)=HNST_F(I_F,1:ND)/HNST_S(I_S,1:ND)
	       END IF
	    END DO
!
! Identify upper level of 2-photon transition
!
	    DO I=1,N_F
	       IF( LEVEL_NAME(I) .EQ. UP_NAME_TWO(J) .OR. 
	1                 LEVEL_NAME(I) .EQ. A_UP_NAME_TWO(J))THEN
	          I_F=I
	          I_S=F_TO_S(I)
	          UP_LEV_TWO(J)=EQSPEC+I_S-1
	          G_UP_TWO(J)=G_F(I_F)
	          FREQ_TWO(J)=FREQ_TWO(J)-EDGE_F(I_F)
	          FS_RAT_UP(1:ND,J)=HNST_F(I_F,1:ND)/HNST_S(I_S,1:ND)
	       END IF
	    END DO
!
! Verify that ordering of level names in 2-photon data file is correct.
!
	    IF(UP_LEV_TWO(J) .LE. LOW_LEV_TWO(J))THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error in SET_TWO_PHOT --- invalid level ordering'
	      WRITE(LUER,*)SPEC_ID_TWO(J)
	      STOP
	    END IF
	    TWO_PHOT_AVAILABLE(J)=.TRUE.
	  END IF
!
	END DO		!J: Which transition
!
	RETURN
	END
