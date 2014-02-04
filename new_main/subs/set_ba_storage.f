!
! Routine to create most of the storage for STEQ_DATA_MOD.
! The LNK_TO_IV and LNK_TO_F vectors were allocated in
! CREATE_IV_LINKS_V2.
!
      SUBROUTINE SET_BA_STORAGE(NT,NUM_BNDS,ND,NION)
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      IMPLICIT NONE
!
! Altered 04-Apr-2013 : MEMORY made REAL to avoid integer overflow.
! Altered 01-Apr-2004 : BA matrices allocated after all STEQ matrices.
! Created 05-Apr-2001
!
      INTEGER NT
      INTEGER NUM_BNDS
      INTEGER ND
      INTEGER NION
!
      INTEGER ID
      INTEGER ISPEC
      INTEGER IOS
      INTEGER NSUM
      INTEGER LU_ER,ERROR_LU
      EXTERNAL ERROR_LU
!
      REAL*8 MEMORY
      INTEGER I,NX,NY
!
      LU_ER=ERROR_LU()
!
! Check consistency of parameters between two modules.
!
      IF(  BA_NUM_SPECIES          .NE. NUM_SPECIES              .OR.
     &     BA_MAX_IONS_PER_SPECIES .NE. MAX_IONS_PER_SPECIES     .OR.
     &     BA_MAX_NUM_IONS         .NE. MAX_NUM_IONS             .OR.
     &     BA_NPHOT_MAX            .NE. NPHOT_MAX)               THEN
        WRITE(LU_ER,*)'Inconsistency in parameters SET_BA_STORAGE'
        WRITE(LU_ER,*)'Check MOD_CMFGEN and STEQ_DATA_MOD for consistency of'
        WRITE(LU_ER,*)'NUM_SPECIES, NPHOT_MAX etc'
      STOP
      END IF
!
! Will need to add an extra equation for X-rays.
!
      MEMORY=0.0D0
      DO ID=1,NION
	IF(SE(ID)%XzV_PRES)THEN
	  NX=SE(ID)%N_SE
	  NY=SE(ID)%N_IV
	ELSE
	  NX=1; NY=1
	END IF
	SE(ID)%IMPURITY_SPECIES=.FALSE.
	IF(ATM(ID)%NXzV_IV .EQ. 0)SE(ID)%IMPURITY_SPECIES=.TRUE.
!
                      ALLOCATE (SE(ID)%STEQ(NX,ND),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (SE(ID)%QFV_R(NX,ND),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (SE(ID)%QFV_P(NX,ND),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (SE(ID)%EQ_IN_BA(NX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SE(ID)%STEQ_ADV(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SE(ID)%STRT_ADV_ID(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( SE(ID)%END_ADV_ID(ND),STAT=IOS)
        MEMORY=MEMORY+NX*NY*ND*(NUM_BNDS+1)
!
        IF(IOS .NE. 0)THEN
          WRITE(LU_ER,*)'Unable to allocate SE(ID)%STEQ in SET_BA_STORAGE'
          WRITE(LU_ER,*)'STAT=',IOS,'ID=',ID
          STOP
        END IF
!
	DO I=1,ATM(ID)%NxZV
	  SE(ID)%EQ_IN_BA(I)=ATM(ID)%EQXzV+I-1
	END DO
	DO I=ATM(ID)%NXzV+1,NX-1
	  SE(ID)%EQ_IN_BA(I)=ATM(ID)%EQXzV+ATM(ID)%NXzV+(SE(ID)%EQ_TO_ION_LEV_PNT(I)-1)
	END DO
        SE(ID)%EQ_IN_BA(NX)=EQ_SPECIES(SPECIES_LNK(ID))
	SE(ID)%NUMBER_BAL_EQ=SE(ID)%N_SE
      END DO
!
	            ALLOCATE (STEQ_ED(ND),STAT=IOS)
      IF(IOS .EQ. 0)ALLOCATE (STEQ_T(ND),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(LU_ER,*)'Unable to allocate STEQ_ED etc in SET_BA_STORAGE'
        WRITE(LU_ER,*)'STAT=',IOS,'ID=',ID
        STOP
      END IF
!
! We try allocate the BA matrices separately, as this might 
! improve memory management.
!
      DO ID=1,NION
	IF(SE(ID)%XzV_PRES)THEN
	  NX=SE(ID)%N_SE
	  NY=SE(ID)%N_IV
	ELSE
	  NX=1; NY=1
	END IF
	IOS=0
        IF(IOS .EQ. 0)ALLOCATE (SE(ID)%BA_PAR(NX,NY,ND),STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE (SE(ID)%BA(NX,NY,NUM_BNDS,ND),STAT=IOS)
!
        IF(IOS .NE. 0)THEN
          WRITE(LU_ER,*)'Unable to allocate memory in SET_BA_STORAGE'
          WRITE(LU_ER,*)'STAT=',IOS,'ID=',ID
          STOP
        END IF
      END DO
!
! Use to store Charge Equilibrium, and Radiative Equilibrium equations.
!
      IF(IOS .EQ. 0)ALLOCATE (BA_ED(NT,NUM_BNDS,ND),STAT=IOS)
      IF(IOS .EQ. 0)ALLOCATE (BA_T(NT,NUM_BNDS,ND),STAT=IOS)
      IF(IOS .EQ. 0)ALLOCATE (BA_T_PAR(NT,ND),STAT=IOS)
      IF(IOS .EQ. 0)ALLOCATE (BA_ADV_TERM(NUM_BNDS,ND),STAT=IOS)
      IF(IOS .NE. 0)THEN
        WRITE(LU_ER,*)'Unable to allocate BA_ED etc in SET_BA_STORAGE'
        WRITE(LU_ER,*)'STAT=',IOS,'ID=',ID
        STOP
      END IF
!
      MEMORY=MEMORY+2*NT*ND*NUM_BNDS+NT*ND
      WRITE(LU_ER,*)' '
      WRITE(LU_ER,'(A,ES17.10,A)')' Amount of memory allocated for BA is:  ',MEMORY,' words'
      MEMORY=DFLOAT(NT)*NT*(NUM_BNDS+1)*ND
      WRITE(LU_ER,'(A,ES17.10,A)')' Memory needed with full dependence is: ',MEMORY,' words'
!
      RETURN
      END
