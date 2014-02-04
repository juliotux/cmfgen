!
! Routine designed to convert the SE%STEQ and SE%BA structures into the
! full STEQ and BA matrices that are used to find the corrections.
! The conversion is done one depth at a time. At each depth, the DIAGONAL
! matrix should be passed first.
!
	SUBROUTINE GENERATE_FULL_MATRIX(C_MAT,STEQ_VEC,POPS,NT,ND,NION,
	1                    NUM_BNDS,BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1                    FIRST_MATRIX,LAST_MATRIX,FIX_IMPURITY)
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 10-Sep-2001 : Using DIAG_INDX as logical variable in an IF statement
! Created 05-Apr-2001
!
!******************************************************************************
!
	INTEGER NT
	INTEGER ND
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER BAND_INDX 
	INTEGER DIAG_INDX
	INTEGER DEPTH_INDX
!
        REAL*8 POPS(NT,ND)
        REAL*8 C_MAT(NT,NT)
	REAL*8 STEQ_VEC(NT)
	LOGICAL FIRST_MATRIX
	LOGICAL LAST_MATRIX
	LOGICAL FIX_IMPURITY
!
! Local matrices and vectors to facilitate construction of C_MAT and STEQ_VEC.
!
! We use ?_ION to create the ionization balance equations for each
! ionization stage.
!
	REAL*8 C_ION(NION,NT)
	REAL*8 STEQ_ION(NION)
!
! We use ?_NC for the Number conservation equation for each species.
! We don't replace them directly into C_MAT and STEQ_VEC for ease of
! programming.
!
	REAL*8 C_NC(NUM_SPECIES,NT)
	REAL*8 STEQ_NC(NUM_SPECIES)
!
! FAC is used as a scale factor to determine at what depth the ground-state
! equilibrium equation is replaced by the ionization equation.
!
!	REAL*8, PARAMETER :: FAC=1.0D+05
	REAL*8, PARAMETER :: FAC=1.0D+02
!
	LOGICAL DIAG_BAND
!
! For consistency with the old version of CMFGEN, we only replace the
! ground-state equations over a continuous set of depths. REPLACE
! is used to indicate whether the current depth is to be rplaced
! (as determined from the DIAGONAL band).
!
	INTEGER, SAVE, ALLOCATABLE ::  REP_CNT(:)
	LOGICAL,   SAVE, ALLOCATABLE ::  REPLACE(:)
!
	INTEGER, SAVE :: LST_DEPTH_INDX=0
!
	INTEGER I,J,K,L,JJ,N
	INTEGER EQ
	INTEGER ID
	INTEGER ISPEC
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! To save typing.
!
	LUER=ERROR_LU()
        K=DEPTH_INDX
	DIAG_BAND=.FALSE.
	IF(DIAG_INDX .EQ. BAND_INDX)DIAG_BAND=.TRUE.
	IF(DEPTH_INDX .NE. LST_DEPTH_INDX)THEN
	  LST_DEPTH_INDX=DEPTH_INDX
	  IF(.NOT. DIAG_BAND)THEN
	    WRITE(LUER,*)'Error in GENERATE_FULL_MATRIX'
	    WRITE(LUER,*)'Diagonal band must be passed first at each depth'
	    WRITE(LUER,*)'DEPTH_INDX=',DEPTH_INDX
	    STOP
	  END IF
	END IF
!
	C_MAT(:,:)=0.0D0
	C_ION(:,:)=0.0D0
	C_NC(:,:)=0.0D0
!
	IF(DIAG_BAND)THEN
	  STEQ_VEC(:)=0.0D0
	  STEQ_ION(:)=0.0D0
	  STEQ_NC(:)=0.0D0
!
	  IF(FIRST_MATRIX)THEN
            I=SE(4)%LNK_TO_IV(1101)
	    WRITE(99,*)I
            WRITE(99,*)SE(1)%BA(1,I,DIAG_INDX,1),SE(1)%BA(2,I,DIAG_INDX,1)
            WRITE(99,*)SE(3)%BA(1,I,DIAG_INDX,1),SE(3)%BA(2,I,DIAG_INDX,1)
            WRITE(99,*)SE(4)%BA(1,I,DIAG_INDX,1),SE(4)%BA(2,I,DIAG_INDX,1)
	  END IF
!
	END IF
!
	IF( .NOT. ALLOCATED(REPLACE))THEN
	  ALLOCATE (REPLACE(NION)); REPLACE(:)=.FALSE.
	  ALLOCATE (REP_CNT(NION)); REP_CNT(:)=0
	END IF
!
! Map the small BA array onto the full BA array
!
! NB: C_ION(1,:) refers to to the ionization/recombination equation
!                for ion 1 (e.g. CI in the carbon sequence). It is
!                dN(CI)/dt. Since the eqaution (i.e. BA(I,:,:,:) with 
!                I > ATMD(ID)%NXzV refers to dN/dt for the recombining 
!                level (i.e. C2) we need a - sign when we evaluate 
!                C_ION and STEQ_ION.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    DO J=1,SE(ID)%N_IV
	      JJ=SE(ID)%LNK_TO_F(J)
	      DO I=1,SE(ID)%N_SE-1
	        EQ=SE(ID)%EQ_IN_BA(I)
                C_MAT(EQ,JJ)=C_MAT(EQ,JJ)+SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	      DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	        C_ION(ID,JJ)=C_ION(ID,JJ)-SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	      C_NC(ISPEC,JJ)=C_NC(ISPEC,JJ)+SE(ID)%BA(SE(ID)%N_SE,J,BAND_INDX,K)
	    END DO
	  END DO
	END DO
!
! Now do STEQ
!
	IF(DIAG_BAND)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,SE(ID)%N_SE-1
	        EQ=SE(ID)%EQ_IN_BA(I)
                STEQ_VEC(EQ)=STEQ_VEC(EQ)+SE(ID)%STEQ(I,K)
	      END DO
	      DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	        STEQ_ION(ID)=STEQ_ION(ID)-SE(ID)%STEQ(I,K)
	      END DO
	      STEQ_NC(ISPEC)=STEQ_NC(ISPEC)+SE(ID)%STEQ(SE(ID)%N_SE,K)
	    END DO
	  END DO
	  STEQ_VEC(NT-1)=STEQ_ED(K)
	  STEQ_VEC(NT)=STEQ_T(K)
	END IF
!
! Update the ionization equations for when X-rays are included.
! The following is photoionizations/recombinations which change z by 2.
! Only other process allowed are /\z=1. The corrections to the ionization 
! equations follow from simple algebraic manipulations.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-2
            I=SE(ID)%XRAY_EQ
	    IF(I .NE. 0)THEN
	      DO J=1,SE(ID)%N_IV
	        JJ=SE(ID)%LNK_TO_F(J)
	        C_ION(ID+1,JJ)=C_ION(ID+1,JJ)-SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	      IF(DIAG_BAND)STEQ_ION(ID+1)=STEQ_ION(ID+1)-SE(ID)%STEQ(I,K)
	    END IF
	  END DO
	END DO
!
! We now replace any equations that need replacing.
!
! Insert the number conservation equation for each species.
!
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    I=EQ_SPECIES(ISPEC)
	    C_MAT(I,:)=C_NC(ISPEC,:)
	  END IF
	END DO
!
	DO ISPEC=1,NUM_SPECIES
	  IF(DIAG_BAND .AND. SPECIES_PRES(ISPEC))THEN
	    I=EQ_SPECIES(ISPEC)
	    STEQ_VEC(I)=STEQ_NC(ISPEC)
	  END IF
	END DO
!
! Charge conservation and radiative equilibrium equations.
!
	C_MAT(NT-1,:)=BA_ED(:,BAND_INDX,DEPTH_INDX)
	C_MAT(NT,:)  =BA_T(:,BAND_INDX,DEPTH_INDX)
!
! Check if we wish to use the ionization conservation equation.
! We only use the diagonal band to do this. We use REPLACE
! array so that we only use the conservation equation over a
! sequence of optical depths.
!
	IF(DIAG_BAND .AND. DEPTH_INDX .EQ. 1)REP_CNT(:)=0
        DO ID=1,NION
	  IF(ATM(ID)%XzV_PRES .AND. DIAG_BAND)THEN
	    EQ=ATM(ID)%EQXzV
	    N=ATM(ID)%NXzV
            IF( ABS(C_MAT(EQ,EQ))*ATM(ID)%XzV(1,K) .GT.
	1          FAC*ABS(C_MAT(EQ,EQ+N))*ATM(ID)%DXzV(K) )THEN
	      REPLACE(ID)=.TRUE.
	      REP_CNT(ID)=REP_CNT(ID)+1
	    ELSE
	      REPLACE(ID)=.FALSE.
	    END IF
	  END IF
	END DO
!
! In all cases, we replace the ground state equation.
! We need the check on XzV_PRES since we initially set
! all the REPLACE to true.
!
	DO ID=1,NION
          IF(REPLACE(ID) .AND. ATM(ID)%XzV_PRES)THEN
	    C_MAT(ATM(ID)%EQXzV,:)=C_ION(ID,:)
	    IF(DIAG_BAND)STEQ_VEC(ATM(ID)%EQXzV)=STEQ_ION(ID)
	  END IF
	END DO
!
	IF(DIAG_BAND .AND. DEPTH_INDX .EQ. ND)THEN
	  DO ID=1,NION
	    IF(REP_CNT(ID) .GT. 0)THEN
              WRITE(LUER,'(X,A,T9,A,I3,A)')
	1         TRIM(ION_ID(ID)),' g.s. eq. replaced by ionization eq. at ',
	1         REP_CNT(ID),' depths.'
	      END IF
	  END DO
	END IF
!
	IF(K .EQ. 1 .AND. DIAG_BAND)THEN
	  CALL WR2D(STEQ_VEC,NT,1,'STEQ_VEC_D1',94)
	  OPEN(UNIT=96,FILE='BA_ASCI_D1',STATUS='UNKNOWN')
	    CALL WR2D(C_MAT,NT,NT,'C_MAT_VEC_D1',96)
	  CLOSE(UNIT=96)
	END IF
	IF(K .EQ. 10 .AND. DIAG_BAND)THEN
	  CALL WR2D(STEQ_VEC,NT,1,'STEQ_VEC_D10',94)
	  OPEN(UNIT=96,FILE='BA_ASCI_D1',STATUS='UNKNOWN',POSITION='APPEND')
	    CALL WR2D(C_MAT,NT,NT,'C_MAT_VEC_D10',96)
	  CLOSE(UNIT=96)
	END IF
!
! Fix any populations.
!        
	CALL FIXPOP_IN_BA_V2(C_MAT,STEQ_VEC,NT,ND,NION,DIAG_BAND,K,
	1       FIRST_MATRIX,LAST_MATRIX,FIX_IMPURITY)
!
! Scale the BA matrix so that we solve for the fractional corrections to
! the populations. This seems to yield better solutions for large matrices,
! especially when we subsequently condition the matrices.
!                
	L=DEPTH_INDX-DIAG_INDX+BAND_INDX
	DO J=1,NT
	  DO I=1,NT
	    C_MAT(I,J)=C_MAT(I,J)*POPS(J,L)
	  END DO
	END DO
!
	RETURN
	END
