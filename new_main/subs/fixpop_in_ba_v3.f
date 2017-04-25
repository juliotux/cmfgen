!
! Routine to alter the Statistical Equilibrium equations so that
! a particular population is held fixed. 
!
	SUBROUTINE FIXPOP_IN_BA_V3(BA,STEQ,ZERO_STEQ,
	1             NT,ND,NION,DIAG_BAND,DEPTH_INDX,
	1             FIRST_MATRIX,LAST_MATRIX)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 15-Sep-2011 : Fixed bug: IS was not being set [ID=SPECIES_END_ID(ISPEC)] for last ionization stage.
! Altered 30-May-2010 : Fixed bug with output describing number of levels held fixed.
! Created 23-Apr-2001 : Based on FIXPOP
!                       Designed to use MOD_SMGEN nd to operate on a single depth.
!                       To be call by GENERATE_FULL_MATRIX.
!
	INTEGER NT
	INTEGER ND
	INTEGER NION
	INTEGER DEPTH_INDX
	REAL*8 BA(NT,NT)
	REAL*8 STEQ(NT)
	LOGICAL ZERO_STEQ(NT)
	LOGICAL FIRST_MATRIX
	LOGICAL LAST_MATRIX
	LOGICAL DIAG_BAND
!
! Local variables.
!
	INTEGER I,J
	INTEGER ID,ISPEC
	INTEGER FIX_N
	INTEGER LOC_EQ
	INTEGER ERROR_LU,LUER
	LOGICAL LOC_IMP
	EXTERNAL ERROR_LU
!
! Variables to allow information to be output regarding the number
! of levels and depths where a population was held fixed.
!
	REAL*8 T1
	INTEGER, SAVE, ALLOCATABLE :: CNT(:)
!
	LUER=ERROR_LU()
	IF(ALLOCATED(CNT)) THEN
	  IF(SIZE(CNT) .NE. NT)THEN
	    WRITE(LUER,*)'Iconsistent dynamic allocation of CNT in FIXPOP_IN_BA_v2'
	    STOP
	  END IF
	ELSE
	  ALLOCATE(CNT(NT))
	END IF
	IF(FIRST_MATRIX .AND. DIAG_BAND)CNT(:)=0
	IF(DIAG_BAND)ZERO_STEQ(:)=.FALSE.
!
! Fix the electron density, if requested.
!
	IF(MOD_FIXED_NE)THEN
	  BA(NT-1,:)=0.0D0
	  IF(DIAG_BAND)THEN
	    BA(NT-1,NT-1)=1.0D0
	    STEQ(NT-1)=0.0D0
	    ZERO_STEQ(NT-1)=.TRUE.
	  END IF
	  IF(LAST_MATRIX)THEN
            WRITE(LUER,'(A)')' Electron density held fixed at all depths.'
	  END IF
	END IF
!
! Fix the Temperature, if requested.
!
	IF(MOD_FIXED_T .AND. DEPTH_INDX .GE. MOD_FIX_T_D_ST .AND.
	1                    DEPTH_INDX .LE. MOD_FIX_T_D_END)THEN
	  BA(NT,:)=0.0D0
	  IF(DIAG_BAND)THEN
	    BA(NT,NT)=1.0D0
	    STEQ(NT)=0.0D0
	    ZERO_STEQ(NT)=.TRUE.
	  END IF
	END IF
	IF(MOD_FIXED_T .AND. LAST_MATRIX)THEN
	   I=MOD_FIX_T_D_END-MOD_FIX_T_D_ST+1
	   IF(I .EQ. ND)THEN
              WRITE(LUER,'(A)')' Temperature held fixed at all depths.'
	   ELSE
              WRITE(LUER,'(A,ES12.3)')
	1              ' Temperature held partially fixed: TAU_SCL_T=',MOD_TAU_SCL_T
              WRITE(LUER,'(A,I3,A)')' Temperature held fixed at ',I,' depths.'
	   END IF
	 END IF
!
! We only write the information mesage out if STEQ has not been zeroed already.
!
	IF( FIX_IN_BOUND_T .AND. DEPTH_INDX .GE. (ND+1-FIX_LST_X_DPTHS) )THEN
	  BA(NT,:)=0.0D0
	  IF(DIAG_BAND)THEN
	    IF(DEPTH_INDX .EQ. ND .AND. .NOT. ZERO_STEQ(NT))THEN
              WRITE(LUER,'(A,I4,A)')' Temperature held fixed at',
	1                 FIX_LST_X_DPTHS,'depths at inner boundary'
	    END IF
	    BA(NT,NT)=1.0D0
	    STEQ(NT)=0.0D0
	    ZERO_STEQ(NT)=.TRUE.
	  END IF
	END IF
!
! This section handles the final ionization stage of each species.
!
	DO ISPEC=1,NUM_SPECIES
!
! Determine whether this depth is to be held fixed, and if so
! update depth counter.
!
	  IF( EQ_SPECIES(ISPEC) .NE. 0 .AND. 
	1              (FIX_SPECIES(ISPEC) .NE. 0 .OR. MOD_FIX_IMPURITY) )THEN
	    LOC_EQ=EQ_SPECIES(ISPEC)
	    ID=SPECIES_END_ID(ISPEC)
	    T1=ATM(ID-1)%DXzV(DEPTH_INDX)/POP_SPECIES(DEPTH_INDX,SPECIES_LNK(ID)) 
	    IF( (FIX_SPECIES(ISPEC) .NE. 0) .OR. T1 .LT. 1.0D-15)THEN
	      BA(LOC_EQ,:)=0.0D0
	      IF(DIAG_BAND)THEN
	        BA(LOC_EQ,LOC_EQ)=1.0D0
	        STEQ(LOC_EQ)=0.0D0
	        ZERO_STEQ(LOC_EQ)=.TRUE.
	        CNT(LOC_EQ)=CNT(LOC_EQ)+1
	      END IF
	    END IF
	  END IF
	END DO
!
! Handle those ionization stages with more than 1 level present.
!
	DO ID=1,NION
	  IF( ATM(ID)%XzV_PRES .AND. 
	1           (ATM(ID)%FIX_NXzV .NE. 0 .OR. MOD_FIX_IMPURITY) )THEN
!
	    IF(ATM(ID)%FIX_NXzV .NE. 0)THEN
	      FIX_N=MIN( ABS(ATM(ID)%FIX_NXzV),ATM(ID)%NXzV )
	      LOC_IMP=.FALSE.
	    ELSE IF(MOD_FIX_IMPURITY)THEN
	      FIX_N=ATM(ID)%NXzV
	      LOC_IMP=.TRUE.
	    ELSE
	      FIX_N=-10
	      LOC_IMP=.FALSE.
	    END IF
!
! Determine whether this depth is to be held fixed, and if so
! update depth counter.
!
	    IF(LOC_IMP)THEN
	      T1=0
	      DO J=1,ATM(ID)%NXzV
	        T1=T1+ATM(ID)%XzV(J,DEPTH_INDX)
	      END DO        
	      IF( T1/POP_SPECIES(DEPTH_INDX,SPECIES_LNK(ID)) .GT. 1.0D-15 )FIX_N=-10
	    END IF
	    LOC_EQ=ATM(ID)%EQXzV
	    IF(FIRST_MATRIX .AND. DIAG_BAND)CNT(LOC_EQ)=0
	    IF(DIAG_BAND .AND. FIX_N .GT. 0)CNT(LOC_EQ)=CNT(LOC_EQ)+1
!
! Zero requested equations.
!
	    IF(FIX_N .GT. 0)THEN
	      DO J=1,NT
	        DO I=ATM(ID)%EQXZV,ATM(ID)%EQXZV+FIX_N-1
	          BA(I,J)=0.0D0
	        END DO
	      END DO
!
	      IF(DIAG_BAND)THEN
	        DO I=ATM(ID)%EQXZV,ATM(ID)%EQXZV+FIX_N-1
	          BA(I,I)=1.0D0
	          STEQ(I)=0.0D0
	          ZERO_STEQ(I)=.TRUE.
 	        END DO
	      END IF
	    END IF
!	
	  END IF
!
	  IF(LAST_MATRIX .AND. CNT(ATM(ID)%EQXZV) .NE. 0)THEN
	    IF( ATM(ID)%XzV_PRES)THEN 
	      WRITE(LUER,100)FIX_N,ATM(ID)%NXzV,TRIM(ION_ID(ID)),CNT(ATM(ID)%EQXZV)
100	      FORMAT(1X,I3,' levels of',I4,' fixed for ',A,' at',I4,' depths')
	    ELSE
	      WRITE(LUER,110)TRIM(SPECIES(SPECIES_LNK(ID))),CNT(ATM(ID)%EQXZV)
110	      FORMAT(1X,'Last level of ',A,' held fixed at ',I4,' depths')
	    END IF
	  END IF
C
	END DO
C
	RETURN
	END
