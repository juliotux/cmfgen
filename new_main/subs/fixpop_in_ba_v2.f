!
! Routine to alter the Statistical Equilibrium equations so that
! a particular population is held fixed. 
!
	SUBROUTINE FIXPOP_IN_BA_V2(BA,STEQ,NT,ND,NION,DIAG_BAND,DEPTH_INDX,
	1             FIRST_MATRIX,LAST_MATRIX,FIX_IMPURITY)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
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
	LOGICAL FIRST_MATRIX
	LOGICAL LAST_MATRIX
	LOGICAL DIAG_BAND
	LOGICAL FIX_IMPURITY
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
!
! FIX_NSPEC takes priority in determining the number of levels
! to be fixed. 
!
	DO ID=1,NUM_IONS
	  ISPEC=SPECIES_LNK(ID)
	  IF(.NOT. ATM(ID)%XZV_PRES)THEN
!
! This section handles the final ionization stage of each species.
!
	    IF(FIX_SPECIES(ISPEC) .EQ. 0 .AND. .NOT. FIX_IMPURITY)THEN
	    ELSE
	      LOC_EQ=ATM(ID-1)%EQXzV+ATM(ID-1)%NXzV
	      LOC_IMP=.TRUE.
	      IF(FIX_SPECIES(ISPEC) .NE. 0)LOC_IMP=.FALSE.
!
! Determine whether this depth is to be held fixed, and if so
! update depth counter.
!
	      IF(LOC_IMP)THEN
	        T1=ATM(ID-1)%DXzV(DEPTH_INDX)/POP_SPECIES(DEPTH_INDX,SPECIES_LNK(ID)) 
	        IF(T1 .LT. 1.0D-15)THEN
	          BA(LOC_EQ,:)=0.0D0
	          BA(LOC_EQ,LOC_EQ)=1.0D0
	          STEQ(LOC_EQ)=0.0D0
	          IF(DIAG_BAND)CNT(LOC_EQ)=CNT(LOC_EQ)+1
	        END IF
	      ELSE
	        BA(LOC_EQ,:)=0.0D0
	        BA(LOC_EQ,LOC_EQ)=1.0D0
	        STEQ(LOC_EQ)=0.0D0
	        IF(DIAG_BAND)CNT(LOC_EQ)=CNT(LOC_EQ)+1
	      END IF
	    END IF
	  ELSE IF(ATM(ID)%FIX_NXzV .EQ. 0 .AND. .NOT. FIX_IMPURITY)THEN
	  ELSE
!
! We now handle the ions with more than 1 level present/
!
	    IF(ATM(ID)%FIX_NXzV .NE. 0)THEN
	      FIX_N=MIN( ABS(ATM(ID)%FIX_NXzV),ATM(ID)%NXzV )
	      LOC_IMP=.FALSE.
	    ELSE
	      FIX_N=ATM(ID)%NXzV
	      LOC_IMP=.TRUE.
	    END IF
	    IF(FIRST_MATRIX .AND. DIAG_BAND)CNT(LOC_EQ)=0
!
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
	    IF(DIAG_BAND .AND. FIX_N .GT. 0)CNT(LOC_EQ)=CNT(LOC_EQ)+1
!
	    DO J=1,NT
	      DO I=ATM(ID)%EQXZV,ATM(ID)%EQXZV+FIX_N-1
	        BA(I,J)=0.0D0
	      END DO
	    END DO
!
	    DO I=ATM(ID)%EQXZV,ATM(ID)%EQXZV+FIX_N-1
	      BA(I,I)=1.0D0
	      STEQ(I)=0.0D0
	    END DO
	  END IF
!
	  IF(LAST_MATRIX .AND. CNT(ATM(ID)%EQXZV) .NE. 0)THEN
	    WRITE(LUER,100)FIX_N,ATM(ID)%NXzV,TRIM(ION_ID(ID)),CNT(ATM(ID)%EQXZV)
100	    FORMAT(1X,I3,' levels of',I4,' fixed for ',A,' at',I4,' depths')
	  END IF
C
	END DO
C
	RETURN
	END
