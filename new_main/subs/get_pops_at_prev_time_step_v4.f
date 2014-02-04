!
! Subroutine to read in a model containing the populations from a
! previously converged model at an earlier time step. The new and current
! models must have an identical number of grid points, and super levels.
! Further, the velocity at the inner boundary must be identical. The older
! model may extend to larger velocities. 
!
! NB: TIME_SEQ_NO refers to the current model. This routine reads the populations
! corresponding to model TIME_SEQ_NO-1
!
        SUBROUTINE GET_POPS_AT_PREV_TIME_STEP_V4(OLD_POPS,OLD_R,DO_ADVECT,DO_RAD_DECAYS,
	1                      NORMALIZE_POPS,TIME_SEQ_NO,ND,NT,LU)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered : 21-Mar-2007 : READ_TIME_MODEL_V1 
! Altered : 01-Jul-2006 : Now allow outer velocity of older model to be larger than current
!                             model, although number of grid points still must be identical.
! Altered : 23-Jun-2006 : Removed R, V from call as in MOD_CMFGEN
! Created : 12-Dec-2005
!
	INTEGER ND
	INTEGER NT
	INTEGER LU
	INTEGER TIME_SEQ_NO
	LOGICAL DO_ADVECT
	LOGICAL DO_RAD_DECAYS
	LOGICAL NORMALIZE_POPS
!
	REAL*8 OLD_R(ND)
	REAL*8 OLD_POPS(NT,ND)
!
! Local arrays, vectors, and variables.
!
	REAL*8, ALLOCATABLE :: TMP_R(:)
	REAL*8, ALLOCATABLE :: TMP_V(:)
	REAL*8, ALLOCATABLE :: TMP_SIGMA(:)
	REAL*8, ALLOCATABLE :: TMP_POPS(:,:)
	REAL*8, ALLOCATABLE :: TMP_VEC(:)
	REAL*8, ALLOCATABLE :: TMP_DENSITY(:)
	REAL*8, ALLOCATABLE :: TMP_POP_ATOM(:)
	REAL*8, ALLOCATABLE :: LOG_TMP_V(:)
        LOGICAL OLD_ION_STAGE_PRES(NUM_IONS)
!
	REAL*8 LOG_V(ND)
	REAL*8 OLD_ED(ND)
	REAL*8 NEW_VEC(ND)
!
	REAL*8 T1,T2
	REAL*8 OLD_SN_AGE
!
	INTEGER ND_OLD
	INTEGER IOS
	INTEGER IREC_RD
	INTEGER I,J,K,ISPEC,ID
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU,EQUAL
	LOGICAL EQUAL
!
	INTEGER, PARAMETER :: IONE=1
	LOGICAL, PARAMETER :: RVSIG_WRITTEN=.TRUE.
!
!	ND_OLD=ND
!
	CALL GET_ND_SEQ_MODEL_FILE(ND_OLD,LU)
	ALLOCATE (TMP_R(ND_OLD),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TMP_V(ND_OLD),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TMP_SIGMA(ND_OLD),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TMP_VEC(ND_OLD),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TMP_DENSITY(ND_OLD),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TMP_POP_ATOM(ND_OLD),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (LOG_TMP_V(ND_OLD),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (TMP_POPS(NT,ND_OLD),STAT=IOS)
!
! Get model from the last time step. TIME_SEQ_NO refers to the CURRENT 
! time model. Therefore we must subtract 1.
!
!	IREC_RD=TIME_SEQ_NO-1
!	CALL  READ_TIME_MODEL_V1(TMP_R,TMP_V,TMP_SIGMA,TMP_POPS,
!	1             IREC_RD,RVSIG_WRITTEN,NT,ND,LU)
        CALL READ_SEQ_TIME_FILE_V1(TMP_R,TMP_V,TMP_SIGMA,TMP_POP_ATOM,
	1            TMP_DENSITY,TMP_POPS,
	1            OLD_ION_STAGE_PRES,OLD_SN_AGE,ND_OLD,NT,LU)
!
! As a Hubble law, we can use V to interpolate. Note that 
! V is a comoving variable.
!
	T1=1.0D-06
	IF(EQUAL(TMP_V(ND_OLD),V(ND),T1))THEN
	  TMP_V(ND_OLD)=V(ND)
	ELSE 
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in GET_POPS_AT_PREV_TIME_STEP_V4' 
	  WRITE(LUER,*)'Velocities at inner boundary are unequal'
	  WRITE(LUER,*)'V(ND)=',V(ND)
	  WRITE(LUER,*)'OLD_V(ND_OLD)=',TMP_V(ND_OLD)
	  WRITE(LUER,*)'V(ND)/OLD_V(ND)=',V(ND)/TMP_V(ND)
	  WRITE(LUER,*)'Accuracy indicator=',T1
	  STOP
	END IF
	IF(EQUAL(TMP_V(1),V(1),T1))THEN
	  TMP_V(1)=V(1)
	ELSE IF(TMP_V(1) .LT. V(1))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in GET_POPS_AT_PREV_TIME_STEP_V4' 
	  WRITE(LUER,*)'Old velocity at outer boundary is too small'
	  WRITE(LUER,*)'V(1)=',V(1)
	  WRITE(LUER,*)'OLD_V(1)=',TMP_V(1)
	  WRITE(LUER,*)'V(1)/OLD_V(1)=',V(1)/TMP_V(1)
	  WRITE(LUER,*)'Accuracy indicator=',T1
	END IF
!
! Perform interpolations in Log-Log plane, with V as the independent variable.
! Since V is proportional to r, this is equivalent to assuming r is the
! independent variable.
!
	LOG_V=LOG(V)
	LOG_TMP_V=LOG(TMP_V)
	DO I=1,NT
	  TMP_VEC=LOG(TMP_POPS(I,:))
	  CALL MON_INTERP(NEW_VEC,ND,IONE,LOG_V,ND,TMP_VEC,ND_OLD,LOG_TMP_V,ND_OLD)
	  OLD_POPS(I,:)=EXP(NEW_VEC)
	END DO
!
! Get the interpolated radius scale in the old model. This radius scale will 
! have the same velocity coordinates as the current model. Since we have a 
! Hubble law, linear interpolation is accurate.
!
	CALL MON_INTERP(OLD_R,ND,IONE,V,ND,TMP_R,ND_OLD,TMP_V,ND_OLD)
!
! From now on the OLD models has the same length as the new model: ND
!
! The following works for an arbitrary expansion factor.
!
	IF(DO_ADVECT)THEN
	  DO I=1,ND
	    OLD_POPS(1:NT-1,I)=OLD_POPS(1:NT-1,I)/VOL_EXP_FAC(I)
	  END DO
	END IF
!
! Adjust the populations for radioactive decays.
!
	IF(DO_RAD_DECAYS)THEN
	  CALL DO_LEV_POP_DECAYS(OLD_POPS,ND,NT)
	END IF
!
! Normalize the populations to ensure continuity equation is satisfied.
! This removes any slight variations introduced by the non-monotonic
! interpolation, rounding errors, and the numerical accuracy of the
! data files used for input. The normalization should not be done if
! radioactive decays are important and DO_RAD_DECAYS is false.
!
	IF(NORMALIZE_POPS)THEN
	  IF(.NOT. DO_RAD_DECAYS)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error in GET_POPS_AT_PREV_TIME_STEP_V4' 
	    WRITE(LUER,*)'You should not normalize pops when DO_RAD_DECAYS is FALSE.'
	    STOP
	  END IF
	  DO ISPEC=1,NUM_SPECIES
	    DO J=1,ND
	      T2=0.0D0
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        DO I=1,ATM(ID)%NXzV
	          K=ATM(ID)%EQXzV+I-1
	          T2=T2+OLD_POPS(K,J)
	        END DO
	      END DO
	      IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	        ID=SPECIES_END_ID(ISPEC)-1
	        K=ATM(ID)%EQXzV+ATM(ID)%NXzV
	        T2=T2+OLD_POPS(K,J)
!
! Can now do the normalization. Valid for all expansion laws.
!
	        IF(DO_ADVECT)THEN
	          T2=POP_SPECIES(J,ISPEC)/T2
	        ELSE
	          T2=POP_SPECIES(J,ISPEC)*VOL_EXP_FAC(I)/T2
	        END IF
	        DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	          DO I=1,ATM(ID)%NXzV
	            K=ATM(ID)%EQXzV+I-1
	            OLD_POPS(K,J)=OLD_POPS(K,J)*T2
	          END DO
	        END DO
	        ID=SPECIES_END_ID(ISPEC)-1
	        K=ATM(ID)%EQXzV+ATM(ID)%NXzV
	        OLD_POPS(K,J)=OLD_POPS(K,J)*T2
	      END IF
	    END DO
	  END DO
	END IF
!
! Now determine the electron density. Since we determine Ne from the scaled
! populations, no scaling is necessary.
!
	OLD_ED=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  DO J=1,ND
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV
	        K=ATM(ID)%EQXzV+I-1
	        OLD_ED(J)=OLD_ED(J)+(ATM(ID)%ZXzV-1.0D0)*OLD_POPS(K,J)
	      END DO
	    END DO
	    IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	      ID=SPECIES_END_ID(ISPEC)-1
	      K=ATM(ID)%EQXzV+ATM(ID)%NXzV
	      OLD_ED(J)=OLD_ED(J)+ATM(ID)%ZXzV*OLD_POPS(K,J)
	    END IF
	  END DO
	END DO
	DO J=1,ND
	  OLD_POPS(NT-1,J)=OLD_ED(J)
	END DO
!
	RETURN
	END
