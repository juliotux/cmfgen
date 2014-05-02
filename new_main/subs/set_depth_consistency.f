!
! Subroutine designed to adjust populations so that there is some sort of
! consistency from one depth to the next.
!
	SUBROUTINE SET_DEPTH_CONSISTENCY(STEQ_VALS,POPS,ND,NT,
	1              BAD_DECREASE_LIMIT,BAD_INCREASE_LIMIT,OPTION)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created 14-Dec-2013
!
	INTEGER ND		!Number of depth points
	INTEGER NT		!Number of equations
!
	REAL*8 STEQ_VALS(NT,ND)         !Suggested corrections to current populations
	REAL*8 POPS(NT,ND)
	REAL*8 BAD_INCREASE_LIMIT	!Should be large and neagtive.
	REAL*8 BAD_DECREASE_LIMIT	!Should be negative (typically just less than 1).
	CHARACTER(LEN=*) OPTION		!Not used currently (installed for later changes)
!
! Local variables.
!
	REAL*8 MAX_COR_VEC(ND)
	REAL*8 MIN_COR_VEC(ND)
	REAL*8 MAX_COR
	REAL*8 MIN_COR
	REAL*8 T1
!
	INTEGER ID
	INTEGER ISPEC
	INTEGER LOW_LEV
	INTEGER HIGH_LEV
	INTEGER L,K
	INTEGER CORRECTION_CNT
!
	LOGICAL BAD_SOLUTION_VEC(ND)
	LOGICAL PREV_DEPTH_OKAY(NUM_SPECIES)
!
	IF(BAD_INCREASE_LIMIT .GT. -1.0D+05)THEN
	  WRITE(6,*)'Error in SET_DEPTH_CONSISTENCY -- BAD_INCREASE_LIMIT is too small'
	  WRITE(6,*)'BAD_INCREASE_LIMIT should be less than -1.0D+05'
	  WRITE(6,*)'BAD_INCREASE_LIMIT is set in FIDDLE_POP_CORRECTIONS.F'
	  WRITE(6,*)'No adjustments made to population'
	  RETURN
	END IF
	IF(BAD_DECREASE_LIMIT .LT. 0.9999)THEN
	  WRITE(6,*)'Error in SET_DEPTH_CONSISTENCY -- BAD_DECREASE_LIMIT is too small'
	  WRITE(6,*)'BAD_DECREASE_LIMIT should be greater than 0.9999'
	  WRITE(6,*)'BAD_DECREASE_LIMIT is set in FIDDLE_POP_CORRECTIONS.F'
	  WRITE(6,*)'No adjustments made to population'
	  RETURN
	END IF
!
! Determine maximum/minum correction at each depth, and whether solution was
! obtained.
!
	DO K=1,ND
	  BAD_SOLUTION_VEC(K)=.FALSE.
	  MAX_COR_VEC(K)=MAXVAL(STEQ_VALS(1:NT,K)) 
	  MIN_COR_VEC(K)=MINVAL(STEQ_VALS(1:NT,K))
	  IF(MAX_COR_VEC(K) .EQ. 0.0D0 .AND. MIN_COR_VEC(K) .EQ. 0.0D0)BAD_SOLUTION_VEC(K)=.TRUE. 
	  WRITE(6,*)K,MIN_COR_VEC(K),MAX_COR_VEC(K),BAD_SOLUTION_VEC(K)
	END DO
!
! At present we loop outwards, and only make checks to a depth if the
! corrections to the current depth for the current species are relatively
! small. "Relative small" is currently hardwired".
!
! NB:  Correction limit should be large and positive (like 10^5) meaning that
!        the poulation is far from its equilibrium value. The current procedure
!        insures that changes will not be made if the populations have nearly
!        converged.
!
	PREV_DEPTH_OKAY(1:NUM_SPECIES)=.TRUE.
	DO K=ND-1,1,-1
	  CORRECTION_CNT=0
	  IF(MAX_COR_VEC(K) .LT. 0.5 .AND. MIN_COR_VEC(K) .GT. -2.0D0 .AND.
	1       .NOT. BAD_SOLUTION_VEC(K))THEN
	    PREV_DEPTH_OKAY(1:NUM_SPECIES)=.TRUE.
	    WRITE(6,*)'Depth ',K,' has good solutions'
	  ELSE
	    DO ISPEC=1,NUM_SPECIES
	      IF(SPECIES_PRES(ISPEC))THEN
	        ID=SPECIES_BEG_ID(ISPEC)
	        LOW_LEV=ATM(ID)%EQXzV
	        ID=SPECIES_END_ID(ISPEC)-1
	        HIGH_LEV=ATM(ID)%EQXzV+ATM(ID)%NXzV          !Points to last ion.
	        MAX_COR=MAXVAL(STEQ_VALS(LOW_LEV:HIGH_LEV,K)) 
	        MIN_COR=MINVAL(STEQ_VALS(LOW_LEV:HIGH_LEV,K)) 
	        IF(MAX_COR .LT. 0.5 .AND. MIN_COR .GT. -2.0D0 .AND.
	1              .NOT. BAD_SOLUTION_VEC(K))THEN
	          PREV_DEPTH_OKAY(ISPEC)=.TRUE.
	        ELSE
!
! At present we simply scale the populations allowing for density changes. We
! make no correction for rapid ionization changes. Thus procedure will work best
! with a well defined grid.
!
	          IF(PREV_DEPTH_OKAY(ISPEC))THEN
	             DO L=LOW_LEV,HIGH_LEV
	               IF(STEQ_VALS(L,K) .LT. BAD_INCREASE_LIMIT .OR. 
	1                   STEQ_VALS(L,K) .GT. BAD_DECREASE_LIMIT .OR. BAD_SOLUTION_VEC(K))THEN
	                  POPS(L,K)=POPS(L,K+1)*(POP_SPECIES(K,ISPEC)/POP_SPECIES(K+1,ISPEC))
	                  CORRECTION_CNT=CORRECTION_CNT+1
	               END IF
	             END DO
	          ELSE IF(K .GT. 1)THEN
	            MAX_COR=MAXVAL(STEQ_VALS(LOW_LEV:HIGH_LEV,K-1)) 
	            MIN_COR=MINVAL(STEQ_VALS(LOW_LEV:HIGH_LEV,K-1)) 
	            IF(MAX_COR .LT. 0.5 .AND. MIN_COR .GT. -2.0D0 .AND. .NOT.
	1               BAD_SOLUTION_VEC(K-1))THEN
	              DO L=LOW_LEV,HIGH_LEV
	                IF(STEQ_VALS(L,K) .LT. BAD_INCREASE_LIMIT .OR. 
	1                   STEQ_VALS(L,K) .GT. BAD_DECREASE_LIMIT .OR. BAD_SOLUTION_VEC(K))THEN
	                  POPS(L,K)=POPS(L,K-1)*(POP_SPECIES(K,ISPEC)/POP_SPECIES(K-1,ISPEC))
	                  CORRECTION_CNT=CORRECTION_CNT+1
	                END IF
	             END DO
	            END IF
	          END IF
	          PREV_DEPTH_OKAY(ISPEC)=.FALSE.
	        END IF
	      END IF
	    END DO                   !Loop over species
	  END IF                     !Depth check
	  IF(CORRECTION_CNT .NE. 0)WRITE(6,*)'Made corrections to depth ',K,CORRECTION_CNT
!
	END DO                       !Loop over depth
!
! Replace widely discrepant values with average value computed from surrounding
! data points. We do the double SQRT to avoid floating point undeflow.
!
	DO K=2,ND-1
	  DO L=1,NT-2
	    IF( (STEQ_VALS(L,K)   .LT. -1000.0D0 .OR. STEQ_VALS(L,K)   .GT. 0.99D0) .AND.
	1        STEQ_VALS(L,K-1) .GT. -2.0D0   .AND. STEQ_VALS(L,K-1) .LT. 0.5D0 .AND.
	1        STEQ_VALS(L,K+1) .GT. -2.0D0   .AND. STEQ_VALS(L,K+1) .LT. 0.5D0 )THEN
	      T1=SQRT(POPS(L,K-1))*SQRT(POPS(L,K+1))
	      IF( POPS(L,K) .GT. MAX(POPS(L,K-1),POPS(L,K+1)) .OR.
	1         POPS(L,K) .LT. MIN(POPS(L,K-1),POPS(L,K+1)))POPS(L,K)=T1
	    END IF
	  END DO
	END DO
!
	RETURN
	END
