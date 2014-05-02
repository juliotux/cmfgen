!
! This routine has been designed to use with SOLVEBA_V9. At present only the
! MAJOR otpion is treated. This section was made a subroutine to faciliate easy modification
! in the hope that we could improve convergence. With the use of the file 'ADJUST_CORRECTIONS'
! we can adjust the relaxation parameter at specific depths, and on the fly.
!
	SUBROUTINE FIDDLE_POP_CORRECTIONS_V2(POPS,STEQ,T_MIN,CHANGE_LIM,MAX_dT_COR,
	1              SCALE_OPT,LAMBDA_IT,LU_SUM,NT,ND)
	IMPLICIT NONE
!
! Altered: 14-Feb-2014 -- Changed to V2 -- added MAX_dT_COR to call.
! Altered: 01-Jan-2014 -- Modified to use standard read routines. Some cleaning and
!                           meaning of POP_LIM altered.  
! Altered: 14-DEc-2013 -- Subroutine calls SET_DEPTH_CONSISTENCY if LAMBDA iteration.
! Altered: 30-Sep-2011 -- T LIMIT and POP LIMIT introduced -- Based on a similar modification
!                            by Chendong. Some cleaning, and better error messages.
! Altered: 01-Feb-2011 -- Can use a : in list to specify a consecutive range.
! Created: 15-Nov-2010
!
	INTEGER NT
	INTEGER ND
	INTEGER LU_SUM
	REAL*8 POPS(NT,ND)
	REAL*8 STEQ(NT,ND)
!
	REAL*8 T_MIN
	REAL*8 CHANGE_LIM
	REAL*8 MAX_dT_COR
	LOGICAL LAMBDA_IT
	CHARACTER(LEN=*) SCALE_OPT
!
! Local variables
!
	REAL*8 RELAX_PARAM(ND)
	REAL*8 POP_LIM(ND)
	REAL*8 T_LIM(ND)
!
	REAL*8 RELAX_VARIABLE
	REAL*8 T_LIM_VARIABLE
	REAL*8 POP_LIM_VARIABLE
!
	REAL*8 LIT_LIM,BIG_LIM
	REAL*8 DPTH_LIT_LIM,DPTH_BIG_LIM
	REAL*8 T1,T2,T3
	REAL*8 MIN_SCALE
	REAL*8 SCALE
	REAL*8 BAD_DECREASE_LIMIT
	REAL*8 BAD_INCREASE_LIMIT
!
	INTEGER, PARAMETER :: LU_SCR=26
	INTEGER, SAVE:: COUNTER=0
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	INTEGER CONSISTENCY_CNT
	INTEGER IOS
	INTEGER I,J,IC
	INTEGER L,L_ST,L_END
	LOGICAL FILE_OPEN
	CHARACTER(LEN=20) DC_OPTION
	CHARACTER(LEN=80) STRING
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
	BIG_LIM=(CHANGE_LIM-1.0D0)/CHANGE_LIM
        LIT_LIM=1.0D0-CHANGE_LIM
	MIN_SCALE=1.0D+20
!
! Set default parameters.
!
	L_ST=1; L_END=ND
	RELAX_VARIABLE=1.0D0
	T_LIM_VARIABLE=MAX_dT_COR
	POP_LIM_VARIABLE=100.0D0*CHANGE_LIM		!=>implies no effect
	DO I=1,ND
	  RELAX_PARAM(I)=RELAX_VARIABLE
	  T_LIM(I)=T_LIM_VARIABLE
	  POP_LIM(I)=POP_LIM_VARIABLE
	END DO
	BAD_DECREASE_LIMIT=1.0D0-1.0D-10
	BAD_INCREASE_LIMIT=-1.0D+10
	CONSISTENCY_CNT=0
!
! Valid ranges:
!              0 < RELAX_PARAM < 2
!              0 < T_LIM < 0.2 
!              POP_LIM > 1
!
! The file, ADJUST_CORRECTIONS, need not be present.
!
! NB: The relaxation parameter and the T LIMIT effect all corrections (if used) at the
!     specified depth. POP LIMIT only effects corrections bigger than POP LIMIT.
!
	
	OPEN(UNIT=LU_SUM,FILE='ADJUST_CORRECTIONS',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .EQ. 0)THEN
	  CALL RD_OPTIONS_INTO_STORE(LU_SUM,LU_SCR)
	  CALL RD_STORE_INT(L_ST,  'L_ST',L_FALSE,'Beginning depth')
	  CALL RD_STORE_INT(L_END,'L_END',L_FALSE,'Final depth')
	  CALL RD_STORE_DBLE(T_LIM_VARIABLE,'T_LIM',L_FALSE,'Maximum fracton correcton to T')
	  CALL RD_STORE_DBLE(RELAX_VARIABLE,'RELAX',L_FALSE,'Relaxation variable')
	  CALL RD_STORE_DBLE(POP_LIM_VARIABLE,'MAX_CHNG',L_FALSE,'1 + Maximum fractional increase (> 1)')
	  CALL RD_STORE_INT(CONSISTENCY_CNT,'CONSIS_CNT',L_FALSE,'Check whether adjacent pops consistent every ? iterations')
	  CALL CLEAN_RD_STORE()
!
	  IF(L_ST .LT. 1 .OR. L_END .GT. ND)THEN
	    WRITE(6,*)'Error in depth indices in ADJUST_CORRECTIONS -- invalid range'
	    WRITE(6,*)'Depth indices=',L_ST,L_END
	    L_ST=1; L_END=ND
	    RELAX_VARIABLE=1.0D0; T_LIM_VARIABLE=0.2D0
	    POP_LIM_VARIABLE=100.0D0*CHANGE_LIM		!=>implies no effect
	  END IF
!
	  IF(T_LIM_VARIABLE .LT. 0.0D0 .OR. T_LIM_VARIABLE .GT. 0.20D0)THEN
	    WRITE(6,*)'Error for T LIMIT in ADJUST_CORRECTIONS -- invalid value'
	    WRITE(6,*)'Valid range is 0.0 to 0.2'
	    WRITE(6,*)'T LIMIT parameter read in is',T_LIM_VARIABLE
	    T_LIM_VARIABLE=0.2D0
	  END IF
	  DO L=L_ST,L_END
	    T_LIM(L)=T_LIM_VARIABLE
	  END DO
!
	  IF(POP_LIM_VARIABLE .LE. 1.0D0)THEN
	    WRITE(6,*)'Error for POP LIMIT in ADJUST_CORRECTIONS -- invalid value'
	    WRITE(6,*)'Valid range is > 1.0'
	    WRITE(6,*)'POP LIMIT parameter read in',POP_LIM_VARIABLE
	    POP_LIM_VARIABLE=100.0D0*CHANGE_LIM		!=>implies no effect
	  END IF
	  DO L=L_ST,L_END
	    POP_LIM(L)=POP_LIM_VARIABLE
	  END DO
!
	  IF(RELAX_VARIABLE .LT. 0.0D0 .OR. RELAX_VARIABLE .GT. 2.0D0)THEN
	    WRITE(6,*)'Error for relaxation parameter in ADJUST_CORRECTIONS -- invalid value'
	    WRITE(6,*)'Valid range is 0.0 to 2.0 '
	    WRITE(6,*)'Relaxation parameter read in is',RELAX_VARIABLE
	    RELAX_VARIABLE=1.0D0
	  END IF
	  DO L=L_ST,L_END
	    RELAX_PARAM(L)=RELAX_VARIABLE
	  END DO
	END IF
!
	INQUIRE(UNIT=LU_SUM,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(LU_SUM)
	IF(LAMBDA_IT)RELAX_PARAM(1:ND)=1.0D0
	IF(LAMBDA_IT)POP_LIM(1:ND)=100.0D0*CHANGE_LIM
!
	IF(SCALE_OPT(1:5) .EQ. 'MAJOR')THEN
	  DO I=1,ND
!
	    T1=BIG_LIM                  !Prevents division by zero and insures
	    T2=LIT_LIM                  !SCALE=1 if small changes.
	    DO J=1,NT-1
	      IF(POPS(J,I) .GT. 1.0D-10*POPS(NT-1,I))THEN
	        T1=MAX(T1,STEQ(J,I))            !Note + means decrease
	        T2=MIN(T2,STEQ(J,I))            !Note - means increase
	      END IF
	    END DO
	    SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2 )
!
! Limit the change in T to a maximum of 20%, and ensure T > T_MIN.
!
	    T3=MAX( T_LIM(I),ABS(STEQ(NT,I)) )
	    SCALE=MIN( T_LIM(I)/T3,SCALE )
	    MIN_SCALE=MIN(SCALE,MIN_SCALE)
	    IF(STEQ(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0D0-STEQ(NT,I)*SCALE) .LT. T_MIN)THEN
	      SCALE=(1.0D0-T_MIN/POPS(NT,I))/STEQ(NT,I)
	    END IF
	    IF(SCALE .GT. 1.0D0)SCALE=1.0D0             !i.e., will not force T to T_MIN
!
! RELAX_PARAM allows for the use of successive over or under relaxation.
! When RELAX_PARAM > 1, BIG_LIM and LIT_LIM ensure that we don't get
! negatve populations.
!
	    IF(SCALE .EQ. 1.0D0)SCALE=RELAX_PARAM(I)
!
! Ensure population change doesn't change population by too large an amount.
! POP_LIM allows us to force smaller corrections at some depths even while
! CMFGEN is running.
!
	    T2=(POP_LIM(I)-1.0D0)/POP_LIM(I)
	    DPTH_BIG_LIM=MIN(BIG_LIM,T2)
	    DPTH_LIT_LIM=MAX(LIT_LIM,1.0D0-POP_LIM(I))
	    DO J=1,NT
	      T1=STEQ(J,I)*SCALE
	      IF(T1 .GT. DPTH_BIG_LIM)T1=DPTH_BIG_LIM
	      IF(T1 .LT. DPTH_LIT_LIM)T1=DPTH_LIT_LIM
	      POPS(J,I)=POPS(J,I)*(1.0D0-T1)
	    END DO
	  END DO
	  WRITE(6,'(A,ES12.4)')' The minimum value of scale for Major species is:',MIN_SCALE
	END IF
!
! Only adjust populatons at adjacent depths when the corrections are rediculously large
! (e.g. -T1 or 1-1/T1.). If CONSISTENCY_CNT .LE. 0, this is never done.
!
	IF(LAMBDA_IT)COUNTER=COUNTER+1
	IF(LAMBDA_IT .AND. COUNTER .EQ. CONSISTENCY_CNT)THEN
	  WRITE(6,*)'Calling SET_DEPTH_CONSISTENCY'
	  DC_OPTION=' ' 		!Not in use at the present time
	  CALL SET_DEPTH_CONSISTENCY(STEQ,POPS,ND,NT,
	1           BAD_INCREASE_LIMIT,BAD_DECREASE_LIMIT,DC_OPTION)
	  COUNTER=0
	END IF
!
	RETURN
	END
