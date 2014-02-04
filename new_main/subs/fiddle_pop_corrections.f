!
! This routine has been designed to use with SOLVEBA_V9. At present only the
! MAJOR otpion is treated. This section was made a subroutine to faciliate easy modification
! in the hope that we could improve convergence. With the use of the file 'ADJUST_CORRECTIONS'
! we can adjust the relaxation parameter at specific depths, and on the fly.
!
	SUBROUTINE FIDDLE_POP_CORRECTIONS(POPS,STEQ,T_MIN,CHANGE_LIM,
	1              SCALE_OPT,LAMBDA_IT,LU_SUM,NT,ND)
	IMPLICIT NONE
!
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
	LOGICAL LAMBDA_IT
	CHARACTER(LEN=*) SCALE_OPT
!
! Local variables
!
	REAL*8 RELAX_PARAM(ND)
	REAL*8 POP_LIM(ND)
	REAL*8 T_LIM(ND)
!
	REAL*8 LIT_LIM,BIG_LIM
	REAL*8 T1,T2,T3
	REAL*8 MIN_SCALE
	REAL*8 SCALE
!
	INTEGER IOS
	INTEGER I,J,IC
	INTEGER L,L_ST,L_END
	LOGICAL FILE_OPEN
	CHARACTER(LEN=80)STRING
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
	BIG_LIM=(CHANGE_LIM-1.0D0)/CHANGE_LIM
        LIT_LIM=1.0D0-CHANGE_LIM
!
! Set default parameters.
!
	DO I=1,ND
	  RELAX_PARAM(I)=1.0D0
	  T_LIM(I)=0.2D0
	  POP_LIM(I)=1.0D+10
	END DO
!
! ADJUST_CORRECTIONS should contain a list of depths with a parameter: The parameter is
! taken to be the relaxation parameter unless the key strings "T LIMIT" or "POP LIMIT" are
! present. The validity of the depth range and parameters are checked.
!
! Valid ranges:
!              Relaxations      0 < T1 < 2
!              T LIMIT          0 < T1 < 0.2 
!              POP LIMIT        0 < T1 < 0.8 (really only valid if T1 < 0.3)
!
! The file need not be present, and read errors will cuass a skip of additional reads.
! Lines beginning with ! and blank lines are ignored.
!
! NB: The relaxation parameter and the T LIMIT effect all corrections (if used) at the
!     specified depth. POP LIMIT only effects corrections bigger than POP LIMIT.
!
	OPEN(UNIT=LU_SUM,FILE='ADJUST_CORRECTIONS',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .EQ. 0)THEN
	  DO WHILE(1 .EQ. 1)
	    READ(LU_SUM,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)EXIT
	    IF(STRING .NE. ' ' .AND. STRING(1:1) .NE. '!')THEN
	      IC=INDEX(STRING,':')
	      IF(IC .EQ. 0)THEN
	        READ(STRING,*,IOSTAT=IOS)L_ST,T1
	        L_END=L_ST
	      ELSE
	        STRING(IC:IC)=' '
	        READ(STRING,*,IOSTAT=IOS)L_ST,L_END,T1
	      END IF
!
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Error reading ADJUST_CORRECTIONS in fiddle_pop_corrections'
	        WRITE(6,*)'IOS=',IOS
	        EXIT
	      END IF
	      IF(L_ST .LT. 1 .OR. L_END .GT. ND)THEN
	        WRITE(6,*)'Error in depth indices in ADJUST_CORRECTIONS -- invalid range'
	        WRITE(6,*)'Depth indices=',L_ST,L_END
	        EXIT
	      END IF
!
	      CALL SET_CASE_UP(STRING,IONE,IZERO)
	      IF(INDEX(STRING,'T LIMIT') .NE. 0)THEN
	        IF(T1 .LT. 0.0D0 .OR. T1 .GT. 0.20D0)THEN
	          WRITE(6,*)'Error for T LIMIT in ADJUST_CORRECTIONS -- invalid value'
	          WRITE(6,*)'Valid range is 0.0 to 0.2'
	          WRITE(6,*)'T LIMIT parameter read in is',T1
	          EXIT
	        END IF
	        DO L=L_ST,L_END
	           T_LIM(L)=T1
	        END DO
	      ELSE IF(INDEX(STRING,'POP LIMIT') .NE. 0)THEN
	        IF(T1 .LT. 0.0D0 .OR. T1 .GT. 0.8D0)THEN
	          WRITE(6,*)'Error for POP LIMIT in ADJUST_CORRECTIONS -- invalid value'
	          WRITE(6,*)'Valid range is 0.0 to 0.8'
	          WRITE(6,*)'POP LIMIT parameter read in',T1
	          EXIT
	        END IF
	        DO L=L_ST,L_END
	          POP_LIM(L)=T1
	        END DO
	      ELSE 
	        IF(T1 .LT. 0.0D0 .OR. T1 .GT. 2.0D0)THEN
	          WRITE(6,*)'Error for relaxation parameter in ADJUST_CORRECTIONS -- invalid value'
	          WRITE(6,*)'Valid range is 0.0 to 2.0 '
	          WRITE(6,*)'Relaxation parameter read in is',T1
	          EXIT
	        END IF
	        DO L=L_ST,L_END
	           RELAX_PARAM(L)=T1
	        END DO
	      END IF
	    END IF
	  END DO
	END IF
!
	INQUIRE(UNIT=LU_SUM,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(LU_SUM)
	IF(LAMBDA_IT)RELAX_PARAM(1:ND)=1.0D0
	IF(LAMBDA_IT)POP_LIM(1:ND)=1.0D+10
!
	IF(SCALE_OPT(1:5) .EQ. 'MAJOR')THEN
	  DO I=1,ND
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
	    IF(RELAX_PARAM(I) .LT. SCALE)SCALE=RELAX_PARAM(I)
	    DO J=1,NT
	      T1=STEQ(J,I)*SCALE
	      IF(T1 .GT. BIG_LIM)T1=BIG_LIM
	      IF(T1 .LT. LIT_LIM)T1=LIT_LIM
	      IF(T1 .GT. POP_LIM(I))T1=POP_LIM(I)
	      IF(T1 .LT. -POP_LIM(I))T1=-POP_LIM(I)
	      POPS(J,I)=POPS(J,I)*(1.0D0-T1)
	    END DO
	  END DO
	  WRITE(6,'(A,ES12.4)')' The minimum value of scale for Major species is:',MIN_SCALE
	END IF
!
	RETURN
	END
