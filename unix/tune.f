C
C Routine to replace VMS TUNE routine for collecting tuning statistics for 
C a program section.
C
C Usage:
C          CALL TUNE(1,'Unique ID')
C              .......
C              Section of code.
C              .......
C          CALL TUNE(2,'Same Unique ID')
C
C To print results:
C
C           CALL TUNE(3,' ')
C
C Output from print:
C
C           CPUTOT() : CPU time accumulated in each IDENT
C           WALLTOT() : Total wall time bewteen each IDENT.
C
C NB: Total time for each code section is accumulated on successive calls.
C i.e. TUNE(1,'Unique ID') does not initialize the counters.
C
        SUBROUTINE TUNE(LRUN,IDENT)
        IMPLICIT NONE
!
! Altered 18-Feb-2013 : MAX_IDS increased. STACK introduced.
!                       Routine should now be much more efficient, with less instructions per call.
! Altered 08-Mar-2010 : Change variable for system clock to 8 bytes.
!                          This prevents loss of elapsed time due to clock rollover.
! Altered 11-Nov-2000 : Call to F90 SYSTEM_CLOCK routine implemented.
!                       Wall time now returened as in original VMS routine.
!                       Counters now initialized if LRUN=0 is passed.
!
	INTEGER LRUN
	CHARACTER*(*) IDENT
!
	INTEGER, PARAMETER :: MAX_IDS=80
	INTEGER, PARAMETER :: LUOUT=55
!
        REAL*8, SAVE :: T0,OVERHEAD
        REAL*8, SAVE :: ST_CPU(MAX_IDS)
	REAL*8, SAVE :: CPUTOT(MAX_IDS)
        REAL*8, SAVE :: WALLTOT(MAX_IDS)
	INTEGER, SAVE :: STACK(MAX_IDS)
        CHARACTER(LEN=20), SAVE ::  IDLIST(MAX_IDS)
!
	REAL*8, SAVE :: RR0
	INTEGER*8, SAVE :: IEND_WALL
	INTEGER*8, SAVE :: IC0,IR0,IM0,IT1
	INTEGER*8, SAVE :: IST_WALL(MAX_IDS)
	INTEGER, SAVE :: NUM_IDS
        INTEGER, SAVE :: NUM_STACK
!
	INTEGER CURRENT_ID 
	INTEGER ACTIVE_ID
        INTEGER I,LU,TERM_OUT
	EXTERNAL TERM_OUT
!
	REAL*4 ETIME,TARRY(2)
!
	LOGICAL, SAVE :: FIRST_TOO_MANY
	LOGICAL, SAVE :: FIRST_UNMATCHED
	LOGICAL, SAVE :: FIRSTTIME
	DATA FIRSTTIME/.TRUE./
	DATA FIRST_TOO_MANY/.TRUE./
	DATA FIRST_UNMATCHED/.TRUE./
!
        IF (FIRSTTIME)THEN
          FIRSTTIME=.FALSE.
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.D0
            IST_WALL(I)=0.D0
            CPUTOT(I)=0.D0
            WALLTOT(I)=0.D0
	    IDLIST(I)=' '
	    STACK(I)=1
          END DO
	  ACTIVE_ID=0
	  NUM_IDS=0
	  NUM_STACK=0
	  CALL SYSTEM_CLOCK(IC0,IR0,IM0);    RR0=IR0
          T0=ETIME(TARRY)
          OVERHEAD=2.0D0*(ETIME(TARRY)-T0)
	  OPEN(UNIT=LUOUT,STATUS='REPLACE',FILE='TIMING')
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,*)'Overhead is ',OVERHEAD
	  WRITE(LUOUT,*)'   Count rate for wall clock is',IR0
	  WRITE(LUOUT,*)'Maximum count for wall clock is',IM0
	  WRITE(LUOUT,*)' '
        ENDIF
!
! If LRUN =1, we are beginning the TIME bracket. Therefore we find the
! correct storage location first. 
!
	IF (LRUN .EQ. 1) THEN
	  DO I=NUM_IDS,1,-1
            IF (IDENT .EQ. IDLIST(I))THEN
	      NUM_STACK=NUM_STACK+1
	      STACK(NUM_STACK)=I
	      CALL SYSTEM_CLOCK(IST_WALL(I))
	      ST_CPU(I)=ETIME(TARRY)
	      RETURN	
	    END IF
	  END DO
	  NUM_IDS=NUM_IDS+1
	  IF(NUM_IDS .LE. MAX_IDS)THEN 
	    I=NUM_IDS
	    IDLIST(I)=IDENT
	    NUM_STACK=NUM_STACK+1
	    STACK(NUM_STACK)=I
	    CALL SYSTEM_CLOCK(IST_WALL(I))
	    ST_CPU(I)=ETIME(TARRY)
	    RETURN	
	  END IF
	  IF(FIRST_TOO_MANY)THEN
	    WRITE (LUOUT,'(A)')' ***** TOO MANY TUNING POINTS '
	    FIRST_TOO_MANY=.FALSE.
	  END IF
	  RETURN
!
	ELSE IF (LRUN .EQ. 2) THEN
!
! If LRUN=2, we are ending the TIME bracket. Therefore we call the timing 
! routine first.
!
! If TUNE has called been called correctly, then STACK should always be set correctly. 
!
          T0=ETIME(TARRY)
	  CALL SYSTEM_CLOCK(IEND_WALL)
	  ACTIVE_ID=STACK(NUM_STACK)
	  IF (IDENT .EQ. IDLIST(ACTIVE_ID))THEN
	    CURRENT_ID=ACTIVE_ID
	    NUM_STACK=NUM_STACK-1
	  ELSE
	    LU=TERM_OUT()
	    WRITE(LU,*)'Error in TUNE: STACK not correct'
	    WRITE(LU,*)'Make sure inner TUNE section is fully contained in outer TUNE section'
            WRITE(LU,*)LRUN,TRIM(IDENT)
	    WRITE(LU,*)' '
	    WRITE(LU,*)'A printout of the STACK follows'
	    WRITE(LU,*)' '
	    DO I=1,NUM_STACK
              WRITE(LU,'(I7,I7,T20,A)')I,STACK(I),TRIM(IDLIST(STACK(I)))
	    END DO
	    STOP
	    CURRENT_ID=0
            DO I=1,MAX_IDS
	      IF (IDENT.EQ.IDLIST(I))THEN
	        CURRENT_ID=I
	        EXIT
	      END IF
	    END DO
	  END IF
!
	  IF(CURRENT_ID .NE. 0)THEN
	    I=CURRENT_ID
	    CPUTOT(I)=CPUTOT(I)+(T0-ST_CPU(I)-OVERHEAD)
	    IT1=IEND_WALL-IST_WALL(I)
	    IF(IT1 .LT. 0)IT1=IT1+IM0
	    WALLTOT(I)=WALLTOT(I)+IT1/RR0
	  ELSE IF(FIRST_UNMATCHED)THEN
	    FIRST_UNMATCHED=.FALSE.
	    WRITE(LUOUT,*)' ***** UNMATCHED TUNING POINT: ',TRIM(IDENT)
	  END IF
	  RETURN
C
	ELSE IF (LRUN .EQ. 3) THEN
	  WRITE(LUOUT,'(8X,''Identifier'',11x,''Elapsed'',11x,''  CPU'')')
	  WRITE(LUOUT,'(29x,''  Time '',11x,''  Time'')')
	  DO I=1,MAX_IDS
	    IF (IDLIST(I).EQ.' ') EXIT
	    WRITE(LUOUT,'(1X,a24,f15.6,2X,f15.6)')
	1   IDLIST(I),WALLTOT(I),CPUTOT(I)
          END DO
	  CLOSE(LUOUT)
	  OPEN(UNIT=LUOUT,STATUS='OLD',ACTION='WRITE',POSITION='APPEND',FILE='TIMING')
	ELSE IF(LRUN .EQ. 0) THEN
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.D0
            IST_WALL(I)=0.D0
            CPUTOT(I)=0.D0
            WALLTOT(I)=0.D0
	    IDLIST(I)=' '
	    STACK(I)=1
          END DO
	  ACTIVE_ID=0
	  NUM_IDS=0
	  NUM_STACK=0
	ELSE
	  WRITE (LUOUT,'(A)')' ***** ILLEGAL VALUE OF LRUN IN CALL TO TUNE '
	  WRITE(LUOUT,*)' LRUN=',LRUN
	  STOP
	ENDIF
        
	RETURN
	END
