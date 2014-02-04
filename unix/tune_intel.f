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
! Altered 11-Nov-2000 : Call to F90 SYSTEM_CLOCK routine implemented.
!                       Wall time now returened as in original VMS routine.
!                       Counters now initialized if LRUN=0 is passed.
!
	INTEGER LRUN
	CHARACTER*(*) IDENT
!
	INTEGER, PARAMETER :: MAX_IDS=50
	INTEGER, PARAMETER :: LUOUT=55
!
        REAL*8 OVERHEAD, StartTime, StopTime
        REAL*8 ST_CPU(MAX_IDS)
	REAL*8 CPUTOT(MAX_IDS)
        REAL*8 WALLTOT(MAX_IDS)
        CHARACTER*30 IDLIST(MAX_IDS)
        INTEGER I
!
	INTEGER IEND_WALL
	INTEGER IC0,IR0,IM0,IT1
	INTEGER IST_WALL(MAX_IDS)
	REAL*8 CLK_PERIOD,RR0
!
!
	LOGICAL FIRSTTIME
	DATA FIRSTTIME/.TRUE./
        SAVE FIRSTTIME,OVERHEAD
        SAVE ST_CPU,IST_WALL,CPUTOT,WALLTOT
	SAVE IC0,IR0,IM0,RR0
        SAVE IDLIST

        IF (FIRSTTIME)THEN
          FIRSTTIME=.FALSE.
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.D0
            IST_WALL(I)=0.D0
            CPUTOT(I)=0.D0
            WALLTOT(I)=0.D0
	    IDLIST(I)=' '
          END DO
	  CALL SYSTEM_CLOCK(IC0,IR0,IM0);    RR0=IR0
!          T0=ETIME(TARRY)
          CALL CPU_TIME(StartTime)
          CALL CPU_TIME(StopTime)
          OVERHEAD=2.0D0*(StopTime - StartTime)
	  OPEN(UNIT=LUOUT,STATUS='REPLACE',FILE='TIMING')
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,*)'Overhead is ',OVERHEAD
	  WRITE(LUOUT,*)' '
        ENDIF
!
! If LRUN =1, we are beginning the TIME bracket. Therefore we find the
! correct storage location first. 
!
	IF (LRUN .EQ. 1) THEN
	  DO I=1,MAX_IDS
            IF (IDENT .EQ. IDLIST(I))THEN
	      CALL SYSTEM_CLOCK(IST_WALL(I))
              CALL CPU_TIME(StartTime)
	      ST_CPU(I)=StartTime
	      RETURN	
	    END IF
	    IF (IDLIST(I) .EQ. ' ') THEN
	      IDLIST(I)=IDENT
	      CALL SYSTEM_CLOCK(IST_WALL(I))
              CALL CPU_TIME(StartTime)
              ST_CPU(I)=StartTime
	      RETURN	
	    END IF
	  END DO
	  WRITE (LUOUT,'(A)')' ***** TOO MANY TUNING POINTS '
	  RETURN
	  
!
	ELSE IF (LRUN .EQ. 2) THEN
!
! If LRUN=2, we are ending the TIME bracket. Therefore we call the timing 
! routine first.
!
          CALL CPU_TIME(StopTime)
	  CALL SYSTEM_CLOCK(IEND_WALL)
          DO I=1,MAX_IDS
	    IF (IDENT.EQ.IDLIST(I))THEN
	      CPUTOT(I)=CPUTOT(I)+(StopTime-ST_CPU(I)-OVERHEAD)
	      IT1=IEND_WALL-IST_WALL(I)
	      IF(IT1 .LT. 0)IT1=IT1+IM0
	      WALLTOT(I)=WALLTOT(I)+IT1/RR0
	      RETURN
	    END IF 
	    IF (IDLIST(I).EQ.' ')EXIT
	  END DO
	  WRITE(LUOUT,*)' ***** UNMATCHED TUNING POINT: ',TRIM(IDENT)
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
	ELSE IF(LRUN .EQ. 0) THEN
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.D0
            IST_WALL(I)=0.D0
            CPUTOT(I)=0.D0
            WALLTOT(I)=0.D0
	    IDLIST(I)=' '
          END DO
	ELSE
	  WRITE (LUOUT,'(A)')' ***** ILLEGAL VALUE OF LRUN IN CALL TO TUNE '
	  WRITE(LUOUT,*)' LRUN=',LRUN
	  STOP
	ENDIF
        
	RETURN
	END
