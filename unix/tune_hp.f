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
C           RUNTOT() : Total CPU elapsed time of the program at the point where
C                            IDENT has finished (LRUN=2)
C
C NB: Total time for each code section is accumulated on successive calls.
C i.e. TUNE(1,'Unique ID') does not initialize the counters.
C
        SUBROUTINE TUNE(LRUN,IDENT)
        IMPLICIT NONE
	INTEGER LRUN
	CHARACTER*(*) IDENT
!
	INTEGER, PARAMETER :: MAX_IDS=50
	INTEGER, PARAMETER :: LUOUT=55
!
        REAL*8 T0,OVERHEAD
        REAL*8 BEFORE(MAX_IDS),AFTER(MAX_IDS),CPUTOT(MAX_IDS)
        REAL*8 RUNTOT(MAX_IDS)
        CHARACTER*30 IDLIST(MAX_IDS)
        INTEGER I
C
	LOGICAL*4 FIRSTTIME
	DATA FIRSTTIME/.TRUE./
        SAVE FIRSTTIME,OVERHEAD
        SAVE BEFORE,AFTER,CPUTOT,RUNTOT
        SAVE IDLIST

        IF (FIRSTTIME) THEN
          FIRSTTIME=.FALSE.
          DO  I=1,MAX_IDS
            BEFORE(I)=0.D0
            AFTER(I)=0.D0
            CPUTOT(I)=0.D0
            RUNTOT(I)=0.D0
	    IDLIST(I)=' '
          END DO
          T0=MCLOCK()
          OVERHEAD=2.0D-06*(MCLOCK()-T0)
	  OPEN(UNIT=LUOUT,STATUS='REPLACE',FILE='TIMING')
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,*)'Overhead is ',OVERHEAD
	  WRITE(LUOUT,*)' '
        ENDIF
!
! If LRUN =1, we are beginning the TIME bracket. Therefore we find the
! correct storage location first. 
!
	IF (LRUN.EQ.1) THEN
	  DO I=1,MAX_IDS
            IF (IDENT .EQ. IDLIST(I))THEN
	      BEFORE(I)=1.0D-06*MCLOCK()
	      RETURN	
	    END IF
	    IF (IDLIST(I) .EQ. ' ') THEN
	      IDLIST(I)=IDENT
	      BEFORE(I)=1.0D-06*MCLOCK()
	      RETURN	
	    END IF
	  END DO
	  WRITE (LUOUT,'(A)')' ***** TOO MANY TUNING POINTS '
	  RETURN
	ELSE IF (LRUN.EQ.2) THEN
!
! If LRUN=2, we are ending the TIME bracket. Therefore we call the timing 
! routine first.
!
          T0=1.0D-06*MCLOCK()
          DO I=1,MAX_IDS
	    IF (IDENT.EQ.IDLIST(I))THEN
              AFTER(I)=T0
	      CPUTOT(I)=CPUTOT(I)+(AFTER(I)-BEFORE(I)-OVERHEAD)
	      RUNTOT(I)=T0
	      RETURN
	    END IF 
	    IF (IDLIST(I).EQ.' ')EXIT
	  END DO
	  WRITE (LUOUT,*)' ***** UNMATCHED TUNING POINT ',TRIM(IDENT)
	  RETURN
C
	ELSE IF (LRUN.EQ.3) THEN
	  WRITE(LUOUT,'(8X,''Identifier'',11x,''Elapsed'',11x,''  CPU'')')
	  WRITE(LUOUT,'(29x,''  Time '',11x,''  Time'')')
	  DO I=1,MAX_IDS
	    IF (IDLIST(I).EQ.' ') EXIT
	    WRITE(LUOUT,'(1X,a24,f15.6,2X,f15.6)')
	1   IDLIST(I),RUNTOT(I),CPUTOT(I)
          END DO
	ELSE
	  WRITE (LUOUT,'(A)')' ***** ILLEGAL VALUE OF LRUN IN CALL TO TUNE '
	  STOP
	ENDIF
        
	RETURN
	END
