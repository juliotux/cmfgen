C
C Dummy routine to replave VMS TUNE rutine for collecting
C tuning statistics for a program.
C
C  CPUTOT() : CPU time accumulated in each IDENT
C  RUNTOT() : Total CPU elapsed time of the program at the point where
C             IDENT has finished (LRUN=2)
C
C  
C  we use unit LU=55 to write out data 
C
        SUBROUTINE TUNE(LRUN,IDENT)
        IMPLICIT NONE
	INTEGER LRUN
	CHARACTER*(*) IDENT
C
        REAL*8 T0,OVERHEAD
        INTEGER LCALL
        INTEGER MAX_IDS
        PARAMETER (MAX_IDS=50)
        REAL*8 BEFORE(MAX_IDS),AFTER(MAX_IDS),CPUTOT(MAX_IDS)
        REAL*8 RUNTOT(MAX_IDS)
        CHARACTER*30 IDLIST(MAX_IDS)
        INTEGER I,J

        REAL*4 ETIME,TARRY(2)
        EXTERNAL ETIME

	LOGICAL*4 FIRSTTIME
	DATA FIRSTTIME/.TRUE./
        SAVE FIRSTTIME
        SAVE BEFORE,AFTER,CPUTOT,RUNTOT
        SAVE IDLIST

        IF (FIRSTTIME) THEN
          OVERHEAD=0.D0
          FIRSTTIME=.FALSE.
          DO  I=1,MAX_IDS
              BEFORE(I)=0.D0
              AFTER(I)=0.D0
              CPUTOT(I)=0.D0
              RUNTOT(I)=0.D0
	      IDLIST(I)=' '
          END DO
        ENDIF
	IF (LRUN.EQ.1) THEN
		DO 300 I=1,MAX_IDS
			IF (IDENT .EQ. IDLIST(I)) GO TO 310
			IF (IDLIST(I).EQ.' ') THEN
				IDLIST(I)=IDENT
				GO TO 310
			ENDIF
300		CONTINUE
		WRITE (55,210)
210		FORMAT(' ***** TOO MANY TUNING POINTS ')
310		CONTINUE
                T0 =ETIME(TARRY)
                OVERHEAD = ETIME(TARRY) - T0
                BEFORE(I) = ETIME(TARRY)
	ELSE IF (LRUN.EQ.2) THEN
		DO 400 I=1,MAX_IDS
			IF (IDENT.EQ.IDLIST(I)) GO TO 410
			IF (IDLIST(I).EQ.' ') GO TO 400
400		CONTINUE
		WRITE (55,211)
211		FORMAT(' ***** UNMATCHED TUNING POINT ')
		RETURN
410		CONTINUE
                AFTER(I)=ETIME(TARRY)
                CPUTOT(I)=CPUTOT(I)+(AFTER(I)-BEFORE(I))-OVERHEAD
		RUNTOT(I)=AFTER(I)
C
	ELSE IF (LRUN.EQ.3) THEN
		WRITE(55,204)
204		FORMAT(8X,'Identifier',11x,'Elapsed',11x,'  CPU')
		WRITE(55,205)
205		FORMAT(29x,'  Time ',11x,'  Time')
		DO 500 I=1,MAX_IDS
			IF (IDLIST(I).EQ.' ') GO TO 501
500		CONTINUE
501		I=I-1
	        WRITE(55,200)(IDLIST(J),RUNTOT(J),CPUTOT(J),J=1,I)
200		FORMAT((1X,a24,f15.6,2X,f15.6))
C
	ELSE
		WRITE (55,201)
201		FORMAT(' ***** ILLEGAL VALUE OF LRUN IN CALL TO TUNE ')
		STOP
	ENDIF
        
	RETURN
	END
