!
! Routine to merge 6 color PGPLOTS onto a single 
! page. The plots should be done with:
!      EXPAND_CHAR=1.2; EXPAND_TICK=1.2; ASR=0.89; Plot Size=10.0 cm
! The plots are located in 2 columns, 2 per row.
!
! Landscape format, with CPS as printer.
!
	PROGRAM N_MULT_MERGE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Finalized 29-May-1997
!
	CHARACTER*132 STRING
	CHARACTER*80 FILE1
	CHARACTER*80 OUTF
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: LU_IN=11
	INTEGER, PARAMETER :: LU_OUT=30
	INTEGER, PARAMETER :: LU_TERM=6
!
	INTEGER IOS,I,J,K,IREC,M
	INTEGER NR,NC
	INTEGER NPLTS,NCOLS,NROWS
!
	CHARACTER*10 TMP_STR
	CHARACTER*4 HOR_OFFSET
	CHARACTER*4 VER_OFFSET
	LOGICAL TEST
!
	WRITE(LU_TERM,*)' '
	WRITE(LU_TERM,*)'Program to merge N plots in C columns'
	WRITE(LU_TERM,*)'2x2 plots (cr): ',
	1    'EXPAND_CHAR=1.5; EXPAND_TICK=1.5; ASR=0.35; Plot Size=10 cm'
	WRITE(LU_TERM,*)' '
!
	NCOLS=2;  CALL GEN_IN(NCOLS,'Number of columns on page')
	NROWS=2;  CALL GEN_IN(NROWS,'Number of plot rows')
	NPLTS=NCOLS*NROWS
!
	I=10000/NCOLS
	WRITE(HOR_OFFSET,'(I4)')I
	I=7000/NROWS
	WRITE(VER_OFFSET,'(I4)')I
	IF(NCOLS .EQ. 2 .AND. NROWS .EQ. 2)THEN
	  HOR_OFFSET='5000'
	  VER_OFFSET='3500'
	END IF
	WRITE(6,*)'Default values only tested for 2,2'
	CALL GEN_IN(HOR_OFFSET,'Horizontal offset')
	CALL GEN_IN(VER_OFFSET,'Vertical offset')
!
	IREC=0
	IOS=0
	CALL GEN_IN(OUTF,'Output file')
	CALL GEN_ASCI_OPEN(LU_OUT,OUTF,'NEW',' ',' ',IREC,IOS)
!
	TEST=.FALSE.
	FILE1='pgplot.ps'
10	CALL GEN_IN(FILE1,'Top PGPLOT  file')
	I=INDEX(FILE1,'(test)')
	IF(I .NE. 0)THEN
	  TEST=.TRUE.
	  FILE1=FILE1(1:I-1)
	END IF
	CALL GEN_ASCI_OPEN(LU_IN,FILE1,'OLD',' ','READ',IREC,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_TERM,*)'Error opening 1st input FILE:',FILE1
	    GOTO 10
	  END IF
	  STRING=' '
	  DO WHILE( INDEX(STRING,'0.072 0.072 scale') .EQ. 0)
	    READ(LU_IN,'(A)')STRING
!	    IF(INDEX(STRING,'%%Orientation: Landscape').NE. 0)THEN
!	      STRING= '%%Orientation: Portrait'
!	    END IF
	    WRITE(LU_OUT,'(A)')TRIM(STRING)
	  END DO
	  READ(LU_IN,'(A)')STRING   !Skip Translate/rotation string
	  WRITE(LU_OUT,'(A)')STRING
	  WRITE(LU_OUT,'(A)')' -500 3000 translate'
!
	  DO WHILE(1 .EQ. 1)
	    READ(LU_IN,'(A)')STRING
	    IF(INDEX(STRING,'PGPLOT restore showpage') .NE. 0)THEN
	      GOTO 100
	    END IF
	    WRITE(LU_OUT,'(A)')TRIM(STRING)
	  END DO
100	  CONTINUE
	CLOSE(UNIT=LU_IN)
C
!
	DO NR=1,NROWS
	  DO NC=1,NCOLS
	    IF(NR .EQ. 1 .AND. NC .EQ. 1)THEN
	      GOTO 500			!Done before as special case
	    ELSE IF(NC .EQ. 1)THEN
	      WRITE(LU_OUT,'(5A)')'  -',TRIM(HOR_OFFSET),' -',TRIM(VER_OFFSET),' translate'
	      DO I=1,NCOLS-2
	        WRITE(LU_OUT,'(4A)')'  -',TRIM(HOR_OFFSET),'    0',' translate'
	      END DO
	    ELSE
	      WRITE(LU_OUT,'(4A)')'  ',TRIM(HOR_OFFSET),'    0',' translate'
	    END IF
!
! Update file name in a systematic way to save typing.
!
	  IF(.NOT. TEST)THEN
	    K=INDEX(FILE1,'_')
	    IF(K .NE. 0)K=INDEX(FILE1(K+1:),'_')+K
	    IF(K .NE. 0)THEN
	      J=K+1
	      DO WHILE(FILE1(J:J) .GE. '0' .AND. FILE1(J:J) .LE. '9')
	        J=J+1
	      END DO
	      IF(J .NE. K+1)THEN
	        READ(FILE1(K+1:J-1),*)I
	        I=I+1
	        WRITE(TMP_STR,'(I6)')I
	        TMP_STR=ADJUSTL(TMP_STR)
	        FILE1=FILE1(1:K)//TRIM(TMP_STR)//FILE1(J:)
	      END IF
	    END IF
	    CALL SET_CASE_LOW(FILE1,IZERO,IZERO)
	    IF(FILE1 .EQ. 'pgplot.ps')FILE1='pgplot_2.ps'
	  END IF
!
150	  CALL GEN_IN(FILE1,'Next PGPLOT file')
	  IF(FILE1 .EQ. ' ')GOTO 1000
	  CALL GEN_ASCI_OPEN(LU_IN,FILE1,'OLD',' ','READ',IREC,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_TERM,*)'Error opening input FILE:',FILE1
	    GOTO 150
	  END IF
	  STRING=' '
	  DO WHILE(INDEX(STRING,'%%Page: 1 1') .EQ. 0)
	    READ(LU_IN,'(A)')STRING
	  END DO
	  DO WHILE(INDEX(STRING,' EP ') .EQ. 0)
	    READ(LU_IN,'(A)')STRING
	  END DO
	  I=INDEX(STRING,' EP ')
	  WRITE(LU_OUT,'(A)')TRIM(STRING(I+4:))
	  DO WHILE(1 .EQ. 1)
	    READ(LU_IN,'(A)')STRING
	    IF(INDEX(STRING,'PGPLOT restore showpage') .NE. 0)THEN
	      GOTO 200
	    END IF
	    WRITE(LU_OUT,'(A)')TRIM(STRING)
	  END DO
200	  CONTINUE
	  CLOSE(UNIT=LU_IN)
500	  CONTINUE
	  END DO
	END DO
C
1000	CONTINUE
	WRITE(LU_OUT,'(A)')STRING
	WRITE(LU_OUT,'(A)')'%%EOF'
C
	STOP
	END
