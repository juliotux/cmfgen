!
! Routine to merge N COLOR PGPLOTS onto a single page. The plots should be 
! done in LANDSCAPE MODE with CPS as printer: The final plot is in
! Portrait mode.
!
!      For N=2: EXPAND_CHAR=2.0; EXPAND_TICK=2.0; ASR=0.35; Plot Size=20 cm
!      For N=3: EXPAND_CHAR=3.0; EXPAND_TICK=3.0; ASR=0.25; Plot Size=20 cm
!
	PROGRAM LAND_COL_MERGE
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 18-Apr-2014 : Fixed so that top plot input first.
!                         Error fixed when input file does not exist.
! Altered 06-Jan-2001 : Automatic file naming. Top plot input first.
! Finalized 29-May-1997
!
	CHARACTER*132 STRING
	CHARACTER*80 FILE1
	CHARACTER*80 OUTF
	CHARACTER*10 TMP_STR
!
	INTEGER, PARAMETER :: IZERO=0 
	INTEGER, PARAMETER :: LU_IN=11
	INTEGER, PARAMETER :: LU_OUT=30
	INTEGER, PARAMETER :: LU_TERM=6
!
	INTEGER IOS,I,J,K,IREC,M
	INTEGER N_PLTS
!
	LOGICAL FILE_EXISTS
	LOGICAL OVER_WRITE
!
	WRITE(LU_TERM,*)' '
	WRITE(LU_TERM,*)'For N=2: EXPAND_CHAR=1.6; EXPAND_TICK=1.6; ASR=0.35; Plot Size=20 cm'
	WRITE(LU_TERM,*)'For N=3: EXPAND_CHAR=3.0; EXPAND_TICK=3.0; ASR=0.25; Plot Size=20 cm'
	WRITE(LU_TERM,*)'For N=3: EXPAND_CHAR=3.0; EXPAND_TICK=3.0; ASR=0.23; Plot Size=20 cm (xlab on all plots)'
!
	FILE1='pgplot_1.ps'
	DO WHILE(1 .EQ. 1)
	  WRITE(LU_TERM,*)' '
	  N_PLTS=2; OUTF='merged' 
	  CALL GEN_IN(N_PLTS,'Number of plots to merge (0 to exit)')
	  IF(N_PLTS .LE. 0)STOP
!
10	  CALL GEN_IN(OUTF,'Output file')
	  INQUIRE(FILE=OUTF,EXIST=FILE_EXISTS)
	  IF(FILE_EXISTS)THEN
	    OVER_WRITE=.FALSE.
	    CALL GEN_IN(OVER_WRITE,'File already exists -- do you wish to overwrite?')
	    IF(OVER_WRITE)THEN
	      CALL GEN_ASCI_OPEN(LU_OUT,OUTF,'OLD',' ',' ',IREC,IOS)
	    ELSE
	      GOTO 10
	    END IF
	  ELSE
	    CALL GEN_ASCI_OPEN(LU_OUT,OUTF,'NEW',' ',' ',IREC,IOS)
	  END IF
!
	  IREC=0
	  IOS=0
!
20	  CALL GEN_IN(FILE1,'Top PGPLOT file')
	  CALL GEN_ASCI_OPEN(LU_IN,FILE1,'OLD',' ','READ',IREC,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LU_TERM,*)'Error opening 1st input FILE:',FILE1
	      GOTO 20
	    END IF
	    STRING=' '
	    DO WHILE( INDEX(STRING,'0.072 0.072 scale') .EQ. 0)
	      READ(LU_IN,'(A)')STRING
	      WRITE(LU_OUT,'(A)')TRIM(STRING)
	    END DO
	    READ(LU_IN,'(A)')STRING   		!Output Translate/rotation string
	    WRITE(LU_OUT,'(A)')TRIM(STRING)
!
	    IF(N_PLTS .EQ. 2)THEN
	      WRITE(LU_OUT,'(A)')'  0 3100 translate'
	    ELSE IF(N_PLTS .EQ. 3)THEN
	      WRITE(LU_OUT,'(A)')'  0 3100 translate'
	    ELSE IF(N_PLTS .EQ. 4)THEN
	      WRITE(LU_OUT,'(A)')'  0 3600 translate'
	    ELSE
	      WRITE(LU_OUT,'(A)')'  0 3600 translate'
	    END IF
	    DO WHILE(1 .EQ. 1)
	      READ(LU_IN,'(A)')STRING
	      IF(INDEX(STRING,'PGPLOT restore showpage') .NE. 0)THEN
	        GOTO 100
	      END IF
	      WRITE(LU_OUT,'(A)')TRIM(STRING)
	    END DO
100	    CONTINUE
	  CLOSE(UNIT=LU_IN)
C
	  DO M=2,N_PLTS
	    IF( N_PLTS .EQ. 2)THEN
	      WRITE(LU_OUT,'(A)')' 0 -3600 translate'
	    ELSE IF( N_PLTS .EQ. 3)THEN
	      WRITE(LU_OUT,'(A)')' 0 -2500 translate'
	    ELSE IF( N_PLTS .EQ. 4)THEN
	      WRITE(LU_OUT,'(A)')' 0 -1800 translate'
	    ELSE
	      I=-(7200)/N_PLTS
	      WRITE(LU_OUT,'(A,I5,A)')' 0 ',I,' translate'
	    END IF
!
! Update file name in a systematic way to save typing.
! NB: A blank filename means that we have no more files.
!
	    CALL UPDATE_PG_FILENAME(FILE1)
150	    CALL GEN_IN(FILE1,'Next PGPLOT file')
	    CALL GEN_ASCI_OPEN(LU_IN,FILE1,'OLD',' ','READ',IREC,IOS)
	    IF(FILE1 .EQ. ' ')GOTO 1000
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
200	    CONTINUE
	    CLOSE(UNIT=LU_IN)
	  END DO
!
1000	  CONTINUE
	  WRITE(LU_OUT,'(A)')STRING
	  WRITE(LU_OUT,'(A)')'%%EOF'
	  CLOSE(LU_OUT)
!
! Update file name in a systematic way to save typing.
!
	  CALL UPDATE_PG_FILENAME(FILE1)
!
	END DO				!Do another merge
!
	STOP
	END
