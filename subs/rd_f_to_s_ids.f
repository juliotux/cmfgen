C
C Routine to read in the vector describing the matching of actual atomic levels
C to that in the model atom with super levels.
C
	SUBROUTINE RD_F_TO_S_IDS(F_TO_S,INT_SEQ,
	1               LEVNAME_F,N_F,N_S,LUIN,FILENAME)
	IMPLICIT NONE
C
C Altered 30-Sep-2005 : Removed warning about mixed parity.
C Altered 31-Jan-1997 : Only a warning is output if SUPER-LEVELS are defined
C                          with mixed parity.
C                       Improved error reporting for bad super-level ID.
C Altered 29-Dec-1996 : Dynamic allocation used for CNT, PAR
C Altered 29-May-1996 : String increased to 132 to accommodate longer names.
C Altered 26-May-1996 : ACTION installed in OPEN statement.
C Altered  2-Jan-1996 : INT_SEQ variable installed.
C Altered 13-Oct-1995 : If N_F=N_S 1 to 1 mapping of levels is assumed.
C Created 31-Mar-1995
C
	INTEGER N_S,N_F,LUIN
	INTEGER F_TO_S(N_F),INT_SEQ(N_F)
	CHARACTER*(*) LEVNAME_F(N_F),FILENAME
	CHARACTER*132 STRING
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER LUER,I,J,K,IOS,ENTRY_NUM
	INTEGER CNT(N_S),PAR(N_S)
C
	LUER=ERROR_LU()
C
C IF N_F .EQ. N_S there must be 1 to 1 correspondence between the levels.
C Thus we do not need to open a link file.
C
	IF(N_F .EQ. N_S)THEN
	  DO I=1,N_F
	    F_TO_S(I)=I
	    INT_SEQ(I)=0
	  END DO
	  RETURN
	END IF
C
	OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS,
	1       ACTION='READ')
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',FILENAME,' in RD_F_TO_S_IDS'
	  STOP
	END IF
C
	STRING=' '
	DO WHILE( INDEX(STRING,'Number of energy levels') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	BACKSPACE(LUIN)
	READ(LUIN,*)I
	IF(I .LT. N_F)THEN
	  WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	  WRITE(LUER,*)'Currently treating ',FILENAME
	  WRITE(LUER,*)'Insufficient levels in file'
	  STOP
	END IF
C
	STRING=' '
	DO WHILE( INDEX(STRING,'Entry number of link') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	BACKSPACE(LUIN)
	READ(LUIN,*)ENTRY_NUM
	IF(ENTRY_NUM .LT. 2 .OR. ENTRY_NUM .GT. 20)THEN
	  WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	  WRITE(LUER,*)'Currently treating ',FILENAME
	  WRITE(LUER,*)'Bad entry number for level link'
	  STOP
	END IF
C
C NB: All entries must be separated by at LEAST 2 spaces.
C
	READ(LUIN,'(A)')STRING		!Blankline
	DO I=1,N_F
	  F_TO_S(I)=0.
	  READ(LUIN,'(A)')STRING
	  DO WHILE(STRING(1:1) .EQ. ' ')	!Strip leading blanks.
	    STRING(1:)=STRING(2:)
	  END DO
	  J=INDEX(STRING,'  ')
	  IF(STRING(1:J) .NE. LEVNAME_F(I))THEN
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Invalid level id'
	    WRITE(LUER,*)'STRING=',STRING(1:J)
	    STOP
	  END IF
C
C NB: ENTRY_NUM-2 as
C              -1 due to 1st entry being level ID
C              -1 as only need to point to entry.
C
	  DO K=1,ENTRY_NUM-2
	    STRING(1:)=STRING(J:)
	    DO WHILE(STRING(1:1) .EQ. ' ')
	      STRING(1:)=STRING(2:)
	    END DO
	    J=INDEX(STRING,'  ')
	  END DO
	  READ(STRING(J:),*)F_TO_S(I),INT_SEQ(I)
	END DO
C
C Now check that:
C (1) All super-levels have at least 1 level from full atom
C (2) Correct number of super levels.
C (3) Parity of super levels matched (Warning only).
C
	DO I=1,N_F
	  IF(F_TO_S(I) .LE. 0 .OR. F_TO_S(I) .GT. N_S)THEN
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Incorrect super level ID for level ',LEVNAME_F(I)
	    STOP
	  END IF
	END DO
C
C Check that all super levels have at least 1 level.
C
	DO I=1,N_S
	  CNT(I)=0
	END DO
	DO J=1,N_F
	  I=F_TO_S(J)
	  CNT(I)=CNT(I)+1
	END DO
	DO I=1,N_S
	  IF(CNT(I) .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RD_F_TO_S_IDS'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Super level ',I,' has no components'
	    STOP
	  END IF
	END DO
C
C Check parity.
C
	DO I=1,N_S
	  PAR(I)=0
	END DO
	DO J=1,N_F
	  I=F_TO_S(J)
	  IF(INDEX(LEVNAME_F(J),'e') .NE. 0)THEN
	     PAR(I)=PAR(I)+1
	  ELSE
	     PAR(I)=PAR(I)-1
	  END IF
	END DO
!
!	J=0
!	DO I=1,N_S
!	  PAR(I)=ABS(PAR(I))/CNT(I)
!	  IF(PAR(I) .NE. 1)THEN
!	    IF(J .EQ. 0)THEN
!	      WRITE(LUER,*)'Warning in RD_F_TO_S_IDS'
!	      WRITE(LUER,*)'Currently treating ',FILENAME
!	    END IF
!	    J=J+1
!	    WRITE(LUER,'(4X,A,I5,A)')'Super level ',
!	1                 I,' has a mix of odd and even levels'
!	  END IF
!	END DO
!
	RETURN
	END
