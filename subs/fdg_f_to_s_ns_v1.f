!
! Routine to read in the vector describing the matching of actual atomic levels
! to that in the model atom with super levels.
!
	SUBROUTINE FDG_F_TO_S_NS_V1(NF,NS,NV,FL_OPTION,SL_OPTION,IL_OPTION,LUIN,FILENAME)
	IMPLICIT NONE
!
! Altered 5-Feb-2009: Inserted FL_OPTION into call. (Still V1).
!
	INTEGER NF
	INTEGER NS
	INTEGER NV
	INTEGER LUIN
	CHARACTER(LEN=*) FILENAME
	CHARACTER(LEN=*) FL_OPTION
	CHARACTER(LEN=*) SL_OPTION
	CHARACTER(LEN=*) IL_OPTION
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER, ALLOCATABLE :: F_TO_S(:)
	INTEGER, ALLOCATABLE :: INT_SEQ(:)
	INTEGER, ALLOCATABLE :: CNT(:)
!
	INTEGER LUER,I,J,K,IOS,ENTRY_NUM
	INTEGER NS_OLD
	CHARACTER*132 STRING
!
	IF(NF .EQ. NS)RETURN
	IF(FL_OPTION .EQ . ' ' .AND. IL_OPTION .EQ. ' ' .AND. SL_OPTION .EQ. ' ')RETURN
!
	LUER=ERROR_LU()
	OPEN(UNIT=LUIN,FILE=FILENAME,STATUS='OLD',IOSTAT=IOS,
	1       ACTION='READ')
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening ',FILENAME,' in FDG_F_TO_S_NS_V1'
	    STOP
	  END IF
!
	  STRING=' '
	  DO WHILE( INDEX(STRING,'Number of energy levels') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  BACKSPACE(LUIN)
	  READ(LUIN,*)I
	  IF(I .LT. NF)THEN
	    WRITE(LUER,*)'Error in FDG_F_TO_S_NS_V1'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Insufficient levels in file'
	    STOP
	  END IF
!
	  IF(FL_OPTION .EQ. ' ')THEN
	  ELSE IF(FL_OPTION .EQ. 'DO_ALL_LEVS')THEN
	    NF=I
	  ELSE IF(FL_OPTION(1:7) .EQ. 'SET_TO_')THEN
	    READ(FL_OPTION(8:),*)J
	    J=MAX(NF,J)			!Keep using existing level being used.
	   NF=MIN(I,J)
	  ELSE
	    WRITE(6,*)'Error in FDG_F_O_S_NS_V1: FL_OPTION not recognized'
	    WRITE(6,*)'FL_OPTION=',TRIM(FL_OPTION)
	    STOP
	  END IF
	  ALLOCATE (F_TO_S(NF))
	  ALLOCATE (INT_SEQ(NF))
	  ALLOCATE (CNT(NF))
!
	  STRING=' '
	  DO WHILE( INDEX(STRING,'Entry number of link') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  BACKSPACE(LUIN)
	  READ(LUIN,*)ENTRY_NUM
	  IF(ENTRY_NUM .LT. 2 .OR. ENTRY_NUM .GT. 20)THEN
	    WRITE(LUER,*)'Error in FDG_F_TO_S_NS_V1'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Bad entry number for level link'
	    STOP
	  END IF
!
! NB: All entries must be separated by at LEAST 2 spaces.
!
	  READ(LUIN,'(A)')STRING		!Blankline
	  DO I=1,NF
	    F_TO_S(I)=0.
	    READ(LUIN,'(A)')STRING
	    STRING=ADJUSTL(STRING)
!
! NB: ENTRY_NUM-2 as
!              -1 due to 1st entry being level ID
!              -1 as only need to point to entry.
!
	    J=INDEX(STRING,'  ')
	    DO K=1,ENTRY_NUM-2
	      STRING(1:)=STRING(J:)
	      STRING=ADJUSTL(STRING)
	      J=INDEX(STRING,'  ')
	    END DO
	    READ(STRING(J:),*)F_TO_S(I),INT_SEQ(I)
	  END DO
	CLOSE(LUIN)
!
	NS_OLD=NS
	IF(FL_OPTION .NE. ' ')THEN
	  NS=MAXVAL(F_TO_S)
	END IF
!
! Now check that:
! (1) All super-levels have at least 1 level from full atom
! (2) Correct number of super levels.
!
	DO I=1,NF
	  IF(F_TO_S(I) .LE. 0 .OR. F_TO_S(I) .GT. NS)THEN
	    WRITE(LUER,*)'Error in FDG_F_TO_S_NS_V1'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Incorrect super level ID for level ',I
	    STOP
	  END IF
	END DO
!
! Check that all super levels have at least 1 level.
!
	DO I=1,NS
	  CNT(I)=0
	END DO
	DO J=1,NF
	  I=F_TO_S(J)
	  CNT(I)=CNT(I)+1
	END DO
	DO I=1,NS
	  IF(CNT(I) .EQ. 0)THEN
	    WRITE(LUER,*)'Error in FDG_F_TO_S_NS_V1'
	    WRITE(LUER,*)'Currently treating ',FILENAME
	    WRITE(LUER,*)'Super level ',I,' has no components'
	    STOP
	  END IF
	END DO
!
	NS_OLD=NS
	CALL DO_SL_ADJUSTEMENT(F_TO_S,INT_SEQ,NF,NS,SL_OPTION,FILENAME)
	IF(IL_OPTION .EQ. 'USE_ALL_SL')THEN
	   NV=NS
	ELSE IF(IL_OPTION(1:7) .EQ. 'SET_TO_')THEN
	  READ(IL_OPTION(8:),*)NV
	  NV=MIN(NV,NS)
	ELSE IF(IL_OPTION .EQ. 'UPDATE')THEN
	  NV=MIN(NS, NV+(NS-NS_OLD))
	ELSE IF(IL_OPTION .EQ. 'NOCHANGE')THEN
	ELSE
	    WRITE(6,*)'Error in FDG_F_O_S_NS_V1: IL_OPTION not recognized'
	    WRITE(6,*)'IL_OPTION=',TRIM(IL_OPTION)
	    STOP
	END IF
!
	DEALLOCATE (F_TO_S,INT_SEQ,CNT)
!
	RETURN
	END
!
	SUBROUTINE DO_SL_ADJUSTEMENT(F_TO_S,INT_SEQ,NF,NS,SL_OPTION,FILENAME)
	IMPLICIT NONE
!
! Altered  1-FEb-2010 : Added FILENAME to call to improve diagnostic output.
!
	INTEGER NF
	INTEGER NS
	INTEGER F_TO_S(NF)
	INTEGER INT_SEQ(NF)
	CHARACTER(LEN=*) SL_OPTION
	CHARACTER(LEN=*) FILENAME
!
	INTEGER OLD_F_TO_S(NF)
	INTEGER TMP_F_TO_S
!
	INTEGER I, J
	INTEGER ID
	INTEGER COUNT
	INTEGER NS_LOW
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
!***************************************************************************
!***************************************************************************
!
! We now make the correction.
!
	IF(SL_OPTION(1:10) .EQ. ' ')THEN
!
! Nothing to do.
!
	ELSE IF(SL_OPTION(1:10) .EQ. 'SPLIT_LOW_')THEN
	  READ(SL_OPTION(11:),*)NS_LOW
	  COUNT=NF
	  OLD_F_TO_S(1:NF)=F_TO_S(1:NF)
	  DO I=1,NF
	    IF(F_TO_S(I) .LE. NS_LOW)THEN
	      COUNT=COUNT+1
	      F_TO_S(I)=F_TO_S(I)+COUNT
	    END IF
	  END DO
!
	  ID=0
	  F_TO_S(1:NF)=-F_TO_S(1:NF)
	  DO I=1,NF
	    IF(F_TO_S(I) .LT. 0)THEN
	      ID=ID+1
	      TMP_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	      DO J=I+1,NF
	        IF(F_TO_S(J) .EQ. TMP_F_TO_S)F_TO_S(J)=ID
	      END DO
	    ELSE IF(F_TO_S(I) .EQ. 0)THEN
	      ID=ID+1
	      TMP_F_TO_S=F_TO_S(I)
	      F_TO_S(I)=ID
	    END IF
	  END DO
!
! Fix sequence number.
!
	  DO I=1,NF
	    IF(INT_SEQ(I) .NE. 0)THEN
	      DO J=1,NF
	        IF(INT_SEQ(I) .EQ. OLD_F_TO_S(J))THEN
	          INT_SEQ(I)=F_TO_S(J)
	          EXIT
	        END IF
	      END DO
	    END IF
	  END DO
!
	  J=0
	  DO I=1,NF
	    J=MAX(J,F_TO_S(I))
	  END DO
	  NS=J
!
	ELSE IF(SL_OPTION(1:5) .EQ. 'HALF_')THEN
	  READ(SL_OPTION(6:),*)NS_LOW
	  COUNT=NF
	  OLD_F_TO_S(1:NF)=F_TO_S(1:NF)
	  DO I=1,NF
	    IF(F_TO_S(I) .LE. NS_LOW)THEN
	    ELSE
	      F_TO_S(I)=NS_LOW+(F_TO_S(I)+1-NS_LOW)/2
	    END IF
	  END DO
!
! Fix sequence number.
!
	  INT_SEQ=0
!
! Check consistency
!
	  NS=MAXVAL(F_TO_S)
	  J=MINVAL(F_TO_S)
	  IF(J .EQ. 0 .OR. NS .GT. NF)THEN
	    WRITE(6,*)'Error in DO_SL_ADJUSTEMENT - invalid SL assignment for HALF option'
	    WRITE(6,*)'NS=',NS,'NF=',NF
	    IF(J .EQ. 0)WRITE(6,*)'Level not assigned:'
	    STOP
	  END IF
	  OLD_F_TO_S=0
	  DO J=1,NF
	   OLD_F_TO_S(F_TO_S(J))=1
	  END DO
	  DO J=1,NS
	    IF(OLD_F_TO_S(J) .EQ. 0)THEN
	      WRITE(6,*)'Error in DO_SL_ADJUSTEMENT - invalid SL assignment for HALF option'
	      WRITE(6,*)'There is no SL assignment for SL',J
	      STOP
	    END IF
	  END DO
!
	ELSE
	  WRITE(6,*)'Error in DO_SL_ADJUSTEMENT: SL_OPTION not recognized'
	  WRITE(6,*)'SL_OPTION=',TRIM(SL_OPTION)
	  STOP
	END IF
!
	IF(SL_OPTION(1:10) .NE. ' ' .AND. FILENAME .NE. ' ')THEN
	  IF(FIRST_TIME)THEN
	    OPEN(UNIT=101,STATUS='UNKNOWN',FILE='MODIFIED_SL_ASSIGNMENTS')
	  END IF
	  WRITE(101,'(A,A)')' Species is ',FILENAME(1:INDEX(FILENAME,'_')-1)
	  WRITE(101,'(A,I4)')'     Old NS',MAXVAL(OLD_F_TO_S)
	  WRITE(101,'(A,I4)')' Revised NS',NS
	  WRITE(101,'(4(2X,A))')'  Level','NEW_F2S','OLD_F2S','INT_SEQ'
	  DO I=1,NF
	    WRITE(101,'(4I9)')I,F_TO_S(I),OLD_F_TO_S(I),INT_SEQ(I)
	  END DO
	  FLUSH(UNIT=101)
	END IF
!
	RETURN
	END
