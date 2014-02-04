	MODULE NEW_GEN_IN_INTERFACE
	INTERFACE NEW_GEN_IN
	  MODULE PROCEDURE NEW_GEN_IN_LOG,
	1                  NEW_GEN_IN_SP,
	1                  NEW_GEN_IN_DBLE,
	1                  NEW_GEN_IN_INT,
	1                  NEW_GEN_IN_STR,
	1                  NEW_GEN_IN_MULGT_INT,
	1                  NEW_GEN_IN_MULT_SP,
	1                  NEW_GEN_IN_MULT_DP
	END INTERFACE
C
	INTEGER, PARAMETER :: GT_IN=5
	INTEGER, PARAMETER :: GT_OUT=6
	INTEGER  LF_IN
	INTEGER  LF_OUT
	LOGICAL TERM_INPUT
	LOGICAL LOG_FILE
	CHARACTER*3 ADV_OPT
	CHARACTER*80 GSTRING
	CHARACTER*80 KEY_GSTRING
	DATA ADV_OPT/'NO'/
	DATA TERM_INPUT/.TRUE./
	DATA LOG_FILE/.FALSE./
!
	CONTAINS
!
	SUBROUTINE NEW_GEN_IN_OPTS(OPTION,FILENAME,LU)
	IMPLICIT NONE
	INTEGER LU
	INTEGER IOS
	CHARACTER*(*) OPTION
	CHARACTER*(*) FILENAME
!
	IF(OPTION .EQ. 'OPEN_LOG_FILE')THEN
	  OPEN(UNIT=LU,FILE=FILENAME,STATUS='UNKNOWN',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(GT_OUT,*)'Unable to open log-file, IOSTAT=',IOS
	    RETURN
	  END IF
	  LF_OUT=LU
	  LOG_FILE=.TRUE.
	ELSE IF(OPTION .EQ. 'CLOSE_LOG_FILE')THEN
	  CLOSE(LF_OUT)
	  LOG_FILE=.FALSE.
	ELSE IF(OPTION .EQ. 'OPEN_FILE_INPUT')THEN
	  OPEN(UNIT=LU,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(GT_OUT,*)'Unable to open input-file, IOSTAT=',IOS
	    RETURN
	  END IF
	  LF_IN=LU
	  TERM_INPUT=.FALSE.
	  ADV_OPT='YES'
	ELSE IF(OPTION .EQ. 'CLOSE_FILE_INPUT')THEN
	  CLOSE(LF_IN)
	  TERM_INPUT=.TRUE.
	  ADV_OPT='NO'
	ELSE
	  WRITE(GT_OUT,*)'Unrecognized option in NEW_GEN_IN_OPTS'
	  WRITE(GT_OUT,*)'No action taken'
	  RETURN
	END IF
!
	RETURN
	END SUBROUTINE NEW_GEN_IN_OPTS
C
C General input routines ---generally read a single value from user only.
C Default value is output to user.
C
	SUBROUTINE NEW_GEN_IN_LOG(VAL,DESC)
	IMPLICIT NONE
C
	LOGICAL VAL
	INTEGER L,IOS
	CHARACTER DESC*(*)
C
900	GSTRING=DESC
	L=LEN_TRIM(GSTRING)
	IF(GSTRING(L:L) .EQ. '=')L=L-1
	WRITE(GT_OUT,100,ADVANCE=ADV_OPT)GSTRING(1:L),VAL
100	FORMAT(1X,A,' [',L1,']: ')
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE 
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
!
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(GSTRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(GT_OUT,*)'Error reading logical value - try again'
	GOTO 900
	END SUBROUTINE NEW_GEN_IN_LOG
C
C 
C
	SUBROUTINE NEW_GEN_IN_STR(STR,DESC)
C
C Altered 13-Oct-1997 : Ability to input single blank character.
C
	INTEGER L,J
	CHARACTER DESC*(*),STR*(*),LOC_STR*80
C
900	GSTRING=DESC
	L=LEN_TRIM(GSTRING)
	IF(GSTRING(L:L) .EQ. '=')L=L-1
	LOC_STR=STR
	J=LEN(TRIM(LOC_STR))
C
	WRITE(GT_OUT,100,ADVANCE=ADV_OPT)GSTRING(1:L),LOC_STR(1:J)
100	FORMAT(1X,A,' [',A,']: ')
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
C
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)RETURN		!Take default value
C
C The following statement allows a single blank character to be returned.
C It assumed that the progammer will never have need to enter these
C string particularly (An altenative and better way would be to use
C the non-standard Q descriptor)
C
	IF(GSTRING .EQ. '""' .OR. GSTRING .EQ. '" "')THEN
	  STR=' '
	  RETURN
	END IF
C
	J=LEN(STR)
	IF(J .LT. L)THEN
	  WRITE(GT_OUT,*)'Error --- string too small'
	  GOTO 900
	END IF
	STR=GSTRING(1:L)
	RETURN
C
	END SUBROUTINE NEW_GEN_IN_STR
C
C 
	SUBROUTINE NEW_GEN_IN_INT(VAL,STR)
	INTEGER, PARAMETER :: GT_IN=5
	INTEGER, PARAMETER :: GT_OUT=6
C
	INTEGER L,IOS,VAL
	CHARACTER STR*(*)
	CHARACTER*20 FORM
C
900	GSTRING=STR
	L=LEN_TRIM(GSTRING)
	IF(GSTRING(L:L) .EQ. '=')L=L-1
C
	WRITE(FORM,'(I15)')ABS(VAL)
	LEN_FORM=15
	DO WHILE(FORM(1:1) .EQ. ' ')
	  FORM(1:)=FORM(2:)
	  LEN_FORM=LEN_FORM-1
	END DO
	IF(VAL .LT. 0)THEN
	  FORM='-'//FORM(1:LEN_FORM)
	  LEN_FORM=LEN_FORM+1
	END IF
C
	WRITE(GT_OUT,'(1X,A,A,A,A)',ADVANCE=ADV_OPT)
	1                    GSTRING(1:L),' [',FORM(1:LEN_FORM),']: '
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
! 
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(GSTRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(GT_OUT,*)'Error reading integer value - try again'
	GOTO 900
	END SUBROUTINE NEW_GEN_IN_INT
C
C 
C
	SUBROUTINE NEW_GEN_IN_DBLE(VAL,STR)
	INTEGER, PARAMETER :: GT_IN=5
	INTEGER, PARAMETER :: GT_OUT=6
C
C Altered 10-Dec-1990 --- FORM_SP_NUM routine installed.
C
	REAL*8 VAL
	INTEGER L,IOS,LEN_FORM
	CHARACTER STR*(*),FORM*20
C
900	GSTRING=STR
	L=LEN_TRIM(GSTRING)
	IF(GSTRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_DP_NUM(VAL,FORM,LEN_FORM)
	WRITE(GT_OUT,100,ADVANCE=ADV_OPT)GSTRING(1:L),FORM(1:LEN_FORM)
100	FORMAT(1X,A,' [',A,']: ')
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
C
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(GSTRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(GT_OUT,*)'Error reading logical value - try again'
	GOTO 900
	END SUBROUTINE NEW_GEN_IN_DBLE
C 
C
	SUBROUTINE NEW_GEN_IN_SP(VAL,STR)
	INTEGER, PARAMETER :: GT_IN=5
	INTEGER, PARAMETER :: GT_OUT=6
C
C Altered 10-Dec-1990 --- FORM_SP_NUM routine installed.
C
	REAL*4 VAL
	INTEGER L,IOS,LEN_FORM
	CHARACTER STR*(*),FORM*20
C
900	GSTRING=STR
	L=LEN_TRIM(GSTRING)
	IF(GSTRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_SP_NUM(VAL,FORM,LEN_FORM)
	WRITE(GT_OUT,100,ADVANCE=ADV_OPT)GSTRING(1:L),FORM(1:LEN_FORM)
100	FORMAT(1X,A,' [',A,']: ')
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
C
  	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(GSTRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(GT_OUT,*)'Error reading logical value - try again'
	GOTO 900
	END SUBROUTINE NEW_GEN_IN_SP
C 
C
	SUBROUTINE NEW_GEN_IN_MULT_DP(VAL,N,NMAX,STR)
C
C Altered 7-Jul-1997 : NMAX installed. Number of values input returned in N.
C                      Previously all NMAX values had to be input.
C
	INTEGER, PARAMETER :: GT_IN=5
	INTEGER, PARAMETER :: GT_OUT=6
C
	INTEGER N,NMAX
	REAL*8 VAL(NMAX)
	CHARACTER STR*(*)
C
	INTEGER K,L,IOS,LEN_FORM,LEN_TOT
	CHARACTER TOT_FORM*80,FORM*20
C
900	GSTRING=STR
	L=LEN_TRIM(GSTRING)
	IF(GSTRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_DP_NUM(VAL(1),TOT_FORM,LEN_TOT)
	DO K=2,NMAX
	  CALL FORM_DP_NUM(VAL(K),FORM,LEN_FORM)
	  TOT_FORM=TOT_FORM(1:LEN_TOT)//','//FORM(1:LEN_FORM)
          LEN_TOT=LEN_TOT+1+LEN_FORM
	END DO
	WRITE(GT_OUT,100,ADVANCE=ADV_OPT)GSTRING(1:L),TOT_FORM(1:LEN_TOT)
100	FORMAT(1X,A,' [',A,']: ')
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
C
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)THEN
	  N=NMAX
	  RETURN		!Take default values
	END IF
C
	DO I=1,NMAX
C
C Strip leading blanks and commas.
C
	  DO WHILE(GSTRING(1:1) .EQ. ' ' .OR. GSTRING(1:1) .EQ. ',')
	    GSTRING(1:)=GSTRING(2:)
	    L=L-1
	    IF(L .EQ. 0)THEN
	      WRITE(GT_OUT,*)'Error reading multiple DP values - try again'
	      WRITE(GT_OUT,*)'Currently reading variable',I
	      GOTO 900
	    END IF
	  END DO
C
C Read next value
C
	  READ(GSTRING(1:L),*,IOSTAT=IOS)VAL(I)
	  IF(IOS .NE. 0)THEN
	     WRITE(GT_OUT,*)'Error reading multiple DP values - try again'
	     WRITE(GT_OUT,*)'Currently reading variable',I
	     GOTO 900
	  END IF
C
C Strip the value just read from the string.
C
	  DO WHILE(GSTRING(1:1) .NE. ' ' .AND. GSTRING(1:1) .NE. ',')
	    GSTRING(1:)=GSTRING(2:)
	    L=L-1
	  END DO
	  IF(L .LE. 0)THEN
	    N=I
	    RETURN
	  END IF
	END DO
	END SUBROUTINE NEW_GEN_IN_MULT_DP
C 
C
	SUBROUTINE NEW_GEN_IN_MULT_SP(VAL,N,NMAX,STR)
C
C Altered 7-Jul-1997 : NMAX installed. Number of values input returned in N.
C                      Previously all NMAX values had to be input.
C
	INTEGER, PARAMETER :: GT_IN=5
	INTEGER, PARAMETER :: GT_OUT=6
C
	INTEGER N,NMAX
	REAL*4 VAL(NMAX)
	CHARACTER STR*(*)
C
	INTEGER K,L,IOS,LEN_FORM,LEN_TOT
	CHARACTER TOT_FORM*80,FORM*20
C
900	GSTRING=STR
	L=LEN_TRIM(GSTRING)
	IF(GSTRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_SP_NUM(VAL(1),TOT_FORM,LEN_TOT)
	DO K=2,NMAX
	  CALL FORM_SP_NUM(VAL(K),FORM,LEN_FORM)
	  TOT_FORM=TOT_FORM(1:LEN_TOT)//','//FORM(1:LEN_FORM)
          LEN_TOT=LEN_TOT+1+LEN_FORM
	END DO
	WRITE(GT_OUT,100,ADVANCE=ADV_OPT)GSTRING(1:L),TOT_FORM(1:LEN_TOT)
100	FORMAT(1X,A,' [',A,']: ')
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
C
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)THEN
	  N=NMAX
	  RETURN		!Take default values
	END IF
C
	DO I=1,NMAX
C
C Strip leading blanks and commas.
C
	  DO WHILE(GSTRING(1:1) .EQ. ' ' .OR. GSTRING(1:1) .EQ. ',')
	    GSTRING(1:)=GSTRING(2:)
	    L=L-1
	    IF(L .EQ. 0)THEN
	      WRITE(GT_OUT,*)'Error reading multiple SP values - try again'
	      WRITE(GT_OUT,*)'Currently reading variable',I
	      GOTO 900
	    END IF
	  END DO
C
C Read next value
C
	  READ(GSTRING(1:L),*,IOSTAT=IOS)VAL(I)
	  IF(IOS .NE. 0)THEN
	     WRITE(GT_OUT,*)'Error reading multiple SP values - try again'
	     WRITE(GT_OUT,*)'Currently reading variable',I
	     GOTO 900
	  END IF
C
C Strip the value just read from the string.
C
	  DO WHILE(GSTRING(1:1) .NE. ' ' .AND. GSTRING(1:1) .NE. ',')
	    GSTRING(1:)=GSTRING(2:)
	    L=L-1
	  END DO
	  IF(L .LE. 0)THEN
	    N=I
	    RETURN
	  END IF
	END DO
	END SUBROUTINE NEW_GEN_IN_MULT_SP
C
C 
C
	SUBROUTINE NEW_GEN_IN_MULGT_INT(VAL,N,NMAX,STR)
C
C Altered 7-Jul-1997 : NMAX installed. Number of values input returned in N.
C                      Previously all NMAX values had to be input.
C
	INTEGER, PARAMETER :: GT_IN=5
	INTEGER, PARAMETER :: GT_OUT=6
C
	INTEGER N,NMAX,VAL(NMAX)
	CHARACTER STR*(*)
C
	INTEGER L,IOS,LEN_FORM,LEN_TOT
	CHARACTER FORM*20,TOT_FORM*80
C
900	CONTINUE
	GSTRING=STR
	L=LEN_TRIM(GSTRING)
C
C Format the string showing the default values.
C
	DO K=1,NMAX
	  WRITE(FORM,'(I15)')ABS(VAL(K))
	  LEN_FORM=15
	  DO WHILE(FORM(1:1) .EQ. '' )
	    FORM(1:)=FORM(2:)
	    LEN_FORM=LEN_FORM-1
	  END DO
	  IF(VAL(K) .LT. 0)THEN
	    FORM='-'//FORM(1:LEN_FORM)
	    LEN_FORM=LEN_FORM+1
	  END IF
	  IF(K .EQ. 1)THEN
	    TOT_FORM=FORM(1:LEN_FORM)
	    LEN_TOT=LEN_FORM
	  ELSE
	    TOT_FORM=TOT_FORM(1:LEN_TOT)//','//FORM(1:LEN_FORM)
	    LEN_TOT=LEN_TOT+LEN_FORM+1
	  END IF
	END DO
C
	WRITE(GT_OUT,'(1X,A,A,A,A)',ADVANCE=ADV_OPT)
	1             GSTRING(1:L),' [',TOT_FORM(1:LEN_TOT),']: '
100	FORMAT(1X,A,' [',A,']: ')
!
	IF(.NOT. TERM_INPUT)THEN
	  READ(LF_IN,'(A)',IOSTAT=IOS)KEY_GSTRING
	  IF(IOS .EQ. 0 .AND. TRIM(KEY_GSTRING(3:)) .EQ. TRIM(GSTRING))THEN
	    READ(LF_IN,'(A)')GSTRING
	  ELSE
	    IF(IOS .NE. 0)THEN
	      WRITE(GT_OUT,*)'Error reading options from log file'
	    ELSE
	      WRITE(GT_OUT,*)'Inconsistency reading options from log file'
	    END IF
	    WRITE(GT_OUT,*)'Returning to keyboard input'
	    CLOSE(LF_IN)
	    TERM_INPUT=.TRUE.
	    ADV_OPT='NO'
	  END IF
	END IF
!
	IF(TERM_INPUT)THEN
	  IF(LOG_FILE)WRITE(LF_OUT,'(A,A)')'>>',TRIM(GSTRING)
	  READ(GT_IN,'(A)')GSTRING
	  IF(LOG_FILE)WRITE(LF_OUT,'(A)')TRIM(GSTRING)
	END IF
C
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)THEN
	  N=NMAX
	  RETURN		!Take default values
	END IF
C
C Now read in up to NMAX values. The actual number of values returned is
C returned in N.
C
	DO I=1,NMAX
C
C Strip leading blanks and commas.
C
	  DO WHILE(GSTRING(1:1) .EQ. ' ' .OR. GSTRING(1:1) .EQ. ',')
	    GSTRING(1:)=GSTRING(2:)
	    L=L-1
	    IF(L .EQ. 0)THEN
	      WRITE(GT_OUT,*)'Error reading multiple INT values - try again'
	      WRITE(GT_OUT,*)'Currently reading variable',I
	      GOTO 900
	    END IF
	  END DO
C
C Read next value
C
	  READ(GSTRING(1:L),*,IOSTAT=IOS)VAL(I)
	  IF(IOS .NE. 0)THEN
	     WRITE(GT_OUT,*)'Error reading multiple INT values - try again'
	     WRITE(GT_OUT,*)'Currently reading variable',I
	     GOTO 900
	  END IF
C
C Strip the value just read from the string.
C
	  DO WHILE(GSTRING(1:1) .NE. ' ' .AND. GSTRING(1:1) .NE. ',')
	    GSTRING(1:)=GSTRING(2:)
	    L=L-1
	  END DO
	  IF(L .LE. 0)THEN
	    N=I
	    RETURN
	  END IF
	END DO
	READ(GT_IN,'(A)')GSTRING
	L=LEN_TRIM(GSTRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(GSTRING(1:L),*,IOSTAT=IOS)(VAL(I),I=1,N)
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(GT_OUT,*)'Error reading integer value - try again'
	GOTO 900
	END SUBROUTINE NEW_GEN_IN_MULGT_INT
C
	END MODULE NEW_GEN_IN_INTERFACE
