	MODULE GEN_IN_INTERFACE
	INTERFACE GEN_IN
	  MODULE PROCEDURE GEN_IN_LOG,
	1                  GEN_IN_SP,
	1                  GEN_IN_DBLE,
	1                  GEN_IN_INT,
	1                  GEN_IN_STR,
	1                  GEN_IN_MULT_INT,
	1                  GEN_IN_MULT_SP,
	1                  GEN_IN_MULT_DP
	END INTERFACE
C
	CONTAINS
C
C General input routines ---generally read a single value from user only.
C Default value is output to user.
C
	SUBROUTINE GEN_IN_LOG(VAL,DESC)
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
C
	LOGICAL VAL
	INTEGER L,IOS
	CHARACTER DESC*(*),STRING*80
C
900	STRING=DESC
	L=LEN_TRIM(STRING)
	IF(STRING(L:L) .EQ. '=')L=L-1
	WRITE(T_OUT,100,ADVANCE='NO')STRING(1:L),VAL
100	FORMAT(1X,A,' [',L1,']: ')
	READ(T_IN,'(A)')STRING
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(STRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(T_OUT,*)'Error reading logical value - try again'
	GOTO 900
	END SUBROUTINE GEN_IN_LOG
C
C 
C
	SUBROUTINE GEN_IN_STR(STR,DESC)
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
C
C Altered 13-Oct-1997 : Ability to input single blank character.
C
	INTEGER L,J
	CHARACTER DESC*(*),STR*(*),STRING*80,LOC_STR*80
C
900	STRING=DESC
	L=LEN_TRIM(STRING)
	IF(STRING(L:L) .EQ. '=')L=L-1
	LOC_STR=STR
	J=LEN_TRIM(LOC_STR)
C
	WRITE(T_OUT,100,ADVANCE='NO')STRING(1:L),LOC_STR(1:J)
100	FORMAT(1X,A,' [',A,']: ')
	READ(T_IN,'(A)')STRING
C
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)RETURN		!Take default value
C
C The following statement allows a single blank character to be returned.
C It assumed that the progammer will never have need to enter these
C string particularly (An altenative and better way would be to use
C the non-standard Q descriptor)
C
	IF(STRING .EQ. '""' .OR. STRING .EQ. '" "')THEN
	  STR=' '
	  RETURN
	END IF
C
	J=LEN(STR)
	IF(J .LT. L)THEN
	  WRITE(T_OUT,*)'Error --- string too small'
	  GOTO 900
	END IF
	STR=STRING(1:L)
	RETURN
C
	END SUBROUTINE GEN_IN_STR
C
C 
	SUBROUTINE GEN_IN_INT(VAL,STR,LOW_LIM,UP_LIM)
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, OPTIONAL :: LOW_LIM
	INTEGER, OPTIONAL :: UP_LIM
!
	INTEGER L,IOS,VAL
	CHARACTER STR*(*),STRING*80
	CHARACTER*20 FORM
C
900	STRING=STR
	L=LEN_TRIM(STRING)
	IF(STRING(L:L) .EQ. '=')L=L-1
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
	WRITE(T_OUT,'(1X,A,A,A,A)',ADVANCE='NO')
	1                    STRING(1:L),' [',FORM(1:LEN_FORM),']: '
	READ(T_IN,'(A)')STRING
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(STRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)THEN
	  IF(PRESENT(LOW_LIM))THEN
	    IF(VAL .LT. LOW_LIM)THEN
	      WRITE(T_OUT,*)'Invalid value: lower limit is ',LOW_LIM
	      GOTO 900
	    END IF
	  END IF
	  IF(PRESENT(UP_LIM))THEN
	    IF(VAL .GT. UP_LIM)THEN
	      WRITE(T_OUT,*)'Invalid value: upper limit is ',UP_LIM
	      GOTO 900
	    END IF
	  END IF
	  RETURN		!Succesful read
	END IF
	WRITE(T_OUT,*)'Error reading integer value - try again'
	GOTO 900
	END SUBROUTINE GEN_IN_INT
C
C 
C
	SUBROUTINE GEN_IN_DBLE(VAL,STR)
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
C
C Altered 10-Dec-1990 --- FORM_SP_NUM routine installed.
C
	REAL*8 VAL
	INTEGER L,IOS,LEN_FORM
	CHARACTER STR*(*),STRING*80,FORM*20
C
900	STRING=STR
	L=LEN_TRIM(STRING)
	IF(STRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_DP_NUM(VAL,FORM,LEN_FORM)
	WRITE(T_OUT,100,ADVANCE='NO')STRING(1:L),FORM(1:LEN_FORM)
100	FORMAT(1X,A,' [',A,']: ')
	READ(T_IN,'(A)')STRING
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(STRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(T_OUT,*)'Error reading logical value - try again'
	GOTO 900
	END SUBROUTINE GEN_IN_DBLE
C 
C
	SUBROUTINE GEN_IN_SP(VAL,STR)
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
C
C Altered 10-Dec-1990 --- FORM_SP_NUM routine installed.
C
	REAL*4 VAL
	INTEGER L,IOS,LEN_FORM
	CHARACTER STR*(*),STRING*80,FORM*20
C
900	STRING=STR
	L=LEN_TRIM(STRING)
	IF(STRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_SP_NUM(VAL,FORM,LEN_FORM)
	WRITE(T_OUT,100,ADVANCE='NO')STRING(1:L),FORM(1:LEN_FORM)
100	FORMAT(1X,A,' [',A,']: ')
	READ(T_IN,'(A)')STRING
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(STRING(1:L),*,IOSTAT=IOS)VAL
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(T_OUT,*)'Error reading logical value - try again'
	GOTO 900
	END SUBROUTINE GEN_IN_SP
C 
C
	SUBROUTINE GEN_IN_MULT_DP(VAL,N,NMAX,STR)
C
C Altered 7-Jul-1997 : NMAX installed. Number of values input returned in N.
C                      Previously all NMAX values had to be input.
C
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
C
	INTEGER N,NMAX
	REAL*8 VAL(NMAX)
	CHARACTER STR*(*)
C
	INTEGER K,L,IOS,LEN_FORM,LEN_TOT
	CHARACTER STRING*80,TOT_FORM*80,FORM*20
C
900	STRING=STR
	L=LEN_TRIM(STRING)
	IF(STRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_DP_NUM(VAL(1),TOT_FORM,LEN_TOT)
	DO K=2,NMAX
	  CALL FORM_DP_NUM(VAL(K),FORM,LEN_FORM)
	  TOT_FORM=TOT_FORM(1:LEN_TOT)//','//FORM(1:LEN_FORM)
          LEN_TOT=LEN_TOT+1+LEN_FORM
	END DO
	WRITE(T_OUT,100,ADVANCE='NO')STRING(1:L),TOT_FORM(1:LEN_TOT)
100	FORMAT(1X,A,' [',A,']: ')
	READ(T_IN,'(A)')STRING
C
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)THEN
	  N=NMAX
	  RETURN		!Take default values
	END IF
C
	DO I=1,NMAX
C
C Strip leading blanks and commas.
C
	  DO WHILE(STRING(1:1) .EQ. ' ' .OR. STRING(1:1) .EQ. ',')
	    STRING(1:)=STRING(2:)
	    L=L-1
	    IF(L .EQ. 0)THEN
	      WRITE(T_OUT,*)'Error reading multiple DP values - try again'
	      WRITE(T_OUT,*)'Currently reading variable',I
	      GOTO 900
	    END IF
	  END DO
C
C Read next value
C
	  READ(STRING(1:L),*,IOSTAT=IOS)VAL(I)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error reading multiple DP values - try again'
	     WRITE(T_OUT,*)'Currently reading variable',I
	     GOTO 900
	  END IF
C
C Strip the value just read from the string.
C
	  DO WHILE(STRING(1:1) .NE. ' ' .AND. STRING(1:1) .NE. ',')
	    STRING(1:)=STRING(2:)
	    L=L-1
	  END DO
	  IF(L .LE. 0)THEN
	    N=I
	    RETURN
	  END IF
	END DO
	END SUBROUTINE GEN_IN_MULT_DP
C 
C
	SUBROUTINE GEN_IN_MULT_SP(VAL,N,NMAX,STR)
C
C Altered 7-Jul-1997 : NMAX installed. Number of values input returned in N.
C                      Previously all NMAX values had to be input.
C
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
C
	INTEGER N,NMAX
	REAL*4 VAL(NMAX)
	CHARACTER STR*(*)
C
	INTEGER K,L,IOS,LEN_FORM,LEN_TOT
	CHARACTER STRING*80,TOT_FORM*80,FORM*20
C
900	STRING=STR
	L=LEN_TRIM(STRING)
	IF(STRING(L:L) .EQ. '=')L=L-1
C
	CALL FORM_SP_NUM(VAL(1),TOT_FORM,LEN_TOT)
	DO K=2,NMAX
	  CALL FORM_SP_NUM(VAL(K),FORM,LEN_FORM)
	  TOT_FORM=TOT_FORM(1:LEN_TOT)//','//FORM(1:LEN_FORM)
          LEN_TOT=LEN_TOT+1+LEN_FORM
	END DO
	WRITE(T_OUT,100,ADVANCE='NO')STRING(1:L),TOT_FORM(1:LEN_TOT)
100	FORMAT(1X,A,' [',A,']: ')
	READ(T_IN,'(A)')STRING
C
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)THEN
	  N=NMAX
	  RETURN		!Take default values
	END IF
C
	DO I=1,NMAX
C
C Strip leading blanks and commas.
C
	  DO WHILE(STRING(1:1) .EQ. ' ' .OR. STRING(1:1) .EQ. ',')
	    STRING(1:)=STRING(2:)
	    L=L-1
	    IF(L .EQ. 0)THEN
	      WRITE(T_OUT,*)'Error reading multiple SP values - try again'
	      WRITE(T_OUT,*)'Currently reading variable',I
	      GOTO 900
	    END IF
	  END DO
C
C Read next value
C
	  READ(STRING(1:L),*,IOSTAT=IOS)VAL(I)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error reading multiple SP values - try again'
	     WRITE(T_OUT,*)'Currently reading variable',I
	     GOTO 900
	  END IF
C
C Strip the value just read from the string.
C
	  DO WHILE(STRING(1:1) .NE. ' ' .AND. STRING(1:1) .NE. ',')
	    STRING(1:)=STRING(2:)
	    L=L-1
	  END DO
	  IF(L .LE. 0)THEN
	    N=I
	    RETURN
	  END IF
	END DO
	END SUBROUTINE GEN_IN_MULT_SP
C
C 
C
	SUBROUTINE GEN_IN_MULT_INT(VAL,N,NMAX,STR)
C
C Altered 7-Jul-1997 : NMAX installed. Number of values input returned in N.
C                      Previously all NMAX values had to be input.
C
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
C
	INTEGER N,NMAX,VAL(NMAX)
	CHARACTER STR*(*)
C
	INTEGER L,IOS,LEN_FORM,LEN_TOT
	CHARACTER STRING*80
	CHARACTER FORM*20,TOT_FORM*80
C
900	CONTINUE
	STRING=STR
	L=LEN_TRIM(STRING)
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
	WRITE(T_OUT,'(1X,A,A,A,A)',ADVANCE='NO')
	1             STRING(1:L),' [',TOT_FORM(1:LEN_TOT),']: '
100	FORMAT(1X,A,' [',A,']: ')
	READ(T_IN,'(A)')STRING
C
	L=LEN_TRIM(STRING)
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
	  DO WHILE(STRING(1:1) .EQ. ' ' .OR. STRING(1:1) .EQ. ',')
	    STRING(1:)=STRING(2:)
	    L=L-1
	    IF(L .EQ. 0)THEN
	      WRITE(T_OUT,*)'Error reading multiple INT values - try again'
	      WRITE(T_OUT,*)'Currently reading variable',I
	      GOTO 900
	    END IF
	  END DO
C
C Read next value
C
	  READ(STRING(1:L),*,IOSTAT=IOS)VAL(I)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error reading multiple INT values - try again'
	     WRITE(T_OUT,*)'Currently reading variable',I
	     GOTO 900
	  END IF
C
C Strip the value just read from the string.
C
	  DO WHILE(STRING(1:1) .NE. ' ' .AND. STRING(1:1) .NE. ',')
	    STRING(1:)=STRING(2:)
	    L=L-1
	  END DO
	  IF(L .LE. 0)THEN
	    N=I
	    RETURN
	  END IF
	END DO
	READ(T_IN,'(A)')STRING
	L=LEN_TRIM(STRING)
	IF(L .EQ. 0)RETURN		!Take default value
	READ(STRING(1:L),*,IOSTAT=IOS)(VAL(I),I=1,N)
	IF(IOS .EQ. 0)RETURN		!Succesful read
	WRITE(T_OUT,*)'Error reading integer value - try again'
	GOTO 900
	END SUBROUTINE GEN_IN_MULT_INT
C
	END MODULE GEN_IN_INTERFACE
