!
! Simple routine to read in data in column format.
! Based on RXY option if gramon_pgplot. To be used in plt_spec, so
!   simple XY files can be read in with out needing model or observational data file.
! Records beginning with ! and blank lines are ignored.
! The number of data points does not need to be specified.
!    If specified (!Number of data points: ?) it should include ! lines and blank lines
!       mixed within the data.
!
	SUBROUTINE RD_XY_DATA(XVEC,YVEC,NX,NX_MAX,FILNAME,LU,IOS)
	USE NEW_GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created: 12-June-2010
!
	INTEGER IOS
	INTEGER NX
	INTEGER NX_MAX
	INTEGER LU
	REAL*8 XVEC(NX_MAX)
	REAL*8 YVEC(NX_MAX)
	CHARACTER(LEN=*) FILNAME
!
	INTEGER, PARAMETER :: MAX_COL=20
	INTEGER, PARAMETER :: T_OUT=6
!
	REAL*8 TEMP_VAR(MAX_COL)
	INTEGER I,K,L,NX_RD
	INTEGER COLUMN(2)
	CHARACTER(LEN=200) STRING
!
! Procedure to read simple plots in column format from file.
! Blank lines and comments (begin with a !) are ignored.
!
	OPEN(UNIT=LU,FILE=FILNAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error opening file: IOS=',IOS
	  RETURN
	END IF
!
	NX_RD=NX_MAX
	STRING='!'
	DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	  READ(LU,'(A200)',IOSTAT=IOS)STRING
	  WRITE(T_OUT,*)TRIM(STRING)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error reading file'
	    CLOSE(UNIT=LU)
	    NX=0
	    RETURN
	  END IF
	  IF(INDEX(STRING,'Number of data points:') .NE. 0)THEN
	    I=INDEX(STRING,'Number of data points:')
	    READ(STRING(I+22:),*)NX_RD
	    STRING='!'
	  END IF
	END DO
	NX_RD=MIN(NX_RD,NX_MAX)
	BACKSPACE(UNIT=LU)
!
	COLUMN(1)=1; COLUMN(2)=2; K=2
	CALL NEW_GEN_IN(COLUMN,I,K,'Data columns')
	K=MAX(COLUMN(1),COLUMN(2))
	IF(K .GT. MAX_COL)THEN
	   WRITE(T_OUT,*)'Maximum number of columns in RD_XY_DATA is currently',MAX_COL
	   IOS=2
	   RETURN
	END IF
!
	IOS=0
	NX=0
	DO I=1,NX_RD
100	  READ(LU,'(A)',END=500,IOSTAT=IOS)STRING
	  IF(IOS .EQ. 0 .AND. STRING(1:1) .NE. '!' .AND. STRING .NE. ' ')THEN
	    READ(STRING,*,IOSTAT=IOS)(TEMP_VAR(L),L=1,K)
	    IF(IOS .EQ. 0)THEN
	      XVEC(NX+1)=TEMP_VAR(COLUMN(1))
	      YVEC(NX+1)=TEMP_VAR(COLUMN(2))
	      NX=NX+1
	    END IF
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error in reading record from file:',TRIM(FILNAME)
	    WRITE(T_OUT,*)'Data point number is:',I
	    WRITE(T_OUT,*)'Erronous record follows'
	    WRITE(T_OUT,'(A)')TRIM(STRING)
	    IF(NX .LT. 2)THEN
	      WRITE(T_OUT,*)'Skipping to next record'
	      GOTO 100
	    ELSE
	      WRITE(T_OUT,*)NX,' data points passed to plotting package'
	      IOS=0
	    END IF
	    CLOSE(UNIT=LU)
	    RETURN
	  END IF
	END DO
500	CONTINUE
	IOS=0
	IF(NX .EQ. NX_MAX)WRITE(T_OUT,*)'All data points may not have been read'
	WRITE(T_OUT,*)'Number of data points passed to plotting package is:',NX
!
	CLOSE(UNIT=LU)
	RETURN
	END
