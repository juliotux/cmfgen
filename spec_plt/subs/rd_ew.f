	SUBROUTINE RD_EW(LAM,EW,NL_MAX,NL,FILENAME,IOS)
	IMPLICIT NONE
C
C Altered 02-May-1995: IOS installed in call.
C Altered 02-Aug-1995: Bug fix --- OBS_FREQ needed to be zeroed for correct
C                         input operations.
C
	INTEGER NL,NL_MAX
	REAL*4 LAM(NL_MAX)
	REAL*4 EW(NL_MAX)
	CHARACTER*(*) FILENAME
C
	INTEGER IOS
	INTEGER, PARAMETER :: T_OUT=6
C
	CHARACTER*132 STRING
	INTEGER J
	REAL*4 T1
C
	OPEN(UNIT=10,FILE=FILENAME,ACTION='READ',STATUS='OLD',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to open file in RD_EW'
	  RETURN
	END IF
C
C The initialization required for when I find out how many numbers I have
C read in.
C
	DO J=1,NL_MAX
	  LAM(J)=0.0D0
	  EW(J)=0.0D0
	END DO
C
	NL=0
	DO WHILE(1 .EQ. 1)
	  READ(10,'(A)',END=100)STRING
          J=INDEX(STRING,'  ')
	  NL=NL+1
	  IF(NL .GT. NL_MAX)THEN
	    WRITE(T_OUT,*)'Error -- nL_MAX too small in RD_EW'
	    RETURN
	  END IF
	  READ(STRING(J+1:),*)LAM(NL),T1,EW(NL)
	END DO
C
100	CONTINUE
C                          
	RETURN
	END
