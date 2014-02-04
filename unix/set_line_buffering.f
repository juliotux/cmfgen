!
! This routine sets buffering to LINE BUFFERING for unit LU.
! It is designed for the PGI compiler. With the intel compiler this
! routine is not necessarily needed, and can be placed by a dummy
! routine of the same name.
!
	SUBROUTINE SET_LINE_BUFFERING(LU)
	IMPLICIT NONE
!
	INTEGER LU
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LEN_BUF=1000
	CHARACTER*(LEN_BUF) MY_BUFFER
	SAVE MY_BUFFER
!
	EXTERNAL SETVBUF
	INTEGER SETVBUF
	INTEGER IOS
!
! Comment out SETVBUF line for non PGI systems
!	
	IOS=0
	IOS=SETVBUF(LU,IONE,LEN_BUF,MY_BUFFER)
	IF(IOS .NE. 0)THEN
	  WRITE(LU,*)'Error setting bufer in SET_LINE_BUFFERING'
	  WRITE(LU,*)'IOS=',IOS
	  STOP
	END IF
!
	RETURN
	END
