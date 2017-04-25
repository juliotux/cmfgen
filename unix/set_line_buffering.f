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
	CHARACTER*(LEN_BUF), SAVE :: MY_BUFFER(10)
	INTEGER, SAVE ::  LU_BUF(1:10)=0
	INTEGER, SAVE :: IBUF=0
!
	EXTERNAL SETVBUF
	INTEGER SETVBUF
	INTEGER IOS
	INTEGER I
!
! Comment out SETVBUF line for non PGI systems
! or put a RETURN here.
!
!
! Check if the logical unit already has an assciated buffer.
!
	IBUF=0
	DO I=1,10
	  IF(LU_BUF(I) .EQ. LU)THEN
	    IBUF=I
	    EXIT
	  END IF
	END DO
!
! If buffer not assigned, get buffer.
!
	IF(IBUF .EQ. 0)THEN
	  DO I=1,10
	    IF(LU_BUF(I) .EQ. 0)THEN
	      IBUF=I
	      LU_BUF(I)=LU
	      EXIT
	    END IF
	  END DO
	END IF
!
! Confirm buffer has been assigned.
!
	IF(IBUF .EQ. 0)THEN
	  WRITE(6,*)'Error in SET_LINE_BUFFEG - Insufficient line buffers'
	  WRITE(6,*)'List of logical units with bufers follows'
	  WRITE(6,'(5I5)')(LU_BUF(I),I=1,10)
	  STOP
	END IF
!
! Allocate the buffer.
! 
	IOS=0
	IOS=SETVBUF(LU,IONE,LEN_BUF,MY_BUFFER(IBUF))
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error setting bufer in SET_LINE_BUFFERING'
	  WRITE(6,*)'IOS=',IOS
	  STOP
	END IF
!
	RETURN
	END
