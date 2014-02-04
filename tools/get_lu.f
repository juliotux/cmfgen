!
! Simple subroutine to return a free logical unit.
! Only unit above 200 are returned. No track of
! used units is kept.  
!
	SUBROUTINE GET_LU(LU,ID)
	IMPLICIT NONE
	INTEGER LU
!
	INTEGER I
	INTEGER LUER,ERROR_LU
	LOGICAL FILE_OPEN
	CHARACTER(LEN=*), OPTIONAL :: ID
!
	DO I=201,500
	  INQUIRE(UNIT=I,OPENED=FILE_OPEN)
	  IF(.NOT. FILE_OPEN)THEN
	    LU=I
	    RETURN
	  END IF
	END DO
!
	LUER=ERROR_LU()
	WRITE(LUER,*)'Error in GET_LU: no valid LUs available'
	IF(PRESENT(ID))THEN
	  WRITE(LUER,*)'Idenitier assoicated with call is ',TRIM(ID)
	END IF
	STOP
	END
