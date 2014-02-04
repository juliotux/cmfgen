!
! Routine designed to read the header of the collisional data file and to output
! all occurnces of [Reference] to OUTGEN.
!
	SUBROUTINE GET_COL_REF(FILENAME,LUIN,LUOUT)
	IMPLICIT NONE
!
! Created 8-Sep-2012
!
	CHARACTER(LEN=*) FILENAME
	INTEGER LUIN
	INTEGER LUOUT
	INTEGER IOS
	CHARACTER(LEN=200) STRING
!
	OPEN(UNIT=LUIN,FILE=FILENAME,ACTION='READ',STATUS='OLD',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to open collisional data file: ',TRIM(FILENAME)
	  STOP
	END IF
!
	STRING=' '
	DO WHILE(INDEX(STRING,'Transition\T') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'[Reference]') .NE. 0)THEN
	    WRITE(LUOUT,'(A)')TRIM(STRING)	
	  END IF
	END DO
	CLOSE(LUIN)
!
	RETURN
	END
