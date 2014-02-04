!
! Subroutine to write out information needed to reopen a DIRECT access file.
! The RECL is not stored with the file on UNIX systems.
!
	SUBROUTINE WRITE_DIRECT_INFO(ND,RECL,FILENAME,LU_EDD)
	IMPLICIT NONE
!
! Created 16-Jun-2000
!
	INTEGER ND			!Number of depth points written
	INTEGER RECL			!Record length (system dependent)
	INTEGER LU_EDD
	CHARACTER*(*) FILENAME
!
	CHARACTER*80 NEW_FILENAME
!                      
	NEW_FILENAME=TRIM(FILENAME)//'_INFO'
        OPEN(UNIT=LU_EDD,FILE=NEW_FILENAME,STATUS='REPLACE')
          WRITE(LU_EDD,'(2(3X,I6))')ND,RECL
          WRITE(LU_EDD,'(2(5X,A))')'  ND','RECL'
        CLOSE(LU_EDD)
!
	RETURN
	END
!
! Subroutine to read out information needed to reopen a DIRECT access file.
! The RECL is not stored with the file on UNIX systems.
!
	SUBROUTINE READ_DIRECT_INFO(ND,RECL,FILENAME,LU_EDD)
	IMPLICIT NONE
!
! Created 16-Jun-2000
!
	INTEGER ND			!Number of depth points written
	INTEGER RECL			!Record length (system dependent)
	INTEGER LU_EDD
	CHARACTER*(*) FILENAME
!
	CHARACTER*80 NEW_FILENAME
	INTEGER IOS,ERROR_LU
	EXTERNAL ERROR_LU
!                      
	NEW_FILENAME=TRIM(FILENAME)//'_INFO'
        OPEN(UNIT=LU_EDD,FILE=NEW_FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(ERROR_LU(),*)'Unable to open',NEW_FILENAME
	    WRITE(ERROR_LU(),*)'Error occured in READ_DIRECT_INFO'
	  END IF
          READ(LU_EDD,'(2(3X,I6))')ND,RECL
        CLOSE(LU_EDD)
!
	RETURN
	END
