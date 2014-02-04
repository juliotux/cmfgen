!
! Subroutine to write out information needed to reopen a DIRECT access file.
! The RECL is not stored with the file on UNIX systems.
!
! Modifed 20-August-2000 so that a format date can be ouput
!
	SUBROUTINE WRITE_DIRECT_INFO_V2(ND,RECL,FILE_DATE,FILENAME,LU_EDD)
	IMPLICIT NONE
!
! Created 16-Jun-2000
!
	INTEGER ND			!Number of depth points written
	INTEGER RECL			!Record length (system dependent)
	INTEGER LU_EDD
	CHARACTER*(*) FILENAME
	CHARACTER*(*) FILE_DATE
!
	CHARACTER*80 NEW_FILENAME
!                      
	NEW_FILENAME=TRIM(FILENAME)//'_INFO'
        OPEN(UNIT=LU_EDD,FILE=NEW_FILENAME,STATUS='REPLACE')
          WRITE(LU_EDD,'(X,A,39X,A)')'20-Aug-2000','!INFO format date'
          WRITE(LU_EDD,'(X,A,39X,A)')TRIM(FILE_DATE),'!File format date'
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
! Modifed 20-August-2000 so that a format date can be read.
!
	SUBROUTINE READ_DIRECT_INFO_V2(ND,RECL,FILE_DATE,FILENAME,LU_EDD)
	IMPLICIT NONE
!
! Created 16-Jun-2000
!
	INTEGER ND			!Number of depth points written
	INTEGER RECL			!Record length (system dependent)
	INTEGER LU_EDD
	CHARACTER*(*) FILENAME
	CHARACTER*(*) FILE_DATE
!
	CHARACTER*80 NEW_FILENAME
	CHARACTER*80 STRING
	INTEGER IOS,IER,J,ERROR_LU
	EXTERNAL ERROR_LU
!                      
	NEW_FILENAME=TRIM(FILENAME)//'_INFO'
        OPEN(UNIT=LU_EDD,FILE=NEW_FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(ERROR_LU(),*)'Unable to open',NEW_FILENAME
	    WRITE(ERROR_LU(),*)'Error occured in READ_DIRECT_INFO'
	  END IF
	  STRING='!'
	  DO WHILE(STRING(1:1) .EQ. '!')
	    IER=1; READ(LU_EDD,'(A)',ERR=100,IOSTAT=IOS)STRING
	  END DO
	  IF(INDEX(STRING,'INFO format date') .EQ. 0)THEN
            IER=2; READ(STRING,*,ERR=100,IOSTAT=IOS)ND,RECL
	    FILE_DATE='Unavailable'
	  ELSE
	    IER=3; READ(LU_EDD,'(A)',ERR=100,IOSTAT=IOS)STRING
	    IER=4; IF(INDEX(STRING,'File format date') .EQ. 0)GOTO 100
	    FILE_DATE=STRING(1:INDEX(STRING,'  '))
	    IER=5; READ(LU_EDD,*,ERR=100,IOSTAT=IOS)ND,RECL
	  END IF
        CLOSE(LU_EDD)
!
	RETURN
!
100	J=ERROR_LU()
	WRITE(J,*)'Error in READ_DIRECT_INFO_V2'
	WRITE(J,*)'Currently reading: ',NEW_FILENAME
	WRITE(J,*)'Internal subroutine error number is:',IER
	WRITE(J,*)'Frotran error number is:',IOS
	STOP
	END
