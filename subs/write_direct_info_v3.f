!
! Subroutine to write out information needed to reopen a DIRECT access file.
! The RECL is not stored with the file on UNIX systems.
!
! Modified 20-Aug-2000: Format date can be output. Alterations made to minimize
!                         changes to other routines. Changes should make 
!                         the transfer of DIRECT access files between PENTIUM and
!                         ALPHA systems transparent.
! Modified 18-Aug-2003: UNIT_SIZE and WORD_SIZE output to INFO file.
!
	SUBROUTINE WRITE_DIRECT_INFO_V3(ND,RECL,FILE_DATE,FILENAME,LU_EDD)
	IMPLICIT NONE
!
! Created 16-Jun-2000
!
	INTEGER ND			!Number of depth points written
	INTEGER RECL			!Record length (system dependent)
	INTEGER LU_EDD
!
! Keywords to describe data format.
!
	INTEGER REC_SIZE_LIM
	INTEGER WORD_SIZE
	INTEGER UNIT_SIZE
	INTEGER MAX_NUM_REC
!
	CHARACTER*(*) FILENAME
	CHARACTER*(*) FILE_DATE
!
	CHARACTER*80 NEW_FILENAME
!
! Get system dependent direct access parameters.
!
	CALL DIR_ACC_PARS(REC_SIZE_LIM,UNIT_SIZE,WORD_SIZE,MAX_NUM_REC)
!
! Create INFO file.
!
! ND is a USE dependent integer. e.g., # of numbers written out.
! RECL is the length of the file in system units.
!
	NEW_FILENAME=TRIM(FILENAME)//'_INFO'
        OPEN(UNIT=LU_EDD,FILE=NEW_FILENAME,STATUS='REPLACE')
          WRITE(LU_EDD,'(X,A,39X,A)')'18-Aug-2003','!INFO format date'
          WRITE(LU_EDD,'(X,A,39X,A)')TRIM(FILE_DATE),'!File format date'
          WRITE(LU_EDD,'(4(6X,I6))')ND,RECL,WORD_SIZE,UNIT_SIZE
          WRITE(LU_EDD,'(4(3X,A))')'       ND','     RECL','WORD_SIZE','UNIT_SIZE'
        CLOSE(LU_EDD)
!
	RETURN
	END
!
! Subroutine to read out information needed to reopen a DIRECT access file.
! The RECL is not stored with the file on UNIX systems.
!
! Modified 10-Sep-2008: Problem with integer divide when computing RECL.
!                       Error occurred switching from PGF compiler (unit size=1)
!                         to Intel compiler (unit size=4). 1/4 was giving 0
!                         instead of 0.25
! Modified 28-Oct-2001: IOS returned in call 
! Modified 20-August-2000 so that a format date can be read.
!
	SUBROUTINE READ_DIRECT_INFO_V3(ND,RECL,FILE_DATE,FILENAME,LU_EDD,IOS)
	IMPLICIT NONE
!
! Altered 24-Apr-2016 : Increaed allowed filename length to 127 (132 with _INFO).
! Created 16-Jun-2000
!
	INTEGER ND			!Number of depth points written
	INTEGER RECL			!Record length (system dependent)
	INTEGER LU_EDD
	INTEGER IOS
	CHARACTER*(*) FILENAME
	CHARACTER*(*) FILE_DATE
!
! Keywords to describe data format.
!
	INTEGER REC_SIZE_LIM
	INTEGER WORD_SIZE
	INTEGER UNIT_SIZE
	INTEGER MAX_NUM_REC
!
	INTEGER OLD_WORD_SIZE
	INTEGER OLD_UNIT_SIZE
!
	CHARACTER(LEN=132) NEW_FILENAME
	CHARACTER(LEN=132) STRING
	INTEGER IER,J,ERROR_LU
	EXTERNAL ERROR_LU
!
	CALL DIR_ACC_PARS(REC_SIZE_LIM,UNIT_SIZE,WORD_SIZE,MAX_NUM_REC)
!
	IOS=0
	J=LEN_TRIM(FILENAME)
	IF(J .GT. 127)THEN
	  WRITE(ERROR_LU(),'(A,A)')' Temporary string too short - LEN(FILENAME)=',J
	  WRITE(ERROR_LU(),*)'FILENAME=',TRIM(FILENAME)
	  WRITE(ERROR_LU(),*)'Error occurred in READ_DIRECT_INFO: IOS=',IOS
	  STOP
	END IF
!
	NEW_FILENAME=TRIM(FILENAME)//'_INFO'
        OPEN(UNIT=LU_EDD,FILE=NEW_FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(ERROR_LU(),'(A,A)')' Unable to open ',TRIM(NEW_FILENAME)
	    WRITE(ERROR_LU(),*)'Error occurred in READ_DIRECT_INFO: IOS=',IOS
	    RETURN
	  END IF
	  STRING='!'
	  DO WHILE(STRING(1:1) .EQ. '!')
	    IER=1; READ(LU_EDD,'(A)',ERR=100,IOSTAT=IOS)STRING
	  END DO
!
! What we read depends on when the INFO file was written.
!
	  IF(INDEX(STRING,'INFO format date') .EQ. 0)THEN
            IER=2; READ(STRING,*,ERR=100,IOSTAT=IOS)ND,RECL
	    FILE_DATE='Unavailable'
	  ELSE
	    IF(INDEX(STRING,'18-Aug-2003') .NE. 0)THEN
	      IER=3; READ(LU_EDD,'(A)',ERR=100,IOSTAT=IOS)STRING
	      IER=4; IF(INDEX(STRING,'File format date') .EQ. 0)GOTO 100
	      FILE_DATE=STRING(1:INDEX(STRING,'  '))
	      IER=5; READ(LU_EDD,*,ERR=100,IOSTAT=IOS)ND,RECL,OLD_WORD_SIZE,OLD_UNIT_SIZE
!	      RECL=RECL*(WORD_SIZE/OLD_WORD_SIZE)*(OLD_UNIT_SIZE/UNIT_SIZE)
	      RECL=RECL*WORD_SIZE*OLD_UNIT_SIZE/OLD_WORD_SIZE/UNIT_SIZE
	    ELSE IF(INDEX(STRING,'20-Aug-2000') .NE. 0)THEN
	      IER=6; READ(LU_EDD,'(A)',ERR=100,IOSTAT=IOS)STRING
	      IER=7; IF(INDEX(STRING,'File format date') .EQ. 0)GOTO 100
	      FILE_DATE=STRING(1:INDEX(STRING,'  '))
	      IER=8; READ(LU_EDD,*,ERR=100,IOSTAT=IOS)ND,RECL
	    ELSE
	      IER=10
	      GOTO 100
	    END IF
	  END IF
        CLOSE(LU_EDD)
!
	RETURN
!
100	J=ERROR_LU()
	WRITE(J,*)'Error in READ_DIRECT_INFO_V3'
	WRITE(J,*)'Currently reading: ',TRIM(NEW_FILENAME)
	WRITE(J,*)'Internal subroutine error number is:',IER
	IF(IOS .NE. 0)WRITE(J,*)'Fortran error number is:',IOS
	RETURN
	END
