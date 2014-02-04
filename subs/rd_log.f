C
C 21-Nov-2000 : BACKSPACE statement removed from all routines.
C               Now use a free-format read from STRING.
C 24-Mar-1993 : Created - Based on routines in RDLOG.
C	        All calls now RD_...
C               KEY inserted to allow checking of input data, especially
C               useful when input file formats have changed.
C               At present only blank lines are allowed in file.
C               In future versions, by using a buffer, it may be possible
C               to read in only necessary keywords.
C
	SUBROUTINE RD_LOG(VALUE,KEY,LUI,LUO,A)
	IMPLICIT NONE
	LOGICAL VALUE
	INTEGER LUI,LUO
	CHARACTER*(*) KEY,A
	CHARACTER*80 STRING
C
	STRING=' '
	DO WHILE (INDEX(STRING,'[') .EQ. 0)		!Skip blank lines
	  READ(LUI,'(A)')STRING
	  IF (INDEX(STRING,'[') .EQ. 0)WRITE(LUO,'(A)')STRING
	END DO
C
C Check if KEY matches
C
	IF( INDEX(STRING,'['//KEY//']') .EQ. 0)THEN	
	  WRITE(LUO,*)'Error in RD_LOG - KEY deos not match'
	  WRITE(LUO,*)'KEY=',KEY
	  STOP
	END IF
C
C Now read in value
C
	READ(STRING,*)VALUE
C
	WRITE(LUO,10)VALUE,KEY,A
10	FORMAT(12X,L1,5X,'[',A,']',T35,A)
C
	RETURN
	END
C
C
	SUBROUTINE RD_2LOG(VALUE1,VALUE2,KEY,LUI,LUO,A)
	IMPLICIT NONE
	LOGICAL VALUE1,VALUE2
	INTEGER LUI,LUO
	CHARACTER*(*) KEY,A
	CHARACTER*80 STRING
C
	STRING=' '
	DO WHILE (INDEX(STRING,'[') .EQ. 0)		!Skip blank lines
	  READ(LUI,'(A)')STRING
	  IF (INDEX(STRING,'[') .EQ. 0)WRITE(LUO,'(A)')STRING
	END DO
C
C Check if KEY matches
C
	IF( INDEX(STRING,'['//KEY//']') .EQ. 0)THEN	
	  WRITE(LUO,*)'Error in RD_2LOG - KEY deos not match'
	  WRITE(LUO,*)'KEY=',KEY
	  STOP
	END IF
C
C Now read in value
C
	READ(STRING,*)VALUE1,VALUE2
	WRITE(LUO,10)VALUE1,VALUE2,KEY,A
10	FORMAT(10X,L1,',',L1,5X,'[',A,']',T35,A)
C
	RETURN
	END
C
C
	SUBROUTINE RD_INT(VALUE,KEY,LUI,LUO,A)
	IMPLICIT NONE
	INTEGER VALUE,LUI,LUO
	CHARACTER*(*) KEY,A
	CHARACTER*80 STRING
C
	STRING=' '
	DO WHILE (INDEX(STRING,'[') .EQ. 0)		!Skip blank lines
	  READ(LUI,'(A)')STRING
	  IF (INDEX(STRING,'[') .EQ. 0)WRITE(LUO,'(A)')STRING
	END DO
C
C Check if KEY matches
C
	IF( INDEX(STRING,'['//KEY//']') .EQ. 0)THEN	
	  WRITE(LUO,*)'Error in RD_INT - KEY deos not match'
	  WRITE(LUO,*)'KEY=',KEY
	  STOP
	END IF
C
C Now read in value
C
	READ(STRING,*)VALUE
	WRITE(LUO,10)VALUE,KEY,A
10	FORMAT(5X,I8,5X,'[',A,']',T35,A)
C
	RETURN
	END
C
C
	SUBROUTINE RD_DBLE(VALUE,KEY,LUI,LUO,A)
	IMPLICIT NONE
	REAL*8 VALUE
	INTEGER LUI,LUO
	CHARACTER*(*) KEY,A
	CHARACTER*80 STRING
C
	STRING=' '
	DO WHILE (INDEX(STRING,'[') .EQ. 0)		!Skip blank lines
	  READ(LUI,'(A)')STRING
	  IF (INDEX(STRING,'[') .EQ. 0)WRITE(LUO,'(A)')STRING
	END DO
C
C Check if KEY matches
C
	IF( INDEX(STRING,'['//KEY//']') .EQ. 0)THEN	
	  WRITE(LUO,*)'Error in RD_DBLE - KEY deos not match'
	  WRITE(LUO,*)'KEY=',KEY
	  STOP
	END IF
C
C Now read in value
C
	READ(STRING,*)VALUE
	WRITE(LUO,10)VALUE,KEY,A
10	FORMAT(1X,1PE12.5,5X,'[',A,']',T35,A)
C
	END
C
C 
C  Now ignores blank lines in input file.
C
	SUBROUTINE RD_CHAR(VALUE,KEY,LUI,LUO,A)
	IMPLICIT NONE
	INTEGER LUI,LUO
	CHARACTER*(*) KEY,A
	CHARACTER VALUE*6,STRING*80
C
	STRING=' '
	DO WHILE (INDEX(STRING,'[') .EQ. 0)		!Skip blank lines
	  READ(LUI,'(A)')STRING
	  IF (INDEX(STRING,'[') .EQ. 0)WRITE(LUO,'(A)')STRING
	END DO
C
C Check if KEY matches
C
	IF( INDEX(STRING,'['//KEY//']') .EQ. 0)THEN	
	  WRITE(LUO,*)'Error in RD_CHAR - KEY deos not match'
	  WRITE(LUO,*)'KEY=',KEY
	  STOP
	END IF
C
C Now set value.
C
	VALUE=STRING(1:6)
	WRITE(LUO,10)VALUE,KEY,A
10	FORMAT(7X,A6,5X,'[',A,']',T35,A)
C
	RETURN
	END
C
C 
C To replace RDCHAR. The number of characters to be read in is
C now secified in the call. Now skips blank strings.
C
	SUBROUTINE RD_NCHAR(VALUE,KEY,NCHAR,LUI,LUO,A)
	IMPLICIT NONE!
!
! Altered 21-Nov-2000: Changed to allow long strings to be read.
!
	INTEGER LUI,LUO,NCHAR,KEY_LOC,K
	CHARACTER*(*) KEY,A,VALUE
	CHARACTER*80 STRING
C
	STRING=' '
	DO WHILE (INDEX(STRING,'[') .EQ. 0)		!Skip blank lines
	  READ(LUI,'(A)')STRING
	  IF (INDEX(STRING,'[') .EQ. 0)WRITE(LUO,'(A)')STRING
	END DO
C
C Check if KEY matches
C
	KEY_LOC=INDEX(STRING,'['//KEY//']') 
	IF(KEY_LOC .EQ. 0)THEN	
	  WRITE(LUO,*)'Error in RD_NCHAR - KEY deos not match'
	  WRITE(LUO,*)'KEY=',KEY
	  STOP
	END IF
C
C Now set value
C
	VALUE=' '
	VALUE=STRING(1:KEY_LOC-1)
	K=LEN_TRIM(VALUE)
	IF(NCHAR .LE. K)THEN
	  STRING=' '
	  STRING(13-K:12)=VALUE
	END IF
	WRITE(LUO,10)STRING(1:MAX(12,K)),KEY,A
10	FORMAT(1X,A,5X,'[',A,']',T35,A)
C
	RETURN
	END
