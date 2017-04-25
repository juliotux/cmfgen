!
! Routine to read a single double precisions vector, of length ND, from RVTJ
! or a similarly formatted file. KEY is used to locate the required vector.
! If the read fails, a non-zero value of IOS is returned. 
!
	SUBROUTINE RD_SING_VEC_RVTJ(XV,ND,KEY,FILE_NAME,LU,IOS)
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER IOS
	INTEGER LU
	REAL*8 XV(ND)
	CHARACTER(LEN=*) FILE_NAME,KEY
!
	INTEGER I
	CHARACTER(LEN=80)STRING
!
	OPEN(UNIT=LU,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Unable to open ',TRIM(FILE_NAME),' in RD_SING_VEC_RVTJ'
	    WRITE(6,*)'IOS=',IOS
	    RETURN
	  END IF
!
	DO WHILE(1 .EQ. 1)
	  READ(LU,'(A)',END=1000)STRING
	  IF(INDEX(STRING,TRIM(KEY)) .NE. 0)THEN
	    READ(LU,*,IOSTAT=IOS)(XV(I),I=1,ND)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Error reading vector in RD_SING_VEC_RVTJ'
	      WRITE(6,*)'IOS=',IOS
	      CLOSE(LU)
	    END IF
	    CLOSE(LU)
	    RETURN
	  END IF
	END DO
!
1000	CONTINUE
	WRITE(6,*)'Key not found in RD_SING_VEC_RVTJ'
	WRITE(6,*)'KEY is: ',TRIM(KEY)
	IOS=1
	CLOSE(LU)
	RETURN
!
	END
