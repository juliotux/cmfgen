!
! Subroutine to read in T and ED from a DC file. This T and E will be used to
! set profile limits. Accuracy is not crucial.
!
	SUBROUTINE RD_T_ED(T,ED,ND,LU_IN,FILE_NAME)
	IMPLICIT NONE
!
! Created: 04-Nov-2013
!
	INTEGER ND
	INTEGER LU_IN
	REAL*8 T(ND)
	REAL*8 ED(ND)
	REAL*8 T1,T2
!
	INTEGER LUER
	INTEGER ERROR_LU
	INTEGER I,J
	INTEGER IOS
	INTEGER NOLD,NDOLD
	EXTERNAL ERROR_LU
	CHARACTER(LEN=*) FILE_NAME
	CHARACTER(LEN=80) STRING
!
	LUER=ERROR_LU()
!
! Read in values from previous model.
!
	OPEN(UNIT=LU_IN,STATUS='OLD',FILE=FILE_NAME,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',TRIM(FILE_NAME),' in RD_T_ED'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
! Check whether the file has a record containing 'Format date'. Its presence
! effects the way we read the file. 
!
        I=0
        STRING=' '
        DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
          I=I+1
          READ(LU_IN,'(A)')STRING
        END DO
        IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU_IN)
!
	READ(LU_IN,*)T1,T2,NOLD,NDOLD
	IF(ND .NE. NDOLD)THEN
	  WRITE(6,*)'Error ND .NE. NDOLD in RD_T_ED'
	  STOP
	END IF
!
	DO I=1,NDOLD
	  READ(LU_IN,*)T1,T2,ED(I),T(I)
	  READ(LU_IN,*)(T1,J=1,NOLD)
	END DO
	CLOSE(LU_IN)
!
	RETURN
	END
