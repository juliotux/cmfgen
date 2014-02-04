	SUBROUTINE RD_NON_THERM_SPEC(FELEC,FION,FEXC,ND,FILE_NAME,LU_IN,IOS)
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER LU_IN
	INTEGER IOS
	REAL*8 FELEC(ND)
	REAL*8 FION(ND)
	REAL*8 FEXC(ND)
	CHARACTER(LEN=*) FILE_NAME
!
	REAL*8 T1
	INTEGER I
	CHARACTER(LEN=132) STRING
!
	OPEN(UNIT=LU_IN,FILE=FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)RETURN
!
	DO WHILE(INDEX(STRING,'Felec') .EQ. 0)
	  READ(LU_IN,'(A)')STRING
	END DO
!
	DO I=1,ND
	  READ(LU_IN,*)T1,FELEC(I),T1,FION(I),T1,FEXC(I)
	END DO
!
	CLOSE(LU_IN)
	RETURN
	END