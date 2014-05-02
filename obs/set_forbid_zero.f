	SUBROUTINE SET_FORBID_ZERO(LU)
	USE MOD_CMF_OBS
	IMPLICIT NONE
!
	INTEGER LU
	INTEGER ID
	INTEGER I,J
	CHARACTER(LEN=20) STRING
	LOGICAL FILE_PRES
!
	INQUIRE(FILE='FORB_LINE_CONTROL',EXIST=FILE_PRES)
	IF(.NOT. FILE_PRES)RETURN
!
	WRITE(6,'(A)')' '
	WRITE(6,*)('#',I=1,75)
	WRITE(6,*)('#',I=1,75)
	WRITE(6,'(A)')' '
	OPEN(UNIT=LU,FILE='FORB_LINE_CONTROL',STATUS='OLD',ACTION='READ')
!
	DO WHILE(1 .EQ. 1)
	  STRING='!'
	  DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	    READ(LU,'(A)',END=100)STRING
	  END DO
	  STRING=ADJUSTL(STRING)
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (TRIM(STRING) .EQ. ION_ID(ID) .OR.
	1                              TRIM(STRING) .EQ. 'ALL') )THEN
	      WRITE(6,*)'Removing forbidden transitions for species ',TRIM(ION_ID(ID))
	      DO I=1,ATM(ID)%NXzV_F-1
	        IF(INDEX(ATM(ID)%XzVLEVNAME_F(I),'e') .NE. 0)THEN
	          DO J=I+1,ATM(ID)%NXzV_F
	            IF(INDEX(ATM(ID)%XzVLEVNAME_F(J),'e') .NE. 0)THEN
	              ATM(ID)%AXzV_F(I,J)=0.0D0
	              ATM(ID)%AXzV_F(J,I)=0.0D0
	            END IF
	          END DO
	        END IF
	        IF(INDEX(ATM(ID)%XzVLEVNAME_F(I),'o') .NE. 0)THEN
	          DO J=I+1,ATM(ID)%NXzV_F
	            IF(INDEX(ATM(ID)%XzVLEVNAME_F(J),'o') .NE. 0)THEN
	              ATM(ID)%AXzV_F(I,J)=0.0D0
	              ATM(ID)%AXzV_F(J,I)=0.0D0
	            END IF
	          END DO
	        END IF
	      END DO
	    END IF
	  END DO
	END DO
100	CONTINUE
!
	WRITE(6,'(A)')' '
	WRITE(6,*)('#',I=1,75)
	WRITE(6,*)('#',I=1,75)
	WRITE(6,'(A)')' '
!
	RETURN
	END
