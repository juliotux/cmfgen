!
! Simple little routine to check whether the ITS_DONE keyword is included in
! the HYDRO_DEFAULTS file. If not, the keyowrd is added. This is done
! for compatability with earlier versions.
!
	SUBROUTINE CHECK_HYDRO_DEF(STRING,LUIN,LUER)
	IMPLICIT NONE
!
! Created: 04-April-2011
!
	CHARACTER(LEN=*) STRING
	INTEGER LUER
	INTEGER LUIN
!
	INTEGER IOS
!
        OPEN(UNIT=LUIN,FILE='HYDRO_DEFAULTS',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
          IF(IOS .NE. 0)THEN
            WRITE(LUER,*)'ERROR -- DO_HYDRO is set but HYDRO_DEFAULTS cannot be opened'
            STOP
          END IF
	  IOS=0
          DO WHILE(IOS .EQ. 0)
            READ(LUIN,'(A)',IOSTAT=IOS)STRING
            IF(INDEX(STRING,'[ITS_DONE]') .NE. 0)EXIT
          END DO
        CLOSE(LUIN)
!
        IF(IOS .NE. 0)THEN
          WRITE(LUER,*)'Adding ITS_DONE keyword to HYDRO_DEFAULTS'
          OPEN(UNIT=LUIN,FILE='HYDRO_DEFAULTS',STATUS='OLD',POSITION='APPEND')
            WRITE(LUIN,'(A,12X,A)')'0','[ITS_DONE]'
          CLOSE(LUIN)
        END IF
!
	RETURN
	END
