	PROGRAM ED_VADAT
	IMPLICIT NONE
!
! Altered 01-Sep-2015: TIME_SEQ_NO changed from integer to real.
!
! Total number of models 
!
	INTEGER NUM_ITS
!
	REAL*8 RSTAR_MIN,RSTAR_MAX
	REAL*8 LSTAR_MIN,LSTAR_MAX
!
	REAL*8 RSTAR
	REAL*8 LSTAR
!
	INTEGER I,J,K
	REAL*8 TIME_SEQ_NO
	INTEGER NUM_STR
	REAL*8 R_SCALE_FACTOR
	REAL*8 L_SCALE_FACTOR
!
	CHARACTER*132 STRING(1000)
	CHARACTER*132 TMP_STR
!
	WRITE(6,*)'Input RSTAR_MIN & RSTAR_MAX'
	READ(5,*)RSTAR_MIN,RSTAR_MAX
	WRITE(6,*)'Input LSTAR_MIN & LSTAR_MAX'
	READ(5,*)LSTAR_MIN,LSTAR_MAX
	WRITE(6,*)'Input number of iteration, inclusive, to cover RSTAR range'
	READ(5,*)NUM_ITS
!
! Evaluate scale factors.
!
	R_SCALE_FACTOR=EXP(LOG(RSTAR_MAX/RSTAR_MIN)/(NUM_ITS-1))
	L_SCALE_FACTOR=EXP(LOG(LSTAR_MAX/LSTAR_MIN)/(NUM_ITS-1))
!
! Read in VADAT file for editing.
!
	OPEN(UNIT=9,FILE='OLD_VADAT',STATUS='OLD',ACTION='READ')
	  NUM_STR=0
	  DO I=1,1000
	    READ(9,'(A)',END=100)STRING(I)
	    NUM_STR=I
	  END DO
	  WRITE(6,*)'Error --- not enough records available for VADAT file'
	  STOP
100	CONTINUE
	CLOSE(UNIT=9)
!
! Modify RSTAR
!
	J=0
	DO I=1,NUM_STR
	  J=INDEX(STRING(I),'[RSTAR]')
	  IF(J .NE. 0)THEN
	    READ(STRING(I),*)RSTAR
	    RSTAR=RSTAR*R_SCALE_FACTOR
	    WRITE(TMP_STR,'(ES12.6,4X)')RSTAR
	    STRING(I)=TMP_STR(1:16)//STRING(I)(J:)
	    EXIT
	  END IF
	END DO
!
! Modify the luminosity
!
	DO I=1,NUM_STR
	  J=INDEX(STRING(I),'[LSTAR]')
	  IF(J .NE. 0)THEN
	    READ(STRING(I),*)LSTAR
	    LSTAR=LSTAR*L_SCALE_FACTOR
	    WRITE(TMP_STR,'(ES12.6,4X)')LSTAR
	    STRING(I)=TMP_STR(1:16)//STRING(I)(J:)
	    EXIT
	  END IF
	END DO
!
! Ensure D/Dt term switched on.
!
	DO I=1,NUM_STR
	  J=INDEX(STRING(I),'[DO_DDT]')
	  IF(J .NE. 0)THEN
	    TMP_STR=' '; WRITE(TMP_STR,'(A)')'T'
	    STRING(I)=TMP_STR(1:16)//STRING(I)(J:)
	    EXIT
	  END IF
	END DO
	IF(J .EQ. 0)THEN
	  WRITE(6,*)'Error DO_DDT string not present'
	  STOP
	END IF
!
! Update time sequence number
!
	DO I=1,NUM_STR
	  J=INDEX(STRING(I),'[TS_NO]')
	  IF(J .NE. 0)THEN
	    READ(STRING(I),*)TIME_SEQ_NO
	    TIME_SEQ_NO=TIME_SEQ_NO+1
	    TMP_STR=' '; WRITE(TMP_STR,'(F8.3)')TIME_SEQ_NO
	    TMP_STR=ADJUSTL(TMP_STR)
	    K=LEN_TRIM(TMP_STR)
	    DO WHILE(K .GT. 1)
	      IF(TMP_STR(K:K) .EQ. '0')THEN
	        TMP_STR(K:K)=' '
	        K=K-1
	      ELSE
	        EXIT
	      END IF
	    END DO
	    STRING(I)=TMP_STR(1:16)//STRING(I)(J:)
	    EXIT
	  END IF
	END DO
	IF(J .EQ. 0)THEN
	  WRITE(6,*)'Error TS_NO (TIME_SEQ_NO) string not present'
	  STOP
	END IF
!
! Ensure D/Dt adiabatic cooling is switched on.
!
	DO I=1,NUM_STR
	  J=INDEX(STRING(I),'[INC_AD]')
	  IF(J .NE. 0)THEN
	    TMP_STR=' '; WRITE(TMP_STR,'(A)')'T'
	    STRING(I)=TMP_STR(1:16)//STRING(I)(J:)
	    EXIT
	  END IF
	END DO
	IF(J .EQ. 0)THEN
	  WRITE(6,*)'Error INC_AD string not present'
	  STOP
	END IF
!
! Finished string modification. Now outout the modified file.
!
	OPEN(UNIT=10,FILE='VADAT',STATUS='NEW',ACTION='WRITE')
	  DO I=1,NUM_STR
	    WRITE(10,'(A)')TRIM(STRING(I))
	  END DO
	CLOSE(UNIT=10)
!
	STOP
	END
