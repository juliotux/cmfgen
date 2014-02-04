	PROGRAM COUNT_PHOT_VALUES
	IMPLICIT NONE
!
	INTEGER COUNT_DATA
	INTEGER COUNT_LEVS
	INTEGER NC
	INTEGER TYPE
	INTEGER I,K
	REAL*8 FREQ,  OLD_FREQ
	REAL*8 CROSS, OLD_CROSS
!	
	CHARACTER*200 STRING
	CHARACTER*80 FILENAME
	CHARACTER*30 LEVEL_NAME
!
	WRITE(6,*)'Input PHOT file name'
	READ(5,*)FILENAME
!
	COUNT_LEVS=0
	COUNT_DATA=0
!
	OPEN(UNIT=10,FILE=FILENAME,STATUS='OLD',ACTION='READ')
	STRING=' '
	DO WHILE (INDEX(STRING,'!Configuration name') .EQ. 0)
	  READ(10,FMT='(A)',END=100)STRING
	  WRITE(6,*)TRIM(STRING)
	END DO
!
	DO WHILE (1 .EQ. 1)
!
	  WRITE(6,*)TRIM(STRING)
	  STRING=ADJUSTL(STRING)
	  K=INDEX(STRING,'  ')
	  LEVEL_NAME=STRING(1:K)
	  READ(10,FMT='(A)',END=100)STRING
	  IF(INDEX(STRING,'!Type of cross-section') .EQ. 0)THEN
	    WRITE(6,*)'Error- TYPE record does not follow Config record'
	    WRITE(6,*)LEVEL_NAME
	    STOP
	  END IF
!
	  READ(STRING,*)TYPE
	  READ(10,FMT='(A)',END=100)STRING
          IF(INDEX(STRING,'!Number of cross-section points') .EQ. 0)THEN
	    WRITE(6,*)'Error- Number of cross-section points does not follow TYPE record'
	    WRITE(6,*)LEVEL_NAME
	    STOP
	  ELSE
	    COUNT_LEVS=COUNT_LEVS+1
	    READ(STRING,*)NC
	    COUNT_DATA=COUNT_DATA+NC
	    IF(TYPE .EQ. 20 .OR. TYPE .EQ. 21 .OR. TYPE .EQ. 22)THEN
	      READ(10,*)FREQ,CROSS
	      DO I=2,NC
	        OLD_FREQ=FREQ; OLD_CROSS=CROSS
	        READ(10,*)FREQ,CROSS
	        IF(FREQ .LT. OLD_FREQ)THEN
	           WRITE(6,*)'Error in cross-section: frequency grid is non-monotonic'
	           WRITE(6,*)'      NAME=',LEVEL_NAME
	           WRITE(6,*)'      TYPE=',TYPE
	           WRITE(6,*)'        NC=',NC
	           WRITE(6,*)' LOW_FREQ=',OLD_FREQ
	           WRITE(6,*)'HIGH_FREQ=',FREQ
	           WRITE(6,*)'    CROSS=',CROSS
	           STOP
	        END IF
	        IF(FREQ .EQ. OLD_FREQ)THEN
	           WRITE(6,*)'Warning in cross-section: equal frequencis'
	           WRITE(6,*)'      NAME=',LEVEL_NAME
	           WRITE(6,*)'      TYPE=',TYPE,'        NC=',NC
	           WRITE(6,*)' LOW_FREQ=',OLD_FREQ,'HIGH_FREQ=',FREQ
	           WRITE(6,*)'    CROSS=',CROSS
	        END IF
	        IF(CROSS .LT. 0)THEN
	           WRITE(6,*)'Warning in cross-section: negatve cross-section'
	           WRITE(6,*)'      NAME=',LEVEL_NAME
	           WRITE(6,*)'      TYPE=',TYPE,'        NC=',NC
	           WRITE(6,*)'    CROSS=',CROSS
	        END IF
	      END DO
	    ELSE
	      READ(10,*)(FREQ,I=1,NC)
	    END IF
	  END IF
!
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Configuration name') .EQ. 0)
	    IF(STRING .NE. ' ')THEN
	      WRITE(6,*)'Invalid format: unexpected string'
	      WRITE(6,*)TRIM(STRING)
	      STOP
	    END IF
	    READ(10,FMT='(A)',END=100)STRING
	  END DO
	  WRITE(6,*)'2nd: ',TRIM(STRING)
	END DO
!
100	CONTINUE
	WRITE(6,*)'Numer of levels is ',COUNT_LEVS
	WRITE(6,*)'Number of data values is',COUNT_DATA
!
	STOP
	END
