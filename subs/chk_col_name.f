!
! Subroutione to check whether non-matching level names found in a
! collisional data file simply belong to higher levels. In that case,
! no error message is output.
!
	SUBROUTINE CHK_COL_NAME(NO_MATCH_NAME,NUM_NO_MATCH,NLEV,LUIN,LUER,FILE_NAME)
	IMPLICIT NONE
!
! Created: 20-Jan-2014
!
	INTEGER NUM_NO_MATCH
	INTEGER NLEV
	INTEGER LUIN
	INTEGER LUER
	CHARACTER(LEN=*) NO_MATCH_NAME(NUM_NO_MATCH)
	CHARACTER(LEN=*) FILE_NAME
!
	CHARACTER(LEN=40) NAME
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) LOC_FILE
	INTEGER I,J,K
	INTEGER NLOC
	INTEGER CNT
!
	I=INDEX(FILE_NAME,'_')
	LOC_FILE=FILE_NAME(1:I)//'F_OSCDAT'
	OPEN(UNIT=LUIN,FILE=TRIM(LOC_FILE),STATUS='OLD',ACTION='READ')
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Number of energy levels') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
	READ(STRING,*)NLOC
	DO WHILE( INDEX(STRING,'!Number of transitions') .EQ. 0)
	  READ(LUIN,'(A)')STRING
	END DO
!
	STRING=' '
	DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	  READ(LUIN,'(A)')STRING
	END DO
!
! We have alredy checked the first NLEVs for a match, and we have already
! read the first name.
!
	DO I=2,NLEV
	  READ(LUIN,'(A)')STRING
	END DO
!
! Now read in additional levels, not included in the model, to see if they match.
!
	CNT=0
	DO I=NLEV+1,NLOC
	  READ(LUIN,'(A)')STRING
	  STRING=ADJUSTL(STRING)
	  NAME=STRING(1:INDEX(STRING,'  '))
	  DO K=1,NUM_NO_MATCH
	    J=INDEX(NAME,'[')
	    IF(NAME .EQ. NO_MATCH_NAME(K))THEN
	      NO_MATCH_NAME(K)=' '
	      CNT=CNT+1
	    ELSE IF(J .NE. 0)THEN
	      IF(NAME(1:J-1) .EQ. NO_MATCH_NAME(K))THEN
	        NO_MATCH_NAME(K)=' '
	        CNT=CNT+1
	      END IF
	    END IF
	  END DO
	  IF(CNT .EQ. NUM_NO_MATCH)EXIT
	END DO 
!
	CLOSE(LUIN)
	CNT=0
	DO K=1,NUM_NO_MATCH
	  IF(NO_MATCH_NAME(K) .NE. ' ')THEN
	    CNT=CNT+1
	    WRITE(LUER,'(A,A,A,A)')' Error(',TRIM(FILE_NAME),
	1                   '): No match found for level ',TRIM(NO_MATCH_NAME(K))
          END IF
	END DO
	IF(CNT .NE. 0)WRITE(LUER,'(A)')' '
!
	RETURN
	END
