! Altered: 23-Mar-2012:  Improved error reporting.
! Altered: 11-Mar-2008:  For output changed T35 to T40
!           3-Mar-2000:  Created - Based on routines in RD_LOG.
!	                 All calls now RD_STORE_...
!                        Data is first read in from a file and stored.  
!                        Options are then read from store in ANY order.
!                        Only those options requested are checked.
!                        KEYS are checked for uniqueness (7-Jun-2000)
!             
!             Option access should begin with:
!	          CALL RD_OPTIONS_INTO_STORE(LU_IN,LU_OUT)
!             and end with
!	          CALL CLEAN_RD_STORE
!
	MODULE RD_VAR_MOD
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NST_MAX=700
	CHARACTER(LEN=80), ALLOCATABLE :: STORE(:)
	INTEGER, ALLOCATABLE ::  KEY_ST(:)
	INTEGER, ALLOCATABLE :: KEY_END(:)
!
	INTEGER NST
	INTEGER I_UP,I_DWN
	INTEGER LUER
	INTEGER LUO
	INTEGER IOS
!
	END MODULE RD_VAR_MOD
!
	SUBROUTINE RD_OPTIONS_INTO_STORE(LU_IN,LU_OUT)
	USE RD_VAR_MOD
	IMPLICIT NONE
!
! Altered: 10-Apr-2009: Now stop code if find inconsistency in presence of '[''s.
!
	INTEGER LU_IN,LU_OUT
!
	INTEGER I,J
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
	LUO=LU_OUT
	I_UP=1
	I_DWN=1
!
	ALLOCATE (STORE(NST_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (KEY_ST(NST_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (KEY_END(NST_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RD_FILEE_OPTIONS'
	  WRITE(LUER,*)'Unable to allocate memory'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP                                                      
	END IF
	KEY_ST(:)=0
	KEY_END(:)=0
!
! We only store STRINGS containing valid KEYS.
!
	I=1
	DO WHILE(I .LE. NST_MAX)
	  READ(LU_IN,'(A)',END=100)STORE(I)
	  IF(STORE(I)(1:1) .NE. '!' .AND. STORE(I) .NE. ' ')THEN
	    KEY_ST(I)=INDEX(STORE(I),'[')
	    KEY_END(I)=INDEX(STORE(I),']')
	  IF(KEY_ST(I) .NE. 0 .AND. KEY_END(I) .GT. KEY_ST(I)+1)THEN
	      KEY_ST(I)=KEY_ST(I)+1
	      KEY_END(I)=KEY_END(I)-1
	      NST=I
	      I=I+1
	    ELSE IF(KEY_ST(I) .NE. 0 .OR. KEY_END(I) .NE. 0)THEN
	      WRITE(LUER,*)'Possible error in RD_OPTIONS_INTO_STORE'
	      WRITE(LUER,*)'String with inconsistent [ brackets found'
	      WRITE(LUER,*)STORE(I)
	      STOP
	    ELSE	
	      WRITE(LUER,*)'Error in RD_OPTIONS_INTO_STORE'
	      WRITE(LUER,*)'Comments (i.e., strings without [KEY]) must begin with !'
	      WRITE(LUER,*)STORE(I)
	      STOP
	    END IF
	  END IF
	END DO
100	CONTINUE
!
	IF(NST .EQ. NST_MAX)THEN
	  WRITE(LUER,*)'Error in RD_OPTIONS_INTO_STORE'
	  WRITE(LUER,*)'Unable to read in all OPTIONS strings from file.'
	  STOP
	END IF
!
! Check whether all keys are unique.
!
	DO I=1,NST-1
	  DO J=I+1,NST
	    IF( STORE(I)(KEY_ST(I):KEY_END(I)) .EQ.
	1             STORE(J)(KEY_ST(J):KEY_END(J))    )THEN
	      WRITE(LUER,*)'Error in RD_OPTIONS_INTO_STORE'
	      WRITE(LUER,*)'The following key is not unique'
	      WRITE(LUER,*)'Record=',I,'  KEY=',STORE(I)(KEY_ST(I):KEY_END(I)) 
	      WRITE(LUER,*)'Record=',J,'  KEY=',STORE(J)(KEY_ST(J):KEY_END(J)) 
	      STOP
	    END IF
	  END DO
	END DO
!
	RETURN
	END
!
	SUBROUTINE CLEAN_RD_STORE
	USE RD_VAR_MOD
	IMPLICIT NONE
!
	DEALLOCATE (STORE,STAT=IOS)
	IF(IOS .EQ. 0)DEALLOCATE (KEY_ST,STAT=IOS)
	IF(IOS .EQ. 0)DEALLOCATE (KEY_END,STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CLEAN_RD_STORE'
	  WRITE(LUER,*)'Unable to deallocate memory'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
	RETURN
	END
!
	SUBROUTINE GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	USE RD_VAR_MOD
	IMPLICIT NONE
!
! Altered 13-Feb-2011: Bug fixed for second return. Was using J insted of K
!                          to specify string length.
! Altered 27-Apr-2007: STRING that is returned no longer has key and comment
!                         included.
!
	CHARACTER(LEN=*) KEY,STRING
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
!
	INTEGER J,K,ISAVE
!
! We search for the string in the neighborhod of the
! last string.
!
	ISAVE=I_DWN
	KEY_FOUND=.TRUE.
	DO WHILE(I_DWN .GE. 1 .OR. I_UP .LE. NST)
	  J=MIN(I_UP,NST);     K=MAX(I_DWN,1)
	  IF(TRIM(KEY) .EQ. STORE(J)(KEY_ST(J):KEY_END(J)))THEN
	    STRING=STORE(J)(1:KEY_ST(J)-2)
	    I_DWN=I_UP-1
	    I_UP=I_UP+1
	    RETURN
	  ELSE IF(TRIM(KEY) .EQ. STORE(K)(KEY_ST(K):KEY_END(K)))THEN
	    STRING=STORE(K)(1:KEY_ST(K)-2)
	    I_UP=I_DWN+1
	    I_DWN=I_DWN-1
	    RETURN
	  ELSE
	    I_UP=I_UP+1
	    I_DWN=I_DWN-1
	  END IF
	END DO
	I_DWN=ISAVE	!Keep earlier position
	I_UP=ISAVE+1
!
	KEY_FOUND=.FALSE.
	IF(.NOT. MUST_BE_PRES)RETURN
	WRITE(LUER,*)'Error in GET_KEY_STRING'
	WRITE(LUER,*)'Unable to locate string contianing the key ',KEY
	STOP
!
	END
!
	SUBROUTINE RD_STORE_LOG(VALUE,KEY,MUST_BE_PRES,A)
	USE RD_VAR_MOD
	IMPLICIT NONE
	LOGICAL VALUE
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
	CHARACTER(LEN=*) KEY,A
	CHARACTER(LEN=80) STRING
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	READ(STRING,*,IOSTAT=IOS)VALUE
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error reading logical value with RD_STORE_LOG'
	  WRITE(LUER,*)'String with error follows:'
	  WRITE(LUER,*)TRIM(STRING)
	  STOP
	END IF
	STRING=ADJUSTL(STRING)
	IF(STRING(1:1) .NE. 'T' .AND. STRING(1:1) .NE. 'F')THEN
	  WRITE(LUER,*)'Error reading logical value with RD_STORE_LOG'
	  WRITE(LUER,*)'Use F [or FALSE] and T [or TRUE] for logical variables'
	  WRITE(LUER,*)'String with error follows'
	  WRITE(LUER,*)TRIM(STRING)
	  STOP
	END IF
	WRITE(LUO,10)VALUE,TRIM(KEY),TRIM(A)
10	FORMAT(12X,L1,5X,'[',A,']',T40,A)
!
	RETURN
	END
!
!
	SUBROUTINE RD_STORE_2LOG(VALUE1,VALUE2,KEY,MUST_BE_PRES,A)
	USE RD_VAR_MOD
	IMPLICIT NONE
	LOGICAL VALUE1,VALUE2
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
	CHARACTER(LEN=*) KEY,A
	CHARACTER(LEN=80) STRING
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	READ(STRING,*)VALUE1,VALUE2
	WRITE(LUO,10)VALUE1,VALUE2,TRIM(KEY),TRIM(A)
10	FORMAT(10X,L1,',',L1,5X,'[',A,']',T40,A)
!
	RETURN
	END
!
!
	SUBROUTINE RD_STORE_INT(VALUE,KEY,MUST_BE_PRES,A)
	USE RD_VAR_MOD
	IMPLICIT NONE
	INTEGER VALUE
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
	CHARACTER(LEN=*) KEY,A
	CHARACTER(LEN=80) STRING
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	READ(STRING,*,IOSTAT=IOS)VALUE
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error reading integer value with RD_STORE_INT'
	  WRITE(LUER,*)'String with error follows:'
	  WRITE(LUER,*)TRIM(STRING)
	  STOP
	END IF
	WRITE(LUO,10)VALUE,TRIM(KEY),TRIM(A)
10	FORMAT(5X,I8,5X,'[',A,']',T40,A)
!
	RETURN
	END
!
!
!
	SUBROUTINE RD_STORE_2INT(VALUE1,VALUE2,KEY,MUST_BE_PRES,A)
	USE RD_VAR_MOD
	IMPLICIT NONE
	INTEGER VALUE1,VALUE2
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
	INTEGER I
	CHARACTER(LEN=*) KEY,A
	CHARACTER(LEN=80) STRING
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	READ(STRING,*,IOSTAT=IOS)VALUE1,VALUE2
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error reading integer values with RD_STORE_2INT'
	  WRITE(LUER,*)'String with error follows:'
	  WRITE(LUER,*)TRIM(STRING)
	  STOP
	END IF
	WRITE(STRING,'(I8,A,I8)')VALUE1,',',VALUE2
	I=1
	DO WHILE(I.LE. LEN_TRIM(STRING))
	  IF(STRING(I:I) .EQ. ' ')THEN
	    STRING(I:)=STRING(I+1:)
	  ELSE
	     I=I+1
	  END IF
	END DO
	WRITE(LUO,10)TRIM(STRING),TRIM(KEY),TRIM(A)
10	FORMAT(A,T18,'[',A,']',T40,A)
!
	RETURN
	END
!
!
	SUBROUTINE RD_STORE_DBLE(VALUE,KEY,MUST_BE_PRES,A)
	USE RD_VAR_MOD
	IMPLICIT NONE
	REAL*8 VALUE
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
	CHARACTER(LEN=*) KEY,A
	CHARACTER(LEN=80) STRING
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	READ(STRING,*,IOSTAT=IOS)VALUE
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error reading double precision value with RD_STORE_DBLE'
	  WRITE(LUER,*)'String with error follows:'
	  WRITE(LUER,*)TRIM(STRING)
	  STOP
	END IF
	WRITE(LUO,10)VALUE,TRIM(KEY),TRIM(A)
10	FORMAT(1X,1PE12.5,5X,'[',A,']',T40,A)
!
	END
!
! 
!  Now ignores blank lines in input file.
!
	SUBROUTINE RD_STORE_CHAR(VALUE,KEY,MUST_BE_PRES,A)
	USE RD_VAR_MOD
!
! Altered 17-Feb-2009 : VALUE can now be of "arbitratry" length.
!                       (Alteration done much earlier on ROSELLA).
!
	IMPLICIT NONE
	CHARACTER(LEN=*) KEY,A
	CHARACTER(LEN=*) VALUE
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=13) TMP_STR	!Must be 13 for output format
	INTEGER I,J
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	VALUE=ADJUSTL(STRING)
!
! Remove any appended tabs
!
	J=LEN(VALUE)
	DO I=J,1,-1
	  IF(VALUE(I:I) .EQ. CHAR(9))THEN
	    VALUE(I:I)=' '
	  ELSE IF(VALUE(I:I) .NE. ' ')THEN
	    EXIT
	  END IF
	END DO
!
	J=LEN_TRIM(VALUE)
	IF(J .LE. 13)THEN
	  TMP_STR=VALUE
	  TMP_STR=ADJUSTR(TMP_STR)
	  WRITE(LUO,10)TMP_STR,TRIM(KEY),TRIM(A)
10	  FORMAT(A13,5X,'[',A,']',T40,A)
	ELSE IF(J+7+LEN_TRIM(KEY) .LT. 35)THEN
	  WRITE(LUO,11)TRIM(VALUE),TRIM(KEY),TRIM(A)
11	  FORMAT(/,A,5X,'[',A,']',T40,A,/)
	ELSE
	  WRITE(LUO,12)TRIM(VALUE),TRIM(KEY),TRIM(A)
12	  FORMAT(/,A,5X,'[',A,']',10X,A,/)
	END IF
!
	RETURN
	END
!
! 
! To replace RDCHAR. The number of characters to be read in is
! now secified in the call. Now skips blank strings.
!
	SUBROUTINE RD_STORE_NCHAR(VALUE,KEY,NCHAR,MUST_BE_PRES,A)
	USE RD_VAR_MOD
	IMPLICIT NONE
	INTEGER NCHAR
	INTEGER I
	CHARACTER(LEN=*) KEY,A,VALUE
	CHARACTER(LEN=80) STRING
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	VALUE=' '
	VALUE=STRING(1:NCHAR)
	IF(NCHAR .LE. 11)THEN
	  STRING=' '
	  STRING(13-NCHAR:12)=VALUE
	END IF
	I=LEN_TRIM(STRING)
	WRITE(LUO,10)STRING(1:MAX(12,I)),TRIM(KEY),TRIM(A)
10	FORMAT(1X,A,5X,'[',A,']',T40,A)
!
	RETURN
	END
!
!
! Routine designed to count the numer of strings of a ceartain key type.
! For example, this routine can be used to count the total nubmer of ionization 
! stages in the MODEL_SPEC file. To do this, pass SUB_KEY as "_NSF]"
!
! NB: We include the ] in the search string, as it provides a more definitive
!     search.
!
	SUBROUTINE CNT_NUM_KEYS(CNT,SUB_KEY)
	USE RD_VAR_MOD
	IMPLICIT NONE
!
! Altered 22-Jun-2000 : When SUB_KEY='_NSF]' routine checks for possible typo's om
!                         the MODEL_SPEC file.
! Created 07-Jun-2000
!
	INTEGER CNT
	CHARACTER(LEN=*) SUB_KEY
!
	INTEGER I,BEG_CNT,END_CNT
!
! Note that the sub-string defined by KEY_ST, KEY_END do not incorporate 
! the []'s.
!
	CNT=0
	BEG_CNT=0
	END_CNT=-1
	DO I=1,NST
	  IF(INDEX(STORE(I)(KEY_ST(I)-1:KEY_END(I)+1),TRIM(SUB_KEY))
	1                                                         .NE. 0)THEN
	    CNT=CNT+1
	    IF(CNT .EQ. 1)BEG_CNT=I
	    END_CNT=I
	  END IF
	END DO
!
! This provides a method for checking for possible typo's. Particularly important
! for the first or last ioization stages. Note that STORE does not contain
! comment strings.
!
	IF(CNT .NE. END_CNT-BEG_CNT+1 .AND. TRIM(SUB_KEY) .EQ. '_NSF]')THEN
	  WRITE(LUER,*)'Possible error detected'
	  WRITE(LUER,*)'The number of occurences of ',TRIM(SUB_KEY),
	1              ' detected was: ',CNT
	  WRITE(LUER,*)'The number of keys from the first to last occurence was: ',
	1               END_CNT-BEG_CNT+1
          WRITE(LUER,*)'This may indicate a typo in MODEL_SPEC'
          WRITE(LUER,*)'Non NSF Keywords should not be mixed among NSF keyords'
	  STOP
	END IF
!
	RETURN
	END
!
!
	SUBROUTINE RD_STORE_3INT(VALUE1,VALUE2,VALUE3,KEY,MUST_BE_PRES,A)
	USE RD_VAR_MOD
	IMPLICIT NONE
	INTEGER VALUE1,VALUE2,VALUE3
	LOGICAL MUST_BE_PRES
	LOGICAL KEY_FOUND
	INTEGER I
	CHARACTER(LEN=*) KEY,A
	CHARACTER(LEN=80) STRING
!
	CALL GET_KEY_STRING(STRING,KEY,MUST_BE_PRES,KEY_FOUND)
	IF(.NOT. KEY_FOUND)RETURN
	READ(STRING,*)VALUE1,VALUE2,VALUE3
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error reading 3 integer values with RD_STORE_3INT'
	  WRITE(LUER,*)'String with error follows:'
	  WRITE(LUER,*)TRIM(STRING)
	  STOP
	END IF
	WRITE(STRING,'(I8,A,I8,A,I8)')VALUE1,',',VALUE2,',',VALUE3
	I=1
	DO WHILE(I.LE. LEN_TRIM(STRING))
	  IF(STRING(I:I) .EQ. ' ')THEN
	    STRING(I:)=STRING(I+1:)
	  ELSE
	     I=I+1
	  END IF
	END DO
	WRITE(LUO,10)TRIM(STRING),TRIM(KEY),TRIM(A)
10	FORMAT(A,T18,'[',A,']',T40,A)
!
	RETURN
	END
