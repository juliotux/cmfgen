C
C Routine to format single precision numbers to their most compact form.
C Only NUM_SIG_DIG signicant digits are assumed.
C
C (The double precision routine [FORM_DP_NUM] calls this routine).
C
	SUBROUTINE FORM_SP_NUM(VAL,FORM,LENGTH)
	IMPLICIT NONE
C
C Altered 29-Apr-1991 - Variable (i.e. <) formats removed for CRAY 
C                       compatibility.
C Created 10-Dec-1990
C
	INTEGER LENGTH,I,J,POW,NUM_SIG_DIG
	PARAMETER (NUM_SIG_DIG=7)
	REAL*4 VAL
	CHARACTER FORM*(*)
	CHARACTER STR*80
C
	IF(VAL .EQ. 0)THEN
	  FORM='0'
	  LENGTH=1
	  RETURN
	END IF
C
	STR=' '
	LENGTH=6+NUM_SIG_DIG
	J=NUM_SIG_DIG-1
	FORM='(1PE  .  )'
	WRITE(FORM(5:6),'(I2.2)')LENGTH
	WRITE(FORM(8:9),'(I2.2)')J
	WRITE(STR,FORM)VAL
	READ(STR(LENGTH-2:LENGTH),'(I3)')POW  	!Assumes E+XX format
C
C Strip unwanted zero''s (Down to decimal point only)
C
	I=LENGTH-4
	DO WHILE( STR(I:I) .EQ. '0')
	  STR(I:LENGTH-1)=STR(I+1:LENGTH)
	  LENGTH=LENGTH-1
	  I=I-1
	END DO
C
C Now decide whether E or F format.
C
	IF(0 .LE. POW .AND. POW .LT. 4)THEN
	  LENGTH=LENGTH-4
	  DO I=4,3+POW
	    IF(I .GT. LENGTH)THEN
	      STR(I-1:I-1)='0'
	      STR(I:I)='.'
	      LENGTH=LENGTH+1
	    ELSE
	      STR(I-1:I-1)=STR(I:I)
	      STR(I:I)='.'
	    END IF
	  END DO
	  IF(VAL .LT. 0)THEN
	    FORM=STR(1:LENGTH)
	  ELSE
	    FORM=STR(2:LENGTH)
	    LENGTH=LENGTH-1
	  END IF
	  IF(FORM(LENGTH:LENGTH) .EQ. '.')LENGTH=LENGTH-1
	  RETURN
	ELSE IF(POW .LE. -1 .AND. POW .GE. -2)THEN
	  LENGTH=LENGTH-4
	  STR(3:3)=STR(2:2)
	  STR(2:2)='.'
	  DO I=-2,POW,-1
	    STR(3:)='0'//STR(3:LENGTH)
	    LENGTH=LENGTH+1
	  END DO
	  IF(VAL .GT. 0)THEN
	    FORM='0.'//STR(3:LENGTH)
	  ELSE
	    FORM='-0.'//STR(3:LENGTH)
	    LENGTH=LENGTH+1
	  END IF
	  RETURN
	ELSE			!E format
	  IF(VAL .LT. 0)THEN
	    FORM=STR(1:LENGTH)
	  ELSE
	    FORM=STR(2:LENGTH)
	    LENGTH=LENGTH-1
	  END IF
	  RETURN
	END IF
C
	END
	  
C
C Routine to format single precision numbers to their most compact form.
C Only NUM_SIG_DIG signicant digits are assumed.
C
C
	SUBROUTINE FORM_DP_NUM(VAL,FORM,LENGTH)
	IMPLICIT NONE
C
C Altered 08-Jun-2013 - Made identical to FORM_SP_NUM except for VAL declaration.
C                         Can now handle very small numbers.
C
	INTEGER LENGTH,I,J,POW,NUM_SIG_DIG
	PARAMETER (NUM_SIG_DIG=7)
	REAL*8 VAL
	CHARACTER(LEN=*) FORM
	CHARACTER(LEN=80) STR
C
	IF(VAL .EQ. 0)THEN
	  FORM='0'
	  LENGTH=1
	  RETURN
	END IF
C
	STR=' '
	LENGTH=6+NUM_SIG_DIG
	J=NUM_SIG_DIG-1
	FORM='(1PE  .  )'
	WRITE(FORM(5:6),'(I2.2)')LENGTH
	WRITE(FORM(8:9),'(I2.2)')J
	WRITE(STR,FORM)VAL
	READ(STR(LENGTH-2:LENGTH),'(I3)')POW  	!Assumes E+XX format
C
C Strip unwanted zero''s (Down to decimal point only)
C
	I=LENGTH-4
	DO WHILE( STR(I:I) .EQ. '0')
	  STR(I:LENGTH-1)=STR(I+1:LENGTH)
	  LENGTH=LENGTH-1
	  I=I-1
	END DO
C
C Now decide whether E or F format.
C
	IF(0 .LE. POW .AND. POW .LT. 4)THEN
	  LENGTH=LENGTH-4
	  DO I=4,3+POW
	    IF(I .GT. LENGTH)THEN
	      STR(I-1:I-1)='0'
	      STR(I:I)='.'
	      LENGTH=LENGTH+1
	    ELSE
	      STR(I-1:I-1)=STR(I:I)
	      STR(I:I)='.'
	    END IF
	  END DO
	  IF(VAL .LT. 0)THEN
	    FORM=STR(1:LENGTH)
	  ELSE
	    FORM=STR(2:LENGTH)
	    LENGTH=LENGTH-1
	  END IF
	  IF(FORM(LENGTH:LENGTH) .EQ. '.')LENGTH=LENGTH-1
	  RETURN
	ELSE IF(POW .LE. -1 .AND. POW .GE. -2)THEN
	  LENGTH=LENGTH-4
	  STR(3:3)=STR(2:2)
	  STR(2:2)='.'
	  DO I=-2,POW,-1
	    STR(3:)='0'//STR(3:LENGTH)
	    LENGTH=LENGTH+1
	  END DO
	  IF(VAL .GT. 0)THEN
	    FORM='0.'//STR(3:LENGTH)
	  ELSE
	    FORM='-0.'//STR(3:LENGTH)
	    LENGTH=LENGTH+1
	  END IF
	  RETURN
	ELSE			!E format
	  IF(VAL .LT. 0)THEN
	    FORM=STR(1:LENGTH)
	  ELSE
	    FORM=STR(2:LENGTH)
	    LENGTH=LENGTH-1
	  END IF
	  RETURN
	END IF
C
	END
	  

