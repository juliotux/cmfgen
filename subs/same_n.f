	FUNCTION SAME_N(LEV1,LEV2)
	IMPLICIT NONE
	LOGICAL SAME_N
	CHARACTER*(*) LEV1,LEV2
C
	INTEGER STN,ENDN
	CHARACTER N1*3,N2*3,CORE1*30,CORE2*30
C
	CALL POSN_SAME(LEV1,STN,ENDN)
	IF(STN .EQ. 0)THEN
	  SAME_N=.FALSE.
	  RETURN
	ELSE
	  N1=LEV1(STN:ENDN)
	  CORE1=LEV1(1:STN)
	END IF
	CALL POSN_SAME(LEV2,STN,ENDN)
	IF(STN .EQ. 0)THEN
	  SAME_N=.FALSE.
	  RETURN
	ELSE
	  N2=LEV2(STN:ENDN)
	  CORE2=LEV2(1:STN)
	END IF
C
	IF( N1 .EQ. N2 .AND. CORE1 .EQ. CORE2)THEN
	  SAME_N=.TRUE.
	ELSE
	  SAME_N=.FALSE.
	END IF
C
	RETURN
	END
C
C 
C
C The folling subroutine determines the location of the principal quantum
C number. It is meant for names of the general form
C
C 3d2(3Fe)4s4Fe[1/2] where n=4 in this case.
C
C If no lower case letter is found in the name (in this case s), STN and ENDN
C are returned as zero.
C
C Altered 20-Sep-1999 LEV1 was being accessed outside its bounds when determing
C                       position of n for names like 4s3So.
C Altered 26/6/1996   DLM Search only for n (number bewteen 0 and 9) instead of
C                       searching for l (letter between a and z) 3 spaces back
C                       from the end of the configuration or 3 spaces back from
C                       parity (e or o).  This change was needed to find n in Fe2
C                       configurations because of the historic spectroscope degeneracy
C                       letter (a-d,u-z) preceeding the some of its terms.
C
C                      **note** Not able to distinguish between l=d and degeneracy=d
C                               may give incorrect n value.
C
C                      **note** Not able to correct n values for configurations such
C                               as 3d54s2De => returns n=54.  Would have to count
C                               electrons in order to determine the correct n, or
C                               create cnofigurations as 3d45d4s2De...
C
C                      These errors should not be important because this ftn is only
C                      called to determine if core configurations and n's are the
C                      same when calculationg collision strenghts.  The value of n
C                      is not actually used.  These error may cause imperical values
C                      of collision strengths to be `incorrect` by a factor of 3, but
C                      for Fe these approximate strengths may be wrong anyway.
C
	SUBROUTINE POSN_SAME(LEV1,STN,ENDN)
	IMPLICIT NONE
	CHARACTER*(*) LEV1
	INTEGER L1,J,STN,ENDN
C
	INTEGER ERROR_LU,ICHRLEN,LUER
	EXTERNAL ERROR_LU,ICHRLEN
C
C Allow for H and HeI names of the form n___ and nSNG
C
	L1=INDEX(LEV1,'__')		!H, He2 check
	IF(L1 .EQ. 0)THEN		!HeI check
	  L1=MAX(INDEX(LEV1,'SNG'),INDEX(LEV1,'TRP'))
	END IF
	IF(L1 .NE. 0)THEN
	  ENDN=L1-1
	  STN=ENDN
	  IF(STN .GT. 1)THEN
	    IF(LEV1(STN-1:STN-1) .GT. '0' .AND. LEV1(STN-1:STN-1) .LT. '9')STN=STN-1
	  END IF
	  IF(STN .GT. 1)THEN
	    IF(LEV1(STN-1:STN-1) .GT. '0' .AND. LEV1(STN-1:STN-1) .LT. '9')STN=STN-1
	  END IF
	  RETURN
	END IF
C
C Check if j value is indicated in name.
C
	L1=INDEX(LEV1,'[')
	IF(L1 .EQ. 0)THEN
	  L1=ICHRLEN(LEV1)
	ELSE
	  L1=L1-1
	END IF
C
C Check for parity indicator.
C
	IF(LEV1(L1:L1) .EQ. 'e' .or. LEV1(L1:L1) .EQ. 'o')THEN
	  L1=L1-1
	END IF
C
C Start search for number (n) 3 spaces back from end of from parity.
C
	STN=0
	ENDN=0
	J=L1-3
C
C If degeneracy letter is present then must start search back one more space
C
	IF(LEV1(J+1:J+1) .LT. 'd' .OR. LEV1(J+1:J+1) .GT. 's')J=J-1
C
	DO WHILE(STN .EQ. 0 .AND. J .GT. 0)
C
C If core configuration is at J then move back to space preceeding '('
C
	  IF(LEV1(J:J) .EQ. ')')THEN
	    DO WHILE(LEV1(J:J) .NE. '(')
	      J=J-1
	    END DO
            J=J-1
	  END IF
C
	  IF(LEV1(J:J) .GE. '0' .AND. LEV1(J:J) .LE. '9')THEN
	    IF((LEV1(J+1:J+1).LT.'d').OR.(LEV1(J+1:J+1).GT.'s'))goto200
	    ENDN=J
	    STN=ENDN
	    IF(STN .GT. 1)THEN
	      IF(LEV1(STN-1:STN-1) .GE. '0' .AND.
	1           LEV1(STN-1:STN-1) .LE. '9')STN=STN-1
	    END IF
	    IF(STN .GT. 1)THEN
	      IF(LEV1(STN-1:STN-1) .GE. '0' .AND.
	1          LEV1(STN-1:STN-1) .LE. '9')STN=STN-1
	    END IF
            IF(LEV1(ENDN+1:ENDN+1) .LT. 'd' .OR.
	1	 LEV1(ENDN+1:ENDN+1) .GT. 's')THEN
              LUER=ERROR_LU()
              WRITE(LUER,*)'Error getting n from name -',LEV1,STN,ENDN
              STOP
	    ENDIF
	    RETURN
 200	  END IF
	  J=J-1
	END DO
C
	RETURN
C
	END
