C
C Function to return the length of a string. Blanks at end of
C string are not counted. The string is not modified.
C
	FUNCTION ICHRLEN(A)
	IMPLICIT NONE
	INTEGER ICHRLEN,I,L
	CHARACTER*(*) A
C
	L=LEN(A)
	DO I=L,1,-1
	  IF(A(I:I) .NE. ' ')THEN
	    ICHRLEN=I
	    RETURN
	  END IF
	END DO
C
	ICHRLEN=0
	RETURN
	END
