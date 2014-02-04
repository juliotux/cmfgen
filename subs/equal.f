C
C Logical function to determine whether two values are equal to within
C 100Z % . If one of the arguments are zero, EQUAL is set false unless
C both are equal to zero in which case it is set true. Neither X, Y or
C Z are altered.
C
	FUNCTION EQUAL(X,Y,Z)
	IMPLICIT NONE
C
C Altered 24-May-1996 - File now contains DP version only (i.e. not SP_EQUAL)
C Altered 19-Jul-1991 - Rearrangement of LOGICAL descriptor for CRAY.
C Altered 14-Apr-1989 - Now divide by the larger (absolute) of X and Y.
C                       This routine should never give a floating point
C                       overflow.
C Altered  4-NOV-86 (Bug for X or Y=0 fixed).
C
	LOGICAL EQUAL
	REAL*8 X,Y,Z
C
	EQUAL=.FALSE.
	IF(X .EQ. 0.0D0 .AND. Y .EQ. 0.0D0)THEN
	  EQUAL=.TRUE.
	ELSE IF( ABS(Y) .GT. ABS(X) )THEN
	  IF( (1.0D0-X/Y) .LE. Z )EQUAL=.TRUE.
	ELSE
	  IF( (1.0D0-Y/X) .LE. Z )EQUAL=.TRUE.
	END IF
C
	RETURN
	END
