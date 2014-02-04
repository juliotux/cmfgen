!
! Simple subroutine to output a colored string to an XTERM.
! There are basically 7 colors, black, red, green, yellow, magenta, and cyan,
! which are represented by the numbers 0 through 6.
!
! To switch to a color, enter |#. Routine defaults to black after ouput.
!
! If you need |# in string as a non-control, add a space before the  #.
!
	SUBROUTINE WR_COL_STR(STR_FOR_OUTPUT)
	IMPLICIT NONE
!
! Created: 2-Nov-2011
!
	CHARACTER(LEN=*) STR_FOR_OUTPUT
!
	INTEGER I,J
	LOGICAL SET_COL
	CHARACTER(LEN=2), PARAMETER :: BEG_COL=achar(27)//'['
	CHARACTER(LEN=4), PARAMETER :: END_COL=achar(27)//'[0m'
	CHARACTER(LEN=200) STRING
!
! Intensity	0	1	2	3	4	5	6	7
!           Black	Red	Green	Yellow	Blue	Magenta	Cyan	White
!
	SET_COL=.FALSE.
	STRING=STR_FOR_OUTPUT
	I=INDEX(STRING,'|')
	DO WHILE(I .NE. 0)
	  IF(STRING(I+1:I+1) .GE. '0' .AND. STRING(I+1:I+1) .LE. '7')THEN
	    STRING(I:)=BEG_COL//'3'//STRING(I+1:I+1)//'m'//STRING(I+2:)
	    SET_COL=.TRUE.
	    I=INDEX(STRING,'|')
	  ELSE
	    J=INDEX(STRING(I+1:I+1),'|')
	    IF(J .EQ. 0)EXIT
	    I=I+J
	  END IF
	END DO
	IF(SET_COL)STRING=TRIM(STRING)//END_COL
	WRITE(6,'(A)')TRIM(STRING)
!
	RETURN
	END
