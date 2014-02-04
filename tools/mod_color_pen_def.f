        MODULE MOD_COLOR_PEN_DEF
!
! Define some standard color pens. Use the DEF_PEN to return terminal
! to a black pen.
!
	CHARACTER(LEN=5), PARAMETER :: DEF_PEN=achar(27)//'[0m'
	CHARACTER(LEN=5), PARAMETER :: BLACK_PEN=achar(27)//'[0m'
	CHARACTER(LEN=5), PARAMETER :: RED_PEN=achar(27)//'[31m'
	CHARACTER(LEN=5), PARAMETER :: GREEN_PEN=achar(27)//'[32m'
	CHARACTER(LEN=5), PARAMETER :: BLUE_PEN=achar(27)//'[34m'
	CHARACTER(LEN=5), PARAMETER :: YELLOW_PEN=achar(27)//'[33m'
	CHARACTER(LEN=5), PARAMETER :: MAGNETA_PEN=achar(27)//'[35m'
	CHARACTER(LEN=5), PARAMETER :: CYAN_PEN=achar(27)//'[36m'
!
! Define pens to match default pens in pgplot. When the pen < 100, we use the null character
! to ensure that an extra space is not output.
!
	INTEGER, PARAMETER :: NUM_PG_PEN=12
	CHARACTER(LEN=11), SAVE :: PG_PEN(0:12)=(/
	1 achar(27)//'[38;5;243m',
	1 achar(0)//achar(27)//'[38;5;00m',
	1 achar(27)//'[38;5;196m',
	1 achar(0)//achar(27)//'[38;5;21m',
	1 achar(0)//achar(27)//'[38;5;28m',
	1 achar(0)//achar(27)//'[38;5;91m',   !56m',
	1 achar(27)//'[38;5;165m',
	1 achar(27)//'[38;5;220m',
	1 achar(27)//'[38;5;208m',
	1 achar(0)//achar(27)//'[38;5;82m',
	1 achar(0)//achar(27)//'[38;5;46m',      !154m',
	1 achar(0)//achar(27)//'[38;5;33m',
	1 achar(27)//'[38;5;129m'/)
!
	END MODULE MOD_COLOR_PEN_DEF
