C
C Simple VMS routine to determine whether a particular Logical unit
C (LU) is attached to a terminal or not. Uses old SYS$TRNLOG instead
C of SYS$TRNLNM as simpler. Machine dependent, and crude.
C
	LOGICAL FUNCTION TERMINAL_LU(LU)
C
C Altered 20-Dec-1990
C Created 10-Oct-1990
C
	IMPLICIT NONE
C	INCLUDE '($SSDEF)'
!	INTEGER SYS$TRNLOG
C
	CHARACTER NAME*10,OUT_NAME*80
	INTEGER LU,LNAME,ISTAT
C
	NAME='FOR'
	OUT_NAME=' '
	WRITE(NAME(4:6),'(I3.3)')LU
	TERMINAL_LU=.TRUE.
	RETURN
!	
!	ISTAT=SYS$TRNLOG(NAME(1:6),LNAME,OUT_NAME,,,)
C
C In first case, no logical name translation, and the UNIT is 5 or 6.
C 5 and 6 and the default interactive LU's for FORTRAN.
C In the second case, a particular unit has been assign to TT:
C (i.e. terminal).
C
	IF(NAME(1:6) .EQ. OUT_NAME .AND. (LU .EQ. 5 .OR. LU .EQ. 6) )THEN
	  TERMINAL_LU=.TRUE.
	ELSE IF(OUT_NAME(1:3) .EQ. 'TT:')THEN
	  TERMINAL_LU=.TRUE.
	ELSE
	  TERMINAL_LU=.FALSE.
	END IF
C
	RETURN
	END
