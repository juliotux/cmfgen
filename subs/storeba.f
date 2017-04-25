C
C Program to output BA and STEQ matrices to scratch file so program
C can be retsarted in event of a crash.
C
	SUBROUTINE STOREBA(BA,STEQ,LU,NT,NS,NL,COMPUTE_BA,DESC)
	IMPLICIT NONE
C
C Altered 18-Oct-2000 - LU1 was not being set on ENTRY call to INIT_BA_PNT
C Altered 13-Dec-1996 - GEN_ASI_OPEN installed. INQUIRE used befor
C                         file closings.
C Altered 14-Jan-1991 - Final Cray compatible version.
C Altered 13-Sep-1989 - IOSTAT specifier placed in write statements to
C                       assist in bug finding,
C                       Minor bug fix - STEQ can now have more than
C                       4000 elemnts (eg > 60* 66).
C                       IMPLICIT NONE installed.
C
C
	INTEGER NT  			!Dimension of BA matrix
	INTEGER NS  			!Dimension of STEQ matrix
        INTEGER NL			!Indicates lower levels that have
                                        !been computed successfully so far.
	INTEGER LU			!Logical unit for BA output.
        INTEGER LU1			!LU+1 is used for BAPNT1
C
C COMPUTE_BA is used to indicate whether BA been computed. This must be
C included since it changes during program execution, and hence is not
C the same as in VADAT file.
C
	LOGICAL COMPUTE_BA
	REAL*8 BA(NT),STEQ(NS)
	CHARACTER DESC*(*)		!File name for BA &STEQ output.
C
	LOGICAL FILE_OPEN
	INTEGER IOS,LUER,ERROR_LU
	INTEGER, PARAMETER :: IZERO=0
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
C
C Indicate in pointer file that BA and STEQ are currently being written.
C
	LU1=LU+1
	CALL GEN_ASCI_OPEN(LU1,DESC//'PNT','UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)GOTO 300
	WRITE(LU1,*,ERR=400,IOSTAT=IOS).FALSE.,'!Bad Write'
C
C Store BA and STEQ Matrices.
C
	OPEN(UNIT=LU,FORM='UNFORMATTED',FILE=DESC,
	1     ACCESS='SEQUENTIAL',STATUS='UNKNOWN',IOSTAT=IOS,ERR=500)
	  WRITE(LU,ERR=600,IOSTAT=IOS)BA
	  WRITE(LU,ERR=700,IOSTAT=IOS)STEQ
	CLOSE(UNIT=LU)
C
C Output to BAPNT file that write of BA and STEQ was successful.
C
	ENTRY INIT_BA_PNT(LU,NT,NS,NL,COMPUTE_BA,DESC)
!
	  LUER=ERROR_LU()		!Must be defined for this entry point.
	  LU1=LU+1
          INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)THEN
	    BACKSPACE(LU1)
	  ELSE
	    CALL GEN_ASCI_OPEN(LU1,DESC//'PNT','UNKNOWN',' ',' ',IZERO,IOS)
	    IF(IOS .NE. 0)GOTO 400
	  END IF
	  WRITE(LU1,'(10X,L1,T20,A)',ERR=800,IOSTAT=IOS)
	1            .TRUE.,'BA and STEQ successfully output'
	  WRITE(LU1,'(1X,I10,T20,A)',ERR=800,IOSTAT=IOS)
	1            NL,'Lower line level coputed successfully (i.e. NL)'
	  WRITE(LU1,'(10X,L1,T20,A)',ERR=800,IOSTAT=IOS)
	1            COMPUTE_BA,'BA is currently being computed'
	  WRITE(LU1,'(1X,I10,T20,A)',ERR=800,IOSTAT=IOS)
	1            NT,'# of elements in BA array'
	  WRITE(LU1,'(1X,I10,T20,A)',ERR=800,IOSTAT=IOS)
	1            NS,'# of elements in STEQ array'
	CLOSE(UNIT=LU1)
	RETURN
C
300	WRITE(LUER,*)'Error opening BA pnter file : ',DESC,'PNT'
        WRITE(LUER,*)'IOSTAT=',IOS
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU1)
	RETURN
400	WRITE(LUER,*)'Error opening BA pnter file : ',DESC,'PNT'
        WRITE(LUER,*)'IOSTAT=',IOS
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU1)
	RETURN
500	WRITE(LUER,*)'Error outputing .FALSE. to : ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU1)
	RETURN
600	WRITE(LUER,*)'Error on outputing BA : ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU1)
	RETURN
700	WRITE(LUER,*)'Error on outputing STEQ : ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU1)
	RETURN
800	WRITE(LUER,*)'Error on finalizing pointer file : ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU1)
	RETURN
C
	END
