	SUBROUTINE READBA(BA,STEQ,LU,NT,NS,NLBEGIN,
	1                   COMPUTE_BA,STATUS,DESC)
	IMPLICIT NONE
C
C Altered 13-Dec-1996 - GEN_ASCI_OPEN used to open BA pointer files.
C                       INQUIRE used befor closing files.
C Altered 26-May-1996 - ACTION=READ installed in OPEN statements.
C Altered 10-Jan-1991 - Pointer file now closed on error.
C Altered 20-Dec-1991 - Extensive rewrite : Cray compatible and simpler.
C Altered 20-Sep-1989 - Filname descriptor is now passed to routine. This
C                       will allow this routine to write out both BA and
C                       BASOL. DESC should only be the first part of the
C                       filname as PNT is appended.
C
	INTEGER NT   			!Dimension of BA array
	INTEGER NS  			!Dimension of STEQ array
	INTEGER NLBEGIN		!Lower level (i.e. NL) for
                                        !which line calculations are complete
        INTEGER LU                    !Input unit for BA and STEQ
	LOGICAL COMPUTE_BA  		!Indicates whether BA is being computed.
	LOGICAL STATUS                !Indicates whether BA/STEQ read successful
	CHARACTER DESC*(*)              !Used for filename
	REAL*8 BA(NT),STEQ(NS)
C
C Local Variables and external functions.
C
	INTEGER LUER,ERROR_LU,IOS,NT1,NS1
	EXTERNAL ERROR_LU
	LOGICAL FILE_OPEN
	INTEGER, PARAMETER :: IZERO=0
C
	LUER=ERROR_LU()
	CALL GEN_ASCI_OPEN(LU,DESC//'PNT','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)GOTO 300
	  READ(LU,*,ERR=400,IOSTAT=IOS)STATUS
	  IF(.NOT. STATUS)THEN
	    WRITE(LUER,*)'Previous store of '//DESC,
	1                ' was not completed succesfully'
	    CLOSE(UNIT=LU)
	    RETURN
	  END IF
	  READ(LU,*,ERR=400,IOSTAT=IOS)NLBEGIN
	  READ(LU,*,ERR=400,IOSTAT=IOS)COMPUTE_BA
	  READ(LU,*,ERR=400,IOSTAT=IOS)NT1
	  READ(LU,*,ERR=400,IOSTAT=IOS)NS1
	CLOSE(UNIT=LU)
C
	IF(NT1 .NE. NT .OR. NS1 .NE. NS)THEN
	  WRITE(LUER,*)'Error : incompatible dimensions in BAREAD'
	  WRITE(LUER,*)'Skipping READ of BA and STEQ'
	  CLOSE(UNIT=LU)
          STATUS=.FALSE.
	  RETURN
	END IF
C
	OPEN(UNIT=LU,FORM='UNFORMATTED',FILE=DESC,IOSTAT=IOS,ERR=500,
	1             ACCESS='SEQUENTIAL',STATUS='OLD',ACTION='READ')
	   READ(LU,ERR=600,IOSTAT=IOS)BA
	   READ(LU,ERR=700,IOSTAT=IOS)STEQ
	CLOSE(UNIT=LU)
	RETURN
C
300	WRITE(LUER,*)'Error opening '//DESC//'PNT in READBA'
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
C
400	WRITE(LUER,*)'Error reading from '//DESC//'PNT in READBA'
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
C
500	WRITE(LUER,*)'Error opening logical unit to recall :',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
c
600	WRITE(LUER,*)'Error on reading BA matrix : ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
C
700	WRITE(LUER,*)'Error on reading STEQ matrix : ',DESC
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
	END
