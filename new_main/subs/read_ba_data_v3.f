!
! Subroutine to store entie BA variation file on disk. The BA data in
! this version is passed via a data module. Unlike earlier routines,
! STEQ is NOT written out. 
!
	SUBROUTINE READ_BA_DATA_V3(LU,NION,NUM_BNDS,ND,COMPUTE_BA,FIXED_T,STATUS,DESC)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 04-Oct-2016 - Added FIXED_T to call.
! Altered 12-Mar-2014 - COMPUTE_BA is nolonger changed by the routine. The value in CMFGEN 
!                         is used to indicate whether the BA matric will be read in.
! Altered 07-Mar-2004 - BA no longer read in if COMPUTE_BA=.TRUE.
! Created 02-Apr-2001 - Created to handle SE data structure.
!                         See READBA for earlier corrections.
!
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER ND
        INTEGER LU                    !Input unit for BA and STEQ
	LOGICAL COMPUTE_BA  		!Indicates whether BA is being computed.
	LOGICAL FIXED_T                 !Indicates whether T is to be held fixed.
	LOGICAL STATUS                  !Indicates whether BA/STEQ read successful
	CHARACTER DESC*(*)              !Used for filename
!
! Local Variables and external functions.
!
	INTEGER NUM_BNDS_RD,ND_RD,NION_RD
	INTEGER ID
	INTEGER LUER,ERROR_LU,IOS
	EXTERNAL ERROR_LU
	LOGICAL FILE_OPEN
	LOGICAL COMPUTE_BA_RD
	INTEGER, PARAMETER :: IZERO=0
!
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
	  READ(LU,*,ERR=400,IOSTAT=IOS)COMPUTE_BA_RD
	  READ(LU,*,ERR=400,IOSTAT=IOS)FIXED_T
	  READ(LU,*,ERR=400,IOSTAT=IOS)NION_RD
	  READ(LU,*,ERR=400,IOSTAT=IOS)NUM_BNDS_RD
	  READ(LU,*,ERR=400,IOSTAT=IOS)ND_RD
	CLOSE(UNIT=LU)
!
	IF(NION_RD .NE. NION .OR. NUM_BNDS_RD .NE. NUM_BNDS .OR. ND_RD .NE. ND)THEN
	  WRITE(LUER,*)'Error : incompatible dimensions in BAREAD'
	  WRITE(LUER,*)'Skipping READ of BA and STEQ'
	  WRITE(LUER,*)'  Read vales:',NION_RD,NUM_BNDS_RD,ND_RD
	  WRITE(LUER,*)'Actual vales:',NION,NUM_BNDS,ND
	  CLOSE(UNIT=LU)
          STATUS=.FALSE.
	  RETURN
	END IF
!
! If we are still computing the BA matrix, there is no need to read it in.
! 
	IF(COMPUTE_BA)RETURN
!
	OPEN(UNIT=LU,FORM='UNFORMATTED',FILE=DESC,IOSTAT=IOS,ERR=500,
	1             ACCESS='SEQUENTIAL',STATUS='OLD',ACTION='READ')
	  DO ID=1,NION
	    READ(LU,ERR=600,IOSTAT=IOS)SE(ID)%BA
	  END DO
	  READ(LU,ERR=605,IOSTAT=IOS)BA_ED
	  READ(LU,ERR=610,IOSTAT=IOS)BA_T
	CLOSE(UNIT=LU)
	RETURN
!
300	WRITE(LUER,*)'Error opening '//DESC//'PNT in READBA'
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
!
400	WRITE(LUER,*)'Error reading from '//DESC//'PNT in READBA'
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
!
500	WRITE(LUER,*)'Error opening logical unit to recall :',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
!
600	WRITE(LUER,*)'Error on reading SE(ID)%BA in READ_BA_DATA_V2: ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
!
605	WRITE(LUER,*)'Error on reading BA_ED in READ_BA_DATA_V2: ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
!
610	WRITE(LUER,*)'Error on reading BA_T in READ_BA_DATA_V2: ',DESC
        WRITE(LUER,*)'IOSTAT=',IOS
	STATUS=.FALSE.
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	RETURN
!
	RETURN
	END
