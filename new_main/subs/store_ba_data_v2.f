C
C Program to output BA matrices to scratch files. When a CMFGEN model
C is nearing convergence, the BA matrix does not need to be recomputed.
C Instead it is read in from disk, when needed.
C
	SUBROUTINE STORE_BA_DATA_V2(LU,NION,NUM_BANDS,ND,COMPUTE_BA,DESC)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Created 05-Apr-2001 - Based on STOREBA
!                       Designed to use TSEQ_DATAMOD.
!                       See REABA for earlier changes.
!
	INTEGER NION  	!Number of ions treated
	INTEGER NUM_BANDS  	!Number of bands in matrix
	INTEGER ND  		!Number of depth points
	INTEGER LU		!Logical unit for BA output.
        INTEGER LU1		!LU+1 is used for BAPNT1
C
C COMPUTE_BA is used to indicate whether BA been computed. This must be
C included since it changes during program execution, and hence is not
C the same as in VADAT file.
C
	LOGICAL COMPUTE_BA
	CHARACTER DESC*(*)		!File name for BA &STEQ output.
C
	LOGICAL FILE_OPEN
	INTEGER IOS,ID,LUER,ERROR_LU
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
	  DO ID=1,NION
	    WRITE(LU,ERR=600,IOSTAT=IOS)SE(ID)%BA
	  END DO
	  WRITE(LU,ERR=600,IOSTAT=IOS)BA_ED
	  WRITE(LU,ERR=600,IOSTAT=IOS)BA_T
	CLOSE(UNIT=LU)
C
C Output to BAPNT file that write of BA and STEQ was successful.
C
	ENTRY INIT_BA_DATA_PNT_V2(LU,NION,NUM_BANDS,ND,COMPUTE_BA,DESC)
	  LU1=LU+1
	  LUER=ERROR_LU()
          INQUIRE(UNIT=LU1,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)THEN
	    BACKSPACE(LU1)
	  ELSE
	    CALL GEN_ASCI_OPEN(LU1,DESC//'PNT','UNKNOWN',' ',' ',IZERO,IOS)
	    IF(IOS .NE. 0)GOTO 400
	  END IF
	  WRITE(LU1,'(10X,L1,T20,A)',ERR=800,IOSTAT=IOS)
	1            .TRUE.,'BA structure successfully output'
	  WRITE(LU1,'(10X,L1,T20,A)',ERR=800,IOSTAT=IOS)
	1            COMPUTE_BA,'BA is currently being computed'
	  WRITE(LU1,'(1X,I10,T20,A)',ERR=800,IOSTAT=IOS)
	1            NION,'Total # of ionization stages'
	  WRITE(LU1,'(1X,I10,T20,A)',ERR=800,IOSTAT=IOS)
	1            NUM_BANDS,'# of bands'
	  WRITE(LU1,'(1X,I10,T20,A)',ERR=800,IOSTAT=IOS)
	1            ND,'# of depth points'
!
! Ouut dimesnions of each structure.
!
	  DO ID=1,NION
	    WRITE(LU1,*)SE(ID)%N_SE,SE(ID)%N_IV
	  END DO
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
