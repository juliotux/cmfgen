!
! Subroutine to READ (and WRITE) a scratch file containg models at
! diffferent time steps. Format of file is identical to SCRTEMP.
! As R is changing, RVSIG_WRITTEN should be true.
!
	SUBROUTINE READ_TIME_MODEL_V1(R,V,SIGMA,POPS,
	1             IREC_RD,RVSIG_WRITTEN,NT,ND,LU)
	IMPLICIT NONE
!
! Created 12-Dec-2005 : Based on READ_TIME_MODEL_V1_V2
!
	INTEGER ND
	INTEGER NT
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
	REAL*8 POPS(NT*ND)
!
!  On input,  IREC_RD is the time step record to be read. If zero, the last
!                time step is read.
!  On output, IREC_RD is the actual record read.
!
	INTEGER IREC_RD
	INTEGER LU
	LOGICAL RVSIG_WRITTEN
!
! Local variables.
!
! REC_SIZE     is the (maximum) record length in bytes.
! REC_LEN      is the record length in computer units.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! NUMRECS      is the # of records required to output POPS.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER REC_SIZE,REC_LEN
	INTEGER NUMRECS
	INTEGER N_PER_REC
	INTEGER ARRAYSIZE  	!Size of POPS array.
!	
	INTEGER LOC_NUMTIMES,ST_REC_M1,RECS_FOR_RV
	INTEGER I,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER IREC
	INTEGER NITSF
	LOGICAL FILE_OPEN
	INTEGER, PARAMETER :: IZERO=0
	EXTERNAL ERROR_LU
	CHARACTER*80 STRING
	CHARACTER*11 FORMAT_DATE
!
	LUER=ERROR_LU()
	FORMAT_DATE=' '
!
! Determine the record size, and the number of records that
! need to be written out to fully write out the population vector.
! These are computer dependent, hence call to DIR_ACC_PARS. NB.
! REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
! passed to the OPEN statement.
 
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Record length was too small fro R,V and'//
	1              ' SIGMA in READ_TIME_MODEL_V1'
	  WRITE(LUER,*)'3ND=',3*ND
	  WRITE(LUER,*)'N_PER_REC=',N_PER_REC
	  STOP
	END IF
	ARRAYSIZE=NT*ND
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
!
! Read in pointer to data file. If pointer file does not exist, or bad
! data, initialize parameters for a NEW MODEL.
!
	CALL GEN_ASCI_OPEN(LU,'TIME_PNT1','OLD',' ','READ',IZERO,IOS)
	  IREC=0; NITSF=0
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening TIME_PNT1 in READ_TIME_MODEL_V1'
	  ELSE
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
              READ(LU,*,IOSTAT=IOS)IREC,NITSF,I,I,RVSIG_WRITTEN
	      IF(IREC .NE. NITSF)THEN
	        WRITE(LUER,*)'Possible error in READ_TIME_MODEL_V1'
	        WRITE(LUER,*)'IREC & NITSF should match in TIME_PNT1'
	        STOP
	      END IF
	    ELSE
	      WRITE(LUER,*)'Unrecognized formate date in TIME_PNT1'
	      STOP
	    END IF
	    IF(IOS .NE. 0)WRITE(LUER,*)'Error reading TIME_PNT1 in READ_TIME_MODEL_V1'
	  END IF
!
	  IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	    IREC=0; NITSF=0
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    CALL GEN_ASCI_OPEN(LU,'TIME_PNT2','OLD',' ','READ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening TIME_PNT2 in READ_TIME_MODEL_V1'
	    ELSE
	      READ(LU,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	        STRING=ADJUSTL(STRING)
	        FORMAT_DATE=STRING(1:11)
                READ(LU,*,IOSTAT=IOS)IREC,NITSF,I,I,RVSIG_WRITTEN
	        IF(IREC .NE. NITSF)THEN
	          WRITE(LUER,*)'Possible error in READ_TIME_MODEL_V1'
	          WRITE(LUER,*)'IREC & NITSF should match in TIME_PNT1'
	          STOP
	        END IF
	      ELSE 
	        WRITE(LUER,*)'Unrecognized formate date in TIME_PNT2'
	        STOP
	      END IF
	    END IF
	  END IF
!
	IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	  WRITE(LUER,*)'Error on reading Pointers in READ_TIME_MODEL'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  STOP
	END IF
!
	IF(IREC_RD .GT. IREC)THEN
	  WRITE(LUER,*)'Error: record specified for reading is too large'
	  WRITE(LUER,*)'Error on reading Pointers in READ_TIME_MODEL'
	  WRITE(LUER,*)'Maximum IREC=',IREC,'Record to be read=',IREC_RD
	  STOP
	END IF
!
! Open the direct access file with the model data.
! 
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
	OPEN(UNIT=LU,FILE='MODELS_FN_TIME',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',
	1       RECL=REC_LEN,IOSTAT=IOS,ACTION='READ')
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening MODEL_FN_TIME for input'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    RETURN
	  END IF
!
! Note that NITSF= # of successful iterations so far.
!
	  READ(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)READ(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading R,V, SIGMA vectors in READ_TIME_MODEL'
	    STOP
	  END IF
	  RECS_FOR_RV=2
!
! Read in the population data.
!
500	CONTINUE		!Try to read an earlier record.
!
! IREC ignores the number of records that it takes to write each time.
! Hence in TIME_PNT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
	IF(IREC_RD .NE. 0)THEN
	  IREC=IREC_RD
	ELSE
	  IREC_RD=IREC
	END IF
	IF(RVSIG_WRITTEN)THEN
	  ST_REC_M1=(IREC-1)*(NUMRECS+1)+RECS_FOR_RV
	  READ(LU,REC=ST_REC_M1+1,IOSTAT=IOS)R,V,SIGMA
	  ST_REC_M1=ST_REC_M1+1
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading iteratuion specific R,V, SIGMA'//
	1           ' vectors in READ_TIME_MODEL'
	    WRITE(LUER,*)'IREC=',IREC
	    STOP
	  END IF
	ELSE
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	END IF
!
	IF(IOS .EQ. 0)THEN
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    READ(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)EXIT
	  END DO
	END IF
!
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error on reading MODELS_FN_TIME in READ_TIME_MODEL'
	  WRITE(LUER,*)'IREC=',IREC
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  STOP
	END IF
!
! Successfull Read !
!
	CLOSE(UNIT=LU)
	RETURN
	END
!
! 
!
! Routine  to save population data for a given time step. The model time
! is an integer counter begining at 1. There is no limit on the
! size of the POPS array.
!
! Entry:
!   IREC=0 :     Next record to be written is read from TIME_PNT1
!   IREC>0:      Indicates that record to be output is IREC
!
! Exit:
!   IREC:        Last record written.
!
	SUBROUTINE RITE_TIME_MODEL(R,V,SIGMA,POPS,IREC,NEW_FILE,WRITE_RVSIG,
	1               NT,ND,LU,WRITFAIL)
	IMPLICIT NONE
!
! Created 12-Dec-2005 : Based on SCR_RITE.
!
	INTEGER ND
	INTEGER NT
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
	REAL*8 POPS(NT*ND)
	INTEGER IREC
	INTEGER LU
	LOGICAL NEW_FILE 
	LOGICAL WRITFAIL
	LOGICAL WRITE_RVSIG
!
! Local variables.
!
! REC_SIZE     is the (maximum) record length in bytes.
! REC_LEN      is the record length in computer units.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! NUMRECS      is the # of records required to output POPS.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER REC_SIZE
	INTEGER REC_LEN
	INTEGER NUMRECS
	INTEGER N_PER_REC
	INTEGER ARRAYSIZE  	!Size of POPS array.
!	
	INTEGER ST_REC_M1,RECS_FOR_RV
	INTEGER I,K,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL FILE_OPEN
	LOGICAL FILE_PRES
	LOGICAL LOC_WR_RVSIG
	EXTERNAL ERROR_LU
	CHARACTER(LEN=11) FORMAT_DATE
	CHARACTER*80 STRING
!
	LUER=ERROR_LU()
!
! Determine the record size, and the number of records that
! need to be written out to fully write out the population vector.
! These are computer dependent, hence call to DIR_ACC_PARS. NB.
! REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
! passed to the OPEN statement.
!
	WRITE(LUER,*)'Entered RITE_TIME'
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Record length is too small to output R,V and'//
	1              ' sigma in RITE_TIME_MODEL'
	  WRITE(LUER,*)'3ND=',3*ND
	  WRITE(LUER,*)'N_PER_REC=',N_PER_REC
	  STOP
	END IF
	ARRAYSIZE=NT*ND			
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
!
! We check FORMAT associated with SCRATCH file. This will allow us to
! preserve the same format for an existing file.
!
	INQUIRE(FILE='MODELS_FN_TIME',EXIST=FILE_PRES)
	IF(NEW_FILE .OR. .NOT. FILE_PRES)THEN
!
! We will create a new scratch file.
!
	  FORMAT_DATE='28-Feb-2004'
	  IREC=0
	ELSE 
!
! Writing to existing SCRATCH file.
!
	  FORMAT_DATE=' '
	  CALL GEN_ASCI_OPEN(LU,'TIME_PNT1','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)
	1     CALL GEN_ASCI_OPEN(LU,'TIME_PNT2','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
	      IF(IREC .EQ. 0)THEN
	        READ(LU,*)IREC,I,I,I,LOC_WR_RVSIG
	        IREC=IREC+1				!Next output record
	      ELSE
	        READ(LU,*)I,I,I,I,LOC_WR_RVSIG
	      END IF
	      IF(LOC_WR_RVSIG .NEQV. WRITE_RVSIG)THEN
	        WRITE(LUER,*)'Error in RITE_TIME_MODEL -- inconsistent WR_RVSIG option'
	        WRITE(LUER,*)'Restart a fresh model by deleting MODELS_FN_TIME etc.'
	        STOP
	      END IF
	    END IF
	  ELSE
	    WRITE(LUER,*)'Error opening TIME_PNT* in RITE_TIME_MODEL'
	    WRITE(LUER,*)'MODELS_FN_TIME exists: check consistency'
	    STOP
	  END IF
	  IF(IOS .NE. 0)FORMAT_DATE='28-Feb-2004'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	END IF
!
!*************************************************************************
!
	IOS=0
	OPEN(UNIT=LU,FILE='MODELS_FN_TIME',FORM='UNFORMATTED'
	1,       ACCESS='DIRECT',STATUS='UNKNOWN',RECL=REC_LEN,IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODELS_FN_TIME in WRITE_MODELS_FN_TIME'
	  WRITE(LUER,*)'Will try to open a new file'
	END IF
	IF(NEW_FIlE .OR. IOS .NE. 0 .OR. IREC .EQ. 0)THEN
	  OPEN(UNIT=LU,FILE='MODELS_FN_TIME',FORM='UNFORMATTED',
	1     ACCESS='DIRECT',STATUS='NEW',
	1     RECL=REC_LEN,IOSTAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error opening MODELS_FN_TIME for output'
	        WRITE(LUER,*)'IOSTAT=',IOS
	        WRITFAIL=.TRUE.
	        INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	        IF(FILE_OPEN)CLOSE(UNIT=LU)
	        RETURN
	      ELSE
	        FORMAT_DATE='28-Feb-2004'
	        IREC=1				!Since new file.
	      END IF
	  END IF
!
! We can now write out R,V and SIGMA on every iteration. This is also 
! done for SN models where R can change with iteration. R read here will
! only be correct for the last iteration. Fixed with a later read.
! We preserve this format for consistency with SCRTEMP.
!
	WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error writing R,V etc vectors in RITE_TIME_MODEL'
	  WRITFAIL=.TRUE.
	  RETURN
	END IF
	RECS_FOR_RV=2
!
! IREC ignores the number of records that it takes to write each time.
! Hence in TIME_PNT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
! NB: We cannot change the type of file that is being written. If we do
!     MODELS_FN_TIME will be corrupted.
!
	IF(WRITE_RVSIG)THEN
	  ST_REC_M1=(IREC-1)*(NUMRECS+1)+RECS_FOR_RV
	  WRITE(LU,REC=ST_REC_M1+1,IOSTAT=IOS)R,V,SIGMA
	  ST_REC_M1=ST_REC_M1+1
	ELSE
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	END IF
	IF(IOS .EQ. 0)THEN
	  DO L=1,NUMRECS
	    IST=(L-1)*N_PER_REC+1
	    IEND=MIN(IST+N_PER_REC-1,ARRAYSIZE)
	    WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(POPS(I),I=IST,IEND)
	    IF(IOS .NE. 0)EXIT
	  END DO
	END IF
!
	IF(IOS .NE. 0)THEN
	  WRITFAIL=.TRUE.
	  WRITE(LUER,*)'Error writing MODELS_FN_TIME in RITE_TIME_MODEL'
	  WRITE(LUER,*)'IOS=',IOS
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  RETURN
	END IF
!
! Successful write.
!
	WRITFAIL=.FALSE.
1000	CLOSE(UNIT=LU)
!
! Write pointer files. The format is identical to POINT1 etc, but only IREC
! and WRITE_RVSIG are used. NITSF should be the same as IREC.
!
	CALL GEN_ASCI_OPEN(LU,'TIME_PNT1','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in TIME_PNT1 in RITE_TIME_MODEL'
	    WRITE(LUER,'(1X,I6)')IREC
	  END IF
	  WRITE(LU,'(1X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,IREC,1,-1000,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
	CALL GEN_ASCI_OPEN(LU,'TIME_PNT2','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in TIME_PNT1 in RITE_TIME_MODEL'
	    WRITE(LUER,'(1X,I6)')IREC
	  END IF
	  WRITE(LU,'(1X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,IREC,1,-1000,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
	RETURN
	END
