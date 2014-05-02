!
	SUBROUTINE SCR_READ_V2(R,V,SIGMA,POPS,
	1             IREC_RD,NITSF,NUM_TIMES,LST_NG,RVSIG_WRITTEN,
	1             NT,ND,LU,NEWMOD)
	IMPLICIT NONE
!
! Altered 04-Nov-2013 : Bug fix -- incorrect read when R, V and SIGMA written, and not
!                          reading the last record.
! Altered 28-Sep-2009 : Changed to allow an increase in the number of depth points.
!                          Previously the R,V, & SIGMA record was limiting ND.
! Altered 02-May-2004 : Changed to V2.
!                          RVSIG_WRITTEN inserted into call.
!                          Only writes R, V & SIGMA when needed, but all
!                          writes to file must be identical.
! Altered 29-Feb-2004 - Adjusted to read new format file for SCRTEMP.
!                          R, V, & SIGMA can be output for every iteration.
!                          IREC_RD must be zero to read last iteration,
!                          other wise iteration is IREC.
!                          NUM_TIMES in no longer used.
! Altered 05-Dec-1996 - INQUIRE statement used before closing files
!                       ASCI files opened by GEN_ASCI_OPEN.
! Altered 25-Jun-1996 - ACTION='READ' installed on OPEN statements.
! Altered 12-Jan-1991 - By using call to DIR_ACC_PARS this version is now
!                        compatible with both CRAY and VAX fortran.
! Altered  3-Apr-1989 - LST_NG installed. LU now transmitted in call.
!
	LOGICAL NEWMOD
!
!                          On input,  IREC_RD is the iteration to be read.
!                          On output, IREC_RD is the actual iteration read.
!
	INTEGER IREC_RD
	INTEGER NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
	LOGICAL RVSIG_WRITTEN
!
! Used to read R, V & SIGMA (allows for 3*ND to be > N_PER_REC.
!
	REAL*8 RVSIG_VEC(3*ND)
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
	INTEGER NUMRECS			!Number of records needed to output POPS
	INTEGER NUM_RV_RECS		!Number of records needed to output R,V & SIGMA.
	INTEGER N_PER_REC		!Maximum number per record.
	INTEGER ARRAYSIZE  		!Size of POPS array.
!	
	INTEGER LOC_NUMTIMES,ST_REC_M1,RECS_FOR_RV
	INTEGER I,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER IREC
	LOGICAL FILE_OPEN
	INTEGER, PARAMETER :: IZERO=0
	EXTERNAL ERROR_LU
	CHARACTER*80 STRING
	CHARACTER*11 FORMAT_DATE
!
	NEWMOD=.FALSE.
	LUER=ERROR_LU()
	FORMAT_DATE=' '
!
! Determine the record size, and the number of records that
! need to be written out to fully write out the population vector.
! These are computer dependent, hence call to DIR_ACC_PARS. NB.
! REC_SIZE is not the same as REC_LEN --- it is REC_LEN which is
! passed to the OPEN statement.
 
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	ARRAYSIZE=NT*ND
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
!
! Read in pointer to data file. If pointer file does not exist, or bad
! data, initialize parameters for a NEW MODEL.
!
	CALL GEN_ASCI_OPEN(LU,'POINT1','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening POINT1 in SCR_READ'
	  ELSE
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
              READ(LU,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG,RVSIG_WRITTEN
	    ELSE IF(IOS .EQ. 0)THEN
	      RVSIG_WRITTEN=.FALSE.
              READ(STRING,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG
	    END IF
	    IF(IOS .NE. 0)WRITE(LUER,*)'Error reading POINT1 in SCR_READ'
	  END IF
!
	  IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    CALL GEN_ASCI_OPEN(LU,'POINT2','OLD',' ','READ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening POINT2 in SCR_READ'
	    ELSE
	      READ(LU,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	        STRING=ADJUSTL(STRING)
	        FORMAT_DATE=STRING(1:11)
                READ(LU,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG,RVSIG_WRITTEN
	      ELSE IF(IOS .EQ. 0)THEN
	        RVSIG_WRITTEN=.FALSE.
                READ(STRING,*,IOSTAT=IOS)IREC,NITSF,LOC_NUMTIMES,LST_NG
	      END IF
	    END IF
	  END IF
!
	IF(IOS .NE. 0 .OR. IREC .LT. 1)THEN
	  WRITE(LUER,*)'Error on reading Pointers in READ_SCRTEMP'
	  NEWMOD=.TRUE.
	  NITSF=0
	  LST_NG=-1000
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  RETURN
	END IF
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
!		' OLD MODEL '
!
	OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',
	1       RECL=REC_LEN,IOSTAT=IOS,ACTION='READ')
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP for input'
	    WRITE(LUER,*)'IOSTAT=',IOS
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)CLOSE(UNIT=LU)
	    NEWMOD=.TRUE.
	    NITSF=0
	    LST_NG=-1000
	    RETURN
	  END IF
!
! Note that NITSF= # of successful iterations so far.
! Note: R,V and SIGMA will be updated if they are output for every iteration.
!
	  NUM_RV_RECS=INT( (3*ND-1)/N_PER_REC ) + 1
	  IF(NUM_RV_RECS .EQ. 1)THEN
	    READ(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	    IF(IOS .NE. 0)READ(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	    RECS_FOR_RV=2
	  ELSE
	    DO L=1,NUM_RV_RECS
	      IST=(L-1)*N_PER_REC+1
	      IEND=MIN(IST+N_PER_REC-1,3*ND)
	      READ(LU,REC=L,IOSTAT=IOS)(RVSIG_VEC(I),I=IST,IEND)
	    END DO
	    R(1:ND)=RVSIG_VEC(1:ND)
	    V(1:ND)=RVSIG_VEC(ND+1:2*ND)
	    SIGMA(1:ND)=RVSIG_VEC(2*ND+1:3*ND)
	    RECS_FOR_RV=NUM_RV_RECS
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading R,V, SIGMA vectors in READ_SCRTEMP'
	    NEWMOD=.TRUE.
	    NITSF=0
	    LST_NG=-1000
	    IREC_RD=0
	    RETURN
	  END IF
!
! IREC ignores the number of records that it takes to write each time.
! Hence in POINT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
	  IF(IREC_RD .NE. 0)THEN
	    IREC=IREC_RD
	  ELSE
	    IREC_RD=IREC
	  END IF
!
! Read in the population data.
!
500	CONTINUE		!Try to read an earlier record.
	IF(RVSIG_WRITTEN)THEN
	  IF(NUM_RV_RECS .EQ. 1)THEN
	    ST_REC_M1=(IREC-1)*(NUMRECS+1)+RECS_FOR_RV
	    READ(LU,REC=ST_REC_M1+1,IOSTAT=IOS)R,V,SIGMA
	    ST_REC_M1=ST_REC_M1+1
	  ELSE
	    ST_REC_M1=(IREC-1)*(NUMRECS+NUM_RV_RECS)+RECS_FOR_RV
	    DO L=1,NUM_RV_RECS
	      IST=(L-1)*N_PER_REC+1
	      IEND=MIN(IST+N_PER_REC-1,3*ND)
	      READ(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(RVSIG_VEC(I),I=IST,IEND)
	    END DO
	    R(1:ND)=RVSIG_VEC(1:ND)
	    V(1:ND)=RVSIG_VEC(ND+1:2*ND)
	    SIGMA(1:ND)=RVSIG_VEC(2*ND+1:3*ND)
	    ST_REC_M1=ST_REC_M1+NUM_RV_RECS
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
	  WRITE(LUER,*)'Error on Scratch Read'
	  IREC=IREC-1
	  NEWMOD=.TRUE.
	  NITSF=0
	  LST_NG=-1000
	  IREC_RD=0
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	  RETURN
	END IF
!
! Successful Read !
!
	CLOSE(UNIT=LU)
	RETURN
	END
!
! 
!
! Routine  to save population data. There is no limit on the
! size of the POPS array.
!
	SUBROUTINE SCR_RITE_V2(R,V,SIGMA,POPS,
	1               IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG,
	1               NT,ND,LU,WRITFAIL)
	IMPLICIT NONE
!
! Altered 28-Sep-2009 : Changed to allow an increase in the number of depth points.
!                          Previously the R,V, & SIGMA record was limiting ND.
! Altered 02-May-2004 : Changed to V2.
!                          RVSIG_WRITTEN inserted into call.
!                          Only writes R, V & SIGMA when needed, but all
!                          writes to file must be identical.
! Altered 29-Feb-2004 - Adjusted to read new format file for SCRTEMP.
!                          R, V, & SIGMA now output for every iteration.
!                          NUM_TIMES in no longer used.
! Altered 27-Feb-2004 : R,V, and SIGMA rewritten on every iteration.
! Altered 25-Jun-1996 : GEN_ASCI_OPEN installed to OPEN POINT files.
!
	LOGICAL WRITFAIL
	LOGICAL WRITE_RVSIG
	INTEGER IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
!
! Used to read R, V & SIGMA (allows for 3*ND to be > N_PER_REC.
!
	REAL*8 RVSIG_VEC(3*ND)
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
	INTEGER NUMRECS			!Number of records needed to output POPS
	INTEGER NUM_RV_RECS		!Number of records needed to output R,V & SIGMA.
	INTEGER N_PER_REC		!Maximum number per record.
	INTEGER ARRAYSIZE  		!Size of POPS array.
!	
	INTEGER ST_REC_M1,RECS_FOR_RV
	INTEGER I,K,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL FILE_OPEN
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
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
        ARRAYSIZE=NT*ND
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
!
! We check FORMAT associated with SCRATCH file. This will allow us to
! preserve the same format for an existing file.
!
	IF(IREC .NE. 1)THEN
!
! Writing to existing SCRATCH file.
!
	  FORMAT_DATE=' '
	  CALL GEN_ASCI_OPEN(LU,'POINT1','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)
	1     CALL GEN_ASCI_OPEN(LU,'POINT2','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
	      READ(LU,*)I,I,I,I,LOC_WR_RVSIG
	      IF(LOC_WR_RVSIG .NEQV. WRITE_RVSIG)THEN
	        WRITE(LUER,*)'Error in SCR_RITE_V2 -- inconsistent WR_RVSIG option'
	        WRITE(LUER,*)'Restart a fresh model by deleting SCRTEMP etc.'
	        STOP
	      END IF
	    END IF
	  END IF
	  IF(IOS .NE. 0)FORMAT_DATE='28-Feb-2004'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	ELSE
!
! Create a new scratch file.
!
	  FORMAT_DATE='28-Feb-2004'
	END IF
!
!*************************************************************************
!
	OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED'
	1,  ACCESS='DIRECT',STATUS='UNKNOWN'
	1,  RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP in WRITE_SCRTEMP'
	    WRITE(LUER,*)'Will try to open a new file'
	    OPEN(UNIT=LU,FILE='SCRTEMP',FORM='UNFORMATTED',
	1     ACCESS='DIRECT',STATUS='NEW',
	1     RECL=REC_LEN,IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening SCRTEMP for output'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      WRITFAIL=.TRUE.
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(UNIT=LU)
	      RETURN
	    ELSE
	      FORMAT_DATE='28-Feb-2004'
	      IREC=0				!Since new file.
	    END IF
	  END IF
!
! We now write out R,V and SIGMA on every iteration. This is also 
! done for SN models where R can change with iteration. R read here will
! only be correct for the last iteration. Fixed with a later read.
!
! We keep two writes for R,V & SIGMA when a single record to preserve
! compatibility woth older versions.
!
	  NUM_RV_RECS=INT( (3*ND-1)/N_PER_REC ) + 1
	  IF(NUM_RV_RECS .EQ. 1)THEN
	    WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	    IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	    RECS_FOR_RV=2
	  ELSE
            RVSIG_VEC(1:ND)=R(1:ND)
            RVSIG_VEC(ND+1:2*ND)=V(1:ND)
            RVSIG_VEC(2*ND+1:3*ND)=SIGMA(1:ND)
            DO L=1,NUM_RV_RECS
              IST=(L-1)*N_PER_REC+1
              IEND=MIN(IST+N_PER_REC-1,3*ND)
              WRITE(LU,REC=L,IOSTAT=IOS)(RVSIG_VEC(I),I=IST,IEND)
            END DO
            RECS_FOR_RV=NUM_RV_RECS
	END IF
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error writing R,V etc vectors in SCR_RITE'
	  WRITFAIL=.TRUE.
	  RETURN
	END IF
!
! WRITE in the population data.
!
	IREC=IREC+1		!Next record output
!
! IREC ignores the number of records that it takes to write each time.
! Hence in POINT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
! NB: We cannot change the type of file that is being written. If we do
!     SCRTEMP will be corrupted.
!
	IF(WRITE_RVSIG)THEN
	  ST_REC_M1=(IREC-1)*(NUMRECS+NUM_RV_RECS)+RECS_FOR_RV
          RVSIG_VEC(1:ND)=R(1:ND)
          RVSIG_VEC(ND+1:2*ND)=V(1:ND)
          RVSIG_VEC(2*ND+1:3*ND)=SIGMA(1:ND)
          DO L=1,NUM_RV_RECS
            IST=(L-1)*N_PER_REC+1
            IEND=MIN(IST+N_PER_REC-1,3*ND)
            WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(RVSIG_VEC(I),I=IST,IEND)
          END DO
	  ST_REC_M1=ST_REC_M1+NUM_RV_RECS
	ELSE
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	END IF
!
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
	  WRITE(LUER,*)'Error writing SCRTEMP in SCR_RITE'
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
! Write pointer to data files.
!
	CALL GEN_ASCI_OPEN(LU,'POINT1','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
	CALL GEN_ASCI_OPEN(LU,'POINT2','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
	RETURN
	END
!
!
! Routine  to save population data. There is no limit on the
! size of the POPS array. This useses the prefix NEW_ attached
! to POINT and SCRTEMP. Otherwise it is identical to SCR_RITE_V2
!
	SUBROUTINE SCR_RITE_NAM_V2(R,V,SIGMA,POPS,
	1               IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG,
	1               NT,ND,LU,WRITFAIL)
	IMPLICIT NONE
!
! Altered 28-Sep-2009 : Changed to allow an increase in the number of depth points.
!                          Previously the R,V, & SIGMA record was limiting ND.
! Altered 02-May-2004 : Changed to V2.
!                          RVSIG_WRITTEN inserted into call.
!                          Only writes R, V & SIGMA when needed, but all
!                          writes to file must be identical.
! Altered 29-Feb-2004 - Adjusted to read new format file for NEW_SCRTEMP.
!                          R, V, & SIGMA now output for every iteration.
!                          NUM_TIMES in no longer used.
! Altered 27-Feb-2004 : R,V, and SIGMA rewritten on every iteration.
! Altered 25-Jun-1996 : GEN_ASCI_OPEN installed to OPEN POINT files.
!
	LOGICAL WRITFAIL
	LOGICAL WRITE_RVSIG
	INTEGER IREC,NITSF,NT,ND,NUM_TIMES,LST_NG,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),POPS(NT*ND)
!
! Used to read R, V & SIGMA (allows for 3*ND to be > N_PER_REC.
!
	REAL*8 RVSIG_VEC(3*ND)
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
	INTEGER NUMRECS			!Number of records needed to output POPS
	INTEGER NUM_RV_RECS		!Number of records needed to output R,V & SIGMA.
	INTEGER N_PER_REC		!Maximum number per record.
	INTEGER ARRAYSIZE  		!Size of POPS array.
!	
	INTEGER ST_REC_M1,RECS_FOR_RV
	INTEGER I,K,L,LUER,ERROR_LU
	INTEGER IOS,IST,IEND
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL FILE_OPEN
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
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
	ARRAYSIZE=NT*ND
	NUMRECS=INT( (ARRAYSIZE-1)/N_PER_REC )+1
	REC_LEN=REC_SIZE/UNIT_SIZE
	WRITE(6,*)NUMRECS,ARRAYSIZE,N_PER_REC
!
! We check FORMAT associated with SCRATCH file. This will allow us to
! preserve the same format for an existing file.
!
	IF(IREC .NE. 1)THEN
!
! Writing to existing SCRATCH file.
!
	  FORMAT_DATE=' '
	  CALL GEN_ASCI_OPEN(LU,'NEW_POINT1','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)
	1     CALL GEN_ASCI_OPEN(LU,'NEW_POINT2','OLD',' ','READ',IZERO,IOS)
	  IF(IOS .EQ. 0)THEN
	    READ(LU,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .EQ. 0 .AND. INDEX(STRING,'!Format date') .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE=STRING(1:11)
	      READ(LU,*)I,I,I,I,LOC_WR_RVSIG
	      IF(LOC_WR_RVSIG .NEQV. WRITE_RVSIG)THEN
	        WRITE(LUER,*)'Error in SCR_RITE_NAM_V2 -- inconsistent WR_RVSIG option'
	        WRITE(LUER,*)'Restart a fresh model by deleting NEW_SCRTEMP etc.'
	        STOP
	      END IF
	    END IF
	  END IF
	  IF(IOS .NE. 0)FORMAT_DATE='28-Feb-2004'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)CLOSE(UNIT=LU)
	ELSE
!
! Create a new scratch file.
!
	  FORMAT_DATE='28-Feb-2004'
	END IF
!
!*************************************************************************
!
	OPEN(UNIT=LU,FILE='NEW_SCRTEMP',FORM='UNFORMATTED'
	1,  ACCESS='DIRECT',STATUS='UNKNOWN'
	1,  RECL=REC_LEN,IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening NEW_SCRTEMP in WRITE_SCRTEMP'
	    WRITE(LUER,*)'Will try to open a new file'
	    OPEN(UNIT=LU,FILE='NEW_SCRTEMP',FORM='UNFORMATTED',
	1     ACCESS='DIRECT',STATUS='NEW',
	1     RECL=REC_LEN,IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening NEW_SCRTEMP for output'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      WRITFAIL=.TRUE.
	      INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	      IF(FILE_OPEN)CLOSE(UNIT=LU)
	      RETURN
	    ELSE
	      FORMAT_DATE='28-Feb-2004'
	      IREC=0				!Since new file.
	    END IF
	  END IF
!
! We now write out R,V and SIGMA on every iteration. This is also 
! done for SN models where R can change with iteration. R read here will
! only be correct for the last iteration. Fixed with a later read.
!
! We keep two writes for R,V & SIGMA when a single record to preserve
! compatibility woth older versions.
!
	  NUM_RV_RECS=INT( (3*ND-1)/N_PER_REC ) + 1
	  IF(NUM_RV_RECS .EQ. 1)THEN
	    WRITE(LU,REC=1,IOSTAT=IOS)R,V,SIGMA
	    IF(IOS .EQ. 0)WRITE(LU,REC=2,IOSTAT=IOS)R,V,SIGMA
	    RECS_FOR_RV=2
	  ELSE
            RVSIG_VEC(1:ND)=R(1:ND)
            RVSIG_VEC(ND+1:2*ND)=V(1:ND)
            RVSIG_VEC(2*ND+1:3*ND)=SIGMA(1:ND)
            DO L=1,NUM_RV_RECS
              IST=(L-1)*N_PER_REC+1
              IEND=MIN(IST+N_PER_REC-1,3*ND)
              WRITE(LU,REC=L,IOSTAT=IOS)(RVSIG_VEC(I),I=IST,IEND)
            END DO
            RECS_FOR_RV=NUM_RV_RECS
	END IF
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error writing R,V etc vectors in SCR_RITE'
	  WRITFAIL=.TRUE.
	  RETURN
	END IF
	WRITE(6,*)'Read R, V & SIGMA'
!
! WRITE in the population data.
!
	IREC=IREC+1		!Next record output
!
! IREC ignores the number of records that it takes to write each time.
! Hence in NEW_POINT it will correspond to the iteration number.
! ST_REC_M1 + 1 is the first output record.
!
! NB: We cannot change the type of file that is being written. If we do
!     NEW_SCRTEMP will be corrupted.
!
	IF(WRITE_RVSIG)THEN
	  ST_REC_M1=(IREC-1)*(NUMRECS+NUM_RV_RECS)+RECS_FOR_RV
          RVSIG_VEC(1:ND)=R(1:ND)
          RVSIG_VEC(ND+1:2*ND)=V(1:ND)
          RVSIG_VEC(2*ND+1:3*ND)=SIGMA(1:ND)
          DO L=1,NUM_RV_RECS
            IST=(L-1)*N_PER_REC+1
            IEND=MIN(IST+N_PER_REC-1,3*ND)
            WRITE(LU,REC=ST_REC_M1+L,IOSTAT=IOS)(RVSIG_VEC(I),I=IST,IEND)
          END DO
	  ST_REC_M1=ST_REC_M1+NUM_RV_RECS
	ELSE
	  ST_REC_M1=(IREC-1)*NUMRECS+RECS_FOR_RV
	END IF
	WRITE(6,*)'Update record',IOS
	WRITE(6,*)'RECS_FOR_RV: ',RECS_FOR_RV
	WRITE(6,*)'NUM_RV_RECS: ',NUM_RV_RECS
	WRITE(6,*)'NUMRECS: ',NUMRECS
!
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
	  WRITE(LUER,*)'Error writing NEW_SCRTEMP in SCR_RITE'
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
! Write pointer to data files.
!
	CALL GEN_ASCI_OPEN(LU,'NEW_POINT1','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in NEW_POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
	CALL GEN_ASCI_OPEN(LU,'NEW_POINT2','UNKNOWN',' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening in NEW_POINT1 in SCR_RITE'
	    WRITE(LUER,'(1X,4(I6,4X))')IREC,NITSF,NUM_TIMES,LST_NG
	  END IF
	  WRITE(LU,'(X,A,50X,A)')FORMAT_DATE,'!Format date'
	  WRITE(LU,'(1X,4(I6,4X),2X,L1)')IREC,NITSF,NUM_TIMES,LST_NG,WRITE_RVSIG
	  WRITE(LU,'(3X,A,5X,A,3X,A,4X,A,4X,A)')'IREC','NITSF','#_TIMES','LST_NG','WR_RVS'
	  INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
	RETURN
	END
