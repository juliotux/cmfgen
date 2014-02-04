!
! Subroutine to perform an NG accleration for a Comoving-Frame Model. Progam
! uses the last 4 iterations which are stored in the last "4 records" 
! (effectively) of SCRTEMP.
!
! The NG acceleration is perfomed separately on each depth, or over a band
! of depths.
!
! Output is to the last record of SCRTEMP.
!
! Input files required:
!                    POINT1(.DAT)
!                    SCRTEMP(.DAT)
! Output files:
!                    POINT1(.DAT)
!                    POINT2(.DAT)
!                    SCRTEMP(.DAT)
!
	SUBROUTINE DO_NG_BAND_ACCEL_V1(POPS,NT,ND,NG_BAND,NG_DONE,
	1             MAXDEC,MAXINC,LUSCR,LUER)
	IMPLICIT NONE
!
	INTEGER NT
	INTEGER ND
	INTEGER NG_BAND
	INTEGER LUSCR
	INTEGER LUER
	REAL*8 POPS(NT,ND)
!
	REAL*8 MAXDEC
	REAL*8 MAXINC
	LOGICAL NG_DONE
!
! Local arrays
!
	REAL*8 RDPOPS(NT,ND,4)
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
!
! Local variables which are adjusted to match the particular model under
! consideration.
!
	INTEGER IOS
	INTEGER IREC
	INTEGER NITSF
	INTEGER LST_NG
	INTEGER IFLAG
!
	INTEGER, PARAMETER :: RITE_N_TIMES=1
!
	LOGICAL NEWMOD
!
	IFLAG=0
!
! Read POPULATIONS that were output on last iteration. This is primarily
! done to get NITSF etc.
! 
	NEWMOD=.FALSE.
	CALL SCR_READ(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                NT,ND,LUSCR,NEWMOD)
!
! Read in the last 4 estimates of the poplations, as output to SCRTEMP.
!
	CALL RD_4_ITS_V2(RDPOPS,NT,ND,LUER,IFLAG)
	IF(IFLAG .NE. 0)THEN
	  WRITE(LUER,*)'Unable to read scratch file'
	  NG_DONE=.FALSE.
	  RETURN
	END IF
!
! Now perform the NG acceleration. The accelerated estimates for the populations
! are returned in POPS.
!
	CALL NG_MIT_BAND_OPT_V1(POPS,RDPOPS,ND,NT,NG_BAND,NG_DONE,MAXDEC,MAXINC,LUER)
!
	RETURN	
	END
!
!
!
	SUBROUTINE RD_4_ITS_V2(RDPOPS,NT,ND,LUER,IFLAG)
	IMPLICIT NONE
!
	INTEGER NT
	INTEGER ND
	INTEGER LUER
	INTEGER IFLAG
	REAL*8 RDPOPS(NT*ND,4)
!
! Note that REC_SIZE is the size of the output record in bytes.
!           REC_LEN  is the size of the output record in COMPUTER units.
!           N_PER_REC is the number of numbers per record.
!           NUM_TIME is the number of times each iteration has been written
!                to the scratch file.
!           WORD_SIZE is the size of the number to be output in bytes.
!           UNIT_SIZE is the number of bytes per unit used to specify
!           the size of a direct access file.
!
	INTEGER ARRAYSIZE,IST,IEND
	INTEGER REC_SIZE,REC_LEN,NUM_RECS,N_PER_REC,NUM_TIME
	INTEGER UNIT_SIZE,WORD_SIZE
!
	INTEGER I,L,NPREV,NITSF,IT_CNT
	INTEGER SRT_REC_M1
	INTEGER CHK
!
! Determine the record size, and the number of records that
! need to be written out to fully write out the population vector.
! As this is computer and installation dependent, we call a subroutine
! to return the parameters.
!
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,
	1                 WORD_SIZE,N_PER_REC)
	ARRAYSIZE=NT*ND
	NUM_RECS=INT( (ARRAYSIZE-1)/N_PER_REC ) + 1
	REC_LEN=REC_SIZE/UNIT_SIZE
!
! Read in pointer to data file. NB. Format has changed. IT_CNT is now
! the iteration count (but may differ from NITSF is problems occurred
! read/riting to SCRTEMP) where as previously it referred to the actual 
! record.
!
	OPEN(UNIT=26,FILE='POINT1',IOSTAT=CHK,STATUS='OLD')
	  IF(CHK .NE. 0)THEN
	    WRITE(LUER,*)'Error opening POINT1 in GENACCEL: IOSTAT=',CHK
	    IFLAG=3
	    RETURN
	  END IF
	  READ(26,*,IOSTAT=CHK)IT_CNT,NITSF,NUM_TIME
	  IF(CHK .NE. 0)THEN
	    WRITE(LUER,*)'Error reading POINT1 in GENACCEL: IOSTAT=',CHK
	    IFLAG=4
	    CLOSE(UNIT=26)
	    RETURN
	  END IF
	CLOSE(UNIT=26)
!
!		' OLD MODEL '
!
	OPEN(UNIT=20,FILE='SCRTEMP',FORM='UNFORMATTED',
	1    ACCESS='DIRECT',STATUS='OLD',RECL=REC_LEN,IOSTAT=CHK)
	  IF(CHK .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP for input in GENACCEL'
	    WRITE(LUER,*)'IOSTAT=',CHK
	    CLOSE(UNIT=20)
	    IFLAG=5
	    RETURN
	  END IF
!
! SRT_REC_M1 + 1 is the first record for the NPREV previous iteration.
! If output more than once, it refers to the first output.
!
	  DO NPREV=1,4
	    SRT_REC_M1=(IT_CNT-NPREV)*NUM_TIME*NUM_RECS+2
	    DO L=1,NUM_RECS
	      IST=(L-1)*N_PER_REC+1
	      IEND=MIN( IST+N_PER_REC-1,ARRAYSIZE )
	      READ(20,REC=SRT_REC_M1+L,IOSTAT=CHK)
	1          (RDPOPS(I,NPREV),I=IST,IEND)
	      IF(CHK .NE. 0)THEN
	        WRITE(LUER,*)'Error on reading scratch file in GENACCEL'
	        WRITE(LUER,200)SRT_REC_M1,CHK
200	        FORMAT(X,'SRT_REC_M1=',I3,5X,'IOSAT=',I7)
	        IFLAG=6
	        CLOSE(UNIT=20)
	        RETURN
	      END IF
	    END DO
	  END DO
	  CLOSE(UNIT=20)
!
	WRITE(6,*)'SCRTEMP file was successfully read'
	RETURN
	END
!
!
!
! Routine is deisgned to allow NG acceleration over NBAND depths
! simultaneously. For use with CMFGEN.
!
	SUBROUTINE NG_MIT_BAND_OPT_V1(POPS,RDPOPS,ND,NT,NBAND,
	1                NG_DONE,MAXDEC,MAXINC,LUER)
	IMPLICIT NONE
!
! Altered 01-Jun-2003: NUM_BAD_NG now correctly initialized to zero.
!
	INTEGER ND
	INTEGER NT
	INTEGER NBAND
	INTEGER LUER
	REAL*8 POPS(NT,ND)              !Initial/corrected populations (I/O) 
	REAL*8 RDPOPS(NT,ND,4)		!Populations to be accelerated
!
! Return maximum % increase and decrease, and whether the NG acceleration
! was succesful.
!
	REAL*8 MAXINC
	REAL*8 MAXDEC
	LOGICAL NG_DONE
!
! Local arrays
!
	REAL*8 NEWPOP(NT,ND)
	REAL*8 VEC_INC(ND)
	REAL*8 VEC_DEC(ND)
	REAL*8 INT_ARRAY(ND)
!
! Local variables
!
	REAL*8 T1
	REAL*8 LOCINC,LOCDEC
!
	INTEGER NS
	INTEGER NUM_BAD_NG
	INTEGER DEC_LOC,INC_LOC
	INTEGER I,K,L
	INTEGER LST,LEND
	INTEGER LOC_NBAND
	INTEGER, PARAMETER :: IONE=1
!
	NUM_BAD_NG=0
	IF(NBAND .LE. 0 .OR. NBAND .GE. ND)THEN
          LOC_NBAND=ND
	  NS=NT*ND
	  CALL DO_TRANS_AND_NG_V1(NEWPOP,RDPOPS,IONE,ND,NT,ND,NS)
	ELSE
	  DO K=ND,1,-NBAND
	    LST=MAX(K-NBAND+1,1)
	    LEND=LST+NBAND-1
	    NS=(LEND-LST+1)*NT
	    CALL DO_TRANS_AND_NG_V1(NEWPOP(1,LST),RDPOPS,LST,LEND,NT,ND,NS)
	  END DO
          LOC_NBAND=NBAND
	END IF
!
! Now check whether the NG acceleation has been reasonable.
! If it has, store the new estimates in POPS.
!
	MAXINC=-1000.0
	MAXDEC=1000.0
	DO L=1,ND 
	  LOCINC=-1000.0
	  LOCDEC=1000.0
	  DO K=1,NT
	    T1=NEWPOP(K,L)/POPS(K,L)
	    LOCINC=MAX(LOCINC,T1)
	    LOCDEC=MIN(LOCDEC,T1)
	  END DO
	  VEC_INC(L)=100.0D0*(LOCINC-1.0D0)
	  VEC_DEC(L)=100.0D0*(1.0D0/LOCDEC-1.0D0)
!
! Before storing the NG acceleration at this depth, we check to
! see whether the predicted corrections are "reasonable".
!
	  IF(LOCINC .GT. 10.1D0 .OR. LOCDEC .LT. 0.09D0)THEN
	    NUM_BAD_NG=NUM_BAD_NG+1
	    WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	    WRITE(LUER,9000)L,LOCINC,LOCDEC
9000	    FORMAT(X,'No NG acceleration performed at depth ',I3,/,
	1          X,'Biggest increase was ',1PE10.2,/,
	1          X,'Biggest decrease was ',E10.2)
	
	  ELSE
	    DO K=1,NT
	      POPS(K,L)=NEWPOP(K,L)
	    END DO
	    IF(LOCINC .GT. MAXINC)THEN
	      MAXINC=LOCINC
	      INC_LOC=L
	    END IF
	    IF(LOCDEC .LT. MAXDEC)THEN
	      MAXDEC=LOCDEC
	      DEC_LOC=L
	    END IF
	  END IF
	END DO
!
	WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	IF(NUM_BAD_NG .GT. 3)THEN
	  WRITE(LUER,*)'Too many bad NG accelerations - '//
	1              'NG acceleration cancelled'
	  NG_DONE=.FALSE. 
!
! Restore old population estimates for all depths.
!
	  DO L=1,ND
	    DO K=1,NT
	      POPS(K,L)=RDPOPS(K,L,1)
	    END DO
	  END DO
	  RETURN
	END IF
!
! By definition MAXINC and MININC are both positive. These are only
! defined for successful NG accelerations (does not include any depths
! at which NG acceleration did not work).
!
	MAXINC=100.0D0*(MAXINC-1.0D0)
	MAXDEC=100.0D0*(1.0D0/MAXDEC-1.0D0)
	WRITE(LUER,'(A,I3)')' NG Acceleration performed using NG_BAND=',LOC_NBAND
	WRITE(LUER,'(A,I3,A,ES10.2)')
	1  ' Max NG % increase at depth ',INC_LOC,' is',MAXINC
	WRITE(LUER,'(A,I3,A,ES10.2)')
	1  ' Max NG % decrease at depth ',DEC_LOC,' is',MAXDEC
!
	NG_DONE=.TRUE.  		!Successful acceleration
!
	RETURN
	END
!
! 
!
	SUBROUTINE DO_TRANS_AND_NG_V1(NEWPOP,RDPOPS,LST,LEND,NT,ND,NS)
	IMPLICIT NONE
!
	INTEGER LST            !Start depth
	INTEGER LEND           !Final depth
	INTEGER NT
	INTEGER ND
!
! Total number of populations to be accelerated.
! This is simplye (LEND-LST+1)*NT
!
	INTEGER NS
!
	REAL*8 NEWPOP(NS)             !Returned with new estimates.
	REAL*8 RDPOPS(NT,ND,4)
!
	REAL*8 TEMP(4,NS)
	INTEGER I,J,K,L
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
	LOGICAL WEIGHT
!
	IF(NS .NE. NT*(LEND-LST+1))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in DO_TRANS_AND_NG_V1'
	  WRITE(LUER,*)'Invalid NS dimension'
	  WRITE(LUER,*)'      NS=',NS
	  WRITE(LUER,*)'Check NS=',NT*(LEND-LST+1)
	  STOP
	END IF
!
! Use the NG acceleration method to improve the population estimates. 
! We can perform the NG acceleration at each depth individually, or
! over a range of depths. Weight is used to indicate that we
! are to minimize the percentage errors - not the absolute magnitudec
! of the errors. The absolute maximum percentage change is also
! determined.
!
! Rewrite the relevant poulations in a form suitable for NGACCEL.
!
	DO J=1,4
	  DO L=LST,LEND
	    DO I=1,NT
	      K=I+NT*(L-LST)
	      TEMP(J,K)=RDPOPS(I,L,J)
	    END DO
	  END DO
	END DO
	WEIGHT=.TRUE.
!
	CALL NGACCEL(NEWPOP,TEMP,NS,WEIGHT)
!
	RETURN
	END
