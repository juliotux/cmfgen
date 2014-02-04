C 
C
C
C Altered 3-April-1989. -Call changed. MAXDEC, MAXINC and IFLAG installed,
C                        ACCELERATE flag removed.
C
C If IFLAG is returned with zero, the NG acceleration was successful.
C Otherwise an error condition occurred.
C
	SUBROUTINE GENACCEL(POPS,MAXINC,MAXDEC,IFLAG,NT,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 - NV removed.
C                       RDPOPS,TEMP and NEWEST now dimensionded using ND and NT.
C
C Altered 10-Jan-1991 - Finay Cray compatible version with bugs fixed.
C Altered 19_dec-1991 - Made Cray compatibility. Parameters for opening
C                       direct access file obtained from subroutine
C                       call.
C Altered 05-Sep-1989 - Bug fix: CUR_REC was not being evaluated correctly
C                       when NUM_RECS was greater than 1. Testing inserted
C                       to check on validity of corrections. Maximum
C                       correction limited to a factor of 2. If correction
C                       larger, NG acceleration is abandoned for this depth.
C                       If more than 3 bad NG accelerations occur, NG
C                       acceleration is cancelled.
C
	INTEGER IFLAG,NT,ND
	REAL*8 POPS(NT*ND),MAXINC,MAXDEC
C
	REAL*8 RDPOPS(ND*NT,4),TEMP(4,NT),NEWEST(NT)
C
C Note that REC_SIZE is the size of the output record in bytes.
C           REC_LEN  is the size of the output record in COMPUTER units.
C           N_PER_REC is the number of numbers per record.
C           NUM_TIME is the number of times each iteration has been written
C                to the scratch file.
C           WORD_SIZE is the size of the number to be output in bytes.
C           UNIT_SIZE is the number of bytes per unit used to specify
C           the size of a direct access file.
C
	REAL*8 T1,LOCINC,LOCDEC
	INTEGER ARRAYSIZE,IST,IEND
	INTEGER REC_SIZE,REC_LEN,NUM_RECS,N_PER_REC,NUM_TIME
	INTEGER UNIT_SIZE,WORD_SIZE
C
	INTEGER I,J,K,L,NPREV,NITSF,IT_CNT
	INTEGER SRT_REC_M1,LUER,ERROR_LU
	INTEGER INDK,INDEXST,CHK,IINC,IDEC,NUM_BAD_NG
	LOGICAL WEIGHT
	EXTERNAL ERROR_LU
C
	NUM_BAD_NG=0
	MAXINC=-100.0D0
	MAXDEC=100.0D0
        LUER=ERROR_LU()
C
C Determine the record size, and the number of records that
C need to be written out to fully write out the population vector.
C As this is computer and installation dependent, we call a subroutine
C to return the parameters.
C
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,
	1                 WORD_SIZE,N_PER_REC)
	ARRAYSIZE=NT*ND			
	NUM_RECS=INT( (ARRAYSIZE-1)/N_PER_REC ) + 1
	REC_LEN=REC_SIZE/UNIT_SIZE
C
C Check we can rite R,V, and SIGMA to one record. This error
C should never occurr, since the limit is >2000 Words. This
C check is not needed here, but will remind me to alter this
C routine if I alter SCR_RITE and SC_READ.
C
	IF(N_PER_REC .LT. 3*ND)THEN
	  WRITE(LUER,*)'Error in GENACCEL --- record size smaller '//
	1              'than length of R,V and sigma',3*ND
	  STOP
	END IF
C
C Read in pointer to data file. NB. Format has changed. IT_CNT is now
C the iteration count (but may differ from NITSF is problems occurred
C read/riting to SCRTEMP) where as previously it referred to the actual
C record.
C
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
C
C		' OLD MODEL '
C
	OPEN(UNIT=20,FILE='SCRTEMP',FORM='UNFORMATTED',
	1    ACCESS='DIRECT',STATUS='OLD',RECL=REC_LEN,IOSTAT=CHK)
	  IF(CHK .NE. 0)THEN
	    WRITE(LUER,*)'Error opening SCRTEMP for input in GENACCEL'
	    WRITE(LUER,*)'IOSTAT=',CHK
	    CLOSE(UNIT=20)
	    IFLAG=5
	    RETURN
	  END IF
C
C SRT_REC_M1 + 1 is the first record for the NPREV previous iteration.
C If output more than once, it refers to the first output.
C
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
200	        FORMAT(1X,'SRT_REC_M1=',I3,5X,'IOSAT=',I7)
	        IFLAG=6
	        CLOSE(UNIT=20)
	        RETURN
	      END IF
	    END DO
	  END DO
	  CLOSE(UNIT=20)
C
C Have now read in the Four previous iterations. We can now use
C the NG acceleration method at each depth to improve the
C population estimates. Weight is used to indicate that we
C are to minimize the percentage errors - not the absolute magnitudec
C of the errors. The absolute maximum percentage change is also
C determined.
C
	DO L=1,ND
	  LOCINC=-100.0D0
	  LOCDEC=100.0D0
	  INDEXST=(L-1)*NT+1
	  DO J=1,4
	    DO K=1,NT
	      INDK=INDEXST+K-1
	      TEMP(J,K)=RDPOPS(INDK,J)
	    END DO
	  END DO
	  WEIGHT=.TRUE.
	  CALL NGACCEL(NEWEST,TEMP,NT,WEIGHT)
	  DO K=1,NT
            INDK=INDEXST+K-1
	    T1=NEWEST(K)/POPS(INDK)
	    LOCINC=MAX(LOCINC,T1)
	    LOCDEC=MIN(LOCDEC,T1)
C
C Following statement is required only for TSTNG.
C
C	    POPS(INDK)=RDPOPS(INDK,1)
	  END DO
C
C Before storing the NG acceleration at this depth, we check to
C see whether the predicted corrections are "reasonable".
C
	  IF(LOCINC .GT. 10.1D0 .OR. LOCDEC .LT. 0.09D0)THEN
	    NUM_BAD_NG=NUM_BAD_NG+1
	    WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	    DO K=1,NT
	      INDK=INDEXST+K-1
	      POPS(INDK)=RDPOPS(INDK,1)
	    END DO
	    WRITE(LUER,9000)L,LOCINC,LOCDEC
9000	    FORMAT(1X,'No NG acceleration performed at depth ',I3,/,
	1          1X,'Biggest increase was ',1P,E10.2,/,
	1          1X,'Biggest decrease was ',E10.2)
	
	  ELSE
	    DO K=1,NT
              INDK=INDEXST+K-1
	      POPS(INDK)=NEWEST(K)
	    END DO
	    IF(LOCINC .GT. MAXINC)THEN
	      MAXINC=LOCINC
	      IINC=L
	    END IF
	    IF(LOCDEC .LT. MAXDEC)THEN
	      MAXDEC=LOCDEC
	      IDEC=L
	    END IF
	  END IF
	END DO
C
	WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	IF(NUM_BAD_NG .GT. 3)THEN
	  WRITE(LUER,*)'Too many bad NG accelerations - '//
	1              'NG acceleration cancelled'
	  IFLAG=1
C
C Restore old population estimates.
C
	  DO K=1,ND*NT
	    POPS(K)=RDPOPS(K,1)
	  END DO
	  RETURN
	END IF
C
C By definition MAXINC and MININC are both positive. These are only
C defined for successful NG accelerations (does not include any depths
C at which NG acceleration did not work).
C
	MAXINC=100.0D0*(MAXINC-1.0D0)
	MAXDEC=100.0D0*(1.0D0/MAXDEC-1.0D0)
	WRITE(LUER,9800)IINC,MAXINC
9800	FORMAT(1X,'Max NG % increase at depth ',I3,' is',1P,E10.2)
	WRITE(LUER,9900)IDEC,MAXDEC
9900	FORMAT(1X,'Max NG % decrease at depth ',I3,' is',1P,E10.2)
C
	IFLAG=0  		!Successful acceleration
C
	RETURN
	END
