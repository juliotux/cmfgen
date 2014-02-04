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
	SUBROUTINE DO_NG_BAND_ACCEL_V2(POPS,R,V,SIGMA,ROLD,
	1             NT,ND,NG_BAND,NG_DONE,MAXDEC,MAXINC,LUSCR,LUER)
	IMPLICIT NONE
!
! Altered 15-Jan-2009 : Minor bug fixed.
!                       Populations were reset to IREC-3 iteration if NG acceleration failed.
! Altered 02-May-2004 : Changed to handle new SCRTEMP file format.
!                       Allows for models in which R, V, and SIGMA are also output to
!                          SCRTEMP on every iteration.
!                       Renamed from V1: R, V, SIGMA, and ROLD added to call.
	INTEGER NT
	INTEGER ND
	INTEGER NG_BAND
	INTEGER LUSCR
	INTEGER LUER
!
	REAL*8 POPS(NT,ND)
	REAL*8 ROLD(ND)
	REAL*8 R(ND)
	REAL*8 V(ND)
	REAL*8 SIGMA(ND)
!
	REAL*8 MAXDEC
	REAL*8 MAXINC
	LOGICAL NG_DONE
!
! Local arrays. We will acelerate R, V, and SIGMA as well. If R is not changing,
! the accleration wiull keep R constant since acceleration is linear.
!
	REAL*8 BIG_POPS(NT+3,ND)
	REAL*8 RDPOPS(NT+3,ND,4)
!
! Local variables which are adjusted to match the particular model under
! consideration.
!
	INTEGER IOS
	INTEGER IREC
	INTEGER NITSF
	INTEGER IT_STEP
	INTEGER LST_NG
	INTEGER IFLAG
!
	INTEGER NT_MOD
	INTEGER I,J
	INTEGER, PARAMETER :: RITE_N_TIMES=1
!
	LOGICAL NEWMOD
	LOGICAL WRITE_RVSIG
!
	IT_STEP=1 		!Use last 4 consecuitive iterations.
	IFLAG=0
!
! Read POPULATIONS that were output on last iteration. This is primarily
! done to get NITSF etc.
! 
	NEWMOD=.FALSE.
	IREC=0
	CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                    WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	IF(NEWMOD)THEN
	  WRITE(LUER,*)'Unable to read last iteration in scratch file in DO_NG_BAND_ACCEL_V2'
	  WRITE(LUER,*)'Unrecoverable error (1) as these should have been written recently'
	  STOP
	  NG_DONE=.FALSE.
	  RETURN
	END IF
	BIG_POPS(1:NT,:)=POPS(1:NT,:)
	BIG_POPS(NT+1,:)=R(:)
	BIG_POPS(NT+2,:)=V(:)
	BIG_POPS(NT+3,:)=SIGMA(:)+1.0D0
!
! Read in the last 4 estimates of the poplations, as output to SCRTEMP.
! We read in the data backwards so that POPS, R, V and SIGMA contain
! the current iteration values on the last read.
!
! If the read fails, we try to recover the current POPS, R, V and SIGMA.
! We should not get a failure here.
!
	NEWMOD=.FALSE.
	DO I=4,1,-1
	  J=IREC-(I-1)*IT_STEP
	  CALL SCR_READ_V2(R,V,SIGMA,POPS,J,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  RDPOPS(1:NT,:,I)=POPS(1:NT,1:ND)
	  RDPOPS(NT+1,:,I)=R(:)
	  RDPOPS(NT+2,:,I)=V(:)
	  RDPOPS(NT+3,:,I)=SIGMA(:)+1.0D0
	  IF(NEWMOD)THEN
	    WRITE(LUER,*)'Unable to read scratch file in DO_NG_BAND_ACCEL_V2'
	    WRITE(LUER,*)'Trying to read IREC ',J
	    CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	    IF(NEWMOD)THEN
	      WRITE(LUER,*)'Error is not recoverable'
	      STOP
	    END IF
            NG_DONE=.FALSE.
	    RETURN
	  END IF
	END DO
!
! Now perform the NG acceleration. The accelerated estimates for the populations
! are returned in BIG_POPS.
!
	NT_MOD=NT+3
	CALL NG_MIT_BAND_OPT_V1(BIG_POPS,RDPOPS,ND,NT_MOD,NG_BAND,NG_DONE,MAXDEC,MAXINC,LUER)
!
! We don't set R(1), R(ND) etc, so as to make sure these remain identical
! to the initial values., These should not change, even if we are changing the
! R grid.
!
	IF(NG_DONE)THEN
	  ROLD(1:ND)=R(1:ND)
          POPS(1:NT,:)=BIG_POPS(1:NT,:)
          R(2:ND-1)=BIG_POPS(NT+1,2:ND-1)
          V(2:ND-1)=BIG_POPS(NT+2,2:ND-1)
          SIGMA(2:ND-1)=BIG_POPS(NT+3,2:ND-1)-1.0D0
	END IF
!
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
	MAXINC=-1000.0D0
	MAXDEC=1000.0D0
	DO L=1,ND 
	  LOCINC=-1000.0D0
	  LOCDEC=1000.0D0
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
9000	    FORMAT(1X,'No NG acceleration performed at depth ',I3,/,
	1          1X,'Biggest increase was ',1PE10.2,/,
	1          1X,'Biggest decrease was ',E10.2)
	
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
