!
! Routine to perform an NG accleration for a Comoving-Frame Model. Progam
! uses the last 4 iterations which are stored in the last "4 records" 
! (effectively) of SCRTEMP.
!
! The NG acceleration is perfomed separately on each depth, or over a band
! of depths.
!
! Output is to the last record of SCRTEMP.
!
! Input files required:
!		     MODEL(.DAT)
!                    POINT1(.DAT)
!                    SCRTEMP(.DAT)
! Output files:
!                    POINT1(.DAT)
!                    POINT2(.DAT)
!                    SCRTEMP(.DAT)
!
	PROGRAM DO_NG_V2
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 22-Feb-2014 : Default for NG acceleration has been reset to use 4 iterations.
! Altered 17-Jan-2014 : Changed I3 to I4 to allow for ND > 99
! Altered 01-Nov-2012 : Bug fix with TG option. Values when r < 1 were not being updated.
! Altered 05-May-2007 : Fixed bug when do NG accleration with band size < ND
! Altered 07-Mar-2006 : Acceleratiion can now start at ND_ST.
! Altered 19-May-2004 : Bug fix: NG_DONE set to true, even when no NG done.
! Altered 28-Mar-2004 : Changed to handle new format SCRTEMP files.
!                       Now choice of 3 options: NG, AV and SOR.
! Altered 01-Jul-2003 : IT_STEP inserted as option.  
! Altered 15-Apr-2003 : IFLAG in RD_4_ITS initialized.
!                       NUM_BAD_NG in NG_MIT_OPST initialized.
! Altered 18-Feb-2001 : Extensive rewrite.
!                       NBAND option installed.
!
! Altered 04-Jan-1998 : Cleaned. ND, NT now read from MODEL file.
!                       Based on a very old version in [JDH.BASOL] and 
!                        REWRITE_SCR.
!
	REAL*8, ALLOCATABLE :: RDPOPS(:,:,:)		!NT+3,ND
	REAL*8, ALLOCATABLE :: BIG_POPS(:,:)		!NT+3,ND
	REAL*8, ALLOCATABLE :: POPS(:,:)		!NT,ND
	REAL*8, ALLOCATABLE :: R(:)			!ND
	REAL*8, ALLOCATABLE :: V(:)			!ND
	REAL*8, ALLOCATABLE :: SIGMA(:)			!ND
	REAL*8 TA(2000)
	REAL*8 TB(2000)
!
! Local variables which are adjusted to match the particular model under
! consideration.
!
	REAL*8 T1,T2,T3
	REAL*8 SCALE_FAC
	REAL*8 BIG_FAC
	REAL*8 SF1,SF2
	REAL*8 REAL_LIMIT
!
	INTEGER ND,NT
	INTEGER NBAND
	INTEGER ND_ST
	INTEGER ND_END
	INTEGER IT1,IT2
!
	INTEGER IOS
	INTEGER IREC
	INTEGER NITSF
	INTEGER LST_NG
	INTEGER IFLAG
	INTEGER N_ITS_TO_RD
	INTEGER N_ITS_TO_AV
	INTEGER IT_STEP
	INTEGER I,J,K,L
	INTEGER LOCATION(1)
	INTEGER IVAR
	INTEGER T_INDEX
        INTEGER ITS_PER_NG
!
	INTEGER, PARAMETER :: RITE_N_TIMES=1
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LUSCR=10
!
	LOGICAL NEWMOD
	LOGICAL NG_DONE
	LOGICAL DO_REGARDLESS
	LOGICAL WRITE_RVSIG
	LOGICAL SCALE_INDIVIDUALLY
	LOGICAL REPLACE_VAL
	CHARACTER*4 OPTION
	CHARACTER*132 STRING
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       MODEL(.DAT)'
	WRITE(T_OUT,*)'                                       POINT1(.DAT)'
	WRITE(T_OUT,*)'                                       SCRTEMP(.DAT)'
	WRITE(T_OUT,*)' '
!
! Get basic Model data (i.e., NT and ND).
!
	OPEN(UNIT=12,FILE='MODEL',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of depth') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)ND
	  DO WHILE(INDEX(STRING,'!Total number of variables') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)NT
!
! Iw we couldn't successfully red the MODEL file, get NT and ND from terminal.
!
100	CONTINUE
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read MODEL file'
	  CALL GEN_IN(NT,'Total number of levels')
	  CALL GEN_IN(ND,'Number of depth points')
	END IF 
	CLOSE(UNIT=12)
!
	ALLOCATE (BIG_POPS(NT+3,ND))
	ALLOCATE (POPS(NT,ND))
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
!
! Set def values.
!
	ND_ST=1
	ND_END=ND
	NBAND=ND
	IT_STEP=1
	N_ITS_TO_RD=4
!
	OPTION='NG'
	WRITE(T_OUT,*)'Options are:'
	WRITE(T_OUT,*)'         NG: NG acceleration'
	WRITE(T_OUT,*)'         TG: NG acceleration based on temperature only'
	WRITE(T_OUT,*)'         AV: Average 2 iterations (last iteration is 1, 2nd past is 2 etc)'
	WRITE(T_OUT,*)'        FID: Linear average populations at a single depth based on SE equations'
	WRITE(T_OUT,*)'       FFID: Linear average populations at multiple depthius based on SE values infile SE'
	WRITE(T_OUT,*)'        SOR: Succesive over relaxation (scale last corrections)'
	WRITE(T_OUT,*)'        NSR: Succesive over relaxation (power scale corrections)'
	WRITE(T_OUT,*)'        REP: Repeat the last N corrections'
	WRITE(T_OUT,*)'       UNDO: Undo part of the last correction'
	CALL GEN_IN(OPTION,'Acceleration method:')
	CALL SET_CASE_UP(OPTION,0,0)
	IF(OPTION(1:2) .EQ.'NG')THEN
	  WRITE(T_OUT,'(A)')' '
	  WRITE(T_OUT,'(A)')' This option was inserted mainly for testing purposes.'
	  WRITE(T_OUT,'(A)')' Code uses LAST, LAST-IT_SEP, LAST-2*IT_STEP and LAST-3*IT_STEP for the acceleration'
	  WRITE(T_OUT,'(A)')' The default value of 1 should generally be used.'
	  CALL GEN_IN(IT_STEP,'Iteration step size for NG acceleration')
	  ITS_PER_NG=0
	  DO WHILE(ITS_PER_NG .NE. 4 .AND. ITS_PER_NG .NE. 8)
	    ITS_PER_NG=4
	    CALL GEN_IN(ITS_PER_NG,'Number of iterations used for acceleration (4 or 8)')
	    IF(ITS_PER_NG .NE. 4 .AND. ITS_PER_NG .NE. 8)WRITE(T_OUT,*)'Error: ITS_PER_NG must be 4 or 8'
	  END DO
          N_ITS_TO_RD=ITS_PER_NG
	  DO_REGARDLESS=.TRUE.
	  WRITE(T_OUT,'(/,A)')' The next parameter indicates the number of depths treated simultaneously'
	  WRITE(T_OUT,'(A)')' The acceleration starts in blocks from the inner boundary'
	  CALL GEN_IN(NBAND,'Band width for NG acceleration')
	  CALL GEN_IN(ND_ST,'Only do NG acceleration for the depth in .GE. ND_ST')
	  CALL GEN_IN(ND_END,'Only do NG acceleration for the depth in .LE. ND_END')
	  CALL GEN_IN(DO_REGARDLESS,'Do acceleration independent of corection size?')
	  SCALE_INDIVIDUALLY=.FALSE.
	  IF(DO_REGARDLESS)THEN
	    CALL GEN_IN(SCALE_INDIVIDUALLY,'Scale each population individually at each depth to limit change')
	  END IF
	ELSE IF(OPTION(1:2) .EQ. 'AV')THEN
	  N_ITS_TO_RD=2
	  CALL GEN_IN(N_ITS_TO_RD,'Number of iterations to read')
	  CALL GEN_IN(ND_ST,'Only do Averaging if depth is .GE. ND_ST')
	  CALL GEN_IN(ND_END,'Only do Averaging if depth is .LE. ND_END')
	ELSE IF(OPTION(1:2) .EQ. 'SM')THEN
	  N_ITS_TO_RD=2
	  CALL GEN_IN(N_ITS_TO_RD,'Number of iterations to read')
	ELSE IF(OPTION(1:3) .EQ. 'FID' .OR. OPTION(1:4) .EQ. 'FFID')THEN
	  N_ITS_TO_RD=3
	  CALL GEN_IN(N_ITS_TO_RD,'Number of iterations to read')
	ELSE IF(OPTION(1:2) .EQ. 'TG')THEN
	  N_ITS_TO_RD=3
	  SCALE_FAC=2.0D0
	  IT_STEP=1
	  N_ITS_TO_AV=1
	  CALL GEN_IN(N_ITS_TO_RD,'Number of iterations to read')
	  CALL GEN_IN(ND_ST,'Only do NG acceleration for the depth in .GE. ND_ST')
	  CALL GEN_IN(ND_END,'Only do NG acceleration for the depth in .LE. ND_END')
	  CALL GEN_IN(N_ITS_TO_AV,'Number of iterations to avergae')
	  IF(N_ITS_TO_AV .EQ. 1)CALL GEN_IN(IT_STEP,'Iteration step size for NG acceleration')
	  N_ITS_TO_RD=3*N_ITS_TO_AV
	ELSE IF(OPTION(1:3) .EQ. 'NSR')THEN
	  N_ITS_TO_RD=2
	  SCALE_FAC=2.0D0
	  CALL GEN_IN(SCALE_FAC,'Exponent (N) to power scale (i.e., (1+T1)**N last correction')
	  CALL GEN_IN(ND_ST,'Only do NG acceleration for the depth in .GE. ND_ST')
	  CALL GEN_IN(ND_END,'Only do NG acceleration for the depth in .LE. ND_END')
	ELSE IF(OPTION(1:3) .EQ. 'SOR')THEN
	  N_ITS_TO_RD=2
	  SCALE_FAC=2.0D0
	  BIG_FAC=5.0D0
	  CALL GEN_IN(SCALE_FAC,'Factor to scale last correction by')
          CALL GEN_IN(BIG_FAC,'Maximum correction to any depth')
	  CALL GEN_IN(ND_ST,'Only do NG acceleration for the depth in .GE. ND_ST')
	  CALL GEN_IN(ND_END,'Only do NG acceleration for the depth in .LE. ND_END')
	ELSE IF(OPTION(1:3) .EQ. 'REP')THEN
	  CALL GEN_IN(N_ITS_TO_RD,'# of iteration to step back')
	  CALL GEN_IN(ND_ST,'Only do REP option if the depth is .GE. ND_ST')
	  CALL GEN_IN(ND_END,'Only do UNDO option if the depth is .LE. ND_END')
	ELSE IF(OPTION(1:4) .EQ. 'UNDO')THEN
	  N_ITS_TO_RD=2
	  CALL GEN_IN(ND_ST,'Only do UNDO option if the depth is .GE. ND_ST')
	  CALL GEN_IN(ND_END,'Only do UNDO option if the depth is .LE. ND_END')
	ELSE
	   WRITE(6,*)'Invalid acceleration option'
	   STOP
	END IF
	ALLOCATE (RDPOPS(NT+3,ND,N_ITS_TO_RD))
!
! Read POPULATIONS that were output on last iteration. This is primarily
! done to get NITSF etc.
! 
	NEWMOD=.FALSE.
	IREC=0			!Get last iteration
	CALL SCR_READ_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	IF(NEWMOD)THEN
	  WRITE(T_OUT,*)'Unable to read last iteration in scratch file'
	  STOP
	END IF
	BIG_POPS(1:NT,:)=POPS(1:NT,:)
	BIG_POPS(NT+1,:)=R(:)
	BIG_POPS(NT+2,:)=V(:)
	BIG_POPS(NT+3,:)=SIGMA(:)+1.0D0
	WRITE(6,*)'Read in last iteration'
!
! Read in the last N_ITS_TO_RD estimates of the populations, as output to SCRTEMP.
!
	NEWMOD=.FALSE.
	DO I=1,N_ITS_TO_RD
	  J=IREC-(I-1)*IT_STEP
	  WRITE(6,*)'Reading record', J
	  CALL SCR_READ_V2(R,V,SIGMA,POPS,J,NITSF,RITE_N_TIMES,LST_NG,
	1                WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  RDPOPS(1:NT,:,I)=POPS(1:NT,1:ND)
	  RDPOPS(NT+1,:,I)=R(:)
	  RDPOPS(NT+2,:,I)=V(:)
	  RDPOPS(NT+3,:,I)=SIGMA(:)+1.0D0
	  IF(NEWMOD)THEN
	    WRITE(T_OUT,*)'Unable to read scratch file'
	    STOP
	  END IF
	END DO

	WRITE(6,*)'Starting the accleration'
	IF(OPTION(1:2) .EQ. 'NG')THEN
!
! Perform the NG acceleration.
!
	  J=NT+3
	  T_INDEX=NT
	  CALL NG_MIT_OPTS_V3(BIG_POPS,RDPOPS,ND,J,NBAND,ND_ST,ND_END,T_INDEX,
	1                     ITS_PER_NG,DO_REGARDLESS,SCALE_INDIVIDUALLY,NG_DONE,T_OUT)
!
	  WRITE(6,*)'Finished NG accleration'
	ELSE IF(OPTION(1:2) .EQ. 'AV')THEN
	  IT1=1; CALL GEN_IN(IT1,'Iteration')
	  IT2=2; CALL GEN_IN(IT2,'Iteration 2')
	  BIG_POPS=RDPOPS(:,:,1)
	  DO J=ND_ST,ND_END
	    DO I=1,NT+3
	      BIG_POPS(I,J)=0.5D0*(RDPOPS(I,J,IT1)+RDPOPS(I,J,IT2))
	    END DO
	  END DO
	  NG_DONE=.TRUE.
!
	ELSE IF(OPTION(1:2) .EQ. 'SM')THEN
	  BIG_POPS=RDPOPS(:,:,1)
	  REAL_LIMIT=1.0D-08
	  CALL GEN_IN(REAL_LIMIT,'Replace value when ratio is larger than this value.')
	  DO I=1,NT
	    TA(1:ND)=BIG_POPS(I,1:ND)
	    DO J=2,ND-1
	      TB(J)=BIG_POPS(I,J)/SQRT(BIG_POPS(I,J-1)*BIG_POPS(I,J+1))
	    END DO
	    TB(1)=TA(1)/TA(2); TB(ND)=TA(ND)/TA(ND-1)
	    LOCATION=MINLOC(TB(1:ND)); K=LOCATION(1)
	    IF(TB(K) .LT. REAL_LIMIT)THEN
	      WRITE(6,'(2I4,3ES14.4)',ADVANCE='NO')I,K,TB(K),TA(MAX(1,K-2):MAX(1,K-1))
	      WRITE(6,'(A,ES14.4,A,2ES14.4)')RED_PEN,TA(K:K),DEF_PEN,TA(MIN(ND,K+1):MIN(ND,K+2))
	      CALL GEN_IN(REPLACE_VAL,'Replace value with harmonic average?')
	      J=K-1; IF(K .EQ. 1)J=2
	      L=K+1; IF(K .EQ. ND)L=ND-1
	      IF(REPLACE_VAL)THEN
	        BIG_POPS(I,K)=SQRT(TA(J)*TA(L))
	      ELSE
	        T1=BIG_POPS(I,K)
	        CALL GEN_IN(T1,'New value (default is old value)')
	        BIG_POPS(I,K)=T1
	      END IF
	    END IF
	  END DO
	  NG_DONE=.TRUE.
!
	ELSE IF(OPTION(1:2) .EQ. 'TG')THEN
	  BIG_POPS=RDPOPS(:,:,1)
	  IVAR=NT
	  CALL GEN_IN(IVAR,'Variable?')
	  DO J=ND_ST,ND_END
	    IF(N_ITS_TO_AV .NE. 1)THEN
	      DO K=1,N_ITS_TO_RD/N_ITS_TO_AV
	        RDPOPS(:,J,K)=0.5D0*(RDPOPS(:,J,2*K-1)+RDPOPS(:,J,2*K))
	      END DO
	    END IF
	    T1=(RDPOPS(IVAR,J,2)-RDPOPS(IVAR,J,3))/(RDPOPS(IVAR,J,1)-RDPOPS(IVAR,J,2))
	    WRITE(6,*)'r = Old Cor./ New Cor =',T1,'Depth=',J
	  END DO
	  CALL GEN_IN(SCALE_FAC,'Exponent (N) to power scale (i.e., (1+T1)**N last correction if r<1')
	  DO J=ND_ST,ND_END
	    T1=(RDPOPS(IVAR,J,2)-RDPOPS(IVAR,J,3))/(RDPOPS(IVAR,J,1)-RDPOPS(IVAR,J,2))
	    IF(T1 .LT. -1.0)THEN
	      DO I=1,NT
	        T2=(RDPOPS(I,J,1)-RDPOPS(I,J,2))/(T1-1.0D0)
	        IF(T2 .GT. RDPOPS(I,J,1))T2=RDPOPS(I,J,1)
	        IF(T2 .LT. -0.5D0*RDPOPS(I,J,1))T2=-0.5D0*RDPOPS(I,J,1)
	        BIG_POPS(I,J)=RDPOPS(I,J,1)+T2
	      END DO
	    ELSE IF(T1 .GT. 1.0)THEN
	      IF(T1 .LT. 1.010D0)T1=1.010D0
	      DO I=1,NT
	        T2=(RDPOPS(I,J,1)-RDPOPS(I,J,2))/(T1-1.0D0)
	        IF(T2 .GT. RDPOPS(I,J,1))T2=RDPOPS(I,J,1)
	        IF(T2 .LT. -0.5D0*RDPOPS(I,J,1))T2=-0.5D0*RDPOPS(I,J,1)
	        BIG_POPS(I,J)=RDPOPS(I,J,1)+T2
	      END DO
	    ELSE
	      K=NINT(SCALE_FAC)
	      DO I=1,NT
	        T1=(RDPOPS(I,J,1)-RDPOPS(I,J,2))/RDPOPS(I,J,2)
	        BIG_POPS(I,J)=RDPOPS(I,J,1)*(1.0D0+T1)**K
	      END DO
	    END IF
	  END DO
	  NG_DONE=.TRUE.
	ELSE IF(OPTION(1:3) .EQ. 'FID')THEN
	  BIG_POPS=RDPOPS(:,:,1)
	  CALL GEN_IN(ND_ST,'Depth to change')
	  IT1=1; IT2=2
	  DO WHILE(ND_ST .GT. 0 .AND. ND_ST .LE. ND)
	    J=ND_ST
	    CALL GEN_IN(IT1,'Highest iteration')
	    CALL GEN_IN(IT2,'Lowest iteration')
	    CALL GEN_IN(SF1,'SE equation IT1') 
	    CALL GEN_IN(SF2,'SE equation IT2') 
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable?')
	    WRITE(6,*)RDPOPS(IVAR,J,IT1),RDPOPS(IVAR,J,IT2)
	    T1=-SF2/(SF1-SF2); T2=SF1/(SF1-SF2)
	    SF1=T1; SF2=T2
	    WRITE(6,*)RDPOPS(IVAR,J,IT1),RDPOPS(IVAR,J,IT2),SF1*RDPOPS(IVAR,J,IT1)+SF2*RDPOPS(IVAR,J,IT2)
	    DO I=1,NT
	      BIG_POPS(I,J)=SF1*RDPOPS(I,J,IT1)+SF2*RDPOPS(I,J,IT2)
	    END DO
	    NG_DONE=.TRUE.
	    CALL GEN_IN(ND_ST,'Depth to change')
	  END DO
	ELSE IF(OPTION(1:4) .EQ. 'FFID')THEN
	  BIG_POPS=RDPOPS(:,:,1)
	  OPEN(UNIT=22,STATUS='OLD',ACTION='READ',FILE='SE')
	  DO WHILE(1 .EQ. 1)
	    READ(22,*,END=200)ND_ST,ND_END
	    READ(22,*,END=200)IT1,IT2
	    READ(22,*,END=200)IVAR
	    READ(22,*)TB(ND_ST:ND_END)
	    READ(22,*)TA(ND_ST:ND_END)
	    DO J=ND_ST,ND_END
	      T1=-TB(J)/(TA(J)-TB(J)); T2=TA(J)/(TA(J)-TB(J))
	      TA(J)=T1; TB(J)=T2
	      WRITE(6,*)RDPOPS(IVAR,J,IT1),RDPOPS(IVAR,J,IT2),TA(J)*RDPOPS(IVAR,J,IT1)+TB(J)*RDPOPS(IVAR,J,IT2)
	      WRITE(21,*)RDPOPS(IVAR,J,IT1),RDPOPS(IVAR,J,IT2),TA(J)*RDPOPS(IVAR,J,IT1)+TB(J)*RDPOPS(IVAR,J,IT2)
	      DO I=1,NT
	        BIG_POPS(I,J)=TA(J)*RDPOPS(I,J,IT1)+TB(J)*RDPOPS(I,J,IT2)
	      END DO
	    END DO
	    NG_DONE=.TRUE.
	  END DO
200	  CONTINUE
	ELSE IF(OPTION(1:3) .EQ. 'NSR')THEN
	  K=NINT(SCALE_FAC)
	  DO J=ND_ST,ND_END
	    DO I=1,NT+3
	      T1=(RDPOPS(I,J,1)-RDPOPS(I,J,2))/RDPOPS(I,J,2)
	      BIG_POPS(I,J)=RDPOPS(I,J,1)*(1.0D0+T1)**K
	    END DO
	  END DO
	  NG_DONE=.TRUE.
	ELSE IF(OPTION(1:3) .EQ. 'SOR')THEN
	  DO J=ND_ST,ND_END
	    T1=-1000.0D0
	    T2=1000.0D0
	    DO I=1,NT+3
	      T1=MAX(T1,(RDPOPS(I,J,1)-RDPOPS(I,J,2))/RDPOPS(I,J,1))
	      T2=MIN(T2,(RDPOPS(I,J,1)-RDPOPS(I,J,2))/RDPOPS(I,J,1))
	    END DO
	    T2=ABS(T2)
	    T3=SCALE_FAC
	    IF(T3*T1 .GT. BIG_FAC)T3=BIG_FAC/T1
	    IF(T3*T2 .GT. (1.0D0-1.0D0/BIG_FAC))T3=(1.0D0-1.0D0/BIG_FAC)/T2
	    WRITE(6,'(I4,3ES14.4)')J,T1,T2,T3
	    DO I=1,NT+3
	      BIG_POPS(I,J)=RDPOPS(I,J,1)+T3*(RDPOPS(I,J,1)-RDPOPS(I,J,2))
	    END DO
	  END DO
	  NG_DONE=.TRUE.
	ELSE IF(OPTION(1:3) .EQ. 'REP')THEN
	  BIG_POPS=RDPOPS(:,:,1)
	  DO J=ND_ST,ND_END
	    DO I=1,NT+3
	      T1=(RDPOPS(I,J,1)-RDPOPS(I,J,N_ITS_TO_RD))/BIG_POPS(I,J)
	      IF(T1 .LT. -0.9)T1=-0.9
	      IF(T1 .GT. 10.0)T1=10.0
	      BIG_POPS(I,J)=BIG_POPS(I,J)*(1.0D0+T1)
	    END DO
	  END DO
	  NG_DONE=.TRUE.
	ELSE IF(OPTION(1:4) .EQ. 'UNDO')THEN
	  BIG_POPS=RDPOPS(:,:,1)
	  DO J=ND_ST,ND_END
	    DO I=1,NT+3
	      BIG_POPS(I,J)=RDPOPS(I,J,2)
	    END DO
	  END DO
	  NG_DONE=.TRUE.
	END IF
!
	POPS(1:NT,:)=BIG_POPS(1:NT,:)
	R(2:ND-1)=BIG_POPS(NT+1,2:ND-1)
	V(2:ND-1)=BIG_POPS(NT+2,2:ND-1)
	SIGMA(2:ND-1)=BIG_POPS(NT+3,2:ND-1)-1.0D0
!
	IF(NG_DONE)THEN
	  NITSF=NITSF+1
	  LST_NG=NITSF
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS,IREC,NITSF,RITE_N_TIMES,
	1                 LST_NG,WRITE_RVSIG,NT,ND,LUSCR,NEWMOD)
	  IF(OPTION(1:2) .EQ.'NG')THEN
	    WRITE(T_OUT,*)'Results of successful NG acceleration ',
	1                 'output to SCRTEMP.'
	  ELSE IF(OPTION(1:3) .EQ. 'NSR')THEN
	    WRITE(T_OUT,*)'Results of successful power SOR (NSR) output to SCRTEMP.'
	  ELSE IF(OPTION(1:3) .EQ.'SOR')THEN
	    WRITE(T_OUT,*)'Results of successful SOR output to SCRTEMP.'
	  ELSE IF(OPTION(1:3) .EQ.'REP')THEN
	    WRITE(T_OUT,*)'Successful extrapolation (REP) output to SCRTEMP.'
	  ELSE IF(OPTION(1:4) .EQ.'UNDO')THEN
	    WRITE(T_OUT,*)'Successfully undid part of the last correction.'
	  ELSE
	    WRITE(T_OUT,*)'Last 2 iterations averaged and output to SCRTEMP.'
	  END IF
	END IF
!
	STOP
	END
!
!
!
! Altered 3-April-1989. -Call changed. MAXDEC, MAXINC and IFLAG installed,
!                        ACCELERATE flag removed.
!
! If IFLAG is returned with zero, the NG acceleration was successful.
! Otherwise an error condition occurred.
!
	SUBROUTINE GENACCEL_V3(NEWPOP,RDPOPS,ITS_PER_NG,LST,LEND,NT,ND,NS)
	IMPLICIT NONE
!
	INTEGER LST,LEND
	INTEGER NT
	INTEGER ND
	INTEGER NS
        INTEGER ITS_PER_NG
	REAL*8 NEWPOP(NS)
	REAL*8 RDPOPS(NT,ND,ITS_PER_NG)
!
	REAL*8 TEMP(ITS_PER_NG,NS)
	INTEGER I,J,K,L
	LOGICAL WEIGHT
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
	DO J=1,ITS_PER_NG
	  DO L=LST,LEND
	    DO I=1,NT
	      K=I+NT*(L-LST)
	      TEMP(J,K)=RDPOPS(I,L,J)
	    END DO
	  END DO
	END DO
	WEIGHT=.TRUE.
!
	WRITE(6,*)'Calling NGACCEL'
!	CALL NGACCEL(NEWPOP,TEMP,NS,WEIGHT)
	I=ITS_PER_NG-2
	CALL NGACCEL_ARB_ORD(NEWPOP,TEMP,NS,I,WEIGHT)
!
	RETURN
	END
!
!
	SUBROUTINE NG_MIT_OPTS_V3(POPS,RDPOPS,ND,NT,NBAND,ND_ST,ND_END,T_INDEX,
	1                           ITS_PER_NG,DO_REGARDLESS,SCALE_INDIVIDUALLY,NG_DONE,LUER)
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 07-MAr-2006: ND_ST inserted in call (Name not changed: alrady have _V1).
! Altered 30-May-2003: Bug fix: When SCALE_INDIVIDUALLY was true, not all
!                         populations were being accelerated.
!
	INTEGER ND
	INTEGER NT
	INTEGER NBAND
	INTEGER ND_ST
	INTEGER ND_END
	INTEGER T_INDEX
	INTEGER LUER
	REAL*8 POPS(NT,ND)
	REAL*8 RDPOPS(NT,ND,ITS_PER_NG)
	LOGICAL DO_REGARDLESS
	LOGICAL SCALE_INDIVIDUALLY
	LOGICAL NG_DONE
!
	REAL*8 NEWPOP(NT,ND)
	REAL*8 VEC_INC(ND)
	REAL*8 VEC_DEC(ND)
	REAL*8 INT_ARRAY(ND)
!
	REAL*8 RINDX(NT)
	REAL*8 RAT(NT)
	REAL*8 TA(MAX(NT,ND))
	REAL*8 TB(MAX(NT,ND))
!
	REAL*8 T1
	REAL*8 LOCINC,LOCDEC
	REAL*8 MAXINC,MAXDEC
!
	INTEGER NS
	INTEGER NUM_BAD_NG
	INTEGER DEC_LOC,INC_LOC
	INTEGER I,K,L
	INTEGER LST,LEND
        INTEGER ITS_PER_NG
	INTEGER IOS
	INTEGER, PARAMETER :: IONE=1
	LOGICAL DO_PLTS
	CHARACTER(LEN=80) STRING
!
	NUM_BAD_NG=0
	VEC_INC(1:ND)=0.0D0
	VEC_DEC(1:ND)=0.0D0
	IF(NBAND .GE. ND .AND. ND_ST .EQ. 1 .AND. ND_END .EQ. ND)THEN
	  NS=NT*ND
	  CALL GENACCEL_V3(NEWPOP,RDPOPS,ITS_PER_NG,IONE,ND,NT,ND,NS)
	ELSE
	  DO K=ND_END,ND_ST,-NBAND
	    LST=MAX(K-NBAND+1,ND_ST)
	    LEND=K				!LST+NBAND-1
	    NS=(LEND-LST+1)*NT
	    CALL GENACCEL_V3(NEWPOP(1,LST),RDPOPS,ITS_PER_NG,LST,LEND,NT,ND,NS)
	  END DO
	  IF(ND_ST .GT. 1)NEWPOP(:,1:ND_ST-1)=POPS(:,1:ND_ST-1)
	  IF(ND_END .LT. ND)NEWPOP(:,ND_END+1:ND)=POPS(:,ND_END+1:ND)
	END IF
!
	WRITE(6,*)' '//RED_PEN
	WRITE(6,'(1X,90A)')('*',I=1,90)
	WRITE(6,*)' '//DEF_PEN
	WRITE(6,*)'Plotting '//RED_PEN//'T correction due to acceleration (red),'//
	1            DEF_PEN//'and '//BLUE_PEN//'the last iterative correction (blue).'//DEF_PEN
	WRITE(6,*)'If the two corrections are not (generally) of the same sign, the acceleration'
	WRITE(6,*)'should be CANCELLED.'
	WRITE(6,*)' '//RED_PEN
	WRITE(6,'(1X,90A)')('*',I=1,90)
	WRITE(6,*)' '//DEF_PEN
!
	T1=0.0D0
	DO L=1,ND
	  INT_ARRAY(L)=L
	  TA(L)=100.0D0*(1.0D0-RDPOPS(T_INDEX,L,1)/NEWPOP(T_INDEX,L))
	  TB(L)=100.0D0*(1.0D0-RDPOPS(T_INDEX,L,3)/RDPOPS(T_INDEX,L,2))
	  T1=MAX(ABS(TA(L)),T1);  T1=MAX(ABS(TB(L)),T1)
	END DO
	IF(T1 .EQ. 0.0D0)THEN
	  OPEN(UNIT=12,FILE='CORRECTION_LINK',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    K=1
	    IF(IOS .NE. 0)THEN
	      DO I=1,27
	        READ(10,'(A)')STRING
	        IF(STRING .NE. ' ')WRITE(6,'(A)')TRIM(STRING)
	      END DO
	      CLOSE(UNIT=12)
	    END IF
	  CALL GEN_IN(K,'Input new variable to be plotted as T correction is zero')
	  DO L=1,ND
	    TA(L)=100.0D0*(1.0D0-RDPOPS(K,L,1)/NEWPOP(K,L))
	    TB(L)=100.0D0*(1.0D0-RDPOPS(K,L,3)/RDPOPS(K,L,2))
	  END DO
	  CALL DP_CURVE(ND,INT_ARRAY,TA)
	  CALL DP_CURVE(ND,INT_ARRAY,TB)
	  CALL GRAMON_PGPLOT('Depth Index','100[X(new)-X(old)]/X(new)',' ',' ')
	ELSE
	  CALL DP_CURVE(ND,INT_ARRAY,TA)
	  CALL DP_CURVE(ND,INT_ARRAY,TB)
	  CALL GRAMON_PGPLOT('Depth Index','100[T(new)-T(old)]/T(new)',' ',' ')
	END IF
!
	WRITE(6,*)' '
	WRITE(6,'(1X,70A)')('*',I=1,70)
	WRITE(6,*)' '
        WRITE(6,*)'Plotting Y(NG accel)/Y(old) for each depth'
        WRITE(6,*)'Ten depths plotted each time.'
        WRITE(6,*)'These plots are for diagnostic purposes only'
        WRITE(6,*)' '
	WRITE(6,'(1X,70A)')('*',I=1,70)
        WRITE(6,*)' '
	DO_PLTS=.FALSE.
        CALL GEN_IN(DO_PLTS,'Plot the above diagnostic information?')
	IF(DO_PLTS)THEN
	  DO L=1,ND
	    DO K=1,NT
	      RINDX(K)=K
	      RAT(K)=NEWPOP(K,L)/POPS(K,L)
	    END DO
	    CALL DP_CURVE(NT,RINDX,RAT)
	    IF(MOD(L,10) .EQ. 0 .OR. L .EQ. ND)THEN
	      WRITE(6,*)' ** Plotting depths',L-9+MOD(L,10),' to ',L
	      CALL GRAMON_PGPLOT('Index','Ratio',' ',' ')
	    END IF
	  END DO
	END IF
	WRITE(6,*)' '
!
! Now check whether the NG acceleation has been reasonable.
! NB: POPS(K,L) is only likely to be zero for sigma.
!
	MAXINC=-1000.0
	MAXDEC=1000.0
	DO L=1,ND 
	  LOCINC=-1000.0
	  LOCDEC=1000.0	
	  T1=LOCDEC
	  DO K=1,NT
	    IF(POPS(K,L) .NE. 0.0D0)T1=NEWPOP(K,L)/POPS(K,L)
	    LOCINC=MAX(LOCINC,T1)
	    LOCDEC=MIN(LOCDEC,T1)
	  END DO
	  VEC_INC(L)=100.0D0*(LOCINC-1.0D0)
	  VEC_DEC(L)=100.0D0*(LOCDEC-1.0D0)
!
	  IF(DO_REGARDLESS)THEN
	    IF(SCALE_INDIVIDUALLY)THEN
	      DO K=1,NT
	        IF(NEWPOP(K,L) .LT. 0.1D0*POPS(K,L))THEN
	          POPS(K,L)=0.1D0*POPS(K,L)
	        ELSE IF(NEWPOP(K,L) .GT. 10.0D0*POPS(K,L))THEN
	          POPS(K,L)=10.0D0*POPS(K,L)
	        ELSE
	          POPS(K,L)=NEWPOP(K,L)
	        END IF
	      END DO
	    ELSE
	      T1=1.0D0
	      DO K=1,NT
                IF(POPS(K,L) .EQ. NEWPOP(K,L))THEN
                ELSE IF(NEWPOP(K,L) .LT. 0.1D0*POPS(K,L))THEN
                  T1=MIN( T1, 0.9D0*POPS(K,L)/(POPS(K,L)-NEWPOP(K,L)) )
                ELSE IF(NEWPOP(K,L) .GT. 10.0D0*POPS(K,L))THEN
                  T1=MIN( T1, 9.0D0*POPS(K,L)/(NEWPOP(K,L)-POPS(K,L)) )
                END IF
	      END DO
	      DO K=1,NT
                POPS(K,L)=POPS(K,L)+T1*(NEWPOP(K,L)-POPS(K,L))
              END DO
	      IF(T1 .NE. 1.0D0)WRITE(LUER,8000)L,LOCINC,LOCDEC
8000	      FORMAT(1X,'Scaling performed at depth ',I4,':',
	1          1X,'Biggest increase/decrease was ',2ES12.2)
	    END IF
!
	    IF(LOCINC .GT. MAXINC)THEN
	      MAXINC=LOCINC
	      INC_LOC=L
	    END IF
	    IF(LOCDEC .LT. MAXDEC)THEN
	      MAXDEC=LOCDEC
	      DEC_LOC=L
	    END IF
!
! Before storing the NG acceleration at this depth, we check to
! see whether the predicted corrections are "reasonable".
!
	  ELSE IF(LOCINC .GT. 10.1D0 .OR. LOCDEC .LT. 0.09D0)THEN
	    NUM_BAD_NG=NUM_BAD_NG+1
	    WRITE(LUER,*)'NUM_BAD_NG=',NUM_BAD_NG
	    WRITE(LUER,9000)L,LOCINC,LOCDEC
9000	    FORMAT(1X,'No NG acceleration performed at depth ',I4,/,
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
	WRITE(6,*)' '
	WRITE(6,'(1X,70A)')('*',I=1,70)
	WRITE(6,*)' '
	WRITE(6,*)'Plotting largest % decrease and increase at each depth.'
	WRITE(6,*)' '
	WRITE(6,'(1X,70A)')('*',I=1,70)
	WRITE(6,*)' '
	CALL DP_CURVE(ND,INT_ARRAY,VEC_DEC)
	CALL DP_CURVE(ND,INT_ARRAY,VEC_INC)
	CALL GRAMON_PGPLOT('Depth Index','100(N/O-1)',' ',' ')
!
! By definition MAXINC and MININC are both positive. These are only
! defined for successful NG accelerations (does not include any depths
! at which NG acceleration did not work).
!
	MAXINC=100.0D0*(MAXINC-1.0D0)
	MAXDEC=100.0D0*(1.0D0/MAXDEC-1.0D0)
	WRITE(LUER,9800)INC_LOC,MAXINC
9800	FORMAT(1X,'Max NG % increase at depth ',I4,' is',1PE10.2)
	WRITE(LUER,9900)DEC_LOC,MAXDEC
9900	FORMAT(1X,'Max NG % decrease at depth ',I4,' is',1PE10.2)
C
	NG_DONE=.TRUE.  		!Successful acceleration
C
	RETURN
	END
