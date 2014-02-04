C
C Routine to compute the convolution of J with the electron scattering
C redistribution function.
C
C Method uses a "Power law [in v] approximation" for the electron scattering
C function (Exponential in ln v].
C
C Basic method used is due to Rybcki and Hummer (A\&A,  ). This has
C been modified so that electron scattering preserves the Planck-function at
C depth. This scaling is necessary to recover LTE depth, and arises since
C the RH formalism neglects both Compton scattering, and stimulated electron 
C scattering. The preservation of the Planck function is achieved by a
C wavelength shift on the Wien side, and a scaling everywhere else.
C
	SUBROUTINE COMP_J_CONV_V2(J_STORE,J_SIZE,NU,TEMP,ND,NCF,
	1             LU_IN,FILE_IN,CONT_REC,RD_NU,ALLOW_UNEQUAL_FREQ,
	1             LU_OUT,FILE_OUT)
	IMPLICIT NONE
C
C Altered 20-Oct-2003 : No scaling peformed when PLANCK_FN has zero''s.
C Altered 18-Aug-2003 : FILE_DATE inserted into INFO call.
C Altered 17-Oct-2002 : Computation of BB changed to avoid floating overflow.
C Altered 16-Jun-2000 : Record length can no be obtained form INFO files.
C Altered 21-Dec-1999 : Allocation of additional memory checked.
C Altered 03-Jan-1999 : J can now be passed in J_STORE in the CALL.
C                       If passed, FILE_IN must contain
C                                        'J PASSED VIA CALL'
C                      Similarly, the convolved J can be returned in
C                      call.  To do this, FILE_OUT must contain
C                                       'RETURN J VIA CALL'
C
C Altered 03-Dec-1998 : Call altered: Changed to version V2
C                         NDEXT can now be passes as ND. 
C                         Better handling of dimensions of C_STORE.
C                         ALLOW_UNEQ_FREQ inserted in call. This will
C                           allow a smaller EDDFACTOR file to be used.
C			  NT, NUM_BNDS removed from call
C
C Altered 10-Feb-1998 : Section introduced to preserve Planck Function at
C                          depth.
C Altered 05-Feb-1998 : Bug fix. Routine was erroneously using T instead of
c                          SQRT(T) to compute BETA. Effect would have been
C                          for P Cygni type models, but large at depth where
C                          T is large.
C
	INTEGER J_SIZE	!Total number of elements in J_STORE
	INTEGER ND		!Number of depth points.
	INTEGER NCF		!Number of continuum frequencies
	INTEGER LU_IN
	INTEGER LU_OUT
	INTEGER CONT_REC	!Record in file that points to beginning of Jc
	LOGICAL RD_NU
	LOGICAL ALLOW_UNEQUAL_FREQ
C
C It is assumed that BA can be passed for J_STORE. J_SIZE is the maximum
C size of the passed vector. 
C
	REAL*8 J_STORE(ND,J_SIZE/ND)		!BA matrix in main program.
	REAL*8 NU(NCF)				!Frequency (10^15 Hz)
	REAL*8 TEMP(ND)				!Temperature (10^4 K)
C
	CHARACTER*(*) FILE_IN
	CHARACTER*(*) FILE_OUT
C
C Simultaneous equations are assumed to have the form
C
C    A(i).X(i-1) - [H(i)+A(i)+C(i)].X(i) - C(i).X(i+1) = D(i)
C
	REAL*8 A(NCF)
	REAL*8 H(NCF)
	REAL*8 C(NCF)
	REAL*8 D(NCF)
C
	REAL*8 A_STORE(NCF)
	REAL*8 C_STORE(NCF)
	REAL*8 J_ES(NCF)
C
	REAL*8 PLANCK_FN(NCF)
	REAL*8 PLANCK_ES(NCF)
	REAL*8 PLANCK_NU(NCF)
C
C Provide some extra storage in case J_STORE dimensions are not large
C enough.
C
	REAL*8, ALLOCATABLE :: EXTRA_J_ST(:,:)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL FILE_OPEN
	LOGICAL UNEQUAL_FREQ
C
	INTEGER IONE
	PARAMETER (IONE=1)
C
	REAL*8 BETA
	REAL*8 T1,T2,T3
	REAL*8 D1,D2,DH
	INTEGER I,K,L,ML,INIT_REC,IOS
	INTEGER IREC_LEN
	INTEGER J_DIM
	INTEGER NCF_RD
	INTEGER ND_RD
	CHARACTER*20 FILE_DATE
C
C Constants for opacity etc [Set in CMFGEN].
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
C
C Parameters for fit to Electrons Scattering redistribution function
C (dipole form). From Rybicki and Hummer (A&A, 290,553)
C
	INTEGER NCOEF
	PARAMETER (NCOEF=2)
	REAL*8 ACOEF(2),BCOEF(2)
	DATA ACOEF/1.690703717290D0,-0.690703717290D0/
	DATA BCOEF/1.614249968779D0,2.154326524957D0/
C
	IF(FILE_IN .EQ. 'J PASSED VIA CALL')THEN
	  IF(J_SIZE .NE. NCF*ND)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in COMP_J_COV_V2: When J is passed in call'
	    WRITE(I,*)'J_SIZE must be the same as NCF*ND'
	    WRITE(I,*)'J_SIZE=',J_SIZE,'NCF*ND=',NCF*ND
	    STOP
	  END IF
	  J_DIM=NCF		!Must be set
	  GOTO 1000		!Jump the section reading in the data
	END IF
!
	IF(FILE_OUT .EQ. 'RETURN J VIA CALL')THEN
	  IF(J_SIZE .NE. NCF*ND)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Error in COMP_J_COV_V2: When J is returned in call'
	    WRITE(I,*)'J_SIZE must be the same as NCF*ND'
	    WRITE(I,*)'J_SIZE=',J_SIZE,'NCF*ND=',NCF*ND
	    STOP
	  END IF
	END IF
C
C It is assumed that we can use BA to store J. If this is not the case we
C allocate some extra storage.
C
	J_DIM=J_SIZE/ND
	IF(J_DIM .LT. NCF)THEN
	  ALLOCATE (EXTRA_J_ST(ND,J_DIM+1:NCF),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    I=ERROR_LU() 
	    WRITE(I,'(A)')'Error in COMP_J_CONV_V2'
	    WRITE(I,'(A)')'Unable to allocate memeory'
	    STOP
	  END IF
	ELSE
	  J_DIM=NCF
	END IF
C
C Open file with mean intensities. Read in J for all depths at all frequencies.
C
	INQUIRE(UNIT=LU_IN,OPENED=FILE_OPEN)
        IF(FILE_OPEN)INQUIRE(FILE=FILE_IN,OPENED=FILE_OPEN)
	IF(.NOT. FILE_OPEN)THEN
	  CALL READ_DIRECT_INFO_V3(I,IREC_LEN,FILE_DATE,FILE_IN,LU_IN,IOS)
	  IF(IOS .NE. 0)STOP
	  OPEN(UNIT=LU_IN,FILE=FILE_IN,STATUS='OLD',RECL=IREC_LEN,
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ')
	END IF
C
	INQUIRE(UNIT=LU_IN,RECL=IREC_LEN)
	READ(LU_IN,REC=CONT_REC)INIT_REC,NCF_RD,ND_RD
	IF(ND_RD .NE. ND)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in COMP_J_COV_V2'
	  WRITE(I,*)'ND and ND_RD must agree'
	  WRITE(I,*)'ND=',ND,'ND_RD=',ND_RD
	  STOP
	END IF
	IF(NCF_RD .GT. NCF)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Error in COMP_J_COV_V2'
	  WRITE(I,*)'NCF must be greater than NCF_RD'
	  WRITE(I,*)'NCF=',NCF,'NCF_RD=',NCF_RD
	  STOP
	END IF
C
	IF(RD_NU)THEN
	  DO ML=1,J_DIM
	    READ(LU_IN,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),NU(ML)
	  END DO
	  DO ML=J_DIM+1,NCF
	    READ(LU_IN,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),NU(ML)
	  END DO
	ELSE IF(ALLOW_UNEQUAL_FREQ)THEN
	  UNEQUAL_FREQ=.FALSE.
	  DO ML=1,NCF_RD
	    IF(ML .LE. J_DIM)THEN
	      READ(LU_IN,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),A(ML)
	    ELSE
	      READ(LU_IN,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),A(ML)
	    END IF
	    IF(A(ML) .NE. NU(ML))UNEQUAL_FREQ=.TRUE.
	  END DO
	  IF(NCF .NE. NCF_RD)UNEQUAL_FREQ=.TRUE.
	  IF(UNEQUAL_FREQ)THEN
	    IF(NU(1) .GT. A(1) .OR. NU(NCF) .LT. A(NCF_RD))THEN
	      I=ERROR_LU()
	      WRITE(I,*)'Insufficient range of frequencies in EDFFACTOR'
	      WRITE(I,*)'Routine is COMP_J_CONV'
	      STOP
	    END IF
	    DO K=1,ND
	      H(1:J_DIM)=J_STORE(K,1:J_DIM)
	      IF(NCF_RD .GT. J_DIM)H(J_DIM+1:NCF_RD)=J_STORE(K,J_DIM+1:NCF_RD)
	      CALL MON_INTERP(C,NCF,IONE,NU,NCF,H,NCF_RD,A,NCF_RD)
	      J_STORE(K,1:J_DIM)=C(1:J_DIM)
	      IF(NCF .GT. J_DIM)EXTRA_J_ST(K,J_DIM+1:NCF)=C(J_DIM+1:NCF)
	    END DO
	  END IF
	ELSE
	  DO ML=1,NCF
	    IF(ML .LE. J_DIM)THEN
	      READ(LU_IN,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),T1
	    ELSE
	      READ(LU_IN,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),T1
	    END IF
	    IF(T1 .NE. NU(ML))THEN
	      I=ERROR_LU()
	      WRITE(I,*)'Invalid frequency in ',FILE_OUT
	      WRITE(I,*)'Routine is COMP_J_CONV'
	      STOP
	    END IF
	  END DO
	END IF
	CLOSE(UNIT=LU_IN)
C
C*****************************************************************************
C*****************************************************************************
C
C All data has now been read in. We can now perform the convolutions.
c
1000	CONTINUE
C
C Compute those parts of the TRIDIAGONAL vectors which are independent of
C depth, and the fitting parameters.
C
	A_STORE(1)=0.0D0
	C_STORE(1)=-2.0D0/( LOG(NU(1)/NU(2)) )**2
	DO ML=2,NCF-1
	  D1=LOG(NU(ML-1)/NU(ML))
	  D2=LOG(NU(ML)/NU(ML+1))
	  DH=0.5D0*(D1+D2)
	  A_STORE(ML)=-1.0D0/D1/DH
	  C_STORE(ML)=-1.0D0/D2/DH
	END DO
	A_STORE(NCF)=-2.0D0/( LOG(NU(NCF-1)/NU(NCF)) )**2
	C_STORE(NCF)=0.0D0
C
C Compute the triadiagonal quantities for performing the convolution, and 
C perform the convolution. Due to the depth dependence of BETA, the vectors 
C are depth dependent. We convolve both J and the Planck Function.
C
!
! We can scale the Placnk function to avoid overflow.
!
	DO K=1,ND
C
	  BETA=1.84D-03*SQRT(TEMP(K))
	  T3=0.5D0*HDKT*NU(1)/TEMP(K)
	  IF(T3 .LT. 1)T3=0.0D0
	  IF(T3 .GT. 700.0D0)T3=700.0D0
	  DO ML=1,NCF
	    T1=EXP( -HDKT*NU(ML)/TEMP(K) )
	    T2=EXP( T3-HDKT*NU(ML)/TEMP(K) )
	    PLANCK_FN(ML)=TWOHCSQ*(NU(ML)**3)*T2/ (1.0D0-T1)
	  END DO
C
	  J_ES(:)=0.0D0
	  PLANCK_ES(:)=0.0D0
C
	  DO L=1,NCOEF
	    T1=BETA*BETA/BCOEF(L)/BCOEF(L)
	    A(:)=T1*A_STORE(:)			!Over frequency
	    H(:)=-1.0D0
	    C(:)=T1*C_STORE(:)
	    D(1:J_DIM)=J_STORE(K,1:J_DIM)
	    IF(J_DIM .LT. NCF)D(J_DIM+1:NCF)=EXTRA_J_ST(K,J_DIM+1:NCF)
	    CALL THOMAS_RH(A,H,C,D,NCF,IONE)
	    J_ES(:)=J_ES(:)+ACOEF(L)*D(:)
	    D(:)=PLANCK_FN(:)
	    CALL SIMPTH_RH(A,H,C,D,NCF,IONE)
	    PLANCK_ES(:)=PLANCK_ES(:)+ACOEF(L)*D(:)
	  END DO				!Fit parameter loop.
C
C Derive the wavelength shift so that the electron scattered planck function
C gives the Planck function. This is only done on the Wien side of the
C BB curve. At longer wavelengths we will use a simple scaling.
C
C Because of the steep variation of B on the Wien side, we operate on Log(B).
C
	  IF( MINVAL(PLANCK_FN) .GT. 0.0D0 .AND. MINVAL(PLANCK_ES) .GT. 0.0D0)THEN
	    PLANCK_ES(:)=LOG(PLANCK_ES(:))
	    PLANCK_FN(:)=LOG(PLANCK_FN(:))
C
	    PLANCK_NU(1)=NU(1)
	    I=2
	    ML=1
	    DO WHILE(NU(I) .GT. 1.5D0*TEMP(K))
	      DO WHILE(PLANCK_ES(I) .GT. PLANCK_FN(ML))
	        ML=ML+1
	      END DO
	      DO WHILE(PLANCK_ES(I) .LT. PLANCK_FN(ML-1))
	         ML=ML+1
	      END DO
	      T1=(PLANCK_ES(I)-PLANCK_FN(ML))/(PLANCK_FN(ML-1)-PLANCK_FN(ML))
              PLANCK_NU(I)=T1*NU(ML-1)+(1.0D0-T1)*NU(ML)
	      I=I+1
	    END DO
C
C Apply same wavelength shift as for our last (i.e. lowest)
C frequency as determined by matching the Planck function.
C
	    T1=1.0D0/NU(I-1)-1.0D0/PLANCK_NU(I-1)		!Wavelength shift
	    DO ML=I,NCF
	      PLANCK_NU(ML)=1.0D0/( 1.0D0/NU(ML) - T1)
	    END DO
C
C We now perform a simple linear interpolation of the electron scattered 
C Planck function back onto the old frequency grid. We apply the same
C interpolation to J_ES also. We use "A" as a temporary store for J_ES,
C and "C" as a temporary store for PLANCK_ES.
C
C No interpolation is done for the end points (which have coherent scattering).
C The immediate interior points, if necessary, are handled by extrapolation.
C
	    A(:)=DLOG(J_ES(:))
	    J_ES(1)=A(1)		!Must be in LOG form: not modified by interp.
	    J_ES(NCF)=A(NCF)
	    C(:)=PLANCK_ES(:)	!Already taken LOG
	    I=1
	    DO ML=2,NCF-1
	      DO WHILE (NU(ML) .LT. PLANCK_NU(I+1) .AND. I .LT. NCF-1)
	         I=I+1
	      END DO
	      T1=(NU(ML)-PLANCK_NU(I))/(PLANCK_NU(I+1)-PLANCK_NU(I))
	      PLANCK_ES(ML)=T1*C(I+1)+(1.0D0-T1)*C(I)
	      J_ES(ML)=T1*A(I+1)+(1.0D0-T1)*A(I)
	    END DO
C
C To remove any residual variations (primarily at low frequencies) we now
C scale J according the the ratio of the PLANCK functions.
C
	     J_ES(:)=EXP( J_ES(:)+ (PLANCK_FN(:)-PLANCK_ES(:)) )
!
	END IF
C                       
C Store the computed  J to be used to treat the electron scattering. We 
C overwrite J_STORE since the mean intensity is no longer required in this 
C routine.
C
	  J_STORE(K,1:J_DIM)=J_ES(1:J_DIM)
	  IF(J_DIM .LT. NCF)EXTRA_J_ST(K,J_DIM+1:NCF)=J_ES(J_DIM+1:NCF)
C
	END DO					!Depth
C
C 
C ****************************************************************************
C ****************************************************************************
C
	IF(FILE_OUT .NE. 'RETURN J VIA CALL')THEN
C
C Now output J convolution to data file. Data file will have exactly the same
C format as EDDFACTOR file (file that contains RJ).
C
	  INQUIRE(UNIT=LU_OUT,OPENED=FILE_OPEN)
          IF(FILE_OPEN)INQUIRE(FILE=FILE_OUT,OPENED=FILE_OPEN)
	  IF(.NOT. FILE_OPEN)THEN
	    CALL WRITE_DIRECT_INFO_V3(ND,IREC_LEN,'20-Aug-2000',FILE_OUT,LU_OUT)
	    OPEN(UNIT=LU_OUT,FILE=FILE_OUT,FORM='UNFORMATTED',
	1      ACCESS='DIRECT',STATUS='UNKNOWN',RECL=IREC_LEN)
	  END IF
	  WRITE(LU_OUT,REC=CONT_REC)INIT_REC,NCF,ND
	  DO ML=1,J_DIM
	    WRITE(LU_OUT,REC=INIT_REC+ML-1)(J_STORE(K,ML),K=1,ND),NU(ML)
	  END DO
	  DO ML=J_DIM+1,NCF
	    WRITE(LU_OUT,REC=INIT_REC+ML-1)(EXTRA_J_ST(K,ML),K=1,ND),NU(ML)
	  END DO
	  CLOSE(UNIT=LU_OUT)
C
	  J_STORE(:,1:J_DIM)=0.0D0
	  IF(J_DIM .NE. NCF)DEALLOCATE (EXTRA_J_ST)
	END IF
C 
	RETURN
	END
