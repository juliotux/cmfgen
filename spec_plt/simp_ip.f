!
! Simple routine for outputing, in asci format, I(p) as a function of p.
! Based in IP_DATA.
!
! I(p) can either be output at a single frequency, of averaged over a
! band of freqencies.
!
!Output:
!       IP_TXT_DATA
!       SPECTRUM (obtained by integrating I(p) over p).
!
	PROGRAM SIMP_IP
	IMPLICIT NONE
!
	INTEGER NCF
	INTEGER ND
	INTEGER NC
	INTEGER NP
	REAL*8, ALLOCATABLE :: IP(:,:)
	REAL*8, ALLOCATABLE :: NU(:)
	REAL*8, ALLOCATABLE :: P(:)
!
! Vectors for passing data to plot package via calls to CURVE.
!
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
	REAL*8, ALLOCATABLE :: ZV(:)
!
	REAL*8 ANG_TO_HZ
	CHARACTER*80 FILENAME
	CHARACTER*80 FILE_DATE
!
! Miscellaneous variables.
!
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,K,L,ML,LS,NX
	INTEGER ST_REC
	INTEGER REC_LENGTH
	REAL*8 T1,T2
	REAL*8 DISTANCE,PI,PARSEC
	REAL*8 LAM_ST,LAM_END
	INTEGER INDX_ST,INDX_END
!
	INTEGER, PARAMETER :: T_IN=5		!For file I/O
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LU_IN=10	!For file I/O
	INTEGER, PARAMETER :: LU_OUT=11
!
	INTEGER GET_INDX_DP
	CHARACTER STRING*80
!
	EXTERNAL GET_INDX_DP
!
	ANG_TO_HZ=2.99792458D+10*1.0D-07  	!10^8/10^15
!
!  Read in model with I(p). P is in units of 10^10 cm, NU in units
! of 10^15 Hz, and I(p) is in cgs units.
!
	FILENAME='IP_DATA'
	WRITE(T_OUT,100,ADVANCE='NO')'File with IP Data',TRIM(FILENAME)
100     FORMAT(1X,A,' [',A,']: ')
        READ(T_IN,'(A)')FILENAME
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read IP_DATA_INFO'
	  STOP
	END IF
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED')
	  READ(LU_IN,REC=3)ST_REC,NCF,NP
	  ALLOCATE (IP(NP,NCF))
	  ALLOCATE (P(NP))
	  ALLOCATE (NU(NCF))
	  IF( INDEX(FILE_DATE,'20-Aug-2000') .NE. 0)THEN
	    READ(LU_IN,REC=ST_REC)(P(I),I=1,NP)
	    ST_REC=ST_REC+1
	  ELSE
	    WRITE(T_OUT,*)'Unrecognized date when reading IP_DATA' 
	    WRITE(T_OUT,*)'Date=',FILE_DATE
	    STOP
	  END IF
	  DO ML=1,NCF
	    READ(LU_IN,REC=ST_REC+ML-1)(IP(I,ML),I=1,NP),NU(ML)
	  END DO
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in IP_DATA file as MODEL A (default)'
!
	WRITE(6,*)' '	
	WRITE(6,*)'Can either enter a single wavelngth (second wavelength 0)  or a band.'	
	WRITE(6,*)'Wavelength is in vacuum'
	WRITE(6,*)' '	
!
	WRITE(T_OUT,'(1X,A)',ADVANCE='NO')'Lam start(A):'
        READ(T_IN,*)LAM_ST
	WRITE(T_OUT,'(1X,A)',ADVANCE='NO')'Lam end(A):'
        READ(T_IN,*)LAM_END
!
	T1=ANG_TO_HZ/LAM_ST
        I=GET_INDX_DP(T1,NU,NCF)
	IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
	INDX_ST=I
!
	IF(LAM_END .GT. LAM_ST)THEN
	  T1=ANG_TO_HZ/LAM_END
          J=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(J)-T1 .GT. T1-NU(J+1))J=J+1
	  INDX_END=J
	END IF
!
	ALLOCATE (XV(NP))
	ALLOCATE (YV(NP))
!
! In this case we return linear axes --- usefule for SN.
!
	 XV(1:NP)=P(1:NP)*1.0D+10          !in cm
!
	IF(LAM_END .GT. LAM_ST)THEN
!
! Average I(p) over frequnecy.
!
	  YV(:)=0.0D0
	  K=MIN(INDX_ST,INDX_END); J=MAX(INDX_ST,INDX_END); I=K
	  IF(J .EQ. I)J=I+1
	  DO K=I,J-1
	    YV(1:NP-1)=YV(1:NP-1)+0.5D0*(IP(1:NP-1,K)+IP(1:NP-1,K+1))*
	1                                 (NU(K)-NU(K+1))
	  END DO
	  T1=ABS(NU(I)-NU(J))
	  YV(1:NP-1)=YV(1:NP-1)/T1
!
	ELSE
	  LAM_END=LAM_ST
	  YV(1:NP)=IP(1:NP,I)
	END IF
!
! Ouput I as a function of p. If LAM_ST=LAM_END, this is output at a single
! frequency. Otherwise I(p,nu) has been averaged over the frequency band.
!
	OPEN(UNIT=12,FILE='IP_TXT_DATA',ACTION='WRITE',STATUS='UNKNOWN')
	  WRITE(12,*)'Lambda start:',LAM_ST
	  WRITE(12,*)'Lambda end:',LAM_END
	  WRITE(12,*)'Number of data values:',NP
	  WRITE(12,*)'Units of p are cm'
	  WRITE(12,*)'Units of I(p) are I(ergs cm\u-1\d s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	  DO I=1,NP
	    WRITE(12,*)I,P(I),YV(I)
	  END DO
	CLOSE(UNIT=12)
!
! Compute spectrum over the interval LAM_ST to LAM_END.
!
	IF(LAM_END .GT. LAM_ST)THEN
	  NX=INDX_END-INDX_ST+1
	  ALLOCATE (YV(NX))
          ZV(1:NX)=0.0D0
          DO I=1,NCF
	    XV(I)=ANG_TO_HZ/NU(INDX_ST+I-1)
	    DO J=1,NP-1
              ZV(I)=ZV(I)+0.5D0*(IP(J,INDX_ST+I-1)*P(J)+IP(J+1,INDX_ST+I-1)*
	1             P(J+1))*(P(J+1)-P(J))
            END DO
          END DO
!
	  DISTANCE=2.3                           !kpc
	  PARSEC=3.0856D+18                      !cm
	  PI=4.0D0*ATAN(1.0D0)
	  T1=DISTANCE*1.0E+03*PARSEC
          T1=2.0D0*PI*1.0D+23*(1.0E+10/T1)**2
          ZV=ZV*T1
!
	  OPEN(UNIT=12,FILE='SPECTRUM',ACTION='WRITE',STATUS='UNKNOWN')
	  DO I=1,NX
	    WRITE(12,*)I,P(I),YV(I)
	  END DO
	END IF
!
	STOP
	END
