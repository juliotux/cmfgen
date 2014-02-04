C
C Routine to plot J from EDDFACTOR file. This J is convolved with the 
C electron redistribution function using the 1 or 2-parameter formulation of 
C Hummer and Rybicki. For comparison, RJ_ES from the ES_J_CONV file, 
C may also be plotted. B may also be convolved with the electron 
C redistribution function.
C
	PROGRAM PLT_RJ
!
! Altered 16-Jun-2000 : Now use DIRECT_INFO files to get file record lengths.
!                       Using V2 of CNVLV routines/
!                       
! Interface routines for IO routines.
!
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE GEN_IN_INTERFACE
!
	IMPLICIT NONE
C
	REAL*8, ALLOCATABLE :: NU(:)
	REAL*8, ALLOCATABLE :: RJ(:,:)
C
	REAL*8, ALLOCATABLE :: RJ_ES_RD(:,:)
C
	REAL*8, ALLOCATABLE :: FLUX_RJ(:)
	REAL*8, ALLOCATABLE :: FLUX_ES(:)
C
	REAL*8, ALLOCATABLE :: PLANCK_FN(:)
C
	REAL*8, ALLOCATABLE :: A(:)
	REAL*8, ALLOCATABLE :: B(:)
	REAL*8, ALLOCATABLE :: C(:)
	REAL*8, ALLOCATABLE :: D(:)
	REAL*8, ALLOCATABLE :: J_ES(:)
	REAL*8, ALLOCATABLE :: RJ_TMP(:)
C
	REAL*4, ALLOCATABLE :: XV(:)
	REAL*4, ALLOCATABLE :: YV(:)
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_OUT=6
C
	REAL*8 RJ_FLUX,ES_FLUX,ES2_FLUX
C
	REAL*8 T_ELEC
	REAL*8 T1
	REAL*8 C_KMS,NU_0
	LOGICAL ONE_PAR
	LOGICAL SET_CONST
C
	INTEGER I,J,K,ML
	INTEGER DEPTH_VAR
	INTEGER ND,NCF
	INTEGER ND2,NCF2
	INTEGER LU_IN
	INTEGER IOS
	INTEGER REC_LENGTH
!
	CHARACTER*80 NAME		!Default title for plot
	CHARACTER*80 XAXIS,XAXSAV	!Label for Absisca
	CHARACTER*80 YAXIS		!Label for Ordinate
C
C USR_OPTION variables
C
	CHARACTER MAIN_OPT_STR*80	!Used for input of the main otion
	CHARACTER X*10			!Used for the idividual option
	CHARACTER*120 DEFAULT
	CHARACTER*120 DESCRIPTION
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	REAL*8 FUN_PI,SPEED_OF_LIGHT
	LOGICAL EQUAL
	CHARACTER*30 UC
	EXTERNAL FUN_PI,SPEED_OF_LIGHT,EQUAL,UC
C                       
	CHARACTER*80 FILENAME
	CHARACTER*80 FILE_DATE
C
C Set constants.
C
	CHIBF=2.815D-06
	CHIFF=3.69D-29
	HDKT=4.7994145D0
	TWOHCSQ=0.0147452575D0
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
	C_KMS=SPEED_OF_LIGHT()/1.0D+05
C
	WRITE(T_OUT,*)
	1  ' Routine to plot J from EDDFACTOR file. This J is convolved'
	WRITE(T_OUT,*)
	1  ' with the electron redistrbution function using the 1-parameter'
	WRITE(T_OUT,*)
	1  ' formulation of Hummer and Rybicki. For comparison, RJ_ES from'
	WRITE(T_OUT,*)
	1  ' the ES_J_CONV file, may also be plotted.'
!
!	CALL GEN_IN(ND,'Number of model depth points')
!	CALL GEN_IN(NCF,'Number of model frequency points')
!
	LU_IN=35
	T_elec=1.0D0		!10^4 K
C
C Open EDDFACTOR file.
C
5	FILENAME='EDDFACTOR'
	CALL GEN_IN(FILENAME,'Filename with RJ')
        CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',RECL=REC_LENGTH,
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ',ERR=5)
	  READ(LU_IN,REC=3)I,NCF,ND
	  ALLOCATE ( NU(NCF) )
	  ALLOCATE ( RJ(ND,NCF) )
	  DO ML=1,NCF
	    READ(LU_IN,REC=I+ML-1)(RJ(K,ML),K=1,ND),NU(ML)
	  END DO
	CLOSE(LU_IN)
!
! Allocate remaining arrays.
!
	ALLOCATE ( RJ_ES_RD(ND,NCF) )
C
	ALLOCATE ( FLUX_RJ(ND) )
	ALLOCATE ( FLUX_ES(ND) )
!
	ALLOCATE ( RJ_TMP(NCF) )
	ALLOCATE ( J_ES(NCF) )
C
	ALLOCATE ( PLANCK_FN(NCF) )
	ALLOCATE ( A(NCF) )
	ALLOCATE ( B(NCF) )
	ALLOCATE ( C(NCF) )
	ALLOCATE ( D(NCF) )
C
	ALLOCATE ( XV(NCF) )
	ALLOCATE ( YV(NCF) )
!
!
10	FILENAME='ES_J_CONV'
	CALL GEN_IN(FILENAME,'Filename with RJ_ES')
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',RECL=REC_LENGTH,
	1      ACCESS='DIRECT',FORM='UNFORMATTED',ACTION='READ',ERR=10)
       	  READ(LU_IN,REC=3)I,NCF2,ND2
          IF(ND2 .NE. ND .OR. NCF2 .NE. NCF)THEN
            WRITE(T_OUT,*)'Error: ND and NCF not the same'
            WRITE(T_OUT,*)'ND=',ND,'NCF=',NCF
            WRITE(T_OUT,*)'ND2=',ND2,'NCF2=',NCF2
            STOP
          END IF
	  DO ML=1,NCF
	    READ(LU_IN,REC=I+ML-1)(RJ_ES_RD(K,ML),K=1,ND),T1
	  END DO
	CLOSE(LU_IN)
C
C Compute integrals as a function of depth to check flux conservation.
C
	FLUX_RJ(1:ND)=+0.0D0
	FLUX_ES(1:ND)=+0.0D0
	DO ML=2,NCF-1
	  DO I=1,ND
	    FLUX_RJ(I)=FLUX_RJ(I)+(NU(ML)-NU(ML+1))*(RJ(I,ML)+
	1                  RJ(I,ML+1))
	    FLUX_ES(I)=FLUX_ES(I)+(NU(ML)-NU(ML+1))*(RJ_ES_RD(I,ML)+
	1                  RJ_ES_RD(I,ML+1))
	  END DO
	END DO
	WRITE(9,*)'% Flux error '
	WRITE(9,'(1P,5E12.5)')200.0D0*(FLUX_RJ(1:ND)-FLUX_ES(1:ND))/
	1                    (FLUX_RJ(1:ND)+FLUX_ES(1:ND))
C
	FLUX_RJ(1:ND)=+0.0D0
	FLUX_ES(1:ND)=+0.0D0
	DO ML=2,NCF-1
	  DO I=1,ND
	    T1=LOG(NU(ML)/NU(ML+1))
	    FLUX_RJ(I)=FLUX_RJ(I)+T1*(RJ(I,ML)+RJ(I,ML+1))
	    FLUX_ES(I)=FLUX_ES(I)+T1*(RJ_ES_RD(I,ML)+RJ_ES_RD(I,ML+1))
	  END DO
	END DO
	WRITE(9,*)' '
	WRITE(9,*)' '
	WRITE(9,*)'% Photon error '
	WRITE(9,'(1P,5E12.5)')200.0D0*(FLUX_RJ(1:ND)-FLUX_ES(1:ND))/
	1                    (FLUX_RJ(1:ND)+FLUX_ES(1:ND))
C
	DEPTH_VAR=ND/2
	XAXIS='\gn(10\u15\dHz)'
	YAXIS=' '
	XAXSAV=XAXIS
	XV(1:NCF)=NU(1:NCF)
	NAME=' '
C
C 
C
C This message will only be printed once
C
	WRITE(T_OUT,*)
	WRITE(T_OUT,"(8X,A)")'(default is to write file '//
	1    'main_option.sve)'
	WRITE(T_OUT,"(8X,A)")'(append sve=filename to '//
	1    'write a new .sve file)'
	WRITE(T_OUT,"(8X,A)")'(box=filename to '//
	1    'write a .box file containing several .sve files)'
	WRITE(T_OUT,"(8X,A)")'(.filename to read .sve file)'
	WRITE(T_OUT,"(8X,A)")'(#filename to read .box file)'
	WRITE(T_OUT,*)
C
C This call resets the .sve algorithm.  Specifically it sets the next
C input answer to be a main option, and all subsequent inputs to be
C sub-options.  
C
 3	CALL SVE_FILE('RESET')
C
	MAIN_OPT_STR='  '
	DEFAULT='GR'
	DESCRIPTION=' '					!Obvious main option
	CALL USR_OPTION(MAIN_OPT_STR,'OPTION',DEFAULT,DESCRIPTION)
C
C   If the main option begins with a '.', a previously
C   written .sve file is read.
C
C   If the main option begins with a '#', a previously
C   written .box file is read.
C
C   If sve= is apended to the end of this main option, a new .sve file
C   is opened with the given name and the main option and all subsequent
C   sub-options are written to this file.
C
C   If box= is input then a .box file is created, which contains the name
C   of several .sve files to process.
C
C   If only a main option is given, the option and subsequent sub-options
C   are saved in a file called 'main option.sve'.  All following main
C   options are saved in separate files.
C
	X=UC(MAIN_OPT_STR)
	I=INDEX(X,'(')
	IF(I .NE. 0)X=X(1:I-1)		!Remove line variables
C
C 
C
	IF(X(1:3) .EQ. 'TIT')THEN
	  CALL USR_OPTION(NAME,'Title',' ',' ')
	ELSE IF(X(1:2) .EQ. 'CJ')THEN
	  DEFAULT=WR_STRING(DEPTH_VAR)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Depth for plotting')
	  DEPTH_VAR=K
	  CALL USR_OPTION(ONE_PAR,'ONE_PAR','.TRUE.',
	1              'Do a one parameter convolution?')
C
	  NU_0=0
	  CALL USR_HIDDEN(NU_0,'NU0','0','Frequency for velocity plot')
	  DO ML=1,NCF
	    XV(ML)=NU(ML)
	    YV(ML)=LOG10(RJ(K,ML))
	  END DO
	  IF(NU_0 .NE. 0)XV(1:NCF)=C_KMS*(NU(1:NCF)-NU_0)/NU_0
	  CALL CURVE(NCF,XV,YV)
C
C Convolve RJ with electron-scattering redistribution function.
C
	  DEFAULT=WR_STRING(T_ELEC)
	  CALL USR_OPTION(T_ELEC,'TEMP',DEFAULT,
	1              'Electron temperature (in units 10^4 K)')
	  RJ_TMP(1:NCF)=RJ(K,1:NCF)
	  IF(ONE_PAR)THEN
	    CALL CNVLV_ES_ONE_PAR_V2(NU,RJ_TMP,J_ES,T_elec,T_OUT,NCF)
	  ELSE
	    CALL CNVLV_ES_TWO_PAR_V2(NU,RJ_TMP,J_ES,T_elec,T_OUT,NCF)
	  END IF
	  YV(1:NCF)=DLOG10(J_ES(1:NCF))
	  CALL CURVE(NCF,XV,YV)
!
! Compute integrals as a function of depth to check flux conservation.
!
	  RJ_FLUX=+0.0D0
	  ES_FLUX=+0.0D0
	  ES2_FLUX=+0.0D0
	  DO ML=2,NCF-1
	    RJ_FLUX=RJ_FLUX+(NU(ML)-NU(ML+1))*(RJ(K,ML)+RJ(K,ML+1))
	    ES_FLUX=ES_FLUX+(NU(ML)-NU(ML+1))*(RJ_ES_RD(K,ML)+RJ_ES_RD(K,ML+1))
	    ES2_FLUX=ES2_FLUX+(NU(ML)-NU(ML+1))*(J_ES(ML)+J_ES(ML+1))
	  END DO
	  ES_FLUX=200.0D0*(RJ_FLUX-ES_FLUX)/(RJ_FLUX+ES_FLUX)
	  ES2_FLUX=200.0D0*(RJ_FLUX-ES2_FLUX)/(RJ_FLUX+ES2_FLUX)
	  WRITE(5,'(A,1X,1P,2E13.5)')'  %Flux errors:  ',ES_FLUX,ES2_FLUX
C
	  RJ_FLUX=+0.0D0
	  ES_FLUX=+0.0D0
	  ES2_FLUX=+0.0D0
	  DO ML=2,NCF-1
	    T1=LOG(NU(ML)/NU(ML+1))
	    RJ_FLUX=RJ_FLUX+T1*(RJ(K,ML)+RJ(K,ML+1))
	    ES_FLUX=ES_FLUX+T1*(RJ_ES_RD(K,ML)+RJ_ES_RD(K,ML+1))
	    ES2_FLUX=ES2_FLUX+T1*(J_ES(ML)+J_ES(ML+1))
	  END DO
	  ES_FLUX=200.0D0*(RJ_FLUX-ES_FLUX)/(RJ_FLUX+ES_FLUX)
	  ES2_FLUX=200.0D0*(RJ_FLUX-ES2_FLUX)/(RJ_FLUX+ES2_FLUX)
	  WRITE(5,'(A,1X,1P,2E13.5)')'  %Photon errors:',ES_FLUX,ES2_FLUX
!
	ELSE IF(X(1:2) .EQ. 'C%')THEN
	  DEFAULT=WR_STRING(DEPTH_VAR)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Depth for plotting')
	  DEPTH_VAR=K
	  CALL USR_OPTION(ONE_PAR,'ONE_PAR','.TRUE.',
	1              'Do a one parameter convolution?')
C
	  XV(1:NCF)=NU(1:NCF)
C
C Convolve RJ with electron-scattering redistribution function.
C
	  DEFAULT=WR_STRING(T_ELEC)
	  CALL USR_OPTION(T_ELEC,'TEMP',DEFAULT,
	1              'Electron temperature (in units 10^4 K)')
	  RJ_TMP(1:NCF)=RJ(K,1:NCF)
	  IF(ONE_PAR)THEN
	    CALL CNVLV_ES_ONE_PAR_V2(NU,RJ_TMP,J_ES,T_elec,T_OUT,NCF)
	  ELSE
	    CALL CNVLV_ES_TWO_PAR_V2(NU,RJ_TMP,J_ES,T_elec,T_OUT,NCF)
	  END IF
!
	  A(:)=J_ES(:)
	  CALL ADJUST_JES_V2(NU,J_ES,NCF,T_ELEC,ONE_PAR,' ',T_OUT)
!
! Compute integrals as a function of depth to check flux conservation.
!
	  RJ_FLUX=+0.0D0
	  ES_FLUX=+0.0D0
	  ES2_FLUX=+0.0D0
	  DO ML=2,NCF-1
	    RJ_FLUX=RJ_FLUX+(NU(ML)-NU(ML+1))*(RJ_TMP(ML)+RJ_TMP(ML+1))
	    ES_FLUX=ES_FLUX+(NU(ML)-NU(ML+1))*(A(ML)+A(ML+1))
	    ES2_FLUX=ES2_FLUX+(NU(ML)-NU(ML+1))*(J_ES(ML)+J_ES(ML+1))
	  END DO
	  ES_FLUX=200.0D0*(RJ_FLUX-ES_FLUX)/(RJ_FLUX+ES_FLUX)
	  ES2_FLUX=200.0D0*(RJ_FLUX-ES2_FLUX)/(RJ_FLUX+ES2_FLUX)
	  WRITE(5,'(A,1X,1P,2E13.5)')'  %Flux errors:  ',ES_FLUX,ES2_FLUX
C
	  RJ_FLUX=+0.0D0
	  ES_FLUX=+0.0D0
	  ES2_FLUX=+0.0D0
	  DO ML=2,NCF-1
	      T1=LOG(NU(ML)/NU(ML+1))
	      RJ_FLUX=RJ_FLUX+T1*(RJ_TMP(ML)+RJ_TMP(ML+1))
	      ES_FLUX=ES_FLUX+T1*(A(ML)+A(ML+1))
	      ES2_FLUX=ES2_FLUX+T1*(J_ES(ML)+J_ES(ML+1))
	  END DO
	  ES_FLUX=200.0D0*(RJ_FLUX-ES_FLUX)/(RJ_FLUX+ES_FLUX)
	  ES2_FLUX=200.0D0*(RJ_FLUX-ES2_FLUX)/(RJ_FLUX+ES2_FLUX)
	  WRITE(5,'(A,1X,1P,2E13.5)')'  %Photon errors:',ES_FLUX,ES2_FLUX
!
	  CALL USR_OPTION(ONE_PAR,'FLUX_PLOT','.TRUE.','FLux plot?')
!
	  IF(ONE_PAR)THEN
	    C(:)=0.0D0
	    D(:)=0.0D0
	    RJ_FLUX=0.0D0
	    DO ML=2,NCF-1
	      RJ_FLUX=RJ_FLUX+(NU(ML)-NU(ML+1))*(RJ_TMP(ML)+RJ_TMP(ML+1))
	      C(ML)=C(ML-1)+(NU(ML)-NU(ML+1))*
	1              ( (RJ_TMP(ML)-A(ML)) +(RJ_TMP(ML+1)-A(ML+1)) )
	      D(ML)=D(ML-1)+(NU(ML)-NU(ML+1))*
	1              ( (RJ_TMP(ML)-J_ES(ML)) +(RJ_TMP(ML+1)-J_ES(ML+1)) )
	    END DO
	    YV(:)=100.0D0*C(:)/RJ_FLUX
	    CALL CURVE(NCF,XV,YV)
	    YV(:)=100.0D0*D(:)/RJ_FLUX
	    CALL CURVE(NCF,XV,YV)
	    YAXIS='%Flux Diff'
	  ELSE
	    I=NCF-1
	    YV(:)=100.0D0*(A(:)-RJ_TMP(:))/RJ_TMP(:)
	    CALL CURVE(I,XV,YV)
	    YV(:)=100.0D0*(J_ES(:)-RJ_TMP(:))/RJ_TMP(:)
	    CALL CURVE(I,XV,YV)
	    YAXIS='%Diff'
	  END IF
!
	ELSE IF(X(1:2) .EQ. 'PC')THEN
!
	  DEFAULT=WR_STRING(DEPTH_VAR)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Depth for plotting')
	  DEPTH_VAR=K
	  YV(1:NCF)=DLOG10( RJ_ES_RD(K,1:NCF) )
	  CALL CURVE(NCF,XV,YV)
	  YAXIS='Log (RJ_ES)'
!
	ELSE IF(X(1:2) .EQ. 'FL')THEN
!
	  DEFAULT=WR_STRING(DEPTH_VAR)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Depth for plotting')
	  DEPTH_VAR=K
	  RJ_TMP(:)=RJ(K,:)
	  C(:)=0.0D0
	  DO ML=2,NCF-1
	    C(ML)=C(ML-1)+(NU(ML)-NU(ML+1))*( RJ_TMP(ML)+RJ_TMP(ML+1))
	  END DO
	  I=NCF-1
	  T1=C(I)
	  YV(:)=100.0D0*C(:)/T1
	  XV(:)=NU(:)
	  CALL CURVE(I,XV,YV)
	  YAXIS='%Flux'
!
	ELSE IF(X(1:1) .EQ. '%')THEN
!
	  DEFAULT=WR_STRING(DEPTH_VAR)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Depth for plotting')
	  DEPTH_VAR=K
	  DO ML=1,NCF
	     YV(ML)=100.0D0*(RJ(K,ML)-RJ_ES_RD(K,ML))/RJ(K,ML)
	  END DO
	  CALL CURVE(NCF,XV,YV)
	  YAXIS='%Difference'
!
	ELSE IF(X(1:2) .EQ. 'B%')THEN
!
	  DEFAULT=WR_STRING(DEPTH_VAR)
	  CALL USR_OPTION(K,'DEPTH',DEFAULT,'Depth for plotting')
	  DEFAULT=WR_STRING(T_ELEC)
	  CALL USR_OPTION(T_ELEC,'TEMP',DEFAULT,
	1              'Electron temperature (in units 10^4 K)')
!
	  PLANCK_FN(1:NCF)=TWOHCSQ*NU(1:NCF)**3/
	1                      (EXP(HDKT*NU(1:NCF)/T_ELEC)-1)
	  DEPTH_VAR=K
	  DO ML=1,NCF
	     YV(ML)=100.0D0*(RJ(K,ML)-PLANCK_FN(ML))/RJ(K,ML)
	  END DO
	  CALL CURVE(NCF,XV,YV)
	  YAXIS='%Difference'
!
	ELSE IF(X(1:2) .EQ. 'CB')THEN
!
	  DEFAULT=WR_STRING(T_ELEC)
	  CALL USR_OPTION(T_ELEC,'TEMP',DEFAULT,
	1              'Electron temperature (in units 10^4 K)')
	  CALL USR_OPTION(ONE_PAR,'ONE_PAR','.TRUE.',
	1              'Do a one parameter convolution?')
!
	  CALL USR_HIDDEN(SET_CONST,'CONST','.FALSE.',
	1              'Set BB to unity?')
!
! Convolve Planck function and plot
!
	  IF(SET_CONST)THEN
	    PLANCK_FN(1:NCF)=1.0D0
	  ELSE
	    PLANCK_FN(1:NCF)=NU(1:NCF)**3/(EXP(4.7994*NU(1:NCF)/T_ELEC)-1)
	  END IF
	  IF(ONE_PAR)THEN
	    CALL CNVLV_ES_ONE_PAR_V2(NU,PLANCK_FN,D,T_elec,T_OUT,NCF)
	  ELSE
	    CALL CNVLV_ES_TWO_PAR_V2(NU,PLANCK_FN,D,T_elec,T_OUT,NCF)
	  END IF
	  PLANCK_FN=LOG10(PLANCK_FN)
	  D=LOG10(D)
C
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=2302.5D0*(D(1:NCF)-PLANCK_FN(1:NCF))	!% Error
	  CALL CURVE(NCF,XV,YV)
	  IF(SET_CONST)GOTO 1
!
! Shift and scale e.s. funtion so still have Placnk function
!
	  A(1)=NU(1)
	  I=2
	  ML=1
	  DO WHILE(NU(I) .GT. 1.5D0*T_ELEC)
	    DO WHILE(D(I) .GT. PLANCK_FN(ML))
	      ML=ML+1
	    END DO
	    DO WHILE(D(I) .LT. PLANCK_FN(ML-1))
	      ML=ML+1
	    END DO
	    T1=(D(I)-PLANCK_FN(ML))/(PLANCK_FN(ML-1)-PLANCK_FN(ML))
            A(I)=T1*NU(ML-1)+(1.0D0-T1)*NU(ML)
	    I=I+1
	  END DO
C
C	  T1=1/NU(I-1)-1/A(I-1)		!Wavelength shift
	  DO J=I,NCF
C	    A(J)=1.0D0/( 1.0D0/NU(J) - T1)
	    A(J)=NU(J)
	  END DO
C
C Perfrom a simple linear interpolation back onto the old frequency grid.
C
	  B(1)=D(1)
	  I=1
	  DO ML=2,NCF-2
	    DO WHILE (NU(ML) .LT. A(I+1))
	     I=I+1
	    END DO
	    T1=(NU(ML)-A(I))/(A(I+1)-A(I))
	    B(ML)=T1*D(I+1)+(1.0D0-T1)*D(I)
	  END DO
	  B(NCF-1:NCF)=D(NCF-1:NCF)
C
	  YV(1:NCF)=2302.5D0*(B(1:NCF)-PLANCK_FN(1:NCF))	!% Error
	  CALL CURVE(NCF,XV,YV)
	  YAXIS='%Diff'
C
C Plot section:
C
	ELSE IF(X(1:2) .EQ. 'GR') THEN
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,' ')
	  XAXIS=XAXSAV
C
	ELSE IF(X(1:4) .EQ.'GRNL') THEN
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  XAXIS=XAXSAV
C
C 
C
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:4) .EQ. 'LIST')THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'TIT  : Set default plot title'
	  WRITE(T_OUT,*)'%    : Plot (J_RD-Jes_RED)/J_RD as a %'
	  WRITE(T_OUT,*)'PC   : Plot Jes_RD'
	  WRITE(T_OUT,*)'CJ   : Convolve J with e.s. r. func. and plot'
	  WRITE(T_OUT,*)'C%   : Plot (Jed-J_RD)/J_RD as a %'
	  WRITE(T_OUT,*)'FL   : Plot Int J_RD dv'
	  WRITE(T_OUT,*)'CB   : Convolve B with e.s. r. func.'
	  WRITE(T_OUT,*)'B%   : Plot (B-Bes)/B as a %'
!
	ELSE IF(X(1:2) .EQ. 'HE' .OR. X(1:4) .EQ. 'HELP')THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'TIT     CJ       C%      PC'
	  WRITE(T_OUT,*)'FL       %       B%      CB'
C
	ELSE IF(X(1:2) .EQ. 'EX') THEN
	  CALL CURVE(0,XV,YV)
	  STOP
	ELSE IF(X(1:3) .EQ. 'BOX') THEN
	  CALL WR_BOX_FILE(MAIN_OPT_STR)
	ELSE
	  PRINT*,'OPTION REQUESTED DOES NOT EXIST'
	END IF
C
1	CONTINUE
	GO TO 3
C
	END
C
C Subroutine to solve a tridiagonal system of N1 simultaneous
C Equations which are tridiagonal in nature for N2 R.H.sides.
C
C The equations to be solved are assumed to have the form
C
C    A(i).X(i-1) - [H(i)+A(i)+B(i)].X(i+1) - C(i).X(i+1)
C
	SUBROUTINE THOMAS_RH(A,H,C,D,N1,N2)
	IMPLICIT NONE
C
	INTEGER N1,N2
	REAL*8 A(N1),H(N1),C(N1),D(N1,N2)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
	INTEGER I,J
	REAL*8 DIV(N1)
C
C
C Compute quantities that will be used repeatedly if the same tridiagonal
C system is used for many R.H. Sides.
C
	C(N1)=0				!As used.
	DIV(1)=1.0/(C(1)+H(1))
	C(1)=C(1)*DIV(1)
	H(1)=H(1)*DIV(1)
	DO I=2,N1
	  DIV(I)=1.0D0/(A(I)*H(I-1)+H(I)+C(I))
	  H(I)=(A(I)*H(I-1)+H(I))*DIV(I)
	  C(I)=C(I)*DIV(I)
	END DO
C
C Entry for Thomas algorithim when H,C have been previously modified.
C
	ENTRY SIMPTH_RH(A,H,C,D,N1,N2)
C
	DO J=1,N2
	  D(1,J)=D(1,J)*DIV(1)
	  DO I=2,N1
	    D(I,J)=(D(I,J)+A(I)*D(I-1,J))*DIV(I)
	  END DO
	END DO
C
	DO J=1,N2
	  D(N1,J)=-D(N1,J)
	  DO I=N1-1,1,-1
	    D(I,J)=C(I)*D(I+1,J)-D(I,J)
	  END DO
	END DO
C
	RETURN
	END
