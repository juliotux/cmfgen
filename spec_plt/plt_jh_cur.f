!
! General routine for plotting and comparing J & H  obtained from the 
! JH_AT_CURRENT_TIME data file.
!
! Various options are available to redden and normalize the model spectra. 
! Several different units can be used for the X and Y axes.
!
	PROGRAM PLT_JH_CUR
!
! Altered 17-Feb-2015   : Code plots r^2H at midpoint of R (as defined).
!                         [OSPREY cur/cmf: 17-Jan-2015] 
! Altered  8-Oct-2011   : Improved plot labeling, pass RSQ? to DP_CNVRT routine.
! Created  5-April-2009 : Based on PLT_JH (plots EDDFACTOR).
! Interface routines for IO routines.
!
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE GEN_IN_INTERFACE
	USE MOD_COLOR_PEN_DEF
	USE READ_KEYWORD_INTERFACE
!
	IMPLICIT NONE
!
	TYPE MODEL_INTENSITY
	  INTEGER NCF
	  INTEGER ND
	  REAL*8, POINTER :: RJ(:,:)
	  REAL*8, POINTER :: HFLUX(:,:)
	  REAL*8, POINTER :: JGREY(:)
	  REAL*8, POINTER :: HGREY(:)
	  REAL*8, POINTER :: NU(:)
	  REAL*8, POINTER :: R(:)
	  REAL*8, POINTER :: V(:)
	  REAL*8, POINTER :: H_INBC(:)
	  REAL*8, POINTER :: H_OUTBC(:)
	  CHARACTER*10 DATA_TYPE
	  CHARACTER*40 FILE_DATE
	  CHARACTER*80 FILENAME
	END TYPE MODEL_INTENSITY
	TYPE (MODEL_INTENSITY) ZM(5)
	INTEGER ND,NCF
	INTEGER NUM_FILES
	INTEGER ID
	INTEGER XV_LENGTH
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: TAU_ROSS(:)
	REAL*8, ALLOCATABLE :: TAU_ES(:)
!
	REAL*8 RMDOT
	REAL*8 RLUM
	REAL*8 ABUND_HYD
	INTEGER ND_ATM,NC_ATM,NP_ATM
	CHARACTER*21 TIME
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ROSS_MEAN(:)
	REAL*8, ALLOCATABLE :: FLUX_MEAN(:)
	REAL*8, ALLOCATABLE :: POP_ATOM(:)
	REAL*8, ALLOCATABLE :: MASS_DENSITY(:)
	REAL*8, ALLOCATABLE :: POPION(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC(:)
!
! Vectors for passing data to plot package via calls to CURVE.
!
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: XV_MID(:)
	REAL*8, ALLOCATABLE :: YV(:)
	REAL*8, ALLOCATABLE :: ZV(:)
!
	CHARACTER*80 NAME		!Default title for plot
	CHARACTER*80 XAXIS,XAXSAV	!Label for Abscissa
	CHARACTER*80 YAXIS		!Label for Ordinate
	CHARACTER*80 XAX_OPTION
!
	REAL*8 ANG_TO_HZ
	REAL*8 KEV_TO_HZ
	REAL*8 C_CMS
	REAL*8 C_KMS
	REAL*8 SN_AGE
!
	LOGICAL LOG_X,LOG_Y
	CHARACTER*10 Y_PLT_OPT,X_UNIT
!
	CHARACTER*6 METHOD,TYPE_ATM
	CHARACTER*10 NAME_CONVENTION
!
        INTEGER ACCESS_F
        INTEGER, PARAMETER :: EDD_CONT_REC=3
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
        INTEGER REC_SIZE
        INTEGER UNIT_SIZE
        INTEGER WORD_SIZE
        INTEGER N_PER_REC

! Miscellaneous variables.
!
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,K,L,ML,ISAV
	INTEGER ST_REC
	INTEGER REC_LENGTH
	INTEGER NEW_ND
	REAL*8 SCALE_FAC
	REAL*8 RVAL
	REAL*8 TEMP
	REAL*8 DTDR
	REAL*8 RADIUS
	REAL*8 LAMBDA
	REAL*8 T1,T2,T3
	REAL*8 PI
	REAL*8 LAMC
	REAL*8 T_ELEC
	REAL*8, ALLOCATABLE :: NEW_R(:)
	LOGICAL AIR_LAM
	LOGICAL FILE_PRES
	LOGICAL TMP_LOG
	LOGICAL PLT_H
	LOGICAL NORM
	LOGICAL VADAT_EXISTS
!
	INTEGER LEN_DIR
	CHARACTER(LEN=80) DIR_NAME
	CHARACTER(LEN=80) VADAT_FILE
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_IN=5		!For file I/O
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LU_IN=10	!For file I/O
	INTEGER, PARAMETER :: LU_OUT=11
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER GET_INDX_DP
!
! USR_OPTION variables
!
	CHARACTER MAIN_OPT_STR*80	!Used for input of the main option
	CHARACTER X*10			!Used for the individual option
	CHARACTER STRING*80
	CHARACTER*120 DEFAULT
	CHARACTER*120 DESCRIPTION
	CHARACTER*80 RVTJ_FILE_NAME
!
	REAL*8 SPEED_OF_LIGHT,FAC,LAM_VAC
	LOGICAL EQUAL
	CHARACTER*30 UC
	EXTERNAL SPEED_OF_LIGHT,EQUAL,FAC,UC,LAM_VAC,GET_INDX_DP
!
! 
! Set constants.
!
	CHIBF=2.815D-06
	CHIFF=3.69D-29
	HDKT=4.7994145D0
	TWOHCSQ=0.0147452575D0
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
	PI=ACOS(-1.0D0)
!
	C_CMS=SPEED_OF_LIGHT()
	C_KMS=1.0D-05*C_CMS
	XAX_OPTION='XLOGR'
	LAMBDA=5000.0D0
!
        CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
! Set defaults.
!
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	NAME=' '
	METHOD='LOGMON'
	TYPE_ATM=' '			!i.e. not Exponential at outer boundary
!
	LOG_X=.FALSE.
	LOG_Y=.FALSE.
	X_UNIT='ANG'
	Y_PLT_OPT='NAT'
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	KEV_TO_HZ=0.241838D+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
	RADIUS=1.0D0; DTDR=1.0D0; TEMP=5.0D0 	!Non zero defaults
!
!  Read in default model.
!
	ID=1
	NUM_FILES=1
	ZM(ID)%FILENAME='JH_AT_CURRENT_TIME'
5	CALL GEN_IN(ZM(ID)%FILENAME,'First data file')
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,ZM(ID)%FILE_DATE,ZM(ID)%FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error opening/reading INFO file: check format'
	  WRITE(T_OUT,*)'Also check error file or fort.2'
	  GOTO 5
	END IF
	OPEN(UNIT=LU_IN,FILE=ZM(ID)%FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error opening ',TRIM(ZM(ID)%FILENAME)
	     WRITE(T_OUT,*)'IOS=',IOS
	     GOTO 5
	  END IF
	  READ(LU_IN,REC=3)ST_REC,ZM(ID)%NCF,ZM(ID)%ND
	  WRITE(6,*)ST_REC,ZM(ID)%NCF,ZM(ID)%ND
	  ST_REC=ST_REC+2
	  ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	  ALLOCATE (ZM(ID)%RJ(ND,NCF))
	  ALLOCATE (ZM(ID)%HFLUX(ND,NCF))
	  ALLOCATE (ZM(ID)%NU(NCF))
	  ALLOCATE (ZM(ID)%JGREY(ND))
	  ALLOCATE (ZM(ID)%HGREY(ND))
	  ALLOCATE (ZM(ID)%R(ND))
	  ALLOCATE (ZM(ID)%V(ND))
	  ALLOCATE (ZM(ID)%H_INBC(NCF))
	  ALLOCATE (ZM(ID)%H_OUTBC(NCF))
	  READ(LU_IN,REC=ST_REC-2,IOSTAT=IOS)(ZM(ID)%R(I),I=1,ZM(ID)%ND),(ZM(ID)%V(I),I=1,ZM(ID)%ND)
	  READ(LU_IN,REC=ST_REC-1,IOSTAT=IOS)(ZM(ID)%JGREY(I),I=1,ZM(ID)%ND),
	1            (ZM(ID)%HGREY(I),I=1,ZM(ID)%ND-1)
	  DO ML=1,ZM(ID)%NCF
	    READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND),
	1            (ZM(ID)%HFLUX(I,ML),I=1,ZM(ID)%ND-1),
	1            ZM(ID)%H_INBC(ML),ZM(ID)%H_OUTBC(ML),ZM(ID)%NU(ML)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error reading all frequencies'
	      ZM(ID)%NCF=ML-1
	      EXIT
	    END IF
	  END DO
	  ZM(ID)%NCF=ZM(ID)%NCF-1
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in ',TRIM(ZM(ID)%FILENAME),' file as MODEL A (default)'
	WRITE(T_OUT,*)'Number of depth points is',ZM(ID)%ND
	WRITE(T_OUT,*)'Number of frequencies is ',ZM(ID)%NCF
!
	I=MAX(ZM(ID)%ND,ZM(ID)%NCF)
	ALLOCATE (XV(I),XV_MID(I),YV(I),ZV(I))
	XV_LENGTH=I
!
! Set default data types
!
	STRING=ZM(ID)%FILENAME
	CALL SET_CASE_UP(STRING,IZERO,IZERO)
	IF(INDEX(STRING,'JH_AT') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='RSQJ'
	ELSE IF(INDEX(STRING,'FLUX') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='RSQH'
	ELSE
	   ZM(ID)%DATA_TYPE='UNKNOWN'
	END IF
	ND=ZM(ID)%ND; XAXIS='Log R(10\u10\d cm)'
	XV(1:ND)=LOG10(ZM(ID)%R(1:ND))
!
!
!
! *************************************************************************
!
! Read in basic model [i.e. R, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
! Get default directory.
!
	DIR_NAME=' '            !Valid DIR_NAME if not present.
	LEN_DIR=0
	J=LEN(ZM(1)%FILENAME)
	DO WHILE(J .GT. 0)
	  IF( ZM(1)%FILENAME(J:J) .EQ. ']' .OR.
	1     ZM(1)%FILENAME(J:J) .EQ. ':' .OR.
	1     ZM(1)%FILENAME(J:J) .EQ. '/'        )THEN
	  DIR_NAME=ZM(1)%FILENAME(1:J)
	    LEN_DIR=J
	    J=0
	  END IF
	  J=J-1
	END DO
!
	RVTJ_FILE_NAME=DIR_NAME(1:LEN_DIR)//'RVTJ'
10	CALL GEN_IN(RVTJ_FILE_NAME,'File with R, V, T etc (RVTJ)')
	IF(INDEX(RVTJ_FILE_NAME,'NULL') .EQ. 0)THEN
	  OPEN(UNIT=LU_IN,FILE=RVTJ_FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	      RVTJ_FILE_NAME='../RVTJ'
	      GOTO 10
	    END IF
	  CLOSE(LU_IN)
	  CALL RD_RVTJ_PARAMS_V2(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1             ND_ATM,NC_ATM,NP_ATM,RVTJ_FILE_NAME,LU_IN)
	ELSE
	  ND_ATM=10
	END IF
	ALLOCATE (R(ND_ATM));			R=0.0D0
	ALLOCATE (V(ND_ATM));			V=0.0D0
	ALLOCATE (SIGMA(ND_ATM));		SIGMA=0.0D0
	ALLOCATE (T(ND_ATM));			T=0.0D0
	ALLOCATE (ED(ND_ATM));			ED=0.0D0
	ALLOCATE (ROSS_MEAN(ND_ATM));		ROSS_MEAN=0.0D0
	ALLOCATE (FLUX_MEAN(ND_ATM));		FLUX_MEAN=0.0D0
	ALLOCATE (POP_ATOM(ND_ATM));		POP_ATOM=0.0D0
	ALLOCATE (MASS_DENSITY(ND_ATM));	MASS_DENSITY=0.0D0
	ALLOCATE (POPION(ND_ATM));		POPION=0.0D0
	ALLOCATE (CLUMP_FAC(ND_ATM));		CLUMP_FAC=1.0D0
	IF(ND_ATM .GT. 10)THEN
	  CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,ND_ATM,LU_IN)
	END IF
	CLOSE(LU_IN)
!
! Now compute the important optical depth scales.
!
	 ALLOCATE (TA(ND_ATM))
	 ALLOCATE (TB(ND_ATM))
	 ALLOCATE (TC(ND_ATM))
	 ALLOCATE (TAU_ROSS(ND_ATM))
	 ALLOCATE (TAU_ES(ND_ATM))
	 IF(ND_ATM .NE. 10)THEN
	   IF(ROSS_MEAN(ND_ATM) .NE. 0)THEN
	     CALL TORSCL(TAU_ROSS,ROSS_MEAN,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
	   ELSE
	     TAU_ROSS(1:ND)=0.0D0
	   END IF
	   TA(1:ND_ATM)=6.65D-15*ED(1:ND_ATM)
	   CALL TORSCL(TAU_ES,TA,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
	 ELSE

	 END IF
!
! 
!
! This message will only be printed once
!
	WRITE(T_OUT,*)
	WRITE(T_OUT,"(8X,A)")'(default is to write file main_option.sve)'
	WRITE(T_OUT,"(8X,A)")'(append sve=filename to write a new .sve file)'
	WRITE(T_OUT,"(8X,A)")'(box=filename to write a .box file containing several .sve files)'
	WRITE(T_OUT,"(8X,A)")'(.filename to read .sve file)'
	WRITE(T_OUT,"(8X,A)")'(#filename to read .box file)'
	WRITE(T_OUT,*)
!
! This call resets the .sve algorithm.  Specifically it sets the next
! input answer to be a main option, and all subsequent inputs to be
! sub-options.  
!
 3	CALL SVE_FILE('RESET')
!
	MAIN_OPT_STR='  '
	DEFAULT='GR'
	DESCRIPTION=' '					!Obvious main option
	CALL USR_OPTION(MAIN_OPT_STR,'OPTION',DEFAULT,DESCRIPTION)
!
!   If the main option begins with a '.', a previously
!   written .sve file is read.
!
!   If the main option begins with a '#', a previously
!   written .box file is read.
!
!   If sve= is appended to the end of this main option, a new .sve file
!   is opened with the given name and the main option and all subsequent
!   sub-options are written to this file.
!
!   If box= is input then a .box file is created, which contains the name
!   of several .sve files to process.
!
!   If only a main option is given, the option and subsequent sub-options
!   are saved in a file called 'main option.sve'.  All following main
!   options are saved in separate files.
!
	X=UC(MAIN_OPT_STR)
	I=INDEX(X,'(')
	IF(I .NE. 0)X=X(1:I-1)		!Remove line variables
!
! 
!
	IF(X(1:3) .EQ. 'TIT')THEN
	  CALL USR_OPTION(NAME,'Title',' ',' ')
!                    
! Set X-Axis plotting options.
!
	ELSE IF(X(1:2) .EQ.'LX' .OR. X(1:4) .EQ. 'LOGX' .OR. 
	1                            X(1:4) .EQ. 'LINX')THEN
	  LOG_X=.NOT. LOG_X
	  IF(LOG_X)WRITE(T_OUT,*)'Now using Logarithmic X axis'
	  IF(.NOT. LOG_X)WRITE(T_OUT,*)'Now using Linear X axis'
	ELSE IF(X(1:2) .EQ.'XU' .OR. X(1:6) .EQ. 'XUNITS')THEN
	  CALL USR_OPTION(X_UNIT,'X_UNIT','Ang',
	1                  'Ang, um, eV, keV, Hz, Mm/s, km/s')
	  CALL SET_CASE_UP(X_UNIT,IZERO,IZERO)
	  IF(X_UNIT .NE. 'ANG' .AND.
	1        X_UNIT .NE. 'UM' .AND.
	1        X_UNIT .NE. 'EV' .AND.
	1        X_UNIT .NE. 'KEV' .AND.
	1        X_UNIT .NE. 'HZ' .AND.
	1        X_UNIT .NE. 'MM/S' .AND.
	1        X_UNIT .NE. 'KM/S')THEN
	     WRITE(T_OUT,*)'Invalid X unit: Try again'
	   END IF
!
! NB: We offer the option to use the central frequency to avoid
! air/vacuum confusions. Model data is in vacuum wavelengths, which
! we use in plotting at all wavelengths.
!
	   IF(X_UNIT .EQ. 'MM/S' .OR. X_UNIT .EQ. 'KM/S')THEN
	     CALL USR_OPTION(LAMC,'LAMC','0.0D0',
	1             'Central Lambda(Ang) [-ve for frequency (10^15 Hz)]')
	     IF(LAMC .LT. 0)THEN
	       LAMC=1.0D-07*C_CMS/ABS(LAMC)
	     ELSE
	       IF(LAMC .GT. 2000.0D0)THEN
                 CALL USR_OPTION(AIR_LAM,'AIR','T',
	1                'Air wavelength [only for Lam > 2000A]?')
	       ELSE
	         AIR_LAM=.FALSE.
	       END IF
	     END IF
	     IF(AIR_LAM)LAMC=LAM_VAC(LAMC)
	   END IF
!
! Set Y axis plotting options.
!
	ELSE IF(X(1:2) .EQ. 'LY' .OR. X(1:4) .EQ. 'LOGY' .OR. 
	1                             X(1:4) .EQ. 'LINY')THEN
	  LOG_Y=.NOT. LOG_Y
	  IF(LOG_Y)WRITE(T_OUT,*)'Now using Logarithmic Y axis'
	  IF(.NOT. LOG_Y)WRITE(T_OUT,*)'Now using Linear Y axis'
	ELSE IF(X(1:2) .EQ.'YU' .OR. X(1:6) .EQ. 'YUNITS')THEN
	  CALL USR_OPTION(Y_PLT_OPT,'X_UNIT',' ',
	1          'NAT(URAL), FNU, NU_FNU, FLAM')
	  CALL SET_CASE_UP(Y_PLT_OPT,IZERO,IZERO)
	  IF(Y_PLT_OPT .NE. 'NAT' .AND.
	1        Y_PLT_OPT .NE. 'FNU' .AND.
	1        Y_PLT_OPT .NE. 'NU_FNU' .AND.
	1        Y_PLT_OPT .NE. 'FLAM')THEN
	     WRITE(T_OUT,*)'Invalid Y Plot option: Try again'
	  END IF
!
! 
!
	ELSE IF(X(1:5) .EQ. 'XLINR')THEN
	   XAX_OPTION='XLINR'
	ELSE IF(X(1:5) .EQ. 'XLOGR')THEN
	   XAX_OPTION='XLOGR'
	ELSE IF(X(1:5) .EQ. 'XLOGV')THEN
	   XAX_OPTION='XLOGV'
	ELSE IF(X(1:2) .EQ. 'XV' .OR. X(1:5) .EQ. 'XLINV')THEN
	   XAX_OPTION='XVEL'
	ELSE IF(X(1:2) .EQ. 'XN')THEN
	   XAX_OPTION='XN'
!
	ELSE IF(X(1:4) .EQ. 'WRID')THEN
	  DO ID=1,NUM_FILES
	    WRITE(T_OUT,'(A,I2,A,A)')' ID=',ID,'          ',TRIM(ZM(ID)%FILENAME)
	  END DO
	ELSE IF(X(1:6) .EQ. 'RD_MOD')THEN
	  NUM_FILES=NUM_FILES+1
	  ID=NUM_FILES
	  ZM(ID)%FILENAME='JH_AT_CURRENT_TIME'
50	  CALL GEN_IN(ZM(ID)%FILENAME,'First data file')
	  CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,ZM(ID)%FILE_DATE,ZM(ID)%FILENAME,LU_IN,IOS)
	  IF(IOS .NE. 0)GOTO 50
	  OPEN(UNIT=LU_IN,FILE=ZM(ID)%FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	     IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error opening ',TRIM(ZM(ID)%FILENAME)
	       WRITE(T_OUT,*)'IOS=',IOS
	       GOTO 50
	    END IF
	    READ(LU_IN,REC=3)ST_REC,ZM(ID)%NCF,ZM(ID)%ND
	    ST_REC=ST_REC+2
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    ALLOCATE (ZM(ID)%RJ(ND,NCF))
	    ALLOCATE (ZM(ID)%HFLUX(ND,NCF))
	    ALLOCATE (ZM(ID)%NU(NCF))
	    ALLOCATE (ZM(ID)%JGREY(ND))
	    ALLOCATE (ZM(ID)%HGREY(ND))
	    ALLOCATE (ZM(ID)%R(ND))
	    ALLOCATE (ZM(ID)%V(ND))
	    ALLOCATE (ZM(ID)%H_INBC(NCF))
	    ALLOCATE (ZM(ID)%H_OUTBC(NCF))
	    READ(LU_IN,REC=ST_REC-2,IOSTAT=IOS)(ZM(ID)%R(I),I=1,ZM(ID)%ND),(ZM(ID)%V(I),I=1,ZM(ID)%ND)
	    READ(LU_IN,REC=ST_REC-1,IOSTAT=IOS)(ZM(ID)%JGREY(I),I=1,ZM(ID)%ND),
	1            (ZM(ID)%HGREY(I),I=1,ZM(ID)%ND-1)
	    DO ML=1,ZM(ID)%NCF
	      READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND),
	1              (ZM(ID)%HFLUX(I,ML),I=1,ZM(ID)%ND-1),
	1              ZM(ID)%H_INBC(ML),ZM(ID)%H_OUTBC(ML),ZM(ID)%NU(ML)
	    END DO
	  CLOSE(LU_IN)
	  WRITE(T_OUT,*)'Successfully read in ',TRIM(ZM(ID)%FILENAME),' file'
	  WRITE(T_OUT,*)'Number of depth points is',ZM(ID)%ND
	  WRITE(T_OUT,*)'Number of frequencies is ',ZM(ID)%NCF
!
	I=MAX(ZM(ID)%ND,ZM(ID)%NCF)
	IF(I .GT. XV_LENGTH)THEN
	  DEALLOCATE (XV,YV,ZV)
	  ALLOCATE (XV(I),YV(I),ZV(I))
	  XV_LENGTH=I
	END IF
!
! Set default data types
!
	  STRING=ZM(ID)%FILENAME
	  CALL SET_CASE_UP(STRING,IZERO,IZERO)
	  IF(INDEX(STRING,'JH_AT') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='RSQJ'
	  ELSE
	    ZM(ID)%DATA_TYPE='UNKNOWN'
	  END IF
!
	ELSE IF(X(1:7) .EQ. 'H_OUTBC')THEN
	  ND=ZM(ID)%ND
	  DO ID=1,NUM_FILES
	    ZV(1)=ZM(ID)%NU(1); YV(1)=ZM(ID)%H_OUTBC(1)
	    T1=0.0D0
	    DO ML=1,NCF
	      YV(ML)=ZM(ID)%H_OUTBC(ML)
	      T1=T1+0.5D0*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML+1))*ZM(ID)%H_OUTBC(ML)*ZM(ID)%RJ(1,ML)
	    END DO
	    T2=1.0D+20*16.0D0*PI*PI/3.286D+33
	    WRITE(6,*)'Integrated HFLUX is',1.0D+15*T1,1.0D+15*T1*T2
	    WRITE(6,*)'      Grey HFLUX is',ZM(ID)%HGREY(1),ZM(ID)%HGREY(1)*T2
	    CALL DP_CURVE(NCF,ZM(ID)%NU,YV)
	    ZV(NCF)=ZM(ID)%NU(NCF); YV(NCF)=ZM(ID)%H_OUTBC(NCF)
	  END DO
	ELSE IF(X(1:6) .EQ. 'H_INBC')THEN
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND
	    T1=0.0D0
	    ZV(1)=ZM(ID)%NU(1); YV(1)=ZM(ID)%H_INBC(1)
	    DO ML=2,NCF-1
	      ZV(ML)=ZM(ID)%NU(ML)
	      YV(ML)=ZM(ID)%H_INBC(ML)
	      T1=T1+0.5D0*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML+1))*ZM(ID)%H_INBC(ML)
	    END DO
	    ZV(NCF)=ZM(ID)%NU(NCF); YV(NCF)=ZM(ID)%H_INBC(NCF)
	    T2=1.0D+20*16.0D0*PI*PI/3.286D+33
	    WRITE(6,*)' '
	    WRITE(6,*)'The following are evaluated for the inner boundary.'
	    WRITE(6,*)' '
	    WRITE(6,'(A,ES12.4,5X,A,ES13.4)')'Integrated HFLUX is',1.0D+15*T1,
	1             '      Luminosity is',1.0D+15*T1*T2
	    WRITE(6,'(A,ES12.4,5X,A,ES13.4)')'      Grey HFLUX is',ZM(ID)%HGREY(ND-1),
	1             ' Grey luminosity is',ZM(ID)%HGREY(ND-1)*T2
	    CALL DP_CNVRT_J_V2(ZV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL DP_CURVE(NCF,ZV,YV)
	  END DO
!
! Options to plot J and H computed by grey optiom.
!
	ELSE IF(X(1:2) .EQ. 'JG')THEN
	  NORM=.FALSE.
	  CALL GEN_IN(NORM,'Normalize J(Grey) by outer boundary value?')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND
	    CALL SET_X_AXIS_V2(XV,XV_MID,XAXIS,ZM(ID)%R,ZM(ID)%V,XAX_OPTION,ND)
	    WRITE(6,*)'R^2.J at outer boundary is',1.0D+20*ZM(ID)%JGREY(1),'ergs/s'
	    IF(NORM)THEN
	      YV(1:ND)=ZM(ID)%JGREY(1:ND)/ZM(ID)%JGREY(1)
	    ELSE
	      YV(1:ND)=ZM(ID)%JGREY(1:ND)
	    END IF
	    CALL DP_CURVE(ND,XV,YV)
	  END DO
	  YAXIS='r^2 J(grey)/10\u20\d'
	  IF(NORM)YAXIS='r\u2\dJG/r(1)\u2\dJG(1)'
	ELSE IF(X(1:2) .EQ. 'HG')THEN
	  NORM=.FALSE.
	  CALL GEN_IN(NORM,'Normalize H(Grey) by outer boundary value?')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND
	    CALL SET_X_AXIS_V2(XV,XV_MID,XAXIS,ZM(ID)%R,ZM(ID)%V,XAX_OPTION,ND)
	    WRITE(6,*)'R^2.H at outer boundary is',1.0D+20*ZM(ID)%HGREY(1),'ergs/s'
	    IF(NORM)THEN
	      YV(1:ND)=ZM(ID)%HGREY(1:ND)/ZM(ID)%HGREY(1)
	    ELSE
	      YV(1:ND)=ZM(ID)%HGREY(1:ND)
	    END IF
	    CALL DP_CURVE(ND,XV,YV)
	  END DO
	  YAXIS='r^2 H(grey)/10\u20\d'
	  IF(NORM)YAXIS='r\u2\dHG/r(1)\u2\dHG(1)'
!
	ELSE IF(X(1:2) .EQ. 'JD' .OR. X(1:2) .EQ. 'HD')THEN
!
	  SCALE_FAC=1.0D0
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0','Scale factor to prevent overflow')
	  CALL USR_OPTION(I,'Depth',' ','Depth index')
	  ISAV=I
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    I=ISAV
	    IF(ID .NE. 1 .AND. ND .NE. ZM(1)%ND)CALL USR_OPTION(I,'Depth',' ','Depth index')
	    IF(I .GT. ND .OR. I .LT. 1)THEN
	      WRITE(T_OUT,*)'Invalid depth; minimum value is',1
	      WRITE(T_OUT,*)'Invalid depth; maximum value is',ND
	      GOTO 1
	    END IF
	    IF(ID .EQ. 1 .AND. ND .EQ. ND_ATM)THEN
	      CALL DERIVCHI(TB,T,R,ND,'LOGLOG')
	      TEMP=T(I); DTDR=TB(I); RADIUS=R(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'    R(I)/R*=',R(I)/R(ND)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'       V(I)=',V(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'       T(I)=',T(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'      ED(I)=',ED(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'TAU_ROSS(I)=',TAU_ROSS(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'  TAU_ES(I)=',TAU_ES(I)
	    END IF
!
	    ZV(1:NCF)=ZM(ID)%NU(1:NCF)
	    IF(X(1:2) .EQ. 'JD')YV(1:NCF)=ZM(ID)%RJ(I,1:NCF)*SCALE_FAC
	    IF(X(1:2) .EQ. 'HD')YV(1:NCF)=ZM(ID)%HFLUX(I,1:NCF)*SCALE_FAC
	    CALL DP_CNVRT_J_V2(ZV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL DP_CURVE(NCF,ZV,YV)
	  END DO
!
	ELSE IF(X(1:2) .EQ. 'DJ' .OR. X(1:2) .EQ. 'DH')THEN
!
	  SCALE_FAC=1.0D0
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0','Scale factor to prevent overflow')
	  CALL USR_OPTION(I,'Depth',' ','Depth index')
	  ISAV=I
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    I=ISAV
	    IF(ID .NE. 1 .AND. ND .NE. ZM(1)%ND)CALL USR_OPTION(I,'Depth',' ','Depth index')
	    IF(I .GT. ND .OR. I .LT. 1)THEN
	      WRITE(T_OUT,*)'Invalid depth; minimum value is',1
	      WRITE(T_OUT,*)'Invalid depth; maximum value is',ND
	      GOTO 1
	    END IF
	    IF(ID .EQ. 1 .AND. ND .EQ. ND_ATM)THEN
	      CALL DERIVCHI(TB,T,R,ND,'LOGLOG')
	      TEMP=T(I); DTDR=TB(I); RADIUS=R(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'    R(I)/R*=',R(I)/R(ND)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'       V(I)=',V(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'       T(I)=',T(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'      ED(I)=',ED(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'TAU_ROSS(I)=',TAU_ROSS(I)
	      WRITE(T_OUT,'(X,A,1P,E14.6)')'  TAU_ES(I)=',TAU_ES(I)
	    END IF
!
	    YV(1)=0.0D0
	    ZV(1:NCF)=ZM(ID)%NU(1:NCF)
	    DO J=2,NCF
	      IF(X(1:2) .EQ.  'DJ')YV(J)=ABS(ZM(ID)%RJ(I,J+1)-ZM(ID)%RJ(I,J))/MIN(ZM(ID)%RJ(I,J+1),ZM(ID)%RJ(I,J))
	    END DO
	    CALL DP_CNVRT_J_V2(ZV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_TRUE)
	    CALL DP_CURVE(NCF,ZV,YV)
!
	  END DO
!
	ELSE IF(X(1:3) .EQ. 'JNU' .OR. X(1:3) .EQ. 'HNU')THEN
	  DEFAULT=WR_STRING(LAMBDA)
	  CALL USR_OPTION(LAMBDA,'Lambda',DEFAULT,'Wavelength in Ang (-ve for 10^15 Hz)')
	  IF(LAMBDA .EQ. 0.0D0)THEN
	    WRITE(6,*)'Invalid wavelength -- cannot be zero'
	    GOTO 3						!return to get new option
	  ELSE IF(LAMBDA .LT. 0.0D0)THEN
	    T1=ABS(LAMBDA)
	  ELSE
	     T1=0.299794D+04/LAMBDA
	  END IF
!
	  SCALE_FAC=1.0D0
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0','Scale factor to prevent overflow')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
!
            I=GET_INDX_DP(T1,ZM(ID)%NU,NCF)
	    WRITE(6,*)'Index=',I,'NCF=',NCF
	    IF(ZM(ID)%NU(I)-T1 .GT. T1-ZM(ID)%NU(I+1))I=I+1
	    WRITE(6,*)'Index=',I,'NCF=',NCF
	    WRITE(6,*)'Freq=',ZM(ID)%NU(I)
!
	    CALL SET_X_AXIS_V2(XV,XV_MID,XAXIS,ZM(ID)%R,ZM(ID)%V,XAX_OPTION,ND)
	    DO J=1,ND
	      IF(X(1:3) .EQ. 'JNU')YV(J)=ZM(ID)%RJ(J,I)
	      IF(X(1:3) .EQ. 'HNU')YV(J)=ZM(ID)%HFLUX(J,I)
	    END DO
	    IF(LOG_Y)THEN
	      DO J=1,ND
	        IF(YV(J) .GT. 0.0D0)THEN
	          YV(J)=LOG10(YV(J))
	        ELSE
	          YV(J)=-200.0
	        END IF
	      END DO
	    END IF
	    IF(X(1:3) .EQ. 'JNU')THEN
	      CALL DP_CURVE(ND,XV,YV)
	    ELSE
	      CALL DP_CURVE(ND-1,XV_MID,YV)
	    END IF
          END DO
!
	ELSE IF(X(1:3) .EQ. 'CFD')THEN
!
	  CALL USR_OPTION(PLT_H,'PLT_H','F','Use H instead of J for plots?')
	  CALL USR_OPTION(I,'Depth','30','Depth index')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    IF(ALLOCATED(TA))DEALLOCATE(TA)
	    ALLOCATE (TA(NCF))
!
	    TA(1:NCF)=0.0D0
	    IF(PLT_H)THEN
	      DO ML=2,NCF
	        T1=0.5D0*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML))
	        TA(ML)=TA(ML-1)+T1*(ZM(ID)%HFLUX(I,ML-1)+ZM(ID)%HFLUX(I,ML))
	      END DO
	      WRITE(T_OUT,'(/,A,ES11.4)')' Luminosity(Lsun) is:',TA(NCF)*4.1274D+03
	    ELSE
	      DO ML=2,NCF
	        T1=0.5D0*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML))
	        TA(ML)=TA(ML-1)+T1*(ZM(ID)%RJ(I,ML-1)+ZM(ID)%RJ(I,ML))
	      END DO
	      WRITE(T_OUT,'(/,A,ES11.4)')' J(Lsun) is:',TA(NCF)*4.1274D+03
	      WRITE(T_OUT,'(A,ES11.4)')'      Beta is:',V(I)/2.99792458D+05
	    END IF
!
! The transformiation from comoving-frame H to observing fram H is H(obs) = H + beta (J+K).
! At the outer boudary, we will assume H=J+K (since we don't have K).
!
	    IF(I .EQ. 1 .AND. PLT_H)THEN
	      T1=(1.0D0+2.0D0*V(1)/2.99792458D+05)*TA(NCF)
	      WRITE(T_OUT,'(/,A,ES11.4)')'Observer''s frame luminosity(Lsun) is:',T1*4.1274D+03
	    END IF
!
	    DO ML=1,NCF
	      XV(ML)=ZM(ID)%NU(ML)
	      YV(ML)=TA(ML)/TA(NCF)
	    END DO
	    CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,'NAT',
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL DP_CURVE(NCF,XV,YV)
	    YAXIS='C\dL\u(\gv)/L'
	  END DO
!
	ELSE IF(X(1:2) .EQ. 'CF')THEN
!
	  WRITE(6,*)' '
	  WRITE(6,*)' This option plots R^2.J by integrating over J'
	  WRITE(6,*)' The plot may be normalized by R^2,J at d=1, or by the local Planck function.'
	  WRITE(6,*)' '
!
	  CALL USR_OPTION(TMP_LOG,'DIVB','T','Normalize by B?')
	  CALL USR_OPTION(PLT_H,'PLT_H','F','Use H instead of J for plots?')
	  DO ID=1,NUM_FILES
	    IF(ALLOCATED(TA))DEALLOCATE(TA)
	    ALLOCATE (TA(ND))
!
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    TA(1:ND)=0.0D0
	    IF(PLT_H)THEN
	      DO ML=1,NCF-1
	        T1=0.5D0*(ZM(ID)%NU(ML)-ZM(ID)%NU(ML+1))
	        DO I=1,ND
	          TA(I)=TA(I)+T1*(ZM(ID)%HFLUX(I,ML)+ZM(ID)%HFLUX(I,ML+1))
	        END DO
	      END DO
	      WRITE(T_OUT,*)'    Boundary Luminosity is:',TA(1)*4.1274D+03
	      WRITE(T_OUT,*)'R^2.H at outer boundary is:',1.0D+35*TA(1),'ergs/s'
	      YAXIS='r\u2\dH/R\u2\dH(d=1)'
	      IF(TMP_LOG)YAXIS='H/B'
	    ELSE
	      DO ML=1,NCF-1
	        T1=0.5D0*(ZM(ID)%NU(ML)-ZM(ID)%NU(ML+1))
	        DO I=1,ND
	          TA(I)=TA(I)+T1*(ZM(ID)%RJ(I,ML)+ZM(ID)%RJ(I,ML+1))
	        END DO
	      END DO
	      WRITE(T_OUT,*)'Boundary Luminosity (based on J) is:',TA(1)*4.1274D+03
	      WRITE(T_OUT,*)'         R^2.J at outer boundary is:',1.0D+35*TA(1),'ergs/s'
	      YAXIS='r\u2\dJ/R\u2\dJ(d=1)'
	      IF(TMP_LOG)YAXIS='J/B'
	    END IF
!
	    CALL SET_X_AXIS_V2(XV,XV_MID,XAXIS,ZM(ID)%R,ZM(ID)%V,XAX_OPTION,ND)
	    IF(TMP_LOG)THEN
	      DO I=1,ND
	        YV(I)=0.1D0*PI*TA(I)/R(I)/R(I)/5.670400D-05/(T(I)**4)
	      END DO
	    ELSE
	      DO I=1,ND
	        YV(I)=TA(I)/TA(1)
	      END DO
	    END IF
	    CALL DP_CURVE(ND,XV,YV)
	  END DO
!
! 
!
	ELSE IF(X(1:2) .EQ. 'BB') THEN
!
	  I=1
	  CALL USR_OPTION(I,'Depth',' ','Depth index: for default R, T')
	  DEFAULT=WR_STRING(T(I))
	  RVAL=ZM(1)%R(I)
	  NCF=ZM(1)%NCF
          CALL USR_OPTION(TEMP,'TEMP',DEFAULT,' ')
          DO I=1,ZM(1)%NCF
            T3=HDKT*ZM(1)%NU(I)/TEMP
            IF(T3 .GT. 1.0D0)THEN
              YV(I)=RVAL*RVAL*TWOHCSQ*(ZM(1)%NU(I)**3)*DEXP(-T3)/(1.0D0-DEXP(-T3))
            ELSE
              YV(I)=RVAL*RVAL*TWOHCSQ*(ZM(1)%NU(I)**3)/(DEXP(T3)-1.0D0)
            END IF
            XV(I)=ZM(1)%NU(I)
         END DO
	 CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         'RSQJ',LAMC,XAXIS,YAXIS,L_FALSE)
	 CALL DP_CURVE(NCF,XV,YV)
!
! Plot energy density in the radiaton field.
!
	ELSE IF(X(1:2) .EQ. 'EJ')THEN
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    DO ML=1,NCF-1
	      DO J=1,ND
	        T1=(ZM(ID)%HFLUX(J,ML)+ZM(ID)%HFLUX(J,ML+1))*2.0D0*ZM(ID)%V(J)/C_KMS
	        YV(J)=YV(J)+(ZM(ID)%NU(ML)-ZM(ID)%NU(ML+1))*(ZM(ID)%RJ(J,ML)+ZM(ID)%RJ(J,ML+1)+T1)
	      END DO
	    END DO
	    T1=1.6D+16*ATAN(1.0D0)*1/SPEED_OF_LIGHT()      !4*PI*1.0D+15
	    YV(1:ND)=0.5D0*T1*YV(1:ND)
	    YV(1:ND)=3.280D-03*YV(1:ND)                    !(4*PI*Dex(+30)/L(sun)
	    CALL LUM_FROM_ETA(YV,R,ND)
	    DO I=ND-1,1,-1
	      YV(I)=YV(I+1)+YV(I)
	    END DO
	    T2=R(ND)
	    XV(1:ND)=DLOG10(R(1:ND)/T2)
	    CALL DP_CURVE(ND,XV,YV)
!
	    VADAT_FILE='VADAT'
	    INQUIRE(FILE=VADAT_FILE,EXIST=VADAT_EXISTS)
	    IF(.NOT. VADAT_EXISTS)THEN
	      VADAT_FILE='../VADAT'
	      INQUIRE(FILE=VADAT_FILE,EXIST=VADAT_EXISTS)
	    END IF
	    SN_AGE=0.0D0
	    IF(VADAT_EXISTS)THEN
	     CALL READ_KEYWORD(SN_AGE,'[SN_AGE]',L_FALSE,VADAT_FILE,L_TRUE,L_TRUE,7)
	    END IF
	    YAXIS='E(rad)(s.L\dsun\u)'
	    WRITE(6,'(A)')RED_PEN
	    WRITE(6,'(A,ES12.4,A)')'   Integerated energy is',YV(1),' s.Lsun'
	    WRITE(6,'(A,ES12.4,A)')'   Integerated energy is',YV(1)*3.826D+33,' ergs'
	    IF(SN_AGE .NE. 0.0D0)THEN
	      T1=YV(1)*SN_AGE*24.0D0*3600.0D0*3.826D+33
	      WRITE(6,'(A,ES12.4,A,F10.4,A)')' t.Integerated energy is',T1,
	1                     ' s ergs [SN age =',SN_AGE,' d]'
	    END IF
	    WRITE(6,'(A)')DEF_PEN
	  END DO
	ELSE IF(X(1:3) .EQ. 'INT') THEN
!
! For diagnostic purposes. Designed specifically to computes
! Int J dv & Int dJ/dlnv dv, and plot as a function of depth.
!
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND
            XV(1:ND)=R(1:ND)/R(ND)
            YV(1:ND)=0.0D0; ZV(1:ND)=0.0D0
            DO ML=2,NCF-1
              DO I=1,ND
                YV(I)=YV(I)+ZM(ID)%RJ(I,ML)*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML+1))
                T1=(ZM(ID)%RJ(I,ML-1)-ZM(ID)%RJ(I,ML))/LOG(ZM(ID)%NU(ML)/ZM(ID)%NU(ML-1))
	        ZV(I)=ZV(I)+T1*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML+1))
              END DO
            END DO
	    CALL SET_X_AXIS_V2(XV,XV_MID,XAXIS,ZM(ID)%R,ZM(ID)%V,XAX_OPTION,ND)
            YV=0.5D+15*YV; ZV=0.5D+15*ZV
            CALL DP_CURVE(ND,XV,YV)
            CALL DP_CURVE(ND,XV,ZV)
          END DO
	  XAXIS='R/R\d*\u'
	  YAXIS='r\u2\d Int J dv; r\u2\d Int dJ/dlnv dv'
!
! 
! Plot section:
!
	ELSE IF(X(1:2) .EQ. 'GR') THEN
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,' ')
	  XAXIS=XAXSAV
!
	ELSE IF(X(1:4) .EQ.'GRNL') THEN
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  XAXIS=XAXSAV
! 
!
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:4) .EQ. 'LIST' .OR. 
	1                             X(1:2) .EQ. 'HE' 
	1          .OR. X(1:4) .EQ. 'HELP')THEN
!
	   WRITE(6,*)' '
	   WRITE(6,*)'The following set X axis options when plottings as a function of depth.'
	   WRITE(6,*)' '
	   WRITE(6,*)'XLOGR    XLINR    XN      XLOGV    XV'
	   WRITE(6,*)' '
	   WRITE(6,*)'YAXIS options'
	   WRITE(6,*)' '
	   WRITE(6,*)'     JD: Plot R^2.Jv at a given depth'
	   WRITE(6,*)'     JG: Plot R^2.J (grey) as a function of depth'
	   WRITE(6,*)'    JNU: Plot R^2.J at a given frequency as a function of depth'
	   WRITE(6,*)'     CF: Plot Integral of R^2.J (or J/B) as a function of depth'
	   WRITE(6,*)'     BB: Plot Blackbody for a given temperature'
	   WRITE(6,*)' '
	   WRITE(6,*)'     HD: Plot R^2.Hv at a given depth'
	   WRITE(6,*)'     HG: Plot R^2.H (grey) as a function of depth'
	   WRITE(6,*)'    HNU: Plot R^2.H at a given frequency as a function of depth'
	   WRITE(6,*)' '
	   WRITE(6,*)' H_INBC: Plot R^2.H at inner boundary'
	   WRITE(6,*)'H_OUTBC: Plot R^2.H at inner boundary'
	   WRITE(6,*)' '
!
	ELSE IF(X(1:2) .EQ. 'EX') THEN
	  CALL DP_CURVE(0,XV,YV)
	  STOP
	ELSE IF(X(1:3) .EQ. 'BOX') THEN
	  CALL WR_BOX_FILE(MAIN_OPT_STR)
	ELSE
	  PRINT*,'OPTION REQUESTED DOES NOT EXIST'
	END IF
!
1	CONTINUE
	GO TO 3
!
	END
!
!
	FUNCTION FAC(N)
	REAL*8 FAC
	INTEGER N
	INTEGER, PARAMETER :: T_OUT=5
!
	IF(N .EQ. 0)THEN
	  FAC=1
	ELSE IF(N .LT. 0)THEN
	  WRITE(T_OUT,*)'Error in FAC --- invalid argument'
	  STOP
	ELSE
	  FAC=1
	  DO I=2,N
	    FAC=FAC*I
	  END DO
	END IF
!
	RETURN
   	END
!
	SUBROUTINE SET_X_AXIS_V2(XV,XV_MID,XAXIS,R,V,XAX_OPTION,ND)
	IMPLICIT NONE
!
	INTEGER ND
	REAL*8 XV(ND)
	REAL*8 XV_MID(ND)
	REAL*8 R(ND)
	REAL*8 V(ND)
	CHARACTER(LEN=*) XAXIS,XAX_OPTION
	INTEGER I
!
	IF(XAX_OPTION .EQ. 'XLOGR')THEN
	  DO I=1,ND-1
	    XV(I)=LOG10(R(I))
	    XV_MID(I)=LOG10(0.5D0*(R(I)+R(I+1)))
	  END DO
	  XV(ND)=LOG10(R(ND))
	  XAXIS='Log R(10\u10\d cm)'
	ELSE IF(XAX_OPTION .EQ. 'XLOGV')THEN
	  DO I=1,ND-1
	    XV(I)=LOG10(V(I))
	    XV_MID(I)=LOG10(0.5D0*(V(I)+V(I+1)))
	  END DO
	  XV(ND)=LOG10(V(ND))
	  XAXIS='Log V(km/s)'
	ELSE IF(XAX_OPTION .EQ. 'XLINR')THEN
	  DO I=1,ND-1
	    XV(I)=R(I)
	    XV_MID(I)=0.5D0*(R(I)+R(I+1))
	  END DO
	  XV(ND)=R(ND)
	  XAXIS='R(10\u10\d cm)'
	ELSE IF(XAX_OPTION .EQ. 'XVEL')THEN
	  DO I=1,ND-1
	    XV(I)=V(I)
	    XV_MID(I)=0.5D0*(V(I)+V(I+1))
	  END DO
	ELSE IF(XAX_OPTION .EQ. 'XN')THEN
	  DO I=1,ND-1
	    XV(I)=I
	    XV_MID(I)=I+0.5D0
	  END DO
	  XV(ND)=ND
	  XAXIS='I'
	END IF
!
	RETURN
	END

