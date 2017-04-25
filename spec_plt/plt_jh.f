!
! General routine for plotting and comparing J (or H) obtained from the 
! EDDFACTOR (or similar) data file.
!
! Various options are available to redden and normalize the model spectra. 
! Several different units can be used for the X and Y axes.
!
	PROGRAM PLT_JH
!
! Altered 02-Nov-2013 :  Modified EXTJ J option -- now uses MON_INTERP.
! Altered 09-Jan-2009 :  No longer require RVTJ file.
! Altered 23-Nov-2007 :  New option inserted to allow J (in EDDFACTOR File) to be extended to
!                           larger radii assuming simple dilution. (alteration done 14-Nov-2007).
! Altered 16-Jun-2000 : DIRECT_INFO call inserted.
!
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
	  REAL*8, POINTER :: NU(:)
	  CHARACTER*10 DATA_TYPE
	  CHARACTER*40 FILE_DATE
	  CHARACTER*80 FILENAME
	END TYPE MODEL_INTENSITY
	TYPE (MODEL_INTENSITY) ZM(5)
	INTEGER ND,NCF
	INTEGER NUM_FILES
	INTEGER ID
!
	INTEGER NCF_B
	INTEGER ND_B
	REAL*8, ALLOCATABLE :: RJ_B(:,:)
	REAL*8, ALLOCATABLE :: NU_B(:)
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
	INTEGER ND_ATM,NC_ATM,NP_ATM,NT_ATM
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
	REAL*8, ALLOCATABLE :: POPS(:,:)

!
!
! Vectors for passing data to plot package via calls to CURVE.
!
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
	REAL*8, ALLOCATABLE :: ZV(:)
!
	CHARACTER*80 NAME		!Default title for plot
	CHARACTER*80 XAXIS,XAXSAV	!Label for Abscissa
	CHARACTER*80 YAXIS		!Label for Ordinate
!
	REAL*8 ANG_TO_HZ
	REAL*8 KEV_TO_HZ
	REAL*8 C_CMS
	REAL*8 C_KMS
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
	REAL*8 TEMP
	REAL*8 DTDR
	REAL*8 RADIUS
	REAL*8 T1,T2,T3
	REAL*8 LAMC
	REAL*8 EDGE_FREQ
	REAL*8 T_ELEC
	REAL*8 SN_AGE
	REAL*8, ALLOCATABLE :: NEW_R(:)
	LOGICAL AIR_LAM
	LOGICAL USE_V
	LOGICAL PLOT_RSQJ
	LOGICAL FILE_PRES
	LOGICAL ZEROV
	LOGICAL NEWMOD
	LOGICAL WRITE_RVSIG
	LOGICAL VADAT_EXISTS
!
	INTEGER NITSF
	INTEGER RITE_N_TIMES
	INTEGER LST_NG
	INTEGER LEN_DIR
	CHARACTER(LEN=80) DIR_NAME
	CHARACTER(LEN=80) VADAT_FILE
	CHARACTER(LEN=100) HELP_FILE
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
!
	C_CMS=SPEED_OF_LIGHT()
	C_KMS=1.0D-05*C_CMS
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
	ZM(ID)%FILENAME='EDDFACTOR'
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
	  ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	  ALLOCATE (ZM(ID)%RJ(ND,NCF))
	  ALLOCATE (ZM(ID)%NU(NCF))
	  DO ML=1,ZM(ID)%NCF
	    READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND),ZM(ID)%NU(ML)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error reading all frequencies'
	      ZM(ID)%NCF=ML-1
	      EXIT
	    END IF
	  END DO
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in ',TRIM(ZM(ID)%FILENAME),' file as MODEL A (default)'
	WRITE(T_OUT,*)'    Number of depth points is',ZM(ID)%ND
	WRITE(T_OUT,*)'Number of frequencies read is',ZM(ID)%NCF
!
! Set default data types
!
	STRING=ZM(ID)%FILENAME
	CALL SET_CASE_UP(STRING,IZERO,IZERO)
	IF(INDEX(STRING,'EDDF') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='J'
	ELSE IF(INDEX(STRING,'FLUX') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='H'
	ELSE IF(INDEX(STRING,'FORCE') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='M(t)'
	ELSE IF(INDEX(STRING,'ETA') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='ETA'
	ELSE IF(INDEX(STRING,'CHI') .NE. 0)THEN
	   ZM(ID)%DATA_TYPE='CHI'
	ELSE
	   ZM(ID)%DATA_TYPE='UNKNOWN'
	END IF
100	CALL GEN_IN(ZM(ID)%DATA_TYPE,'Default data type is')
	IF( ZM(ID)%DATA_TYPE .NE. 'J' .AND.
	1   ZM(ID)%DATA_TYPE .NE. 'H' .AND.
	1   ZM(ID)%DATA_TYPE .NE. 'M(t)' .AND.
	1   ZM(ID)%DATA_TYPE .NE. 'ETA' .AND.
	1   ZM(ID)%DATA_TYPE .NE. 'CHI')THEN
	   WRITE(6,*)'Invalid data type: Valid types are J, H, M(t), ETA and CHI'
	   GOTO 100
	END IF
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
	NITSF=0; NT_ATM=0
	RVTJ_FILE_NAME=DIR_NAME(1:LEN_DIR)//'RVTJ'
10	CALL GEN_IN(RVTJ_FILE_NAME,'File with R, V, T etc (RVTJ)')
	IF(INDEX(RVTJ_FILE_NAME,'SCRTEMP') .NE. 0)THEN
	  CALL GET_ND_NT_NIT(RVTJ_FILE_NAME,ND_ATM,NT_ATM,NITSF,LU_IN,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to open or READ SCRTEMP and related files: IOS=',IOS
	    WRITE(T_OUT,*)'Enter NULL if no file available'
	    RVTJ_FILE_NAME='../RVTJ'
	    GOTO 10
	  END IF
	ELSE IF(INDEX(RVTJ_FILE_NAME,'NULL') .EQ. 0)THEN
	  OPEN(UNIT=LU_IN,FILE=RVTJ_FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	      WRITE(T_OUT,*)'Enter NULL if no file available'
	      WRITE(T_OUT,*)'If file name contains SCRTEMP, SCRTEMP will be read'
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
	IF(NITSF .NE. 0 .AND. NT_ATM .NE. 0)THEN
	  K=NITSF
	  RITE_N_TIMES=1
	  NEWMOD=.TRUE.
	  ALLOCATE (POPS(NT_ATM,ND_ATM));               POPS=0
	  CALL SCR_READ_V2(R,V,SIGMA,POPS,K,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT_ATM,ND,LU_IN,NEWMOD)
	  ED(1:ND_ATM)=POPS(NT_ATM-1,:)
	  T(1:ND_ATM)=POPS(NT_ATM,:)
	ELSE IF(ND_ATM .GT. 10)THEN
	  CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,ND_ATM,LU_IN)
	END IF
	CLOSE(LU_IN)
!
	 IF(ND_ATM .NE. ZM(1)%ND)THEN
	   WRITE(6,*)' ' 
	   WRITE(6,*)' WARNING -- ND in RVTJ differs from that assoicated with main input file ' 
	   WRITE(6,*)' ' 
	 END IF
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
	 IF(ZM(1)%DATA_TYPE .EQ. 'H' .AND.  ND .EQ. ZM(1)%ND)THEN
	   DO ML=1,NCF
	     DO I=1,ND
	       ZM(ID)%RJ(I,ML)=ZM(ID)%RJ(I,ML)/R(I)/R(I)
	     END DO
	   END DO
	 ELSE IF(ZM(1)%DATA_TYPE .EQ. 'H')THEN
	   WRITE(6,*)'Unable to comput H from R^2.H since R non-matching R grid'
	   WRITE(6,*)'Be warned -- quantities may contain an extra factor of R^2'
	 END IF
!
! 
!
! This message will only be printed once
!
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
	     CALL USR_OPTION(LAMC,'LAMC','0.0',
	1             'Central Lambda(Ang) [-ve for frequency (10^15 Hz)]')
	     IF(LAMC .LT. 0)THEN
	       LAMC=1.0D-07*C_CMS/ABS(LAMC)
	     ELSE
	       IF(LAMC .GT. 2000)THEN
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
	ELSE IF(X(1:6) .EQ. 'FIXNCF')THEN
	  ID=1
	  OPEN(UNIT=LU_IN,FILE=ZM(ID)%FILENAME,STATUS='OLD',
	1             RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	    READ(LU_IN,REC=3)ST_REC,I,J
	    WRITE(LU_IN,REC=3)ST_REC,ZM(ID)%NCF,ZM(ID)%ND
	  CLOSE(LU_IN)
	ELSE IF(X(1:4) .EQ. 'WRID')THEN
	  DO ID=1,NUM_FILES
	    WRITE(T_OUT,'(A,I2,A,A)')' ID=',ID,'          ',TRIM(ZM(ID)%FILENAME)
	  END DO
	ELSE IF(X(1:6) .EQ. 'RD_MOD')THEN
	  NUM_FILES=NUM_FILES+1
	  ID=NUM_FILES
	  ZM(ID)%FILENAME='EDDFACTOR'
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
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    ALLOCATE (ZM(ID)%RJ(ND,NCF))
	    ALLOCATE (ZM(ID)%NU(NCF))
	    DO ML=1,ZM(ID)%NCF
	      READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND),ZM(ID)%NU(ML)
	      IF(IOS .NE. 0)THEN
	        WRITE(T_OUT,*)'Error reading all frequencies'
	        ZM(ID)%NCF=ML-1
	        EXIT
	      END IF
	    END DO
	  CLOSE(LU_IN)
	  WRITE(T_OUT,*)'Successfully read in ',TRIM(ZM(ID)%FILENAME),' file'
	  WRITE(T_OUT,*)'Number of depth points is',ZM(ID)%ND
	  WRITE(T_OUT,*)'Number of frequencies is ',ZM(ID)%NCF
!
! Set default data types
!
	  STRING=ZM(ID)%FILENAME
	  CALL SET_CASE_UP(STRING,IZERO,IZERO)
	  IF(INDEX(STRING,'EDDF') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='J'
	  ELSE IF(INDEX(STRING,'FLUX') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='H'
	  ELSE IF(INDEX(STRING,'FORCE') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='M(t)'
	  ELSE IF(INDEX(STRING,'ETA') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='ETA'
	  ELSE IF(INDEX(STRING,'CHI') .NE. 0)THEN
	    ZM(ID)%DATA_TYPE='CHI'
	  ELSE
	    ZM(ID)%DATA_TYPE='UNKNOWN'
	  END IF
200	  CALL GEN_IN(ZM(ID)%DATA_TYPE,'Default data type is')
	  IF( ZM(ID)%DATA_TYPE .NE. 'J' .AND.
	1     ZM(ID)%DATA_TYPE .NE. 'H' .AND.
	1     ZM(ID)%DATA_TYPE .NE. 'M(t)' .AND.
	1     ZM(ID)%DATA_TYPE .NE. 'ETA' .AND.
	1     ZM(ID)%DATA_TYPE .NE. 'CHI')THEN
	      WRITE(6,*)'Invalid data type: Valid types are J, H, M(t), ETA and CHI'
	     GOTO 200
	  END IF
	  IF(ZM(ID)%DATA_TYPE .EQ. 'H' .AND.  ND .EQ. ZM(ID)%ND)THEN
	    DO ML=1,NCF
	      DO I=1,ND
	       ZM(ID)%RJ(I,ML)=ZM(ID)%RJ(I,ML)/R(I)/R(I)
	      END DO
	    END DO
	  ELSE IF(ZM(ID)%DATA_TYPE .EQ. 'H')THEN
	    WRITE(6,*)'Unable to comput H from R^2.H since R non-matching R grid'
	    WRITE(6,*)'Be warned -- quantities may contain an extra factor of R^2'
	  END IF
!
	ELSE IF(X(1:4) .EQ. 'EXTJ')THEN
!
! Quick and dirty option to extend J onto a larger grid. For best results, grid in inner
! region should be identical to grid in outer region.
!
	  STRING='NEW_R_GRID'
	  CALL GEN_IN(STRING,'File with new R grid')
	  OPEN(UNIT=LU_IN,FILE=TRIM(STRING),STATUS='OLD',ACTION='READ')
	    READ(LU_IN,'(A)')STRING
	    READ(LU_IN,'(A)')STRING
	    READ(LU_IN,'(A)')STRING
	    READ(LU_IN,*)T1,T1,K,NEW_ND
	    IF(ALLOCATED(NEW_R))DEALLOCATE(NEW_R)
	    ALLOCATE(NEW_R(NEW_ND))
	    DO I=1,NEW_ND
	      READ(LU_IN,*)NEW_R(I)
	      READ(LU_IN,*)(T1,J=1,K)
	    END DO
	   CLOSE(LU_IN)
!
! Make sure the temporary vectors are of sufficient length.
!
	  IF(NEW_ND .GT. ND_ATM)THEN
	    DEALLOCATE (TA,TB,TC)
	    ALLOCATE (TA(NEW_ND),TB(NEW_ND),TC(NEW_ND))
	  END IF
!
! Set up the extension & interpolation vectors.
!
	  K=0
	  DO I=1,NEW_ND
	    IF(NEW_R(I) .GT. R(1))THEN
	      TA(I)=(R(1)/NEW_R(I))**2
	      K=I
	    END IF
	  END DO
	  IF(K .EQ. 1)THEN
	    IF(ABS(NEW_R(1)-R(1))/(R(1)-R(2)) .LT. 0.01D0)THEN
	      K=0
	      NEW_R(1)=R(1)
	    END IF
	  END IF
!
	  IF(ABS(NEW_R(NEW_ND)-R(ND))/(R(ND-1)-R(ND)) .LT. 0.01D0)THEN
	    NEW_R(NEW_ND)=R(ND)
	  END IF
!
	  ACCESS_F=5
	  I=WORD_SIZE*(NEW_ND+1)/UNIT_SIZE; J=83
	  ZM(1)%FILE_DATE='20-Aug-2000'
	  I=WORD_SIZE*(NEW_ND+1)/UNIT_SIZE; J=83
	  INQUIRE(FILE='EDDFACTOR',EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    CALL WRITE_DIRECT_INFO_V3(NEW_ND,I,ZM(1)%FILE_DATE,'J_DATA',J)
	    OPEN(UNIT=83,FILE='J_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	  ELSE
	    CALL WRITE_DIRECT_INFO_V3(NEW_ND,I,ZM(1)%FILE_DATE,'EDDFACTOR',J)
	    OPEN(UNIT=83,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	  END IF
	  WRITE(83,REC=EDD_CONT_REC)ACCESS_F,NCF,NEW_ND
	  WRITE(6,*)NEW_R(K+1),NEW_R(NEW_ND)
	  WRITE(6,*)R(1),R(ZM(1)%ND)
	  DO ML=1,ZM(1)%NCF
	    DO I=1,K
	      TB(I)=TA(I)*ZM(1)%RJ(1,ML)
	    END DO
	    I=NEW_ND-K
            CALL MON_INTERP(TB(K+1),I,IONE,NEW_R(K+1),I,ZM(1)%RJ(1,ML),ZM(1)%ND,R,ZM(1)%ND)
            WRITE(83,REC=ACCESS_F-1+ML)(TB(I),I=1,NEW_ND),ZM(ID)%NU(ML)
	  END DO
	  CLOSE(UNIT=83)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)'New J data written to J_DATA and J_DATA_INFO'
	    WRITE(6,*)'Use these to replace EDDFACTOR and EDDFACTOR_INFO if extending R grid'
	  ELSE
	    WRITE(6,*)'New J data written to EDDFACTOR and EDDFACTOR_INFO in the current directory'
	  END IF
!
	ELSE IF(X(1:4) .EQ. 'SPHJ')THEN
!
! Quick and dirty option to modify J in outer wind for sphericity effects.
! Code assumes same grid.
!
	  T1=R(ND)
	  CALL GEN_IN(T1,'Radius at Tau=1')
	  DO I=1,ND
	    IF(R(I) .GT. T1)THEN
	      T2=(T1/R(I))**2
	      TA(I)=SQRT(1.0D0-SQRT(1.0D0-T2))
	    ELSE
	      TA(I)=1.0D0
	    END IF
	  END DO
!
	  ACCESS_F=5
	  I=WORD_SIZE*(ND+1)/UNIT_SIZE; J=83
	  ZM(1)%FILE_DATE='20-Aug-2000'
	  I=WORD_SIZE*(ND+1)/UNIT_SIZE; J=83
	  INQUIRE(FILE='EDDFACTOR',EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    CALL WRITE_DIRECT_INFO_V3(ND,I,ZM(1)%FILE_DATE,'J_DATA',J)
	    OPEN(UNIT=83,FILE='J_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	  ELSE
	    CALL WRITE_DIRECT_INFO_V3(ND,I,ZM(1)%FILE_DATE,'EDDFACTOR',J)
	    OPEN(UNIT=83,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	  END IF
	  WRITE(83,REC=EDD_CONT_REC)ACCESS_F,NCF,ND
	    DO ML=1,ZM(1)%NCF
	      DO I=1,ND
	        TB(I)=TA(I)*ZM(1)%RJ(I,ML)
	      END DO
	      WRITE(83,REC=ACCESS_F-1+ML)(TB(I),I=1,ND),ZM(ID)%NU(ML)
	    END DO
	  CLOSE(UNIT=83)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)'New J data written to J_DATA and J_DATA_INFO'
	    WRITE(6,*)'Use these to replace EDDFACTOR and EDDFACTOR_INFO if extending R grid'
	  ELSE
	    WRITE(6,*)'New J data written to EDDFACTOR and EDDFACTOR_INFO in the current directory'
	  END IF
!
	ELSE IF(X(1:4) .EQ. 'WSMJ')THEN
!
! Quick and dirty option to write out J on the normal size grid.
! For this to work, points must be inserted equally.
!
	  ID=1
	  K=(ZM(ID)%ND-1)/(ND_ATM-1)	!Number of points inserted/interval
	  IF( MOD(ZM(ID)%ND-1,ND_ATM-1) .EQ. 0)THEN
	    ACCESS_F=5
	    I=WORD_SIZE*(ND_ATM+1)/UNIT_SIZE; J=83
	    ZM(1)%FILE_DATE='20-Aug-2000'
	    CALL WRITE_DIRECT_INFO_V3(ND_ATM,I,ZM(1)%FILE_DATE,'J_DATA',J)
	    OPEN(UNIT=83,FILE='J_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	      WRITE(83,REC=EDD_CONT_REC)ACCESS_F,NCF,ND_ATM
	      DO ML=1,NCF
	        WRITE(83,REC=ACCESS_F-1+ML)(ZM(ID)%RJ(I,ML),I=1,ZM(ID)%ND,K),ZM(ID)%NU(ML)
	      END DO
	    CLOSE(UNIT=83)
	  ELSE
	    WRITE(T_OUT,*)'Error - nonuniform grid extension'
	    WRITE(T_OUT,*)'Unable to write out small J file'
	   STOP
	  END IF
!
	ELSE IF(X(1:3) .EQ. 'R3J')THEN
	  ID=1; ND=ZM(ID)%ND
	  IF(ND .NE. ND_ATM)THEN
	    WRITE(6,*)'Unable to plot r^3.J since ND is not equal to ND_ATM'
	    GOTO 1
	  END IF
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(ND))
	  ALLOCATE (YV(ND))
!
	  CALL USR_OPTION(USE_V,'USE_V','F','Use V for x-axis (otherwise R)')
	  IF(USE_V)THEN
	    XV(1:ND)=V(1:ND)
	    XAXIS='V(km/s)'
	  ELSE
	    XV(1:ND)=R(1:ND)
	    XAXIS='R(10\u10\dcm)'
	  END IF
	  TA=0.0D0
	  DO ML=1,NCF
	    DO I=1,ND
	      TA(I)=TA(I)+ZM(ID)%RJ(I,ML)
	    END DO
	  END DO
	  YV(1:ND)=TA(1:ND)*(R(1:ND)**3)
	  YV(1:ND)=YV(1:ND)
	  CALL DP_CURVE(ND,XV,YV)
	  YAXIS='r\u3\dJ'
	ELSE IF(X(1:4) .EQ. 'PHOT')THEN
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
	  CALL USR_OPTION(I,'Depth',' ','Depth index')
	  CALL USR_OPTION(LAMC,'WAVE',' ','Ionization edge(Ang)')
	  EDGE_FREQ=0.01D0*C_KMS/LAMC
	  WRITE(6,*)'Edge freq is',EDGE_FREQ
	  ND=ZM(1)%ND; NCF=ZM(1)%NCF
	  XV(1)=0.0D0; YV(1)=0.0D0
	  DO ML=2,NCF
	    IF(ZM(1)%NU(ML) .LT. EDGE_FREQ)EXIT
	    XV(ML)=0.01D0*C_KMS/ZM(1)%NU(ML)
	    YV(ML)=YV(ML-1)+0.5D0*(ZM(1)%NU(ML-1)-ZM(1)%NU(ML+1))*ZM(1)%RJ(I,ML)*
	1            (EDGE_FREQ/ZM(1)%NU(ML))**3
	    J=ML
	  END DO
	  T1=YV(J)
	  YV(1:J)=YV(1:J)/T1
	  CALL DP_CURVE(J,XV,YV)
	  YAXIS='Phot'
!
! Plot energy density in the radiaton field.
!
	ELSE IF(X(1:2) .EQ. 'EJ')THEN
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(ND)); XV=0.0D0
	  ALLOCATE (YV(ND)); YV=0.0D0
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    DO ML=1,NCF-1
	      DO J=1,ND
	        YV(J)=YV(J)+(ZM(ID)%NU(ML)-ZM(ID)%NU(ML+1))*(ZM(ID)%RJ(J,ML)+ZM(ID)%RJ(J,ML+1))
	      END DO
	    END DO
	    T1=1.6D+16*ATAN(1.0D0)*1/SPEED_OF_LIGHT()      !4*PI*1.0D+15
	    YV(1:ND)=0.5D0*T1*YV(1:ND)
	    YV(1:ND)=3.280D-03*YV(1:ND)*R(1:ND)*R(1:ND)    !(4*PI*Dex(+30)/L(sun)
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
!
	ELSE IF(X(1:2) .EQ. 'JD' .OR. X(1:5) .EQ. 'RSQJD')THEN
!
	  CALL USR_OPTION(I,'Depth',' ','Depth index')
	  ISAV=I
	  SCALE_FAC=1.0D0
	  PLOT_RSQJ=.FALSE.
	  ZEROV=.FALSE.
	  IF(X(1:5) .EQ. 'RSQJD')THEN
	    IF(ND .NE. ND_ATM)THEN
	      WRITE(6,*)'Unable to plot r^2.J since ND is not equal to ND_ATM'
	      GOTO 1
	    END IF
	    PLOT_RSQJ=.TRUE.
	    CALL USR_OPTION(ZEROV,'ZEROV','T','Correct wavelenghts to zero V')
	  END IF
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0','Scale factor to prevent overflow')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    I=ISAV
	    IF(ID .NE. 1 .AND. ND .NE. ZM(1)%ND)CALL USR_OPTION(I,'Depth',' ','Depth index')
	    IF(I .GT. ND)THEN
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
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(NCF))
	    ALLOCATE (YV(NCF))
!
	    XV(1:NCF)=ZM(ID)%NU(1:NCF)
	    YV(1:NCF)=ZM(ID)%RJ(I,1:NCF)*SCALE_FAC
	    IF(PLOT_RSQJ)YV(1:NCF)=YV(1:NCF)*R(I)*R(I)
	    IF(ZEROV)THEN
	      T1=(1.0D0+V(I)/C_KMS)/(1.0-V(I)/C_KMS) 
	      XV(1:NCF)=XV(1:NCF)*SQRT(T1)
	      YV(1:NCF)=YV(1:NCF)*T1
	    END IF
	    CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL DP_CURVE(NCF,XV,YV)
	  END DO
!
	ELSE IF(X(1:4) .EQ. 'DJNU')THEN
	  CALL USR_OPTION(L,'Depth',' ','Depth index')
	  IF(NUM_FILES .EQ. 1)THEN
	    WRITE(6,*)'Error -- this option only works on two (identical) models'
	  ELSE
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(NCF))
	    ALLOCATE (YV(NCF))
	    XV(1:NCF)=ZM(1)%NU(1:NCF)
	    DO ML=1,NCF
	      T1=ZM(1)%RJ(L,ML)+ZM(2)%RJ(L,ML)
	      IF(T1 .NE. 0)THEN
	        YV(ML)=200.0D0*(ZM(1)%RJ(L,ML)-ZM(2)%RJ(L,ML))/T1
	      ELSE
	        YV(ML)=200.0D0
	      END IF
	    END DO
	    YAXIS='Percentage difference'
	    XAXIS='Frequency (10^15 Hz)'
	    CALL DP_CURVE(NCF,XV,YV)
	  END IF
!
	ELSE IF(X(1:3) .EQ. 'JNU')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Vacuum wavelength in Ang (-ve for Hz)')
	  IF(T1 .LE. 0)THEN
	    T1=ABS(T1)
	  ELSE
	    T1=0.299794D+04/T1
	  END IF
!
	  SCALE_FAC=1.0D0
	  PLOT_RSQJ=.FALSE.
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0','Scale factor to prevent overflow')
	  CALL USR_HIDDEN(PLOT_RSQJ,'RSQJ','F','Plot r^2 Gamma J?')
	  CALL USR_OPTION(USE_V,'USE_V','F','Use V for x-axis (otherwise R)')
!
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
!
            I=GET_INDX_DP(T1,ZM(ID)%NU,NCF)
	    IF(ZM(ID)%NU(I)-T1 .GT. T1-ZM(ID)%NU(I+1))I=I+1
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(ND))
	    ALLOCATE (YV(ND))
!
	    IF(ND .EQ. ND_ATM)THEN
	      IF(USE_V)THEN
	        XV(1:ND)=V(1:ND)
	        XAXIS='V(km\u \ds\u-1\d)'
	      ELSE
	        WRITE(6,*)'R(1)=',R(1)
	        WRITE(6,*)'R(ND)=',R(ND)
	        T2=R(ND)
	        XV(1:ND)=DLOG10(R(1:ND)/T2)
	        XAXIS='R/R(ND)'
	      END IF
	    ELSE
	     WRITE(6,*)'Plotting against depth index since ND is not equal to ND_ATM'
	      DO I=1,ND
	        XV(I)=I
	      END DO
	      XAXIS='Depth index'
	    END IF
	    YAXIS=ZM(ID)%DATA_TYPE
	    IF(PLOT_RSQJ)THEN
	      YAXIS='Log r\u2\d'//YAXIS
	    ELSE
	      YAXIS='Log '//YAXIS
	    END IF
	    DO J=1,ND
	      IF(ZM(ID)%RJ(J,I) .GT. 0)THEN
	        IF(PLOT_RSQJ)THEN
	          YV(J)=DLOG10(ZM(ID)%RJ(J,I)*R(J)*R(J)/SQRT(1.0D0-V(J)*V(J)/C_KMS/C_KMS))
	        ELSE
	          YV(J)=DLOG10(ZM(ID)%RJ(J,I))
	        END IF
	      ELSE
	        YV(J)=-100.0
	      END IF
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
          END DO
!
	ELSE IF(X(1:2) .EQ. 'MT')THEN
	  CALL USR_OPTION(USE_V,'USE_V','T','Use V for x-axis (otherwise R)')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
!
            I=GET_INDX_DP(T1,ZM(ID)%NU,NCF)
	    IF(ZM(ID)%NU(I)-T1 .GT. T1-ZM(ID)%NU(I+1))I=I+1
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(ND))
	    ALLOCATE (YV(ND))
!
	    IF(USE_V)THEN
	      XAXIS='V(km/s)'
	      XV(1:ND)=V(1:ND)
	    ELSE
	      XAXIS='R/R(ND)'
	      T2=R(ND)
	      XV(1:ND)=R(1:ND)/T2
	    END IF
	    YV(1:ND)=ZM(ID)%RJ(1:ND,NCF)
	    YAXIS='M(t)'
	    CALL DP_CURVE(ND,XV,YV)
          END DO

	ELSE IF(X(1:3) .EQ. 'CFD')THEN
!
	  I=ND/2
	  CALL USR_OPTION(I,'Depth',' ','Depth index')
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    IF(ALLOCATED(TA))DEALLOCATE(TA)
	    ALLOCATE (XV(NCF))
	    ALLOCATE (YV(NCF))
	    ALLOCATE (TA(NCF))
!
	    TA(1:NCF)=0.0D0
	    DO ML=2,NCF
	      T1=0.5D0*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML))
	      TA(ML)=TA(ML-1)+T1*(ZM(ID)%RJ(I,ML-1)+ZM(ID)%RJ(I,ML))
	    END DO
	    IF( ZM(ID)%DATA_TYPE .EQ. 'H')THEN
	      IF(ND .EQ. ND_ATM)THEN
	        WRITE(T_OUT,'(/,A,ES11.4)')' Luminosity is:',R(I)*R(I)*TA(NCF)*4.1274D+03
	      ElSE 
	        WRITE(T_OUT,'(/,A,ES11.4)')' Luminosity/R(10^10cm)^2 is:',TA(NCF)*4.1274D+03
	      END IF
	    ELSE
	      IF(ND .EQ. ND_ATM)THEN
	        STRING=' 16 pi^2 r^2 '//TRIM(ZM(ID)%DATA_TYPE)//'/Lsun is:'
	        WRITE(T_OUT,'(/,A,ES11.4)')TRIM(STRING),R(I)*R(I)*TA(NCF)*4.1274D+03
	      ElSE 
	        STRING=' 1.6D+21 pi^2 '//TRIM(ZM(ID)%DATA_TYPE)//'/Lsun is:'
	        WRITE(T_OUT,'(/,A,ES11.4)')TRIM(STRING),TA(NCF)*4.1274D+03
	      END IF
	    END IF
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
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    IF(ALLOCATED(TA))DEALLOCATE(TA)
	    ALLOCATE (XV(ND))
	    ALLOCATE (YV(ND))
	    ALLOCATE (TA(ND))
!
	    TA(1:ND)=0.0D0
	    DO ML=1,NCF-1
	      T1=0.5D0*(ZM(ID)%NU(ML)-ZM(ID)%NU(ML+1))
	      DO I=1,ND
	        TA(I)=TA(I)+T1*(ZM(ID)%RJ(I,ML)+ZM(ID)%RJ(I,ML+1))
	      END DO
	    END DO
	    IF(ND .EQ. ND_ATM)THEN
	      IF(ZM(ID)%DATA_TYPE .EQ. 'J' .OR. ZM(ID)%DATA_TYPE .EQ. 'H')THEN
	        TA(1:ND)=TA(1:ND)*R(1:ND)*R(1:ND)
	        YAXIS='L/L(d=1)'
	        IF(ZM(ID)%DATA_TYPE .EQ. 'J')YAXIS='R\u2\d/J/R\u2\dJ(d=1)'
	      END IF
	      WRITE(T_OUT,*)'Boundary Luminosity is:',TA(1)*4.1274D+03
	    ELSE
	      YAXIS='H/H(d=1)'
	      IF(ZM(ID)%DATA_TYPE .EQ. 'J')YAXIS='J/J(d=1)'
	    END IF
!
	    DO I=1,ND
	      XV(I)=I
	      YV(I)=TA(I)/TA(1)
	    END DO
	    CALL DP_CURVE(ND,XV,YV)
	    XAXIS='I'
	  END DO
!
! 
! Convolve J with Electron scattering redistribution function.
!
	ELSE IF(X(1:2) .EQ. 'ES')THEN
	  CALL USR_OPTION(K,'Depth',' ','Depth index')
	  IF(K .LE. 0 .OR. K .GE. ND)THEN
	    WRITE(T_OUT,*)'Invalid depth index, ND_MAX=',ND
	    GOTO 1
	  END IF
	  ID=1; CALL GEN_IN(ID,'Which data file')
	  ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(ZM(ID)%NCF))
	  ALLOCATE (YV(ZM(ID)%NCF))
!
	  XV(1:NCF)=ZM(ID)%NU(1:NCF)
	  YV(1:NCF)=ZM(ID)%RJ(K,1:NCF)
	  CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
!
! Convolve RJ with electron-scattering redistribution function.
!
!	  DEFAULT=T(K)
!	  CALL USR_OPTION(T_ELEC,'Depth',Default,'Depth index')
	  T_ELEC=T(K)
!
	  IF(ALLOCATED(TB))DEALLOCATE(TB)	!J input
	  ALLOCATE (TB(NCF))
	  IF(ALLOCATED(TC))DEALLOCATE(TC)	!J e.s. output
	  ALLOCATE (TC(NCF))
!
	  TB(1:NCF)=ZM(ID)%RJ(K,1:NCF)
!	  IF(ONE_PAR)THEN
!	    CALL CNVLV_ES_ONE_PAR_V2(NU,TB,TC,T_ELEC,T_OUT,NCF)
!	  ELSE
	    CALL CNVLV_ES_TWO_PAR_V2(ZM(ID)%NU,TB,TC,T_ELEC,T_OUT,NCF)
!	  END IF
!
	  YV(1:NCF)=TC(1:NCF)
	  CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         ZM(ID)%DATA_TYPE,LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL DP_CURVE(NCF,XV,YV)
!
	ELSE IF(X(1:2) .EQ. 'BB') THEN
!
	NCF=ZM(1)%NCF
	IF(ALLOCATED(XV))DEALLOCATE(XV)
	IF(ALLOCATED(YV))DEALLOCATE(YV)
	ALLOCATE (XV(NCF))
	ALLOCATE (YV(NCF))
	I=1
	CALL USR_OPTION(I,'Depth',' ','Depth index: for default R, T')
	DEFAULT=WR_STRING(T(I))
	CALL USR_OPTION(TEMP,'TEMP',DEFAULT,' ')
          DO I=1,ZM(1)%NCF
            T3=HDKT*ZM(1)%NU(I)/TEMP
            IF(T3 .GT. 1.0D0)THEN
              YV(I)=TWOHCSQ*(ZM(1)%NU(I)**3)*DEXP(-T3)/(1.0D0-DEXP(-T3))
            ELSE
              YV(I)=TWOHCSQ*(ZM(1)%NU(I)**3)/(DEXP(T3)-1.0D0)
            END IF
            XV(I)=ZM(1)%NU(I)
         END DO
	 CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         'J',LAMC,XAXIS,YAXIS,L_FALSE)
	 CALL DP_CURVE(NCF,XV,YV)
!
        ELSE IF(X(1:4) .EQ. 'DBDR')THEN
!
! This option is designed to compute 1/3 . dB/dr. This can
! be compared directly with chi.flux. They should be equal if
! the diffusion approximation is valid.
!
! You should call JD option first.
!
	  WRITE(6,*)'Calculating diffusive flux: 1/3. dB/dr'
	  WRITE(6,*)'Call JD option first to set R, T, DTDR'
	  NCF=ZM(1)%NCF
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
	  DEFAULT=WR_STRING(TEMP)
	  CALL USR_OPTION(TEMP,'TEMP',DEFAULT,'(Program units)')
	  DEFAULT=WR_STRING(DTDR)
	  CALL USR_OPTION(DTDR,'DTDR',DEFAULT,'(Program units)')
	  DEFAULT=WR_STRING(RADIUS)
	  CALL USR_OPTION(RADIUS,'RADIUS',DEFAULT,'(Program units)')
	  DO I=1,ZM(1)%NCF
	    T3=HDKT*ZM(1)%NU(I)/TEMP
	    YV(I)=RADIUS*RADIUS*ABS(DTDR)*TWOHCSQ*T3*(ZM(1)%NU(I)**3)/TEMP/3.0D0
	    IF(T3 .GT. 1.0D0)THEN
	      YV(I)=YV(I)*DEXP(-T3)/(1.0D0-DEXP(-T3))/(1.0D0-DEXP(-T3))
	    ELSE
	      YV(I)=YV(I)*DEXP(T3)/(DEXP(T3)-1.0D0)/(DEXP(T3)-1.0D0)
	    END IF
	    XV(I)=ZM(1)%NU(I)
	  END DO
	  CALL DP_CNVRT_J_V2(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         'H',LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL DP_CURVE(NCF,XV,YV)
!
	ELSE IF(X(1:3) .EQ. 'INT') THEN
!
! For diagnostic purposes. Designed specifically to computes
! Int J dv & Int dJ/dlnv dv, and plot as a function of depth.
!
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    IF(ALLOCATED(ZV))DEALLOCATE(ZV)
	    ALLOCATE (XV(ND),YV(ND),ZV(ND))
!
            XV(1:ND)=R(1:ND)/R(ND)
            YV(1:ND)=0.0D0; ZV(1:ND)=0.0D0
            DO ML=2,NCF-1
              DO I=1,ND
                YV(I)=YV(I)+ZM(ID)%RJ(I,ML)*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML+1))
                T1=(ZM(ID)%RJ(I,ML-1)-ZM(ID)%RJ(I,ML))/LOG(ZM(ID)%NU(ML)/ZM(ID)%NU(ML-1))
	        ZV(I)=ZV(I)+T1*(ZM(ID)%NU(ML-1)-ZM(ID)%NU(ML+1))
              END DO
            END DO
            YV=0.5D+15*YV; ZV=0.5D+15*ZV
            CALL DP_CURVE(ND,XV,YV)
            CALL DP_CURVE(ND,XV,ZV)
          END DO
	  XAXIS='R/R\d*\u'
	  YAXIS='Int J dv; Int dJ/dlnv dv'
!
! For testing and diagnostic purposes. Plots c.dNU/NU as a function on NU.
!
	ELSE IF(X(1:3) .EQ. 'DNU') THEN
!
	  DO ID=1,NUM_FILES
	    ND=ZM(ID)%ND; NCF=ZM(ID)%NCF
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(YV))DEALLOCATE(YV)
	    ALLOCATE (XV(NCF),YV(NCF))
!
            XV(1:NCF)=ZM(ID)%NU(1:NCF)
            DO ML=1,NCF-1
               YV(ML)=(ZM(ID)%NU(ML)-ZM(ID)%NU(ML+1))/(ZM(ID)%NU(ML)+ZM(ID)%NU(ML+1))
            END DO
            YV=2.0D0*YV*C_KMS
            CALL DP_CURVE(NCF-1,XV,YV)
          END DO	
	  YAXIS='c.dNU/NU(km/s)'
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
!
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:4) .EQ. 'LIST' .OR. 
	1                             X(1:2) .EQ. 'HE' 
	1          .OR. X(1:4) .EQ. 'HELP')THEN
	  CALL PLT_JH_OPT_DESC
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
