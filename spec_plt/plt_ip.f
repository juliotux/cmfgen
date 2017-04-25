!
! General routine for plotting and comparing I(p,nu) obtained from IP_DATA.
!
! Various options are available to redden and normalize the model spectra. 
! Several different units can be used for the X and Y axes.
!
	PROGRAM PLT_IP
!
! Altered 01-Mar-2016 : Updated to better handle plane-parallel models [5-Feb-2016].
! Altered 24-Jun-2009 : Use INU and INU2 instead of JNU,JNU2
! Altered 18-Aug-2003 : IOS added to DIRECT_INFO call
! Altered 16-Jun-2000 : DIRECT_INFO call inserted.
!
! Interface routines for IO routines.
!
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE GEN_IN_INTERFACE
!
	IMPLICIT NONE
!
	INTEGER NCF
	INTEGER ND
	INTEGER NC
	INTEGER NP
	REAL*8, ALLOCATABLE :: IP(:,:)
	REAL*8, ALLOCATABLE :: NU(:)
	REAL*8, ALLOCATABLE :: P(:)
	REAL*8, ALLOCATABLE :: MU(:)
	REAL*8, ALLOCATABLE :: TMP_MU(:)
	REAL*8, ALLOCATABLE :: IP_NEW(:)
	REAL*8, ALLOCATABLE :: P_NEW(:)
	REAL*8, ALLOCATABLE :: HQW_AT_RMAX(:)
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: TAU_ROSS(:)
	REAL*8, ALLOCATABLE :: TAU_ES(:)
!
	REAL*8, ALLOCATABLE :: TEMP_P(:)
	REAL*8, ALLOCATABLE :: TEMP_IP(:)
!
	REAL*8 RMDOT
	REAL*8 RLUM
	REAL*8 ABUND_HYD
	REAL*8 LAMBDA
	INTEGER NC2,NP2
	CHARACTER*21 TIME
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ROSS_MEAN(:)
	REAL*8, ALLOCATABLE :: FLUX_MEAN(:)
	REAL*8, ALLOCATABLE :: POPTOM(:)
	REAL*8, ALLOCATABLE :: MASS_DENSITY(:)
	REAL*8, ALLOCATABLE :: POPION(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC(:)

!
!
! Vectors for passing data to plot package via calls to CURVE.
!
	REAL*4, ALLOCATABLE :: XV(:)
	REAL*4, ALLOCATABLE :: YV(:)
	REAL*4, ALLOCATABLE :: WV(:)
	REAL*4, ALLOCATABLE :: ZV(:)
!
	CHARACTER*80 NAME		!Default title for plot
	CHARACTER*80 XAXIS,XAXSAV	!Label for Absisca
	CHARACTER*80 YAXIS		!Label for Ordinate
!
	REAL*8 SCALE_FAC
	REAL*8 RAD_VEL
	REAL*8 ADD_FAC
	REAL*8 WT(30)
	INTEGER NCF_MAX
	INTEGER IST,IEND
	INTEGER OBS_COLS(2)
	LOGICAL SMOOTH
	LOGICAL CLEAN
	LOGICAL NON_MONOTONIC
!
	REAL*8 ANG_TO_HZ
	REAL*8 KEV_TO_HZ
	REAL*8 C_CMS
	REAL*8 C_KMS
!
	LOGICAL LOG_X,LOG_Y
	CHARACTER*10 Y_PLT_OPT,X_UNIT
	CHARACTER*80 FILENAME
	CHARACTER*80 FILE_DATE
!
	CHARACTER*6 METHOD,TYPETM
	CHARACTER*10 NAME_CONVENTION
!
	REAL*8 DISTANCE			!kpc
	REAL*8 SLIT_WIDTH		!arcseconds
	REAL*8 PIXEL_LENGTH		!arcseconds
!
	REAL*8 X_CENT
	REAL*8 Y_CENT
	REAL*8 S_WIDTH
	REAL*8 S_LNGTH
	REAL*8 APP_SIZE
	REAL*8 TEL_FWHM
	REAL*8 MOD_RES
!
! Miscellaneous variables.
!
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,K,L,ML,LS
	INTEGER ST_REC
	INTEGER REC_LENGTH
	INTEGER NX
	INTEGER NP_NEW
	INTEGER NINS
	INTEGER NTMP
	INTEGER K_ST,K_END
	REAL*8 T1,T2
	REAL*8 LAMC
	REAL*8 DELV
	REAL*8 FRAC
	REAL*8 PI
	REAL*8 T_ELEC
	LOGICAL AIR_LAM
	LOGICAL COMPUTE_P
	LOGICAL USE_ARCSEC
	LOGICAL MULT_BY_PSQ
	LOGICAL MULT_BY_P
	LOGICAL PLANE_PARALLEL_MOD
!
	REAL*8,  PARAMETER :: RZERO=0.0D0
	REAL*8,  PARAMETER :: RONE=1.0D0
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
	CHARACTER MAIN_OPT_STR*80	!Used for input of the main otion
	CHARACTER X*10			!Used for the idividual option
	CHARACTER STRING*80
	CHARACTER*120 DEFAULT
	CHARACTER*120 DESCRIPTION
	CHARACTER(LEN=20) TMP_STR
!
	REAL*8 SPEED_OF_LIGHT,FAC,LAM_VAC,PARSEC
	LOGICAL EQUAL
	CHARACTER*30 UC
	EXTERNAL SPEED_OF_LIGHT,EQUAL,FAC,UC,LAM_VAC,GET_INDX_DP,PARSEC
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
	PI=ACOS(-1.0D0)
!
! Set defaults.
!
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	NAME=' '
	METHOD='LOGMON'
	TYPETM=' '			!i.e. not Exponential at outer boundary
!
	LOG_X=.FALSE.
	LOG_Y=.FALSE.
	X_UNIT='ANG'
	Y_PLT_OPT='FNU'
!
	DISTANCE=1.0D0		  !kpc
	SLIT_WIDTH=0.1D0          !arcseconds
	PIXEL_LENGTH=0.0254D0     !arcseconds
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	KEV_TO_HZ=0.241838D+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
!
!  Read in model.
!
	FILENAME='IP_DATA'
	CALL GEN_IN(FILENAME,'I(p) file')
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
	  ALLOCATE (MU(NP)); MU=0.0D0
	  ALLOCATE (NU(NCF))
	  IF( INDEX(FILE_DATE,'20-Aug-2000') .NE. 0)THEN
	    READ(LU_IN,REC=ST_REC)(P(I),I=1,NP)
	    ST_REC=ST_REC+1
	    COMPUTE_P=.FALSE.
	  ELSE IF( INDEX(FILE_DATE,'10-Apr-2007') .NE. 0 .OR. INDEX(FILE_DATE,'25-Jan-2016') .NE. 0)THEN
	    READ(LU_IN,REC=ST_REC)(P(I),I=1,NP)
	    ST_REC=ST_REC+1
	    READ(LU_IN,REC=ST_REC)(MU(I),I=1,NP)
	    ST_REC=ST_REC+1
	    COMPUTE_P=.FALSE.
	  ELSE IF( INDEX(FILE_DATE,'Unavailable') .NE. 0)THEN
	    COMPUTE_P=.TRUE.
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
	WRITE(T_OUT,*)'Number of angles is',NP
	WRITE(T_OUT,*)'Number of frequenct points is',NCF
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
10	FILENAME='../RVTJ'
	CALL GEN_IN(FILENAME,'File with R, V, T etc (RVTJ)')
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	    GOTO 10
	  END IF
	CLOSE(LU_IN)
	CALL RD_RVTJ_PARAMS_V2(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1             ND,NC2,NP2,FILENAME,LU_IN)
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
	ALLOCATE (T(ND))
	ALLOCATE (ED(ND))
	ALLOCATE (ROSS_MEAN(ND))
	ALLOCATE (FLUX_MEAN(ND))
	ALLOCATE (POPTOM(ND))
	ALLOCATE (MASS_DENSITY(ND))
	ALLOCATE (POPION(ND))
	ALLOCATE (CLUMP_FAC(ND))
	CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POPTOM,POPION,MASS_DENSITY,CLUMP_FAC,ND,LU_IN)
	CLOSE(LU_IN)
!
! Now compute the important optical depth scales.
!
	 ALLOCATE (TA(ND))
	 ALLOCATE (TB(ND))
	 ALLOCATE (TC(ND))
	 ALLOCATE (TAU_ROSS(ND))
	 ALLOCATE (TAU_ES(ND))
	 IF(ROSS_MEAN(ND) .NE. 0)THEN
	   CALL TORSCL(TAU_ROSS,ROSS_MEAN,R,TB,TC,ND,METHOD,TYPETM)
	 END IF
	 TA(1:ND)=6.65D-15*ED(1:ND)
	 CALL TORSCL(TAU_ES,TA,R,TB,TC,ND,METHOD,TYPETM)
!
	 IF(COMPUTE_P)THEN
	   NC=NP-ND
	   P(NC+1:NP)=R(ND:1:-1)
	   DO I=1,NC
	     P(I)=R(ND)*(I-1)/NC
	   END DO
	   WRITE(T_OUT,*)'P computed internally'
	 END IF
!
	IF(NP .LT. ND)THEN
	  ALLOCATE (TMP_MU(NP))
	  ALLOCATE (HQW_AT_RMAX(NP))
	  CALL GAULEG(RZERO,RONE,TMP_MU,HQW_AT_RMAX,NP)
	  PLANE_PARALLEL_MOD=.TRUE.
	  WRITE(6,*)'Model assumed to be plae-parallel'
	ELSE
	  ALLOCATE (HQW_AT_RMAX(NP))
	  CALL HWEIGHT(MU,HQW_AT_RMAX,NP)
	  PLANE_PARALLEL_MOD=.FALSE.
	END IF
!
	CALL GEN_IN(DISTANCE,'Distance to star (kpc)')
	CALL GEN_IN(SLIT_WIDTH,'Width of spectrograph slit (arcsec)')
	CALL GEN_IN(PIXEL_LENGTH,'Size of pixed in spatial direction (arcsec)')
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
!   If sve= is apended to the end of this main option, a new .sve file
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
! Set X-Ais plotting options.
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
! air/vacuum confusions. Model data is in vacuum wavelngths, which
! we use in plotting at all wavelngths.
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
	1          'FNU, NU_FNU, FLAM')
	  CALL SET_CASE_UP(Y_PLT_OPT,IZERO,IZERO)
	  IF(Y_PLT_OPT .NE. 'FNU' .AND.
	1        Y_PLT_OPT .NE. 'NU_FNU' .AND.
	1        Y_PLT_OPT .NE. 'FLAM')THEN
	     WRITE(T_OUT,*)'Invalid Y Plot option: Try again'
	  END IF
!
! 
!
! Simple option to compute inetgrated spectrum inside impact index I,
! and outside index I.
!
	ELSE IF(X(1:2) .EQ. 'SP')THEN
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
!
	  WRITE(6,*) HQW_AT_RMAX,R(ND)
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=0.0D0
	  DO J=1,NP
	    YV(1:NCF)=YV(1:NCF)+HQW_AT_RMAX(J)*IP(J,1:NCF)
	  END DO
          T1=DISTANCE*1.0D+03*PARSEC()
	  T1=2.0D0*R(ND)*R(ND)*PI*1.0D+23*(1.0E+10/T1)**2
	  YV(1:NCF)=YV(1:NCF)*T1
!
! NB: J and I have the same units, apart from per steradian.
!
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL CURVE(NCF,XV,YV)
!
	ELSE IF(X(1:3) .EQ. 'ISP')THEN
	  CALL USR_OPTION(I,'P',' ','Impact parameter index cuttoff')
	  IF(I .GT. NP)THEN
	    WRITE(T_OUT,*)'Invalid depth; maximum value is',NP
	    GOTO 1
	  END IF
!
	  T1=P(I)*1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	  WRITE(T_OUT,'(1X,A,1P,E14.6,A)')'       P(I)=',P(I)*1.0E+10,' cm'
	  WRITE(T_OUT,'(1X,A,1P,E14.6,A)')'       P(I)=',T1,' arcsec'
	  WRITE(T_OUT,'(1X,A,1P,E14.6)')'    P(I)/R*=',P(I)/R(ND)
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
!
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=0.0D0
	  DO J=1,I-1
	    YV(1:NCF)=YV(1:NCF)+0.5D0*(IP(J,1:NCF)*P(J)+IP(J+1,1:NCF)*P(J+1))*(P(J+1)-P(J))
	  END DO
          T1=DISTANCE*1.0D+03*PARSEC()
	  T1=2.0D0*PI*1.0D+23*(1.0E+10/T1)**2
	  YV(1:NCF)=YV(1:NCF)*T1
!
! NB: J and I have the same units, apart from per steradian.
!
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  CALL CURVE(NCF,XV,YV)
!
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=0.0D0
	  DO J=I,NP-1
	    YV(1:NCF)=YV(1:NCF)+0.5D0*(IP(J,1:NCF)*P(J)+IP(J+1,1:NCF)*P(J+1))*(P(J+1)-P(J))
	  END DO
          T1=DISTANCE*1.0D+03*PARSEC()
	  T1=2.0D0*PI*1.0D+23*(1.0E+10/T1)**2
	  YV(1:NCF)=YV(1:NCF)*T1
!
! NB: J and I have the same units, apart from per steradian/
!
	  WRITE(6,*)XV(1),YV(1)
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL CURVE(NCF,XV,YV)
	  WRITE(6,*)XV(1),YV(1)
!
	  YAXIS(J:J)='I'
	  J=INDEX(YAXIS,')')
	  YAXIS(J:)=' Jy'
!
	ELSE IF(X(1:2) .EQ. 'IP' .OR. X(1:2) .EQ. 'FP' .OR.
	1          X(1:3) .EQ. 'FAP')THEN
	  IF(X(1:2) .EQ. 'IP')THEN
	    WRITE(6,*)'Option to plot the intensity for a given impact parameter'
	  END IF
	  CALL USR_OPTION(I,'P',' ','Impact parameter index')
	  IF(I .GT. NP)THEN
	    WRITE(T_OUT,*)'Invalid depth; maximum value is',NP
	    GOTO 1
	  END IF
!
! If NP < ND, we must have a plane-parallel model.
!
	  IF(NP .LT. ND)THEN
	    WRITE(T_OUT,'(1X,A,1P,E14.6,A)')'       MU(I)=',MU(I)
	    WRITE(T_OUT,'(1X,A,1P,E14.6,A)')'  P(I)/R(ND)=',P(I)/R(ND)
	  ELSE
	    T1=P(I)*1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    WRITE(T_OUT,'(1X,A,1P,E14.6,A)')'       P(I)=',P(I)*1.0E+10,' cm'
	    WRITE(T_OUT,'(1X,A,1P,E14.6,A)')'       P(I)=',T1,' arcsec'
	    WRITE(T_OUT,'(1X,A,1P,E14.6)')'    P(I)/R*=',P(I)/R(ND)
	  END IF
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
!
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=IP(I,1:NCF)
!
! 2.3504D-11 = (AU[cm])^2 / d(cm)^2
!
	  IF(X(1:2) .EQ. 'FP')THEN
!
! Convert to flux per/pixel
!
	    YV=YV*SLIT_WIDTH*PIXEL_LENGTH*2.3504D-11
!
	  ELSE IF(X(1:3) .EQ. 'FAP')THEN
!
! Convert to flux per/arcsecond.
!
	    YV=YV*2.3504D-11
	  END IF
!
! NB: J and I have the same units, apart from per steradian/
!
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  J=INDEX(YAXIS,'J')
	  IF(X(1:2) .EQ. 'FP')THEN
            IF(Y_PLT_OPT .EQ. 'FNU')THEN
	      YAXIS='F(Jy\d \upixel\u-1\d)'
              IF(LOG_Y)YAXIS='Log F(Jy\d \upixel\u-1\d)'
	      YV(1:NCF)=YV(1:NCF)*1.0D+23
	    END IF
	  ELSE IF(X(1:3) .EQ. 'FAP')THEN
            IF(Y_PLT_OPT .EQ. 'FNU')THEN
	      YAXIS='F(Jy\d \uarcsec\u-2\d)'
              IF(LOG_Y)YAXIS='Log F(Jy\d \uarcsec\u-2\d)'
	      YV(1:NCF)=YV(1:NCF)*1.0D+23
	    END IF
	  ELSE
	    YAXIS(J:J)='I'
	    J=INDEX(YAXIS,')')
	    YAXIS(J:)=' \gW\u-1\d)'
	  END IF
!
	  CALL CURVE(NCF,XV,YV)
!
	ELSE IF(X(1:4) .EQ. 'SZ')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Start wavelength in Ang')
	  T1=0.299794D+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  CALL USR_OPTION(T1,'Lambda',' ','End wavelength in Ang')
	  T1=0.299794D+04/T1
          J=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(J)-T1 .GT. T1-NU(J+1))J=J+1
!
	  CALL USR_OPTION(DELV,'/\V','0.0D0','Smoothing size in km/s')
	  CALL USR_OPTION(FRAC,'/\V','0.5D0',
	1            'Fraction used to define radius (0<FRAC<1)')
!
	  K_ST=MIN(I,J)
	  K_END=MAX(I,J)
	  NX=K_END-K_ST+1
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  IF(ALLOCATED(ZV))DEALLOCATE(WV)
	  IF(ALLOCATED(ZV))DEALLOCATE(ZV)
	  ALLOCATE (XV(NX))
	  ALLOCATE (YV(NX))
	  ALLOCATE (ZV(NP))
	  ALLOCATE (WV(NP))
!
          WRITE(6,*)'Allocated arrays'
          WRITE(6,*)K_ST,K_END
          WRITE(6,*)NU(1),NU(NCF)
          WRITE(6,*)DELV
!
	  DO K=K_ST,K_END
	    XV(K-K_ST+1)=0.299794D+04/NU(K)
!
	    IF(DELV .EQ. 0.0D0)THEN
	      ZV(1:NP)=IP(1:NP,K)
	    ELSE
	      ZV(:)=0.0D0
	      DO ML=K-1,1,-1
	        IF( C_KMS*(NU(ML)/NU(K)-1.0D0) .GT. DELV)EXIT
                ZV(1:NP)=ZV(1:NP)+0.5D0*(NU(ML)-NU(ML+1))*(IP(1:NP,ML)+IP(1:NP,ML+1))
	      END DO
	      DO ML=K+1,NCF
	        IF( C_KMS*(NU(K)/NU(ML)-1.0D0) .GT. DELV)EXIT
                ZV(1:NP)=ZV(1:NP)+0.5D0*(NU(ML-1)-NU(ML))*(IP(1:NP,ML)+IP(1:NP,ML-1))
	      END DO
	    END IF
!
	    WV(1:NP)=0.0D0
	    DO I=1,NP-1
	      T1=0.5D0*(P(I)*ZV(I)+P(I+1)*ZV(I+1))*(P(I+1)-P(I))
	      WV(I+1)=WV(I)+T1
	    END DO
!
	    T1=WV(NP)
	    DO I=2,NP
	      WV(I-1)=WV(I)/T1
	    END DO
!
! If FRAC=0.8, we get the radius below which 80% of the light is emitted.
!
	    DO I=2,NP
	      IF(WV(I) .GT. FRAC)THEN
                 T1=(FRAC-WV(I-1))/(WV(I)-WV(I-1))
	         YV(K-K_ST+1)=(1.0D0-T1)*P(I-1)+T1*P(I)
	         EXIT
	      END IF
	    END DO
!
	  END DO 
!
          WRITE(6,*)'Found radii'
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
!
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=1,NX
	      YV(K)=LOG10(YV(K)*T1)
	    END DO
	    YAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    YV(1:NX)=LOG10(YV(1:NX)/T2)
	    YAXIS='Log P/R\d*\u'
	  END IF
	  XAXIS='\gl(\V)'
!
	  CALL CURVE(NX,XV,YV)
!
	ELSE IF(X(1:4) .EQ. 'IF2')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Start wavelength in Ang')
	  T1=0.299794D+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  CALL USR_OPTION(T1,'Lambda',' ','End wavelength in Ang')
	  T1=0.299794D+04/T1
          J=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(J)-T1 .GT. T1-NU(J+1))J=J+1
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  IF(ALLOCATED(ZV))DEALLOCATE(ZV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
	  ALLOCATE (ZV(NP))
!
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=1,NP-1
	      XV(K)=LOG10(P(K+1)*T1)
	    END DO
	    XAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    XV(1:NP-1)=LOG10(P(2:NP)/T2)
	    XAXIS='Log P/R\d*\u'
	  END IF
	  YV(:)=0.0D0
	  DO K=MIN(I,J),MAX(I,J)
	    YV(1:NP)=YV(1:NP)+IP(1:NP,K)*0.5D0*(NU(K-1)-NU(K+1))
	  END DO
	  T1=ABS(NU(I)-NU(J))
	  YV(1:NP)=YV(1:NP)/T1		!Normalize so per Hz
!
	  ZV(1:NP)=0.0D0
	  DO I=1,NP-1
	    T1=0.5D0*(P(I)*YV(I)+P(I+1)*YV(I+1))*(P(I+1)-P(I))
	    ZV(I+1)=ZV(I)+T1
	  END DO
	  T1=ZV(NP)
	  DO I=2,NP
	    ZV(I-1)=ZV(I)/T1
	  END DO
          T2=DISTANCE*1.0E+03*PARSEC()
	  T2=2.0D0*PI*1.0D+23*(1.0E+10/T2)**2
	  WRITE(6,'(A,ES12.3,A)')'The average flux in band is',T1*T2,'Jy'
!
	  YAXIS='F(p)'
	  CALL CURVE(NP-1,XV,ZV)
!
	ELSE IF(X(1:4) .EQ. 'INU2' .OR. X(1:5) .EQ. 'WINU2')THEN
	  CALL USR_OPTION(T1,'lam_st',' ','Start wavelength in Ang')
	  T1=0.299794D+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  CALL USR_OPTION(T1,'lam_end',' ','End wavelength in Ang')
	  T1=0.299794D+04/T1
          J=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(J)-T1 .GT. T1-NU(J+1))J=J+1
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
	  CALL USR_OPTION(MULT_BY_PSQ,'PSQ','F','Multiply by P^2?')
	  MULT_BY_P=.FALSE.
	  IF(.NOT. MULT_BY_PSQ)THEN
	    CALL USR_OPTION(MULT_BY_P,'P','F','Multiply by P?')
	  END IF
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
!
! In this case we return linear axes --- usefule for SN.
!
	  IF(MULT_BY_P)THEN
	    IF(USE_ARCSEC)THEN
	      T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	      XV(1:NP)=P(1:NP)*T1 
	      XAXIS='P(")'
	    ELSE
	      XV(1:NP)=P(1:NP)/R(ND)
	      XAXIS='P/R\d*\u'
	    END IF
!
! Average I(p) over frequnecy.
!
	    YV(:)=0.0D0
	    K=MIN(I,J); J=MAX(I,J); I=K
	    IF(J .EQ. I)J=I+1
	    DO K=I,J-1
	      YV(1:NP-1)=YV(1:NP-1)+0.5D0*(IP(1:NP-1,K)+IP(1:NP-1,K+1))*
	1                                 (NU(K)-NU(K+1))
	    END DO
	    T1=ABS(NU(I)-NU(J))
	    YV(1:NP-1)=YV(1:NP-1)/T1
!
! Now multilply by P. We keep out the factor of 10^10 (arsing from the
! units of P) as it keep the numbers closer to 1.
!
	    YV(1:NP-1)=YV(1:NP-1)*P(1:NP-1)
	    YAXIS='10\u-10 \dpI\d\gn\u(ergs cm\u-1\d s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	    CALL CURVE(NP-1,XV,YV)
!
	  ELSE
	    IF(USE_ARCSEC)THEN
	      T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	      DO K=1,NP-2
	        XV(K)=LOG10(P(K+1)*T1)
	      END DO
	      XAXIS='Log P(")'
	    ELSE
	      T2=R(ND)
	      XV(1:NP-2)=LOG10(P(2:NP-1)/T2)
	      XAXIS='Log P/R\d*\u'
	    END IF
	    YV(:)=0.0D0
	    K=MIN(I,J); J=MAX(I,J); I=K
	    IF(J .EQ. I)J=I+1
!
	    IF(X(1:5) .EQ. 'WINU2')THEN
	      WRITE(30,*)NP-1,J-I+2
	      WRITE(30,'(1000ES14.6)')P(2:NP)
	      WRITE(30,'(1000ES15.7)')XV(1:NP-1)
	      WRITE(30,'(1000ES15.7)')ANG_TO_HZ/NU(I:J+1)
	      DO K=I,J+1
	        WRITE(30,'(500ES12.4)')IP(2:NP,K)
	      END DO
	    ELSE
	      DO K=I,J-1
	        YV(1:NP-2)=YV(1:NP-2)+0.5D0*(IP(2:NP-1,K)+IP(2:NP-1,K+1))*
	1                                 (NU(K)-NU(K+1))
	      END DO
	      T1=ABS(NU(I)-NU(J))
	      YV(1:NP-2)=LOG10(YV(1:NP-2)/T1)
	      IF(MULT_BY_PSQ)THEN
	        DO I=1,NP-2
	          YV(I)=YV(I)+2.0D0*LOG10(P(I+1))+20.0D0
	        END DO
	        YAXIS='Log p\u2\dI\gn\u(ergs s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	      ELSE
	        YAXIS='Log I\d\gn\u(ergs cm\u-2\d s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	      END IF
	      XAXIS='\gm'
	      CALL CURVE(NP-2,XV,YV)
	    END IF
	  END IF
!
	ELSE IF(X(1:3) .EQ. 'IMU')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Wavelength in Ang')
	  T1=0.299794E+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
	  YV(1:NP)=IP(1:NP,I)
!
	  IF(MU(1) .EQ. 0.0D0)THEN
	    CALL USR_OPTION(T1,'Rstar',' ','Radius of star in units of 10^10cm')
	    DO I=1,NP
	      T2=1.0D0-(P(I)/T1)**2
	      IF(T2 .GT. 0.0D0)MU(I)=SQRT(T2)
	      IF(T2 .LT. 0.0D0)MU(I)=-SQRT(-T2)
	    END DO
	    XV(1:NP)=MU(1:NP)
	    MU=0.0D0
	  ELSE
	    XV(1:NP)=MU(1:NP)
	  END IF
	  T1=MAXVAL(YV(1:NP))
	  YAXIS='I(ergs cm\u-2\d s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	  IF(T1 .LT. 1.0D-03 .OR. T1 .GT. 1.0D+04)THEN
	    J=-LOG10(T1)
	    YV(1:NP)=YV(1:NP)*(10**J)
	    WRITE(TMP_STR,*)J; TMP_STR=ADJUSTL(TMP_STR)
	    YAXIS=TRIM(YAXIS)//' x 10\u'//TRIM(TMP_STR)//'\d'
	  END IF
	  CALL CURVE(NP,XV,YV)

	ELSE IF(X(1:3) .EQ. 'INU')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Wavelength in Ang')
	  T1=0.299794E+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','F','Use arcseconds?')
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=2,NP-1
	      XV(K)=LOG10(P(K+1)*T1)
	    END DO
	    XV(1)=XV(2)-3
	    XAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    XV(2:NP-1)=LOG10(P(2:NP-1)/T2)
	    XAXIS='Log P/R\d*\u'
	    XV(1)=XV(2)-2
	  END IF
	  YV(1:NP-1)=IP(1:NP-1,I)
	  DO I=1,NP-1
	    WRITE(6,*)I,XV(I),YV(I)
	  END DO
	  YAXIS='I\d\gn\u(ergs cm\u-2\d s\u-1\d Hz\u-1\d steradian\u-1\d)' 
	  CALL CURVE(NP-1,XV,YV)
!
	ELSE IF(X(1:3) .EQ. 'WIP')THEN
	  CALL USR_OPTION(T1,'Lambda',' ','Start wavelength in Ang')
	  CALL USR_OPTION(T2,'Lambda',' ','End wavelength in Ang')
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
	  LAMBDA=NINT(T1-1.0D0)
	  DO WHILE(LAMBDA .LE. T2)
	    LAMBDA=LAMBDA+1
	    T1=0.299794E+04/LAMBDA
            I=GET_INDX_DP(T1,NU,NCF)
	    IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
	    YV(1:NP-2)=IP(2:NP-1,I)
	    J=NINT(LAMBDA)
	     TMP_STR='J'
	     WRITE(TMP_STR(2:5),'(I4.4)')J
	      OPEN(FILE=TRIM(TMP_STR),UNIT=12,STATUS='UNKNOWN',ACTION='WRITE')
	        DO I=1,NP
	          WRITE(12,'(F12.7)')YV(I)
	        END DO 
	      CLOSE(UNIT=12)
	  END DO
!
	ELSE IF(X(1:2) .EQ. 'CF')THEN
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  IF(ALLOCATED(TA))DEALLOCATE(TA)
	  ALLOCATE (XV(NP))
	  ALLOCATE (YV(NP))
	  ALLOCATE (TA(NP))
!
	  TA(1:ND)=0.0D0
	  DO ML=1,NCF-1
	    T1=0.5D0*(NU(ML)-NU(ML+1))
	    DO I=1,NP
	      TA(I)=TA(I)+T1*(IP(I,ML)+IP(I,ML+1))
	    END DO
	  END DO
!
	  CALL USR_OPTION(USE_ARCSEC,'Arcsec','T','Use arcseconds?')
	  IF(USE_ARCSEC)THEN
	    T1=1.0E+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    DO K=1,NP-2
	      XV(K)=LOG10(P(K+1)*T1)
	    END DO
	    XAXIS='Log P(")'
	  ELSE
	    T2=R(ND)
	    XV(1:NP-2)=LOG10(P(2:NP-1)/T2)
	    XAXIS='Log P/R\d*\u'
	  END IF
	  DO I=1,NP
	    YV(I)=TA(I)/TA(1)
	  END DO
!	  CALL CNVRT_J(XV,YV,ND,LOG_X,LOG_Y,' ',Y_PLT_OPT,
!	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  YAXIS='I(ergs cm\u-2\d s\u-1\d steradian\u-1\d)' 
	  CALL CURVE(ND,XV,YV)
! 
!
! Read in observation data as done in PLT_SPEC
!
	ELSE IF(X(1:6) .EQ. 'RD_OBS')THEN
	  FILENAME=' ' 
	  CALL USR_OPTION(FILENAME,'File',' ',' ')
!
	  SCALE_FAC=1.0D0
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  ADD_FAC=0.0D0
	  CALL USR_HIDDEN(ADD_FAC,'ADD','0.0D0',' ')
!
	  RAD_VEL=0.0D0
	  CALL USR_HIDDEN(RAD_VEL,'RAD_VEL','0.0D0',
	1             'Radial velcoity (+ve if away)')
!
	  CLEAN=.FALSE.
	  CALL USR_HIDDEN(CLEAN,'CLEAN','F',' ')
!
	  SMOOTH=.FALSE.
	  CALL USR_HIDDEN(SMOOTH,'SMOOTH','F',' ')
!
	  IF(SMOOTH)THEN
	    K=5
	    CALL USR_OPTION(K,'HAN','5','Number of points for HAN [ODD]')
	    K=2*(K/2)+1 !Ensures odd.
	  END IF
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  IF(ALLOCATED(ZV))DEALLOCATE(ZV)
	  NCF_MAX=100000
	  ALLOCATE (XV(NCF_MAX))
	  ALLOCATE (YV(NCF_MAX))
!
	  CALL USR_HIDDEN(OBS_COLS,2,2,'COLS','1,2','Columns with data')
	  CALL RD_OBS_DATA_V2(XV,YV,NCF_MAX,J,FILENAME,OBS_COLS,IOS)
	  IF(IOS .NE. 0)GOTO 1          !Get another option
	  DO I=1,J
	    XV(I)=ANG_TO_HZ/XV(I)
	    YV(I)=YV(I)*SCALE_FAC+ADD_FAC
	  END DO
!
	  IF(RAD_VEL .NE. 0)THEN
	    DO I=1,J
	      XV(I)=XV(I)*(1.0D0-1.0D+05*RAD_VEL/C_CMS)
	    END DO
	  END IF
!
! Procdure to remove single pizels that are zero due quirks with IUE.
!
	  IF(CLEAN)THEN
	    DO I=1,J
	      ZV(I)=YV(I)
	    END DO
	    DO I=2,J-1
	      IF(YV(I) .EQ. 0)THEN
	        YV(I)=0.5D0*(ZV(I-1)+ZV(I+1))
	      END IF
	    END DO
	  END IF
!
! We check whether the X axis is monotonic. If not, we smooth each section
! separately. Designed for non-merged overlapping ECHELLE orders.
!
	  NON_MONOTONIC=.FALSE.
	  IF(SMOOTH)THEN
	    T2=XV(2)-XV(1)
	    DO I=1,J-1
	      IF( (XV(I)-XV(I+1))*T2 .LT. 0)THEN
	        NON_MONOTONIC=.TRUE.
	        EXIT
	      END IF
	    END DO
	  END IF
!
          IF(SMOOTH .AND. NON_MONOTONIC)THEN
	    ALLOCATE (ZV(J))
            ZV(1:J)=YV(1:J)
            DO I=1,K
              WT(I)=FAC(K-1)/FAC(I-1)/FAC(K-I)
            END DO
            IST=1
            IEND=0
            T2=XV(2)-XV(1)
            DO WHILE(IEND .LT. J)
              IEND=J
              DO I=IST,J-1
                IF( (XV(I+1)-XV(I))*T2 .LT. 0)THEN
                  IEND=I
                  EXIT
                END IF
              END DO
              WRITE(6,*)IST,IEND
              DO I=IST,IEND
                T1=0.0D0
                YV(I)=0.0D0
                DO L=MAX(IST,I-K/2),MIN(IEND,I+k/2)
                  ML=L-I+K/2+1
                  T1=T1+WT(ML)
                  YV(I)=YV(I)+ZV(L)*WT(ML)
                END DO
                YV(I)=YV(I)/T1
              END DO
              IST=IEND+1
            END DO
          ELSE IF(SMOOTH)THEN
            DO I=1,J
              ZV(I)=YV(I)
            END DO
            DO I=1,K
              WT(I)=FAC(K-1)/FAC(I-1)/FAC(K-I)
            END DO
            DO I=1,J
              T1=0.0D0
              YV(I)=0.0D0
              DO L=MAX(1,I-K/2),MIN(J,I+k/2)
                ML=L-I+K/2+1
                T1=T1+WT(ML)
                YV(I)=YV(I)+ZV(L)*WT(ML)
              END DO
              YV(I)=YV(I)/T1
            END DO
          END IF
!
	  CALL CNVRT(XV,YV,J,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
          CALL CURVE(J,XV,YV)
!
	ELSE IF(X(1:3) .EQ. 'SLT')THEN
	  CALL USR_OPTION(X_CENT,'XC','0.0D0','X center of aperture (across width)')
	  CALL USR_OPTION(Y_CENT,'YC','0.0D0','Y center of aperture (along length)')
	  CALL USR_OPTION(S_WIDTH,'SW','0.05D0','Slit width (arcsec)')
	  CALL USR_OPTION(S_LNGTH,'SL','0.1D0','Slit length (arcsec)')
	  CALL USR_OPTION(TEL_FWHM,'FWHM','0.1D0','Telescope FWHM (arcsec)')
	  APP_SIZE=S_WIDTH*S_LNGTH 		!In square arcseconds
	  CALL USR_OPTION(NINS,'NINS','2','# of points to insert to improve accuracy')
!
	  CALL USR_OPTION(T1,'lam_st',' ','Start wavelength in Ang')
	  T1=0.299794D+04/T1
          I=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(I)-T1 .GT. T1-NU(I+1))I=I+1
!
	  CALL USR_OPTION(T1,'lam_end',' ','End wavelength in Ang')
	  T1=0.299794D+04/T1
          J=GET_INDX_DP(T1,NU,NCF)
	  IF(NU(J)-T1 .GT. T1-NU(J+1))J=J+1
!
! Compute the intensity across the band.
!
	  IF(ALLOCATED(P_NEW))DEALLOCATE(P_NEW)
	  IF(ALLOCATED(IP_NEW))DEALLOCATE(IP_NEW)
!
	  IF(TEL_FWHM .EQ. 0.0D0)THEN
	    ALLOCATE (P_NEW(NP))
	    ALLOCATE (IP_NEW(NP))
	    IP_NEW=0.0D0
	    DO ML=I,J
	      IP_NEW(:)=IP_NEW(:)+IP(:,ML)
	    END DO
	    P_NEW=P
	    NP_NEW=NP
	  ELSE
!
	    NP_NEW=2*NP-1
	    ALLOCATE (P_NEW(NP_NEW))
	    ALLOCATE (IP_NEW(NP_NEW))
	    IP_NEW=0.0D0
	    DO ML=I,J
	      DO LS=1,NP
	        IP_NEW(NP+LS-1)=IP_NEW(NP+LS-1)+IP(LS,ML)
	      END DO
	    END DO
	    DO LS=1,NP
	      IP_NEW(NP-LS+1)=IP_NEW(NP+LS-1)
	      P_NEW(NP+LS-1)=P(LS)
	      P_NEW(NP-LS+1)=-P_NEW(NP+LS-1)
	    END DO
!
	    T1=1.0D+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    MOD_RES=0.1D0*TEL_FWHM
	    WRITE(6,*)'MOD_RES=',MOD_RES
	    NTMP=T1*2*P(NP)/MOD_RES+1
	    IF(NTMP .LT. 2*NP)NTMP=2*NP
	    NP_NEW=2*NP-1
!
	    ALLOCATE (TEMP_P(NTMP))
	    ALLOCATE (TEMP_IP(NTMP))
	    WRITE(6,*)'NTMP=',NTMP
!
	    TEMP_P(1:NP_NEW)=T1*P_NEW(1:NP_NEW)
	    TEMP_IP(1:NP_NEW)=IP_NEW(1:NP_NEW)
	    CALL DP_CURVE(NP_NEW,TEMP_P,TEMP_IP)
	    WRITE(6,*)'Begin linearize'
	    CALL LINEARIZE_V2(TEMP_P,TEMP_IP,NP_NEW,NTMP,MOD_RES)
	    WRITE(6,*)'Done linearize'
	    T1=TEL_FWHM/2.35482
	    CALL CONVOLVE(TEMP_P,TEMP_IP,NTMP,T1,RZERO,RZERO,L_FALSE)
	    CALL DP_CURVE(NTMP,TEMP_P,TEMP_IP)
	    WRITE(6,*)'Done convolution'
	    T1=1.0D+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	    TEMP_P=TEMP_P/T1
	    IP_NEW=0.0D0
	    CALL MAP(TEMP_P,TEMP_IP,NTMP,P,IP_NEW,NP)
	    TEMP_P(1:NP)=P(1:NP)*T1
	    CALL DP_CURVE(NP,TEMP_P,IP_NEW)
	    CALL GRAMON_PGPLOT(' ',' ',' ',' ')
!
	    DEALLOCATE (TEMP_P)
	    DEALLOCATE (TEMP_IP)
	  END IF
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
!
	  WRITE(6,*)'Calling INT_SEQ'
	  WRITE(6,*)X_CENT,Y_CENT,S_WIDTH,S_LNGTH
!
	  K=1.0D0/S_WIDTH+1
	  T1=1.0D+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	  Y_CENT=Y_CENT/T1	
	  S_WIDTH=S_WIDTH/T1	
	  S_LNGTH=S_LNGTH/T1	
!
	  DO I=1,K
	    XV(I)=X_CENT+S_WIDTH*(I-1)*T1
	    T2=(X_CENT+T1*S_WIDTH*(I-1))/T1
	    CALL INT_REC_AP(IP_NEW,P,YV(I),T2,Y_CENT,S_WIDTH,S_LNGTH,
	1                     NP,IONE,NINS)
	  END DO
	  CALL CURVE(K,XV,YV)
!
	ELSE IF(X(1:2) .EQ. 'SQ')THEN
	  CALL USR_OPTION(X_CENT,'XC','0.0D0','X center of aperture (across width)')
	  CALL USR_OPTION(Y_CENT,'YC','0.0D0','Y center of aperture (along length)')
	  CALL USR_OPTION(S_WIDTH,'SW','0.1D0','Slit width (arcsec)')
	  CALL USR_OPTION(S_LNGTH,'SL','0.2D0','Slit length (arcsec)')
	  APP_SIZE=S_WIDTH*S_LNGTH 		!In square arcseconds
	  CALL USR_OPTION(NINS,'NINS','2','# of points to insert to improve accuracy')
!
	  T1=1.0D+10*206265.0D0/(DISTANCE*1.0E+03*PARSEC())
	  X_CENT=X_CENT/T1	
	  Y_CENT=Y_CENT/T1	
	  S_WIDTH=S_WIDTH/T1	
	  S_LNGTH=S_LNGTH/T1	
!
	  IF(ALLOCATED(XV))DEALLOCATE(XV)
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (XV(NCF))
	  ALLOCATE (YV(NCF))
	  XV(1:NCF)=NU(1:NCF)
!
	  WRITE(6,*)'Calling INT_SEQ'
	  WRITE(6,*)X_CENT,Y_CENT,S_WIDTH,S_LNGTH
	  CALL INT_REC_AP(IP,P,YV,X_CENT,Y_CENT,S_WIDTH,S_LNGTH,
	1               NP,NCF,NINS)
!
! NB: J and I have the same units, apart from per steradian/
!
	  CALL CNVRT_J(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  T2=1.0D+03*PARSEC()*DISTANCE
	  T1=1.0D+20/T2/T2/APP_SIZE
	  YV(1:NCF)=YV(1:NCF)*T1
          IF(Y_PLT_OPT .EQ. 'FNU')THEN
	    YV(1:NCF)=YV(1:NCF)*1.0D+23
	    YAXIS='F(Jy\d \uarcsec\u-2\d)'
            IF(LOG_Y)YAXIS='Log F(Jy\d \uarcsec\u-2\d)'
	  END IF
	  CALL CURVE(NCF,XV,YV)
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
	  IF(X(1:2) .EQ. 'LI')THEN
	    CALL GEN_ASCI_OPEN(LU_IN,'PLT_IP_OPT_DESC','OLD',' ','READ',
	1                  IZERO,IOS)
	  ELSE 
	    CALL GEN_ASCI_OPEN(LU_IN,'PLT_IP_OPTIONS','OLD',' ','READ',
	1                  IZERO,IOS)
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening HELP or DESCRIPTER file for PLT_JH'
	    WRITE(6,*)' '
	    WRITE(6,*)' TIT  LX  LY  XU '
	    WRITE(6,*)' CF:   Plot cummulative spectrum as a function of impact parameter'
	    WRITE(6,*)' IP:   Plot spectrum at a given impact parameter'
	    WRITE(6,*)' IMU:  Plot intensity as a function of MU for a given frequency'
	    WRITE(6,*)' SP:   Plot spectrum inside and outside impact parameter p'
	    WRITE(6,*)' INU:  Plot I(p) for a given frequency'
	    WRITE(6,*)' INU2: Plot I(p) for a given frequency band'
	    WRITE(6,*)' IF2:  Plot normalize Flux originating inside p for a given frequency band'
	    WRITE(6,*)' '
	    WRITE(6,*)' GR:   Enter polot packgae to plot passed data'
	    WRITE(6,*)' GRL:  As for GR but no labels passed'
	    WRITE(6,*)' '
	    GOTO 1
	  END IF
	  READ(LU_IN,*)I,K			!For page formating (I=22,K=12)
	  DO WHILE(1.EQ. 1)
	    DO J=1,I
	      READ(LU_IN,'(A)',END=700)STRING
	      L=LEN(TRIM(STRING))
	      IF(L .EQ. 0)THEN
	        WRITE(T_OUT,'(X)')
	      ELSE
	        WRITE(T_OUT,'(1X,A)')STRING(1:L)
	      END IF
	    END DO
	    READ(T_IN,'(A)')STRING
	    IF(STRING(1:1) .EQ. 'E' .OR.
	1       STRING(1:1) .EQ. 'e')GOTO 700		!Exit from listing.
	    I=K
	  END DO
700	  CONTINUE
	  CLOSE(UNIT=LU_IN)
!
	ELSE IF(X(1:3) .EQ. 'WRP') THEN
	  WRITE(6,'(3X,A,11X,A,13X,A,11X,A)')'LS','P(LS)','MU(LS)','HQW(LS)'
	  DO LS=1,NP
	    WRITE(6,'(I5,ES16.8,1X,2F18.12)')LS,P(LS),MU(LS),HQW_AT_RMAX(LS)
	  END DO
	ELSE IF(X(1:2) .EQ. 'EX') THEN
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
