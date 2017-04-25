C
C General routine for plotting and comparing observational and model 
C spectra. Model input spectra are read from OBSLFLUX.
C
C Various options are available to redden and normalize the model spectra. 
C Several different units can be used for the X and Y axes.
C
	PROGRAM PLT_SPEC
C
	USE FILT_PASS_BAND
C
C Interface routines for IO routines.
C
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE GEN_IN_INTERFACE
C
	IMPLICIT NONE
C
C Altered 18-Jan-2015 : Added REM_BP (automatic removal of bad pixel).
C Altered 04-Nov-2015 : Fixed SMC redenning law in optical. Previous formula only
C                         valid for UV. Joined UV smoothly (at 2950A) to CCM law with R=2.74.
C Altered 15-Mar-2011 : SMC reddening law added (done by Kathryn Neugent).
C                       RED option installed -- crude method to get redenning.
C                       Crude procedure to remove cosmic rays (or spikes) form observed
C                         data implemented with rd_obs option.
C                       Fixed comvolution option to make more transparent.
C Altered 17-Jun-1996 : Bug fixed with Wavelength normalization for BB option.
C Altered 16-mar-1997 : Cleaned: USR_OPTION installation finalized.
C Altered 26-Nov-1996 : Norm option fixed so that entire MOD spectrum is 
C                        plotted.
C
C
C Determines largest single plot that can read in.
C
	INTEGER, PARAMETER :: NCF_MAX=3000000
C
C Used to indicate number of data points in BB spectrum
c
	INTEGER, PARAMETER :: NBB=2000
C
	INTEGER NCF		!Number of data points in default data
	INTEGER NCF_MOD		!Used when plotting another model data set
c
	REAL*8 NU(NCF_MAX)
	REAL*8 OBSF(NCF_MAX)
	REAL*8 FQW(NCF_MAX)
	REAL*8 AL_D_EBmV(NCF_MAX)
C
	INTEGER NCF_CONT
	REAL*8 NU_CONT(NCF_MAX)
	REAL*8 OBSF_CONT(NCF_MAX)
!
	INTEGER NOBS
	REAL*8 NU_OBS(NCF_MAX)
	REAL*8 OBSF_OBS(NCF_MAX)
C
C Indicates which columns the observatoinal data is in.
C
	INTEGER OBS_COLS(2)
C
C Vectors for passing data to plot package via calls to CURVE.
C
	REAL*4 XV(NCF_MAX)
	REAL*4 YV(NCF_MAX)
	REAL*4 ZV(NCF_MAX)
	REAL*4 SP_VAL
C
	CHARACTER*80 NAME		!Default title for plot
	CHARACTER*80 XAXIS,XAXSAV	!Label for Absisca
	CHARACTER*80 YAXIS		!Label for Ordinate
C
	REAL*8 LUM
	REAL*8 TOT_LUM
	REAL*8 N_PHOT
	REAL*8 ZAN_FREQ(10)
C
	REAL*8 ANG_TO_HZ
	REAL*8 KEV_TO_HZ
	REAL*8 NORM_WAVE
	REAL*8 NORM_FREQ
	REAL*8 NORM_FLUX
	REAL*8 DNU
	REAL*8 BB_FLUX
	REAL*8 SCALE_FAC
	REAL*8 XFAC
	REAL*8 ADD_FAC
	REAL*8 LAMC
	REAL*8 RAD_VEL			!Radial velcity in km/s
	REAL*8 C_CMS
	REAL*8 WT(30)			!Used when smoothing observational data.
	REAL*8 LAM_RANGE(2)
C
	LOGICAL NON_MONOTONIC
	LOGICAL SMOOTH			!Smooth observational data?
	LOGICAL CLEAN			!Remove IUE bad pixels?
	LOGICAL REMOVE_BAD_PIX
	LOGICAL CLN_CR			!Remove cosmic-ray spikes
	LOGICAL TREAT_AS_MOD		
	LOGICAL READ_OBS
	LOGICAL AIR_LAM
C
	LOGICAL WR_PLT,OVER
C
	LOGICAL LIN_INT
	LOGICAL UNEQUAL
	LOGICAL LOG_X,LOG_Y
	CHARACTER*10 Y_PLT_OPT,X_UNIT
	CHARACTER*80 IS_FILE
	CHARACTER*80 FILENAME
	CHARACTER*80 DIRECTORY 
	CHARACTER*80 XKEY,YKEY
!
! Variable for applying interstellar absorption to model spectrum.
!
	LOGICAL HI_ABS                   ! correct for HI absorption
	LOGICAL H2_ABS                   ! correct for H2 absorption
	REAL*8 T_IN_K                    ! temp in K of intersteallar H&HII
	REAL*8 V_TURB                    ! turbulent velocity of "     "
	REAL*8 V_R                       ! radial v of star w.r.t. ISM
	REAL*8 LOG_NTOT                  ! log of H column density
	REAL*8 LOG_H2_NTOT               ! log of H2 column density
	LOGICAL FFT_CONVOLVE             ! use FFT method to convolve data
	REAL*8 INST_RES                  ! desired instrument resolution (dl) 
	REAL*8 MIN_RES_KMS               ! minimum resolution for model data
	REAL*8 NUM_RES                   ! number of res. elements to cons.
	REAL*8 RESOLUTION                ! desired resolution (R=l/dl)
	REAL*8 WAVE_MAX                  ! max wavelength to convolve over
	REAL*8 WAVE_MIN                  ! min wavelength to convolve over
	REAL*8 VSINI
	REAL*8 EPSILON
C
C Miscellaneous variables.
C
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,K,L,ML,MLST
	INTEGER IST,IEND
	INTEGER CNT
	INTEGER NHAN
	REAL*8 T1,T2,T3
	REAL*8 SUM
	REAL*8 TEMP
	REAL*8 TMP_FREQ
C
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_IN=5		!For file I/O
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: LU_IN=10	!For file I/O
	INTEGER, PARAMETER :: LU_OUT=11
C
	REAL*8, PARAMETER :: EDGE_HYD=3.28808662499619D0
	REAL*8, PARAMETER :: EDGE_HEI=5.94520701882481D0
	REAL*8, PARAMETER :: EDGE_HE2=13.1581564178623D0
C
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
	LOGICAL ADD_NOISE
	INTEGER IRAN
	REAL*8 R_SEED
	REAL*8 COUNTS
	REAL*8 POIDEV
C
C USR_OPTION variables
C
	CHARACTER MAIN_OPT_STR*80	!Used for input of the main otion
	CHARACTER X*10			!Used for the idividual option
	CHARACTER STRING*80
	CHARACTER*120 DEFAULT
	CHARACTER*120 DESCRIPTION
C
	REAL*8 SPEED_OF_LIGHT,FAC,LAM_VAC
	INTEGER GET_INDX_SP
	LOGICAL EQUAL
	CHARACTER*30 UC
	CHARACTER*30 FILTER_SET 
	EXTERNAL SPEED_OF_LIGHT,EQUAL,FAC,UC,LAM_VAC,GET_INDX_DP
C
C 
C Filter and extinction data:
C
	REAL*8 NORM_LUM
	REAL*8 DIST
	REAL*8 R_EXT
	REAL*8 EBMV_GAL
	REAL*8 EBMV_LMC
	REAL*8 EBMV_CCM
	REAL*8 EBMV_SMC         !KN - SMC law
	REAL*8 RAX,RBX
	REAL*8 FILTLAM(8),FILTZP(8),FLAM(24),ZERO_POINT
	CHARACTER*1 FILT(8)
	REAL*8 C1,C2,C3,C4,D,F  !KN For SMC
!
	REAL*8 RESPONSE1
	REAL*8 RESPONSE2
	REAL*8 FILT_INT_BEG
	INTEGER IF
	INTEGER ML_ST 
	INTEGER ML_END

C
	DATA FILT/' ','u','b','v','r','J','H','K'/
	DATA FILTLAM/0.3644,0.3650,0.427,0.5160,0.60,1.25,1.65,2.20/
	DATA FLAM/0.1,.125,0.15,0.175,0.20,0.225,0.25,0.275,0.300,0.325
	1, 0.350,0.4,0.5,0.6,0.7,0.8,0.9,1.25,1.65,2.2,3.6,4.8,10.0,20.0/
C
C Zero point is Based on Flagstaff calibration which gives F_nu(Alpha Lyrae)=3560Jy. Note the
C V magnitude of (Alpha Lyrae) is 0.03.
C
	DATA ZERO_POINT/3560/
	DATA FILTZP/3560.0,3560.0,3560.0,3560.0,3560.0,1564.0,1008.0,628.0/
C
C 
C
C Set constants.
C
	CHIBF=2.815D-06
	CHIFF=3.69D-29
	HDKT=4.7994145D0
	TWOHCSQ=0.0147452575D0
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
C
	C_CMS=SPEED_OF_LIGHT()
C
C Set defaults.
C
	XAXIS=' '
	XAXSAV=' '
	YAXIS=' '
	NAME=' '
!	PI=FUN_PI()
C
	LOG_X=.FALSE.
	LOG_Y=.FALSE.
	X_UNIT='ANG'
	Y_PLT_OPT='FNU'
!
	NCF=0
	NOBS=0
C
C Conversion factor from Kev to units of 10^15 Hz.
C Conversion factor from Angstroms to units of 10^15 Hz.
C
	KEV_TO_HZ=0.241838E+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
C
C  Read in default model.
C
	NCF=0			!In case no model read
	FILENAME='OBSFLUX'
	CALL GEN_IN(FILENAME,'Model file')
	CALL RD_MOD(NU,OBSF,NCF_MAX,NCF,FILENAME,IOS)
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
	IF(MAIN_OPT_STR .EQ. ' ')GOTO 3
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
C                    
C Set X-Ais plotting options.
C
	ELSE IF(X(1:2) .EQ.'LX' .OR. X(1:4) .EQ. 'LOGX' .OR. 
	1                            X(1:4) .EQ. 'LINX')THEN
	  LOG_X=.NOT. LOG_X
	  IF(LOG_X)WRITE(T_OUT,*)'Now using Logarithmic X axis'
	  IF(.NOT. LOG_X)WRITE(T_OUT,*)'Now using Linear X axis'
	ELSE IF(X(1:2) .EQ.'XU' .OR. X(1:6) .EQ. 'XUNITS')THEN
	  CALL USR_OPTION(X_UNIT,'X_UNIT','Ang',
	1                  'Ang, AA[Air Ang], um, eV, keV, Hz, Mm/s, km/s')
	  CALL SET_CASE_UP(X_UNIT,IZERO,IZERO)
	  IF(X_UNIT .NE. 'ANG' .AND.
	1        X_UNIT .NE. 'AA' .AND.
	1        X_UNIT .NE. 'UM' .AND.
	1        X_UNIT .NE. 'EV' .AND.
	1        X_UNIT .NE. 'KEV' .AND.
	1        X_UNIT .NE. 'HZ' .AND.
	1        X_UNIT .NE. 'MM/S' .AND.
	1        X_UNIT .NE. 'KM/S')THEN
	     WRITE(T_OUT,*)'Invalid X unit: Try again'
	   END IF
C
C NB: We offer the option to use the central frequency to avoid
C air/vacuum confusions. Model data is in vacuum wavelngths, which
C we use in plotting at all wavelngths.
C
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
C
C Set Y axis plotting options.
C
	ELSE IF(X(1:2) .EQ. 'LY' .OR. X(1:4) .EQ. 'LOGY' .OR. 
	1                             X(1:4) .EQ. 'LINY')THEN
	  LOG_Y=.NOT. LOG_Y
	  IF(LOG_Y)WRITE(T_OUT,*)'Now using Logarithmic Y axis'
	  IF(.NOT. LOG_Y)WRITE(T_OUT,*)'Now using Linear Y axis'
	ELSE IF(X(1:2) .EQ.'YU' .OR. X(1:6) .EQ. 'YUNITS')THEN
	  CALL USR_OPTION(Y_PLT_OPT,'Y_UNIT',' ',
	1          'FNU, NU_FNU, FLAM, LAM_FLAM')
	  CALL SET_CASE_UP(Y_PLT_OPT,IZERO,IZERO)
	  IF(Y_PLT_OPT .NE. 'FNU' .AND.
	1        Y_PLT_OPT .NE. 'NU_FNU' .AND.
	1        Y_PLT_OPT .NE. 'LAM_FLAM' .AND.
	1        Y_PLT_OPT .NE. 'FLAM')THEN
	     WRITE(T_OUT,*)'Invalid Y Plot option: Try again'
	  END IF
C
C 
C
	ELSE IF(X(1:6) .EQ. 'RD_MOD')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Model file')
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  CALL USR_HIDDEN(XFAC,'XFAC','1.0D0',' ')
	  CALL USR_HIDDEN(RAD_VEL,'RAD_VEL','0.0D0','Radial velocity(km/s) of star')
	  IF(XFAC .NE. 1.0D0 .AND. RAD_VEL .NE. 0.0D0)THEN
	    WRITE(6,*)'Only one of XFAC and RAD_VEL can be changed from their default values of 1 and 0'
	    GOTO 1
	  ELSE IF(RAD_VEL .NE. 0.0D0)THEN
	    XFAC=(1.0D0+1.0D+05*RAD_VEL/C_CMS)
	  END IF
	  IF(OVER)THEN
C
C This option allows all normal model options to be done on the data
C (e.g. redenning).
C
	    CALL RD_MOD(NU,OBSF,NCF_MAX,NCF,FILENAME,IOS)
	    IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error reading model data -- no data read'
	       GOTO 1		!Get another option
	    END IF
	    OVER=.FALSE.
!
	    IF(SCALE_FAC .NE. 1.0D0 .OR. XFAC .NE. 1.0D0)THEN
	      OBSF(1:NCF)=OBSF(1:NCF)*SCALE_FAC
	      NU(1:NCF)=NU(1:NCF)*XFAC
	      WRITE(T_OUT,*)'Model has been scaled!'
	    ELSE
	      WRITE(T_OUT,*)'No scaling done with new model data'
	    END IF
! 
	    WRITE(T_OUT,*)'New model data replaces old data'
	    WRITE(T_OUT,*)'No plots done with new model data'
	  ELSE
C
C This option is now similar to RD_CONT
C 
	    CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_MOD,FILENAME,IOS)
	    IF(IOS .NE. 0)THEN
	       WRITE(T_OUT,*)'Error reading model data -- no data read'
	       GOTO 1		!Get another option
	    END IF
	    DO I=1,NCF_MOD
	      XV(I)=NU_CONT(I)*XFAC
	      YV(I)=OBSF_CONT(I)*SCALE_FAC
	    END DO
	    CALL CNVRT(XV,YV,NCF_MOD,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL CURVE(NCF_MOD,XV,YV)
	  END IF
C
	ELSE IF(X(1:7) .EQ. 'RD_CONT')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Continuum file')
	  CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_CONT,FILENAME,IOS)
	  IF(IOS .NE. 0)GOTO 1		!Get another option
	  SCALE_FAC=1.0D0
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  DO I=1,NCF_CONT
	    XV(I)=NU_CONT(I)
	    YV(I)=OBSF_CONT(I)*SCALE_FAC
	  END DO
	  CALL CNVRT(XV,YV,NCF_CONT,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	  IF(WR_PLT)THEN
	    DO I=1,NCF_CONT
	      WRITE(50,*)XV(I),YV(I)
	    END DO
	  ELSE
	    CALL CURVE(NCF_CONT,XV,YV)
	  END IF
!
! This option simply reads in data in XY format. No conversion is done to the data.
! Comments (bginning with !) and blank lines are ignored in the data file.
!
	ELSE IF(X(1:3) .EQ. 'RXY')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Model file')
	  CALL RD_XY_DATA_USR(NU_CONT,OBSF_CONT,NCF_CONT,NCF_MAX,FILENAME,LU_IN,IOS)
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  IF(OVER .AND. IOS .EQ. 0)THEN
	    NCF=NCF_CONT
	    NU(1:NCF)=NU_CONT(1:NCF)
	    OBSF(1:NCF)=OBSF_CONT(1:NCF)
	    WRITE(6,*)'As read in to PLT_SPEC buffer, X will assumed to be NU'
	  ELSE IF(IOS .EQ. 0)THEN
	    T1=MAXVAL(OBSF_CONT(1:NCF_CONT))
	    IF(T1 .GT. 1.0D+38)THEN
	      WRITE(6,*)'Data exceeds single precision range: Maximum=',T1 
	      WRITE(6,*)'Necessary to scale data for plotting'
	      T1=1.0D0
	      CALL USR_OPTION(T1,'SCL_FAC','1.0D+40','Factor to divide data by')
	      OBSF_CONT(1:NCF_CONT)=OBSF_CONT(1:NCF_CONT)/T1
	    END IF
	    CALL DP_CURVE(NCF_CONT,NU_CONT,OBSF_CONT)
	  END IF
!
	ELSE IF(X(1:7) .EQ. 'RROW')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Model file')
	  CALL USR_OPTION(XKEY,'XKEY',' ','Key with X data')
	  CALL USR_OPTION(YKEY,'YKEY',' ','Key with Y data')
	  CALL USR_OPTION(NCF_CONT,'N',' ','Number of data points')
	  CALL RD_ROW_DATA(NU_CONT,OBSF_CONT,NCF_CONT,XKEY,YKEY,FILENAME,LU_IN,IOS)
	  IF(IOS .EQ. 0)CALL DP_CURVE(NCF_CONT,NU_CONT,OBSF_CONT)
C
	ELSE IF(X(1:4) .EQ. 'NORM')THEN
C
	  READ_OBS=.FALSE.
	  CALL USR_HIDDEN(READ_OBS,'RD_OBS','F',' ')
C
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','Continuum file')
	  IF(READ_OBS)THEN
	    CALL USR_HIDDEN(OBS_COLS,2,2,'COLS','1,2','Columns with data')
	    CALL RD_OBS_DATA_V2(XV,YV,NCF_MAX,J,FILENAME,OBS_COLS,IOS)
	    IF(IOS .NE. 0)GOTO 1		!Get another option
	    DO I=1,J 
	      NU_CONT(I)=ANG_TO_HZ/XV(I)
	      OBSF_CONT(I)=YV(I)
	    END DO
	    NCF_CONT=J
	  ELSE
	    CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_CONT,FILENAME,IOS)
	    IF(IOS .NE. 0)GOTO 1		!Get another option
	  END IF
C
	  T1=1.0D-08
	  I=1
	  UNEQUAL=.FALSE.
	  IF(NCF_CONT .NE. NCF)UNEQUAL=.TRUE.
	  DO WHILE(.NOT. UNEQUAL .AND. I .LE. NCF_CONT)
	    IF( EQUAL(NU_CONT(I),NU(I),T1) )THEN
	      XV(I)=NU_CONT(I)
	      I=I+1
	    ELSE
	      UNEQUAL=.TRUE.
	    END IF
	  END DO
	  IF(UNEQUAL)CALL USR_HIDDEN(LIN_INT,'LIN','F',
	1               'Use linear interpolation')
	  IF(UNEQUAL .AND. LIN_INT)THEN
	    L=1
	    DO I=1,NCF
  	      XV(I)=NU(I)
	      IF(NU(I) .GT. NU_CONT(1))THEN
	        YV(I)=0.0
	      ELSE IF(NU(I) .LT. NU_CONT(NCF_CONT))THEN
	        YV(I)=0.0
	      ELSE 
	        DO WHILE (NU(I) .LT. NU_CONT(L+1))
	          L=L+1           
	        END DO
	        T1=(NU(I)-NU_CONT(L+1))/(NU_CONT(L)-NU_CONT(L+1))
	        T2=(1.0D0-T1)*OBSF_CONT(L+1)+T1*OBSF_CONT(L)
	        YV(I)=0.0
	        IF(T2 .NE. 0)THEN
	          T2=OBSF(I)/T2
	          IF(T2 .LT. 1.0E+020)YV(I)=T2
	        END IF
	      END IF
	    END DO
	  ELSE IF(UNEQUAL)THEN
C
C We will use monotonic cubic interpolation. We first verify the range.
C I & J are temporary variables for the callt o MON_INTERP. I denotes the 
C first element. Initially J denotes the last element, then the numer of
C elements that can be interpolated.
C
	    I=1
	    DO WHILE(NU(I) .GT. NU_CONT(1))
	      I=I+1
	    END DO
	    J=NCF
	    DO WHILE(NU(J) .LE. NU_CONT(NCF_CONT))
	      J=J-1
	    END DO
	    J=J-I+1
C
	    FQW(1:NCF)=0.0D0				!Temporary usage
	    CALL MON_INTERP(FQW(I),J,IONE,NU(I),J,
	1            OBSF_CONT,NCF_CONT,NU_CONT,NCF_CONT)
  	    XV(1:NCF)=NU(1:NCF)
	    DO I=1,NCF
	      IF(FQW(I) .GT. 0)THEN
	         T2=OBSF(I)/FQW(I)
	         IF(T2 .LT. 1.0E+20)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  ELSE
	    DO I=1,NCF
	      YV(I)=0
	      IF(OBSF_CONT(I) .GT. 0)THEN
	         T2=OBSF(I)/OBSF_CONT(I)
	         IF(T2 .LT. 1.0E+20)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  END IF
!
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	  CALL USR_HIDDEN(ADD_NOISE,'ADD','F','Add poisonian noise?')
	  OVER=.FALSE.
	  CALL USR_OPTION(OVER,'OVER','F','Overwrite existing model (buffer) data')
	  IF(ADD_NOISE)THEN
	    CALL USR_OPTION(COUNTS,'CNTS','100','Counts in continuum')
	    CALL USR_OPTION(R_SEED,'R_SEED','0.2',
	1              'Random number seed (0 to 1)')
	    CALL USR_OPTION(LAM_RANGE,2,2,'LAM=','3000.0,7000.0',
	1              'Wavelength Range (A)')
	    IRAN=-R_SEED*1234567
	    DO I=1,NCF
	      T1=YV(I)*COUNTS
	      IF(T1 .GE. 1 .AND. XV(I) .GT. LAM_RANGE(1) .AND. XV(I) 
	1                                       .LT. LAM_RANGE(2))
	1           YV(I)=POIDEV(T1,IRAN)/COUNTS
	    END DO
	  END IF
	  IF(OVER)THEN
	     OBSF(1:NCF)=YV(1:NCF)
	  ELSE IF(WR_PLT)THEN
	    DO I=1,NCF
	      WRITE(50,*)XV(I),YV(I)
	    END DO
	  ELSE
	    CALL CURVE(NCF,XV,YV)
	  END IF
	  YAXIS='F\d\gn\u/F\dc\u'
C
C 
C
	ELSE IF(X(1:6) .EQ. 'RD_EW')THEN
	  FILENAME=' '
	  CALL USR_OPTION(FILENAME,'File',' ','EW file')
C
	    CALL RD_EW(XV,YV,NCF_MAX,J,FILENAME,IOS)
	    IF(IOS .NE. 0)GOTO 1		!Get another option
	    CALL CURVE(J,XV,YV)
C
C 
C
	ELSE IF(X(1:6) .EQ. 'RD_OBS')THEN
	  FILENAME=' '  
	  CALL USR_OPTION(FILENAME,'File',' ',' ')
C
	  SCALE_FAC=1.0D0
	  CALL USR_HIDDEN(SCALE_FAC,'SCALE','1.0D0',' ')
	  ADD_FAC=0.0D0
	  CALL USR_HIDDEN(ADD_FAC,'ADD','0.0D0',' ')
C
	  RAD_VEL=0.0D0
	  CALL USR_HIDDEN(RAD_VEL,'RAD_VEL','0.0D0','Radial velocity of star (+ve if away)')
C
	  CLEAN=.FALSE.
	  CALL USR_HIDDEN(CLEAN,'CLEAN','F',' ')
C
	  REMOVE_BAD_PIX=.FALSE.
	  CALL USR_HIDDEN(REMOVE_BAD_PIX,'REM_BP','F','Remove bad pixels? ')
C
	  CLN_CR=.FALSE.
	  CALL USR_HIDDEN(CLN_CR,'CLN_CR','F','Remove cosmic ray spikes? ')
C
	  SMOOTH=.FALSE.
	  CALL USR_HIDDEN(SMOOTH,'SMOOTH','F',' ')
C
C Option allows the observational data to be treated as though it were
C a model. Thus can reden the data etc. No plot is done, and the model
C is overwritten.
C
	  TREAT_AS_MOD=.FALSE.
	  CALL USR_HIDDEN(TREAT_AS_MOD,'OVER','F',' ')
C
	  IF(SMOOTH)THEN
	    NHAN=5
	    CALL USR_OPTION(NHAN,'HAN','5','Number of points for HAN [ODD]')
	    NHAN=2*(NHAN/2)+1	!Ensures odd. 
	  END IF
C
	  CALL USR_HIDDEN(OBS_COLS,2,2,'COLS','1,2','Columns with data')
	  CALL RD_OBS_DATA_V2(XV,YV,NCF_MAX,J,FILENAME,OBS_COLS,IOS)
	  IF(IOS .NE. 0)GOTO 1		!Get another option
!
! Convert from Ang (Vacuum) to Hz.
!
	  DO I=1,J
	    XV(I)=ANG_TO_HZ/XV(I)
	    YV(I)=YV(I)*SCALE_FAC+ADD_FAC
	  END DO
C
	  IF(RAD_VEL .NE. 0)THEN
	    DO I=1,J
	     XV(I)=XV(I)*(1.0D0+1.0D+05*RAD_VEL/C_CMS)
	    END DO
	  END IF
C
C Procdure to remove single pizels that are zero due quirks with IUE.
C
	  IF(CLEAN)THEN
	    DO I=1,J
	      ZV(I)=YV(I)
	    END DO
	    DO I=2,J-1
	      IF(YV(I) .EQ. 0)THEN
	        YV(I)=0.5D0*(ZV(I-1)+ZV(I+1))
	      ELSE IF(YV(I) .LT. -1.0D+10)THEN
	        YV(I)=0.0D0
	      END IF
	    END DO
	  END IF
!
	  IF(REMOVE_BAD_PIX)THEN
	    DO L=3,J-50,90
	      T1=0.0D0; T2=0.0D0; T3=0.0D0
	      DO K=L,MIN(L+99,J)
	        T1=T1+YV(K)
	        T2=T2+YV(K)*YV(K)
	        T3=T3+1
	      END DO
	      T1=T1/T3
	      T2=SQRT( (T2-T3*T1*T1)/(T3-1) )
	      DO K=L+1,MIN(L+98,J-1)
	        IF( ABS(YV(K)-T1) .GT. 5.0*T2 .AND.
	1           ABS(YV(K-1)-T1) .LT. 3.0*T2 .AND.
	1           ABS(YV(K+1)-T1) .LT. 3.0*T2)THEN
	           YV(K)=0.5D0*(YV(K-1)+YV(K+1))
	        END IF
	      END DO
	    END DO
	  END IF
!
	  IF(CLN_CR)THEN
	    ZV(1:J)=YV(1:J)
	    DO L=3,J-10
	      T1=MAXVAL(YV(L:L+9))
	      DO I=L,MIN(J-2,L+9)
	         IF(YV(I) .EQ. T1)THEN
	           K=I
	           EXIT
	         END IF
	      END DO 
	      T1=0.0D0; T2=0.0D0; CNT=0
	      DO I=MAX(1,K-10),MIN(K+10,J)
	        IF(I .LT. K-1 .OR. I .GT. K+1)THEN
	          T1=T1+YV(I)
	          T2=T2+YV(I)*YV(I)
	          CNT=CNT+1
	        END IF
	      END DO
	      IF(CNT .GE. 4)THEN
	        T1=T1/CNT
	        T2=T2-T1*T1*CNT
	        IF(T2 .GT. 0.0D0)THEN
	          T2=SQRT(T2/(CNT-1.0D0))
	        ELSE
	          T1=YV(K)
	          T2=1.0D+20
	        END IF
	      ELSE
	        T1=YV(K)
	        T2=1.0D+20
	      END IF
	      IF(YV(K) .GT. T1+4.0*T2)THEN
	        YV(K)=YV(K-1)
	        IF(YV(K-1) .GT. T1+4.0*T2)YV(K)=YV(K-2)
		T3=YV(K+1) 
	        IF(T3 .GT. T1+4.0*T2)T3=YV(K+2)
	        YV(K)=0.5D0*(YV(K)+T3)
	      END IF 
	    END DO
	  END IF
C
C We check whether the X axis is monotonic. If not, we smooth each section
C separately. Designed for non-merged overlapping ECHELLE orders.
C
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
C
	  K=NHAN
	  IF(SMOOTH .AND. NON_MONOTONIC)THEN
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
C
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

C
	  IF(TREAT_AS_MOD)THEN
	    NCF=J
	    DO I=1,NCF
	      NU(I)=XV(I)
	      OBSF(I)=YV(I)
	    END DO
	    WRITE(T_OUT,*)'Observational data replaces model data'
	  ELSE
	    NOBS=J
	    NU_OBS(1:NOBS)=XV(1:NOBS)
	    OBSF_OBS(1:NOBS)=YV(1:NOBS)
	    CALL CNVRT(XV,YV,J,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    CALL USR_HIDDEN(WR_PLT,'WR','F','Write data to file')
	    IF(WR_PLT)THEN
	      DO I=1,J
	        WRITE(50,*)XV(I),YV(I)
	      END DO
	    ELSE
	      CALL CURVE(J,XV,YV)
	    END IF
	  END IF
!
! 
!
	ELSE IF(X(1:5) .EQ. 'ISABS') THEN
	  CALL USR_OPTION(T_IN_K,'T_IN_K','100d0','Temp. in Kelvin')
	  CALL USR_OPTION(V_TURB,'V_TURB','10d0','Turbulent Velocity (km/s)')
	  CALL USR_OPTION(LOG_NTOT,'LOG_NTOT','20d0','Log of H column density')
	  CALL USR_OPTION(LOG_H2_NTOT,'LOG_H2_NTOT','20d0','Log H2 column dens')
!
	  CALL USR_HIDDEN(HI_ABS,'HI_ABS','T','Correct for HI absorption')
	  CALL USR_HIDDEN(H2_ABS,'H2_ABS','T','Correct for HII absorption')
	  CALL USR_HIDDEN(V_R,'V_R','0.0D0','Radial Velocity (km/s)')
	  CALL USR_HIDDEN(WAVE_MIN,'WAVE_MIN','900d0','Minimum Wavelength')
	  CALL USR_HIDDEN(WAVE_MAX,'WAVE_MAX','3000d0','Maximum Wavelength')
	  CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','2.0D0','Minimum Model Resolution')
	  CALL USR_HIDDEN(IS_FILE,'IS_FILE','IS_LINE_LIST','File wth list and strengths of IS lines')
!
	  CALL UVABS_V2(NU,OBSF,NCF,NCF_MAX,
     1              T_IN_K,V_TURB,LOG_NTOT,
     1	            LOG_H2_NTOT,V_R,MIN_RES_KMS,WAVE_MAX,WAVE_MIN,
     1	            HI_ABS,H2_ABS,IS_FILE)
!
!
!
	ELSE IF(X(1:5) .EQ. 'CNVLV') THEN
!
	 IF(NCF .EQ. 0)THEN
	    WRITE(T_OUT,*)'Error: no data in buffer to operate on.'
	    WRITE(T_OUT,*)'Use rd_mod(OVER=T) with rd_mod option, or'
	    WRITE(T_OUT,*)'Use rd_obs(MOD=T) wit rd_obs option.'
	    GOTO 1		!Get another option
	 END IF
!
! Instrumental profiles is assumed to be Gaussian.
!
	 WRITE(T_OUT,*)' '
	 WRITE(T_OUT,*)' Two choices are possible: '
	 WRITE(T_OUT,*)'   (1) Fixed resolution (INST_RES) in angstroms'
	 WRITE(T_OUT,*)'   (2) Fixed velocity resoluton (Lam/dLam) '
	 WRITE(T_OUT,*)' The non-zero value is used'
	 WRITE(T_OUT,*)' '
!
! If the RESOLUTION is non zero, it gets used instead of INST_RES.
! Resolution allows a spectrum to be smoothed to a constant velocity
! resolution at all wavelengths.
!
100	 CALL USR_OPTION(INST_RES,'INST_RES','0.0D0','Instrumental Resolution in Angstroms [dLam - FWHM]')
	 CALL USR_OPTION(RESOLUTION,'RES','0.0D0','Resolution [Lam/dLam(FWHM)] (km/s if -ve)')
	 IF(RESOLUTION .LT. 0.0D0)THEN
	   RESOLUTION=1.0D-05*C_CMS/ABS(RESOLUTION)
	 ELSE IF(RESOLUTION .EQ. 0.0D0 .AND. INST_RES .EQ. 0.0D0)THEN
	   WRITE(T_OUT,*)'Only one INST_RES and RES can be zero'
	   GOTO 100
	 ELSE IF(RESOLUTION .NE. 0.0D0 .AND. INST_RES .NE. 0.0D0)THEN
	   WRITE(T_OUT,*)'Only one INST_RES and RES can be non-zero'
	   GOTO 100
	 END IF
! 
! Defaults are those of HUT.
!
	 CALL USR_OPTION(WAVE_MIN,'WAVE_MIN','900d0','Minimum Wavelength')
	 CALL USR_OPTION(WAVE_MAX,'WAVE_MAX','10000d0','Maximum Wavelength')
	 CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','1.0d0',
	1           'Minimum Model Resolution (km/s)')
	 CALL USR_HIDDEN(NUM_RES,'NUM_RES','5.0d0',
     1              'Number of resolution elements condidered outside bandpass')
	 CALL USR_HIDDEN(FFT_CONVOLVE,'FFT','F',
     1              'Use FFT methods for convolution')
!
	 VSINI=0.0D0		!For rotational broadening, so set to zero
	 EPSILON=0.0D0
	 CALL SMEAR_V2(NU,OBSF,NCF,
	1	      WAVE_MAX,WAVE_MIN,
	1             INST_RES,RESOLUTION,VSINI,EPSILON, 
	1             MIN_RES_KMS,NUM_RES,FFT_CONVOLVE)
!
	ELSE IF(X(1:3) .EQ. 'ROT') THEN
!
	  IF(NCF .EQ. 0)THEN
	    WRITE(T_OUT,*)'Error: no data in buffer to operate on.'
	    WRITE(T_OUT,*)'Use rd_mod(OVER=T) with rd_mod option, or'
	    WRITE(T_OUT,*)'Use rd_obs(MOD=T) wit rd_obs option.'
	    GOTO 1		!Get another option
	  END IF
!
! Perform a crude rotational broadening. This is not meant to be rigorous.
!
	 CALL USR_OPTION(VSINI,'VSINI','100.0D0','Vsini')
	 CALL USR_OPTION(EPSILON,'EPS','0.5D0',
	1            'I(mu)/I(mu=1) = 1-eps + eps*mu')
!
! Defaults are observable spectral region.
!
	 CALL USR_OPTION(WAVE_MIN,'WAVE_MIN','900.0d0','Minimum Wavelength')
	 CALL USR_OPTION(WAVE_MAX,'WAVE_MAX','7000.0d0','Maximum Wavelength')
!
! If the RESOLUTION is non zero, it gets used instead of INST_RES.
! Resolution allows a spectrum to be smoothed to a constant velocity
! resolution at all wavelengths.
!
	 CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','1.0d0',
	1                      'Minimum Model Resolution (km/s)')
	 CALL USR_HIDDEN(NUM_RES,'NUM_RES','5.0d0',
     1                   'Number of resolution elements condidered')
!
! Rotational broadening cannot use FFT option at present.
!
	 FFT_CONVOLVE=.FALSE.
	 CALL SMEAR_V2(NU,OBSF,NCF,
	1	         WAVE_MAX,WAVE_MIN,
	1                INST_RES,RESOLUTION,VSINI,EPSILON, 
	1                MIN_RES_KMS,NUM_RES,FFT_CONVOLVE)
!
	ELSE IF(X(1:3) .EQ. 'EXT') THEN
!
! Extract spectrum at a given pixel resolution.
!
	  CALL USR_OPTION(RESOLUTION,'RES','1000.0D0','R / /\R ')
	  I=1
	  T1=NU(1)
	  ML=1
	  DO WHILE(T1 .GT. NU(NCF))
	    DO WHILE( T1 .LT. NU(ML+1) )
	      ML=ML+1
	    END DO
	    T2=(T1-NU(ML))/(NU(ML+1)-NU(ML))
	    YV(I)=(1.0D0-T2)*OBSF(ML)+T2*OBSF(ML+1)
	    XV(I)=T1
	    T1=T1*RESOLUTION/(1.0D0+RESOLUTION)
	    I=I+1
	  END DO
	  NCF=I-1
	  NU(1:NCF)=XV(1:NCF)
	  OBSF(1:NCF)=YV(1:NCF)
C 
C
C This option s a simpliied combination of RD_MOD, ROT, and NORM. It is specifically
C design for examining normalized model spectra. The separate options are more general.
C
	ELSE IF(X(1:3) .EQ. 'GEN')THEN
	  DIRECTORY=' '
	  CALL USR_OPTION(DIRECTORY,'Dir',' ','Model directory')
	  CALL USR_HIDDEN(FILENAME,'File','obs_fin','File (def=obs_fin)')
C
	  IF(INDEX(DIRECTORY,']') .EQ. 0)THEN
	    FILENAME=TRIM(DIRECTORY)//'/obs/'//TRIM(FILENAME)
	  ELSE
	    I=LEN_TRIM(DIRECTORY)
	    CALL SET_CASE_UP(FILENAME,IZERO,IZERO)
	    FILENAME=DIRECTORY(1:I-1)//'.OBS]'//TRIM(FILENAME)
	  END IF
	  CALL RD_MOD(NU,OBSF,NCF_MAX,NCF,FILENAME,IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error reading model data -- no data read'
	     GOTO 1		!Get another option
	  END IF
	  WRITE(T_OUT,*)'New model data replaces old data'
!
! Perform a crude rotational broadening. This is not meant to be rigorous.
!
	  CALL USR_OPTION(VSINI,'VSINI','100.0D0','Vsini')
	  CALL USR_OPTION(EPSILON,'EPS','0.5D0',
	1            'I(mu)/I(mu=1) = 1-eps + eps*mu')
!
! Defaults are observable spectral region.
!
	  CALL USR_OPTION(WAVE_MIN,'WAVE_MIN','900.0d0','Minimum Wavelength')
	  CALL USR_OPTION(WAVE_MAX,'WAVE_MAX','7000.0d0','Maximum Wavelength')
!
! If the RESOLUTION is non zero, it gets used instead of INST_RES.
! Resolution allows a spectrum to be smoothed to a constant velocity
! resolution at all wavelengths.
!
	  CALL USR_HIDDEN(MIN_RES_KMS,'MIN_RES','1.0d0',
	1                     'Minimum Model Resolution (km/s)')
	  CALL USR_HIDDEN(NUM_RES,'NUM_RES','5.0d0',
     1                   'Number of resolution elements condidered')
!
! Rotational broadening cannot use FFT option at present.
!
	  FFT_CONVOLVE=.FALSE.
	  CALL SMEAR_V2(NU,OBSF,NCF,
	1	         WAVE_MAX,WAVE_MIN,
	1                INST_RES,RESOLUTION,VSINI,EPSILON, 
	1                MIN_RES_KMS,NUM_RES,FFT_CONVOLVE)
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Rotational smoothing has been performed'
!
	  CALL USR_OPTION(T_IN_K,'T_IN_K','100d0','Temp. in Kelvin')
	  CALL USR_OPTION(V_TURB,'V_TURB','10d0','Turbulent Velocity (km/s)')
	  CALL USR_OPTION(LOG_NTOT,'LOG_NTOT','20d0','Log of H column density')
	  CALL USR_OPTION(LOG_H2_NTOT,'LOG_H2_NTOT','20d0','Log H2 column dens')
!
	  CALL USR_HIDDEN(HI_ABS,'HI_ABS','T','Correct for HI absorption')
	  CALL USR_HIDDEN(H2_ABS,'H2_ABS','T','Correct for HII absorption')
	  CALL USR_HIDDEN(V_R,'V_R','0.0d0','Radial Velocity (km/s)')
	  CALL USR_HIDDEN(WAVE_MIN,'WAVE_MINI','900d0','Minimum Wavelength')
	  CALL USR_HIDDEN(WAVE_MAX,'WAVE_MAXI','1300d0','Maximum Wavelength')
	  CALL USR_HIDDEN(IS_FILE,'IS_FILE','IS_LINE_LIST','File wth list and strengths of IS lines')
!
	  CALL UVABS_V2(NU,OBSF,NCF,NCF_MAX,
     1              T_IN_K,V_TURB,LOG_NTOT,
     1	            LOG_H2_NTOT,V_R,MIN_RES_KMS,WAVE_MAX,WAVE_MIN,
     1	            HI_ABS,H2_ABS)
!
	  WRITE(T_OUT,*)'ISABS correction has been performed'
!

	  IF(INDEX(DIRECTORY,']') .EQ. 0)THEN
	    FILENAME=TRIM(DIRECTORY)//'/obs/obs_cont'
	  ELSE
	    I=LEN_TRIM(DIRECTORY)
	    FILENAME=DIRECTORY(1:I-1)//'.OBS]OBS_cont'
	  END IF
	  CALL RD_MOD(NU_CONT,OBSF_CONT,NCF_MAX,NCF_CONT,FILENAME,IOS)
	  IF(IOS .NE. 0)GOTO 1		!Get another option
!
	  T1=1.0D-08
	  I=1
	  UNEQUAL=.FALSE.
	  IF(NCF_CONT .NE. NCF)UNEQUAL=.TRUE.
	  DO WHILE(.NOT. UNEQUAL .AND. I .LE. NCF_CONT)
	    IF( EQUAL(NU_CONT(I),NU(I),T1) )THEN
	      XV(I)=NU_CONT(I)
	      I=I+1
	    ELSE
	      UNEQUAL=.TRUE.
	    END IF
	  END DO
	  IF(UNEQUAL)CALL USR_HIDDEN(LIN_INT,'LIN','F',
	1               'Use linear interpolation')
	  IF(UNEQUAL .AND. LIN_INT)THEN
	    L=1
	    DO I=1,NCF
  	      XV(I)=NU(I)
	      IF(NU(I) .GT. NU_CONT(1))THEN
	        YV(I)=0.0
	      ELSE IF(NU(I) .LT. NU_CONT(NCF_CONT))THEN
	        YV(I)=0.0
	      ELSE 
	        DO WHILE (NU(I) .LT. NU_CONT(L+1))
	          L=L+1           
	        END DO
	        T1=(NU(I)-NU_CONT(L+1))/(NU_CONT(L)-NU_CONT(L+1))
	        T2=(1.0D0-T1)*OBSF_CONT(L+1)+T1*OBSF_CONT(L)
	        YV(I)=0.0
	        IF(T2 .NE. 0)THEN
	          T2=OBSF(I)/T2
	          IF(T2 .LT. 1.0E+020)YV(I)=T2
	        END IF
	      END IF
	    END DO
	  ELSE IF(UNEQUAL)THEN
C
C We will use monotonic cubic interpolation. We first verify the range.
C I & J are temporary variables for the callt o MON_INTERP. I denotes the 
C first element. Initially J denotes the last element, then the numer of
C elements that can be interpolated.
C
	    I=1
	    DO WHILE(NU(I) .GT. NU_CONT(1))
	      I=I+1
	    END DO
	    J=NCF
	    DO WHILE(NU(J) .LE. NU_CONT(NCF_CONT))
	      J=J-1
	    END DO
	    J=J-I+1
C
	    FQW(1:NCF)=0.0D0				!Temporary usage
	    CALL MON_INTERP(FQW(I),J,IONE,NU(I),J,
	1            OBSF_CONT,NCF_CONT,NU_CONT,NCF_CONT)
  	    XV(1:NCF)=NU(1:NCF)
	    DO I=1,NCF
	      IF(FQW(I) .GT. 0)THEN
	         T2=OBSF(I)/FQW(I)
	         IF(T2 .LT. 1.0E+20)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  ELSE
	    DO I=1,NCF
	      YV(I)=0
	      IF(OBSF_CONT(I) .GT. 0)THEN
	         T2=OBSF(I)/OBSF_CONT(I)
	         IF(T2 .LT. 1.0E+20)YV(I)=T2
	      ELSE
	        YV(I)=0
	      END IF
	    END DO
	  END IF
!
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL CURVE(NCF,XV,YV)
	  YAXIS='F\d\gn\u/F\dc\u'
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Data has been normalized and plotted'
!
!
!
! Option to do a least fit CCM redenning Law. The observational data
! should have been read in using RD_OBS, and must be contained in one file.
!
	ELSE IF(X(1:3) .EQ. 'RED')THEN
	  IF(NOBS .EQ. 0)THEN
	    WRITE(6,*)'Error -- observational data has not been read in'
	  ELSE IF(NCF .EQ. 0)THEN
	    WRITE(6,*)'Error -- model data has not been strored in the buffer'
	  ELSE
	    T1=2.5; T2=5.5
	    CALL DETERM_REDDENING(OBSF_OBS,NU_OBS,NOBS,OBSF,NU,NCF,T1,T2)
	  END IF
! 
	ELSE IF(X(1:4) .EQ. 'FLAM' .OR. 
	1     X(1:4) .EQ. 'WRFL' .OR.  X(1:3) .EQ. 'FNU' .OR.
	1                X(1:4) .EQ. 'EBMV') THEN
!
! If NCF is defined, we use that freuqency grid when computing the extinction curve.
! Otherwise, we define the grid.
!
	  IF(X(1:4) .EQ. 'EBMV' .AND. NCF .EQ. 0)THEN
	    T1=ANG_TO_HZ/900.0D0; T2=ANG_TO_HZ/5.0E+04
	    NCF=1000
	    T2=EXP(LOG(T1/T2)/(NCF-1))
	    NU(1)=T1
	    DO I=2,NCF
	      NU(I)=NU(I-1)/T2
	    END DO
	  END IF
!
	  CALL USR_OPTION(EBMV_CCM,'EBMV_CCM','0.0D0',
	1             'CCM E(B-V) to correct for I.S. extinction')
	  IF(EBMV_CCM .NE. 0)THEN
	    CALL USR_OPTION(r_ext,'R_EXT','3.1D0',
	1         'R_EXT for Cardelli, Clayton, Mathis extinction law')
	  END IF
!                                                                          
	  CALL USR_OPTION(EBMV_GAL,'EBMV_GAL','0.0',
	1             'Galactic E(B-V) to correct for I.S. extinction')
	  CALL USR_OPTION(EBMV_LMC,'EBMV_LMC','0.0',
	1             'LMC E(B-V) to correct for I.S. extinction')
	  CALL USR_OPTION(EBMV_SMC,'EBMV_SMC','0.0',
	1      'SMC E(B-V) to correct for I.S. extinction') !KN SMC
!
	  DIST=1.0D0
	  CALL USR_OPTION(DIST,'DIST','1.0D0',' (in kpc) ')
!
	  OVER=.FALSE.
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite store')
!
	  DO I=1,NCF
	    XV(I)=NU(I)
	    YV(I)=OBSF(I)/DIST/DIST
	  END DO
C
	  IF(EBMV_CCM .NE. 0.0D0)THEN
	    DO I=1,NCF
	      T1=ANG_TO_HZ/NU(I)
	      T1=(10000.0/T1)				!1/Lambda(um)
	      IF(T1 .LT. 1.1)THEN
	        RAX=0.574*(T1**1.61)
	        RBX=-0.527*(T1**1.61)
	      ELSE IF(T1. LT. 3.3)THEN
	        T2=T1-1.82
	        RAX=1+T2*(0.17699-T2*(0.50447+T2*(0.02427-T2*(0.72085
	1                +T2*(0.01979-T2*(0.77530-0.32999*T2))))))
	        RBX=T2*(1.41338+T2*(2.28305+T2*(1.07233-T2*(5.38434
	1                +T2*(0.62251-T2*(5.30260-2.09002*T2))))))
	      ELSE IF(T1 .lT. 5.9)THEN
	        RAX=1.752-0.316*T1-0.104/((T1-4.67)**2+0.341)
	        RBX=-3.090+1.825*T1+1.206/((T1-4.62)**2+0.263)
	      ELSE IF(T1 .LT. 8.0)THEN
  	        T2=T1-5.9
	        RAX=1.752-0.316*T1-0.104/((T1-4.67)**2+0.341) -
	1                       T2*T2*(0.04773+0.009779*T2)
	        RBX=-3.090+1.825*T1+1.206/((T1-4.62)**2+0.263)+
	1                       T2*T2*(0.2130+0.1207*T2)
	      ELSE IF(T1 .LT. 10)THEN
	        T2=T1-8
	        RAX=-1.073-T2*(0.628-T2*(0.137-0.070*T2))
	        RBX=13.670+T2*(4.257-T2*(0.420-0.374*T2))               
	      ELSE 
	        T1=10
	        T2=T1-8
	        RAX=-1.073-T2*(0.628-T2*(0.137-0.070*T2))
	        RBX=13.670+T2*(4.257-T2*(0.420-0.374*T2))
	      END IF
              AL_D_EBmV(I)=R_EXT*(RAX+RBX/R_EXT)
	    END DO
	    IF(X(1:4) .EQ. 'EBMV')THEN
	      DO I=1,NCF
	        XV(I)=NU(I)
	        YV(I)=EBMV_CCM*AL_D_EBmV(I)
	      END DO
 	      XAXIS='\gl(\A)'
	      YAXIS='A\d\gl\u'
	    ELSE
	      DO I=1,NCF
	        YV(I)=YV(I)*( 10.0**(-0.4D0*EBMV_CCM*AL_D_EBmV(I)) )
	      END DO
	    END IF
	  END IF
C
C Set galactic interstellar extinction curve. Curve is from Howarth
C (1983, MNRAS, 203, 301) and Seaton (1979, MNRAS, 187, 73P).
C
	  IF(EBMV_GAL .NE. 0)THEN
	    R_EXT=3.1
	    DO I=1,NCF
	      T1=ANG_TO_HZ/NU(I)
	      T1=(10000.0/T1)				!1/Lambda(um)
	      IF( T1 .LT. 1.83)THEN
	        AL_D_EBmV(I)=((1.86-0.48*T1)*T1-0.1)*T1
	      ELSE IF( T1 .LT. 2.75)THEN
	        AL_D_EBmV(I)=R_EXT+2.56*(T1-1.83)-0.993*(T1-1.83)**2
	      ELSE IF( T1 .LT. 3.65)THEN
	        AL_D_EBmV(I)=(R_EXT-1.64)+1.048*T1+1.01/( (T1-4.60)**2+0.28 )
	      ELSE IF( T1 .LT. 7.14)THEN
	        AL_D_EBmV(I)=(R_EXT-0.91)+0.848*T1+1.01/( (T1-4.60)**2+0.28 )
	      ELSE IF( T1 .LT. 11)THEN                                      
                AL_D_EBmV(I)=(R_EXT+12.97)-3.20*T1+0.2975*T1*T1
	      ELSE
	        T1=11.0
                AL_D_EBmV(I)=(R_EXT+12.97)-3.20*T1+0.2975*T1*T1
	      END IF
	    END DO
	  END IF
	  IF(X(1:4) .EQ. 'EBMV' .AND. EBMV_GAL .NE. 0)THEN
	    DO I=1,NCF
	      XV(I)=NU(I)
	      YV(I)=EBMV_GAL*AL_D_EBmV(I)
	    END DO
	    XAXIS='\gl(\A)'
	    YAXIS='A\d\gl\u'
	  ELSE IF(X(1:4) .NE. 'EBMV' .AND. EBMV_GAL .NE. 0)THEN
	    DO I=1,NCF
	      YV(I)=YV(I)*( 10.0**(-0.4D0*EBMV_GAL*AL_D_EBmV(I)) )
	    END DO
	  END IF
C
C Set LMC interstellar extinction curve. Curve is from Howarth
C (1983, MNRAS, 203, 301).
C
	  IF(EBMV_LMC .NE. 0)THEN
	    R_EXT=3.1
	    DO I=1,NCF
	      T1=ANG_TO_HZ/NU(I)
	      T1=(10000.0/T1)				!1/Lambda(um)
	      IF( T1 .LT. 1.83)THEN
    	        AL_D_EBmV(I)=((1.86-0.48*T1)*T1-0.1)*T1
	      ELSE IF( T1 .LT. 2.75)THEN
	        AL_D_EBmV(I)=R_EXT+2.04*(T1-1.83)+0.094*(T1-1.83)**2
	      ELSE IF( T1 .LT. 11.0)THEN	!Strictly only below 9
	        AL_D_EBmV(I)=(R_EXT-0.236)+0.462*T1+0.105*T1*T1+
	1                       0.454/( (T1-4.557)**2+0.293 )
	      ELSE
	        T1=11.0
	        AL_D_EBmV(I)=(R_EXT-0.236)+0.462*T1+0.105*T1*T1+
	1                       0.454/( (T1-4.557)**2+0.293 )
	      END IF
	    END DO
	  END IF
C
	  IF(X(1:4) .EQ. 'EBMV' .AND. EBMV_LMC .NE. 0)THEN
	    DO I=1,NCF
	      XV(I)=NU(I)
!	      XV(I)=ANG_TO_HZ/NU(I)
	      YV(I)=EBMV_LMC*AL_D_EBmV(I)
	    END DO
	    XAXIS='\gl(\A)'
	    YAXIS='A\d\gl\u'
	  ELSE IF(X(1:4) .NE. 'EBMV' .AND. EBMV_LMC .NE. 0)THEN
	    DO I=1,NCF
	      YV(I)=YV(I)*( 10.0**(-0.4D0*EBMV_LMC*AL_D_EBmV(I)) )
	    END DO
	  END IF
C
C       KN: Set SMC interstellar extinction curve.
C
	  IF(EBMV_SMC .NE. 0)THEN
	     R_EXT=2.74
	     DO I=1,NCF
	        T1=ANG_TO_HZ/NU(I)
	        T1=(10000.0/T1) !1/Lambda(um)
	        IF(T1 .LT. 1.1)THEN
	          RAX=0.574*(T1**1.61)
	          RBX=-0.527*(T1**1.61)
	          AL_D_EBmV(I)=R_EXT*(RAX+RBX/R_EXT)
	        ELSE IF(T1. LT. 3.39)THEN
	          T2=T1-1.82
	          RAX=1+T2*(0.17699-T2*(0.50447+T2*(0.02427-T2*(0.72085
	1                +T2*(0.01979-T2*(0.77530-0.32999*T2))))))
	          RBX=T2*(1.41338+T2*(2.28305+T2*(1.07233-T2*(5.38434
	1                 +T2*(0.62251-T2*(5.30260-2.09002*T2))))))
	          AL_D_EBmV(I)=R_EXT*(RAX+RBX/R_EXT)
	        ELSE 
	          C1=-4.959
	          C2=2.264*T1
	          D=(T1**2)/(((T1**2-4.6**2)**2)+(T1**2))
	          C3=0.389*D
	          IF(T1 .LT. 5.9)THEN
	             F=0
	          ELSE
	             F=0.5392*((T1-5.9)**2)+0.05644*((T1-5.9)**3)
	          END IF
	          C4=0.461*F
                  AL_D_EBmV(I)=C1+C2+C3+C4+R_EXT
               END IF
	     END DO
          END IF
	  IF(X(1:4) .EQ. 'EBMV' .AND. EBMV_SMC .NE. 0)THEN
	     DO I=1,NCF
	        XV(I)=NU(I)
	        YV(I)=EBMV_SMC*AL_D_EBmV(I)
	     END DO
	     XAXIS='\gl(\A)'
	     YAXIS='A\d\gl\u'
	  ELSE IF(X(1:4) .NE. 'EBMV' .AND. EBMV_SMC .NE. 0)THEN
	     DO I=1,NCF
	        YV(I)=YV(I)*( 10.0**(-0.4D0*EBMV_SMC*AL_D_EBmV(I)) )
	     END DO
	  END IF
C
	  IF(X(1:4) .NE. 'EBMV')THEN
	    CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	  ELSE
	    CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  END IF
 	  IF(X(1:4) .Eq. 'WRFL')THEN
	    DO I=1,NCF
	      WRITE(50,*)XV(I),YV(I)
	    END DO
	  ELSE IF(OVER)THEN
	    OBSF(1:NCF)=YV(1:NCF)
	  ELSE 
	    CALL CURVE(NCF,XV,YV)
	  END IF
	ELSE IF(X(1:2) .EQ.'RV')THEN
	  RAD_VEL=0.0D0
	  CALL USR_OPTION(RAD_VEL,'RAD_VEL','0.0D0','Radial velcoity (+ve if away)(km/s)')
	  RAD_VEL=1.0D+05*RAD_VEL
	  DO I=1,NCF
	    NU(I)=NU(I)*(1.0D0-RAD_VEL/C_CMS)
	  END DO 
C
C The follwing option cumputes the luminosity below the 3 main H/He edges,
C and the number of photons emitted.
C
	ELSE IF(X(1:3) .EQ.'ZAN')THEN
	  CALL TRAPUNEQ(NU,FQW,NCF)
	  CALL USR_OPTION(ZAN_FREQ,5,1,'LEVS','0,0,0,0,0,0',
	1      'Edge frequencies to compute photon flux (not H/He)')
C
	  ZAN_FREQ(4:8)=ZAN_FREQ(1:5)
	  ZAN_FREQ(1)=EDGE_HYD
	  ZAN_FREQ(2)=EDGE_HEI
	  ZAN_FREQ(3)=EDGE_HE2
C
	  TOT_LUM=0.0D0
	  DO I=1,NCF
	    TOT_LUM=TOT_LUM+OBSF(I)*FQW(I)
	  END DO
	  TOT_LUM=TOT_LUM*312.7			!4pi*(1kpc)**2*(1E+15)*(1D-23)/Lsun
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Total luminosity is:               ',TOT_LUM
C
	  ML=1
	  DO WHILE(ZAN_FREQ(ML) .GT. 0.0D0)
	    TMP_FREQ=ZAN_FREQ(ML)
	    N_PHOT=0.0D0
	    LUM=0.0D0
	    J=1
	    DO WHILE(NU(J+1) .GT. TMP_FREQ)
	      J=J+1
	    END DO
	    CALL TRAPUNEQ(NU,FQW,J)
	    DO I=1,J
              N_PHOT=N_PHOT+FQW(I)*OBSF(I)/NU(I)
	      LUM=LUM+FQW(I)*OBSF(I)
	    END DO
	    LUM=LUM*312.7			!4pi*(1kpc)**2*(1E+15)*(1D-23)/Lsun
	    N_PHOT=47.2566+DLOG10(N_PHOT) !DLOG10(4pi*(1kpc)**2*(1D-23)/h)
	    T1=ANG_TO_HZ/TMP_FREQ
	    WRITE(T_OUT,*)' '
	    WRITE(T_OUT,'(1X,A,F4.0,A,1PE11.4)')
	1        'Luminosity shortward of        ',T1,'A is:   ',LUM
	    WRITE(T_OUT,'(1X,A,F4.0,A,1PE11.4)')
	1        'Log(#) of photons shortward of ',T1,'A is:   ',N_PHOT
	    ML=ML+1
	  END DO
C
C
C The follwing option plot the cummulative phopton flux as a function of
C frequency.
C
	ELSE IF(X(1:4) .EQ. 'PZAN')THEN
C
	  XV(1:NCF)=NU(1:NCF)
	  YV(1)=0.0D0
	  DO ML=2,NCF
	    YV(ML)=YV(ML-1)+ 0.5D0*(NU(ML-1)-NU(ML))*
	1              (OBSF(ML-1)/NU(ML-1)+OBSF(ML)/NU(ML))
	  END DO
	  DO ML=2,NCF
	    IF(YV(ML) .LT. 1.0D-30)YV(ML)=1.0D-30
	    YV(ML)=47.2566D0+LOG10(YV(ML)) 		!DLOG10(4pi*(1kpc)**2*(1E-23)/h)
	  END DO
C
	  I=NCF-1
	  CALL CNVRT(XV(2),YV(2),I,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL CURVE(I,XV(2),YV(2))
	  YAXIS='Log(N)'
C
	ELSE IF(X(1:3) .EQ.'CUM')THEN
	  NORM_LUM=0.0D0
	  CALL USR_OPTION(NORM_LUM,'NORM_LUM',' ',' (in Lsun) ')
	  YV(1)=0.0D0
	  DO I=2,NCF
	    YV(I)=YV(I-1)+(OBSF(I)+OBSF(I-1))*(NU(I-1)-NU(I))*0.5D0
	  END DO
	  T2=YV(NCF)*312.7D0		!4pi*(1kpc)**2*(1D+15)*(1D-23)/Lsun
	  IF(NORM_LUM .EQ. 0.0D0)THEN
	    NORM_LUM=YV(NCF)
	  ELSE
	    NORM_LUM=NORM_LUM/312.7D0
	  END IF
	  DO I=1,NCF
	    XV(I)=NU(I)
	    YV(I)=YV(I)/NORM_LUM
	  END DO
	  WRITE(T_OUT,*)'Total Luminosoty is',T2
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL CURVE(NCF,XV,YV)
C
	ELSE IF(X(1:4) .EQ. 'FILT')THEN
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=OBSF(1:NCF)
	  CALL CNVRT(XV,YV,NCF,L_FALSE,L_FALSE,'ANG','FLAM',
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	  DO IF=1,NF
!
	    SP_VAL=FILT_ST_LAM(IF); ML_ST=GET_INDX_SP(SP_VAL,XV,NCF)-1
	    SP_VAL=FILT_ST_LAM(IF)+24.0D0*FILT_DEL_LAM(IF)
	    ML_END=GET_INDX_SP(SP_VAL,XV,NCF)+1
!
	    WRITE(T_OUT,*)'Filter ',IF
	    WRITE(T_OUT,*)FILT_ST_LAM(IF),XV(ML_ST),XV(ML_END),YV(ML_ST),YV(ML_END)
	    WRITE(T_OUT,*)ML_ST,ML_END
!
	    ZV(IF)=0.0D0
	    SUM=0.0D0
	    RESPONSE2=0.0D0
	    I=1
	    DO ML=ML_ST,ML_END
	      IF(XV(ML) .GE. FILT_ST_LAM(IF)+24.0D0*FILT_DEL_LAM(IF))EXIT
	      DO WHILE(XV(ML) .GT. FILT_ST_LAM(IF)+I*FILT_DEL_LAM(IF))
	        I=I+1
	      END DO
	      FILT_INT_BEG =FILT_ST_LAM(IF)+(I-1)*FILT_DEL_LAM(IF)
	      RESPONSE1=RESPONSE2
	      T1=(XV(ML+1)-FILT_INT_BEG)/FILT_DEL_LAM(IF)
	      RESPONSE2=(1.0D0-T1)*NORM_PASS(I,IF)+T1*NORM_PASS(I+1,IF) 
	      IF(ABS(T1) .GT. 1)RESPONSE2=0.0D0
	      ZV(IF)=ZV(IF)+0.5D0*(RESPONSE1*YV(ML)+RESPONSE2*YV(ML+1))*(XV(ML+1)-XV(ML))
	      SUM=SUM+0.5D0*(RESPONSE1+RESPONSE2)*(XV(ML+1)-XV(ML))
	    END DO
	    WRITE(T_OUT,*)ZV(IF),SUM,ZV(IF)/SUM
	    ZV(IF)=-FILT_ZP(IF)-2.5D0*LOG10(ZV(IF)/SUM)
	    WRITE(T_OUT,'(A,1PE10.3)')TRIM(FILT_NAME(IF)),ZV(IF)
	  END DO
	  DO IF=1,NF
            WRITE(T_OUT,'(A,F12.3)')TRIM(FILT_NAME(IF)),ZV(IF)
	  END DO

	ELSE IF(X(1:3) .EQ. 'MAG')THEN
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'The magnitudes computed here are NOT accurate.'
	  WRITE(T_OUT,*)'Filter convolutions need to be implemented.'
	  WRITE(T_OUT,*)'Originally designed for crude estimates using continuum fluxe.s'
	  WRITE(T_OUT,*)' '
	  DIST=1.0
	  CALL USR_OPTION(DIST,'DIST','1.0D0',' (in kpc) ')
	  CALL USR_OPTION(FILTER_SET,'FSET','uvby','Filter set (case sensitive)')
	  CALL GEN_ASCI_OPEN(LU_OUT,'MAG','UNKNOWN',' ',' ',IZERO,IOS)
!
	  IF( UC(FILTER_SET) .EQ. 'ALL')THEN
	    FILTER_SET='ubvy'
	    CALL GET_MAG(NU,OBSF,NCF,DIST,FILTER_SET,LU_OUT)
	    FILTER_SET='UBV'
	    CALL GET_MAG(NU,OBSF,NCF,DIST,FILTER_SET,LU_OUT)
	  ELSE
	    CALL GET_MAG(NU,OBSF,NCF,DIST,FILTER_SET,LU_OUT)
	  END IF 
!
!	  WRITE(LU_OUT,104)DIST	  
!	  DO L=1,8
!	    MLST=1
!	    DO ML=MLST,NCF
!	      IF(NU(ML) .LT. 0.2998/FILTLAM(L))THEN
!	        T1=DLOG10(OBSF(ML-1))-DLOG10(NU(ML-1)*FILTLAM(L)/0.2998)
!	1       *(DLOG10(OBSF(ML)/OBSF(ML-1))/(DLOG10(NU(ML)/NU(ML-1))))
!	        T1=5.0*DLOG10(DIST)-2.5*T1+2.5*DLOG10(FILTZP(L))
!	        MLST=ML         
!	        GOTO 102
!	       END IF
!	    END DO
!102	  CONTINUE         
!	    WRITE(LU_OUT,103)FILTLAM(L),FILTZP(L),FILT(L),T1
!103	    FORMAT(2X,F7.4,5X,F7.1,8X,A1,5X,F6.2)
!104	    FORMAT(2X,' Assumed Distance is',F5.1,' kpc',/
!	1  ,/,2x,'    Lam          ZP     Filt        Mag'//)
!	  END DO
!	  WRITE(LU_OUT,107)
!	  DO L=1,24
!	    MLST=1
!	    DO ML=MLST,NCF
!	      IF(NU(ML) .LT. 0.2998/FLAM(L))THEN
!	        T1=DLOG10(OBSF(ML-1))-DLOG10(NU(ML-1)*FLAM(L)/0.2998)
!	1       *(DLOG10(OBSF(ML)/OBSF(ML-1))/(DLOG10(NU(ML)/NU(ML-1))))
!	        T1=5.0*DLOG10(DIST)-2.5*T1+2.5*DLOG10(ZERO_POINT)
!	        MLST=ML
!	        GOTO 105
!	       END IF
!	    END DO
!105	  CONTINUE
!	    WRITE(LU_OUT,106)FLAM(L),ZERO_POINT,T1
!106	    FORMAT(2X,F7.4,5X,F7.1,5X,F6.2)
!107	    FORMAT(//,'     Lam          ZP     Mag//')
!	  END DO
!	  CLOSE(UNIT=LU_OUT)
!
C
	ELSE IF(X .EQ. 'BB')THEN
	  
	  CALL USR_HIDDEN(OVER,'OVER','F','Overwrite existing model data')
	  CALL USR_OPTION(TEMP,'TEMP','3.0',' ')
!
	  IF(NCF .EQ. 0)THEN
	     T2=EXP(LOG(1.0D+4)/(NBB-1))
	     T1=3.289D0*912.0D0/10.0D0
	     CALL USR_HIDDEN(LIN_INT,'NORM','T','Normalize peak to unity?')
	     DO I=1,NBB
	       NU(I)=T1/(T2**I)
	       T3=HDKT*NU(I)/TEMP
	       IF(T3 .GT. 1.0D0)THEN
	         YV(I)=TWOHCSQ*(NU(I)**3)*DEXP(-T3)/(1.0D0-DEXP(-T3))
	       ELSE
	         YV(I)=TWOHCSQ*(NU(I)**3)/(DEXP(T3)-1.0D0)
	       END IF
	       XV(I)=NU(I)
	     END DO
	     IF(LIN_INT)THEN
	       T1=MAXVAL(YV(1:NBB))
	       YV(1:NBB)=YV(1:NBB)/T1
	     END IF
	  ELSE
	    CALL USR_OPTION(NORM_WAVE,'NW',' ',
	1       'Norm Wave (Angstroms) or Radius (<0 [in Rsun])')
C
C Rather than choose the full NCF points we adopt a set that extends from
C NU_MAX to NU_MIN but with a larger spacing.
C
	    DNU=DLOG10(NU(NCF)/NU(1))/(NBB-1)
	    DO I=1,NBB
	      T1=NU(1)*10.0D0**(DNU*(I-1))
	       T3=HDKT*T1/TEMP
	       IF(T3 .GT. 1.0D0)THEN
	         YV(I)=TWOHCSQ*(T1**3)*DEXP(-T3)/(1.0D0-DEXP(-T3))
	       ELSE
	         YV(I)=TWOHCSQ*(T1**3)/(DEXP(T3)-1.0D0)
	       END IF
	      XV(I)=T1
	    END DO         
	    IF(NORM_WAVE .GT. 0.0D0)THEN
	      NORM_FREQ=2.998D+03/NORM_WAVE
	      DO I=2,NCF
	        IF(NORM_FREQ .LE. NU(I-1) .AND. NORM_FREQ .GT. NU(I))K=I
	      END DO
	      T1=LOG(NU(K-1)/NORM_FREQ)/LOG(NU(K-1)/NU(K))
	      NORM_FLUX=EXP( (1.0D0-T1)*LOG(OBSF(K-1))+T1*LOG(OBSF(K)) )
	      BB_FLUX=TWOHCSQ*(NORM_FREQ**3)/(DEXP(HDKT*NORM_FREQ/TEMP)-1.0D0)
	      SCALE_FAC=NORM_FLUX/BB_FLUX
	      DO I=1,NBB
	        YV(I)=YV(I)*SCALE_FAC
	      END DO
	    ELSE
	      SCALE_FAC=NORM_WAVE*NORM_WAVE*159.8413D0     	 !1E+23*pi*(Rsun/1kpc)**2
	      DO I=1,NBB
	        YV(I)=YV(I)*SCALE_FAC
	      END DO
	    END IF
	  END IF
	  IF(OVER)THEN
	    OBSF(1:NBB)=YV(1:NBB)
	    NU(1:NBB)=XV(1:NBB)
	    NCF=NBB
	    WRITE(T_OUT,*)'New model data replaces old data'
	    WRITE(T_OUT,*)'No plots done with new model data'
	    WRITE(T_OUT,*)'No scaling done with new model data'
	  ELSE
	    CALL CNVRT(XV,YV,NBB,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1                 LAMC,XAXIS,YAXIS,L_FALSE)
	    IF(LIN_INT .AND. NCF .EQ. 0)THEN
	       T1=MAXVAL(YV(1:NBB))
	       YV(1:NBB)=YV(1:NBB)/T1
	    END IF
	    CALL CURVE(NBB,XV,YV)
	  END IF
!
	ELSE IF(X(1:5) .EQ. 'AV_EN')THEN
	  T1=0.0D0
	  T2=0.0D0
	  DO I=1,NCF-1
	    T1=T1+0.5D0*(NU(I)-NU(I+1))*(OBSF(I)+OBSF(I+1))
	    T2=T2+0.5D0*(NU(I)-NU(I+1))*(OBSF(I)/NU(I)+OBSF(I+1)/NU(I+1))
	  END DO
	  WRITE(6,*)'Average photon frequency is (in 10^15 Hz)',T1/T2
	  WRITE(6,*)'Effective photon energy is',6.626075D-12*T1/T2/2.7D0/1.380658D-16,'K'
C 
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
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:4) .EQ. 'LIST' .OR. 
	1                             X(1:2) .EQ. 'HE' 
	1          .OR. X(1:4) .EQ. 'HELP')THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Please see the HTML web pages in $CMFDIST/txt_files for help'
	  WRITE(T_OUT,*)' '
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
C
	FUNCTION FAC(N)
	REAL*8 FAC
	INTEGER N
	INTEGER, PARAMETER :: T_OUT=5
C
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
C
	RETURN
   	END
