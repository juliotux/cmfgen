!
! General routine for plotting and comparing dFR(R,nu) obtained from dFR_DATA.
!
! dFR(I,ML) gives the contribution to the observed flux at frequency NU(ML)
!           for the interval R(I) to R(I+1). The data should be plotted in
!           historgam mode.
!
	PROGRAM PLT_dFR
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
! Altered:  2-Nov-2013 - Remove an excess allocation of XV_SAV.
!                          USe R_RVTJ for R grid from RVTJ file to allow a check for inconsistencies.
! Created: 20-Feb-2012 - Based in PLT_IP
!
	INTEGER NCF
	INTEGER ND
	INTEGER NC
	INTEGER NC2,ND2,NP2
	REAL*8, ALLOCATABLE :: dFR(:,:)
	REAL*8, ALLOCATABLE :: NU(:)
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: TAU_FLUX(:)
	REAL*8, ALLOCATABLE :: TAU_ROSS(:)
	REAL*8, ALLOCATABLE :: TAU_ES(:)
!
	REAL*8 RMDOT
	REAL*8 RLUM
	REAL*8 ABUND_HYD
	REAL*8 LAMBDA
	REAL*8, ALLOCATABLE :: R_RVTJ(:)
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
	CHARACTER(LEN=21) TIME
!
! Vectors for passing data to plot package via calls to CURVE.
!
	REAL*4, ALLOCATABLE :: XV(:)
	REAL*4, ALLOCATABLE :: XV_SAV(:)
	REAL*4, ALLOCATABLE :: YV(:)
	REAL*4, ALLOCATABLE :: WV(:)
	REAL*4, ALLOCATABLE :: ZV(:)
!
	CHARACTER(LEN=80) NAME		!Default title for plot
	CHARACTER(LEN=80) XAXIS
	CHARACTER(LEN=80) XAXIS_SAV
	CHARACTER(LEN=80) YAXIS		!Label for Ordinate
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
!
! Miscellaneous variables.
!
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,K,L,ML,LS
	INTEGER ST_REC
	INTEGER REC_LENGTH
	REAL*8 PI
	REAL*8 FREQ
	REAL*8 T1,T2,T3
	REAL*8 RPHOT,VPHOT
	REAL*8 FRAC
	REAL*8 LAMC
	LOGICAL AIR_LAM
!
	REAL*8, PARAMETER :: RZERO=0.0D0
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
	CHARACTER(LEN=80) MAIN_OPT_STR		!Used for input of the main otion
	CHARACTER(LEN=10) X			!Used for the idividual option
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=120) DEFAULT
	CHARACTER(LEN=120) DESCRIPTION
	CHARACTER(LEN=20) TMP_STR
!
	REAL*8 SPEED_OF_LIGHT,FAC,LAM_VAC,PARSEC,LAMVACAIR
	LOGICAL EQUAL
	CHARACTER(LEN=30) UC
	EXTERNAL SPEED_OF_LIGHT,EQUAL,FAC,UC,LAM_VAC,LAMVACAIR,GET_INDX_DP,PARSEC
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
	DISTANCE=1.0D0		!kpc
!
! Conversion factor from Kev to units of 10^15 Hz.
! Conversion factor from Angstroms to units of 10^15 Hz.
!
	KEV_TO_HZ=0.241838D+03
	ANG_TO_HZ=SPEED_OF_LIGHT()*1.0D-07  	!10^8/10^15
!
!  Read in model.
!
	FILENAME='dFR_DATA'
	CALL GEN_IN(FILENAME,'dF(R) file')
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read dFR_DATA_INFO'
	  STOP
	END IF
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED')
	  READ(LU_IN,REC=3)ST_REC,NCF,ND
	  ALLOCATE (dFR(ND,NCF))
	  ALLOCATE (NU(NCF))
	  ALLOCATE (R(ND))
	  IF( INDEX(FILE_DATE,'20-Aug-2000') .NE. 0)THEN
	    READ(LU_IN,REC=ST_REC)(R(I),I=1,ND)
	    ST_REC=ST_REC+1
	  ELSE
	    WRITE(T_OUT,*)'Unrecognized date when reading IP_DATA' 
	    WRITE(T_OUT,*)'Date=',FILE_DATE
	    STOP
	  END IF
	  DO ML=1,NCF
	    READ(LU_IN,REC=ST_REC+ML-1)(dFR(I,ML),I=1,ND),NU(ML)
	  END DO
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in dFR_DATA file as MODEL A (default)'
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
	1             ND2,NC2,NP2,FILENAME,LU_IN)
	ALLOCATE (R_RVTJ(ND))
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
	CALL RD_RVTJ_VEC(R_RVTJ,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POPTOM,POPION,MASS_DENSITY,CLUMP_FAC,ND,LU_IN)
	CLOSE(LU_IN)
!
	DO I=1,ND
	  IF( ABS(1.0D0-R(I)/R_RVTJ(I)) .GT. 1.0D-08 )THEN
	    WRITE(6,*)' '
	    WRITE(6,*)' ERORR -- R in dFR_DATA not the same as in RVTJ'
	    WRITE(6,*)' The may cause unexpected an wrong bahaviour'
	    WRITE(6,*)' '
	    STOP
	  END IF
	END DO
!
! Now compute the important optical depth scales.
!
	 ALLOCATE (TA(ND))
	 ALLOCATE (TB(ND))
	 ALLOCATE (TC(ND))
	 ALLOCATE (TAU_ROSS(ND)); TAU_ROSS=0.0D0
	 ALLOCATE (TAU_FLUX(ND)); TAU_FLUX=0.0D0
	 ALLOCATE (TAU_ES(ND));   TAU_ES=0.0D0
	 IF(ROSS_MEAN(ND) .NE. 0)THEN
	   TA(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	   CALL TORSCL(TAU_ROSS,TA,R,TB,TC,ND,METHOD,TYPETM)
	 END IF
	 IF(FLUX_MEAN(ND) .NE. 0)THEN
	   TA(1:ND)=FLUX_MEAN(1:ND)*CLUMP_FAC(1:ND)
	   CALL TORSCL(TAU_FLUX,TA,R,TB,TC,ND,METHOD,TYPETM)
	 END IF
	 TA(1:ND)=6.65D-15*ED(1:ND)
	 CALL TORSCL(TAU_ES,TA,R,TB,TC,ND,METHOD,TYPETM)
!
	I=MAX(ND,NCF)
	ALLOCATE (XV(I))
	ALLOCATE (XV_SAV(I))
	ALLOCATE (YV(I))
	ALLOCATE (WV(I))
	ALLOCATE (ZV(I))
	XV_SAV(1:ND)=LOG10(R/R(ND))
	XAXIS_SAV='Log R/R\d*\u)'
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
	ELSE IF(X(1:2) .EQ. 'XN')THEN
	  DO I=1,ND
	    XV_SAV(I)=I
	  END DO
	  XAXIS_SAV='Depth index'
	ELSE IF(X(1:5) .EQ. 'XLOGR')THEN
	  XV_SAV(1:ND)=LOG10(R/R(ND))
	  XAXIS_SAV='Log R/R\d*\u)'
	ELSE IF(X(1:5) .EQ. 'XLINR')THEN
	  XV_SAV(1:ND)=R/R(ND)
	  XAXIS_SAV='R/R\d*\u)'
	ELSE IF(X(1:5) .EQ. 'XLOGV')THEN
	  XV_SAV(1:ND)=LOG10(V)
	  XAXIS_SAV='Log V(km\u \ds\u-1\d)'
	ELSE IF(X(1:5) .EQ. 'XED')THEN
	  XV_SAV(1:ND)=LOG10(ED)
	  XAXIS_SAV='Log Ne(cm\u-3\d)'
	ELSE IF(X(1:2) .EQ. 'XV')THEN
	  IF(V(1) .GT. 10000.0D0)THEN
	    XV_SAV(1:ND)=1.0D-03*V
	    XAXIS_SAV='V(Mm\u \ds\u-1\d)'
	  ELSE
	    XV_SAV(1:ND)=V
	    XAXIS_SAV='V(km\u \ds\u-1\d)'
	  END IF
	ELSE IF(X(1:5) .EQ. 'XROSS')THEN
	  IF(TAU_ROSS(10) .NE. 0.0D0)THEN
	    XV_SAV(1:ND)=LOG10(TAU_ROSS)
	    XAXIS_SAV='Log \gt(Ross)'
	  ELSE
	    WRITE(6,*)'Tau(Ross) is not available'
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
	ELSE IF(X(1:2) .EQ. 'SP')THEN
!
	  XV(1:NCF)=NU(1:NCF)
	  DO ML=1,NCF
	    YV(ML)=SUM(dFR(:,ML))
	  END DO
!
! NB: J and I have the same units, apart from per steradian/
!
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
	  CALL CURVE(NCF,XV,YV)
	  J=INDEX(YAXIS,')')
	  YAXIS(J:)='\d \ukpc\u2\d)'
!
	ELSE IF(X(1:2) .EQ. 'IR')THEN
	  WRITE(6,*)'Option to plot the contribution at a given radius.'
	  CALL USR_OPTION(I,'R',' ','Radial index')
	  IF(I .GT. ND)THEN
	    WRITE(T_OUT,*)'Invalid depth; maximum value is',ND
	    GOTO 1
	  END IF
!
	  XV(1:NCF)=NU(1:NCF)
	  YV(1:NCF)=dFR(I,1:NCF)
!
! NB: J and I have the same units, apart from per steradian/
!
	  CALL CNVRT(XV,YV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_FALSE)
!
	  J=INDEX(YAXIS,'J')
	  YAXIS(J:J)='I'
	  CALL CURVE(NCF,XV,YV)
!
	ELSE IF(X(1:2) .EQ. 'OR' .OR. X(1:2) .EQ. 'OV' .OR. X(1:2) .EQ. 'OD')THEN
	  WRITE(6,*)'Option to plot the radius for a given fractional contribution'
	  CALL USR_OPTION(FRAC,'FRAC',' ','Fractional value')
	  IF(FRAC .LE. 0 .OR. FRAC .GT. 1.0)THEN
	    WRITE(T_OUT,*)'Invalid fractional value'
	    GOTO 1
	  END IF
!
	  DO ML=1,NCF
	    ZV(ML)=SUM(dFR(:,ML))
	  END DO
!
	  XV(1:NCF)=NU(1:NCF)
	  DO ML=1,NCF
	    TA(1)=0
	    IF(ZV(ML) .GT. 0.0D0)THEN
	      DO I=2,ND
	        TA(I)=TA(I-1)+dFR(I,ML)
	        IF(TA(I) .GT. FRAC*ZV(ML))THEN
	          T1=TA(I-1)/ZV(ML)
	          T2=TA(I)/ZV(ML)
	          T3=(FRAC-T1)/(T2-T1)
	          IF(X(1:2) .EQ. 'OR')THEN
	            YV(ML)=(T3*R(I)+(1.0D0-T3)*R(I-1))/R(ND)
	          ELSE IF(X(1:2) .EQ. 'OD')THEN
	            YV(ML)=(T3*I+(1.0D0-T3)*(I-1))
	          ELSE
	            YV(ML)=1.0D-03*(T3*V(I)+(1.0D0-T3)*V(I-1))
	          END IF
	          EXIT
	        END IF
	      END DO
	    ELSE
	      YV(ML)=0.0D0
	    END IF
	  END DO
!
	  IF(X(1:2) .EQ. 'OR')THEN
	    YAXIS='R/R(ND)'
	  ELSE IF(X(1:2) .EQ. 'OD')THEN
	    YAXIS='Depth index'
	  ELSE
	    YAXIS='V(Mm/s)'
	  END IF
!
	  T1=0.0D0; T2=0.0D0
	  DO ML=2,NCF
	    T1=T1+(ZV(ML+1)+ZV(ML))*(XV(ML)-XV(ML+1))
	    T2=T2+(YV(ML+1)*ZV(ML+1)+YV(ML)*ZV(ML))*(XV(ML)-XV(ML+1))
	  END DO
!
	  IF(X(1:2) .EQ. 'OR')THEN
	    RPHOT=R(ND)*T2/T1
	    DO I=2,ND
	      IF(RPHOT .GE. R(I))EXIT
	    END DO
	    T3=(RPHOT-R(I-1))/(R(I)-R(I-1))
	    VPHOT=V(I-1)+T3*(V(I)-V(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric  velocity is',VPHOT,' Mm/s'
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric   radius  is',RPHOT,' 10^10 cm'
	    T1=TAU_FLUX(I-1)+T3*(TAU_FLUX(I)-TAU_FLUX(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric tau(flux) is',T1
	    T1=TAU_ROSS(I-1)+T3*(TAU_ROSS(I)-TAU_ROSS(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric tau(Ross) is',T1
	    T1=TAU_ES(I-1)+T3*(TAU_ES(I)-TAU_ES(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric   tau(es) is',T1
	  ELSE IF(X(1:2) .EQ. 'OV')THEN
	    VPHOT=T2/T1
	    DO I=2,ND
	      IF(VPHOT .GE. 0.001D0*V(I))EXIT
	    END DO
	    T3=(1000.0D0*VPHOT-V(I-1))/(V(I)-V(I-1))
	    RPHOT=R(I-1)+T3*(R(I)-R(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric  velocity is',VPHOT,' Mm/s'
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric   radius  is',RPHOT,' 10^10 cm'
	    T1=TAU_FLUX(I-1)+T3*(TAU_FLUX(I)-TAU_FLUX(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric tau(flux) is',T1
	    T1=TAU_ROSS(I-1)+T3*(TAU_ROSS(I)-TAU_ROSS(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric tau(Ross) is',T1
	    T1=TAU_ES(I-1)+T3*(TAU_ES(I)-TAU_ES(I-1))
	    WRITE(6,'(A,ES14.4,2A)')' The flux weighted photospheric   tau(es) is',T1
	  END IF
!
	  CALL CNVRT(XV,ZV,NCF,LOG_X,LOG_Y,X_UNIT,Y_PLT_OPT,
	1         LAMC,XAXIS,YAXIS,L_TRUE)
	  CALL CURVE(NCF,XV,YV)
!
	ELSE IF(X .EQ. 'WF2')THEN
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
	  WRITE(30,*)ND,J-I+2
	  WRITE(30,'(1000ES14.6)')(R(L),L=1,ND)
	  WRITE(30,'(1000ES15.7)')ANG_TO_HZ/NU(I:J+1)
	  DO K=I,J+1
	    TA(1:ND)=0.0D0
	    DO ML=K-4,K+3
	      TA(1:ND-1)=TA(1:ND-1)+(dFR(1:ND-1,ML)+dFR(1:ND-1,ML+1))*(NU(ML)-NU(ML+1))
	    END DO
	    TA(1:ND-1)=0.5D0*TA(1:ND-1)/(NU(K-4)-NU(K+4))
	    WRITE(30,'(500ES12.4)')(TA(L),L=1,ND-1)
	  END DO
!
	ELSE IF(X .EQ. 'DF2' .OR. X .EQ. 'DDF2')THEN
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
! Average dFR over frequnecy.
!
	  XV(1:ND)=XV_SAV(1:ND)
	  XAXIS=XAXIS_SAV
	  YV(1:ND)=0.0D0
	  K=MIN(I,J); J=MAX(I,J); I=K
	  IF(J .EQ. I)J=I+1
	  DO K=I,J-1
	    YV(1:ND)=YV(1:ND)+0.5D0*(dFR(1:ND,K)+dFR(1:ND,K+1))*(NU(K)-NU(K+1))
	  END DO
	  T1=ABS(NU(I)-NU(J))
	  YV(1:ND)=YV(1:ND)/T1
	  WRITE(6,*)'Start wavelength of interval is (AIR if > 2000A)',LAMVACAIR(NU(I))
	  WRITE(6,*)'  End wavelength of interval is (AIR if > 2000A)',LAMVACAIR(NU(J))
!
	  FREQ=0.5D0*(NU(I)+NU(J))
	  IF(Y_PLT_OPT .EQ. 'NU_FNU')THEN
	    T1=FREQ*1.0D-08
	    YV(1:ND)=YV(1:ND)*T1
	    YAXIS='\gvF\d\gv\u(ergs\d cm\u-2 \ds kpc\u2\d)'
	  ELSE IF(Y_PLT_OPT .EQ. 'FLAM')THEN
	    T1=T1*T1*1.0D-01/C_CMS
	    YV(1:ND)=YV(1:ND)*T1
	    YAXIS='F\d\gl\u(ergs\d cm\u-2 \ds\u-1\da kpc\u2\d\A'
	  ELSE IF(Y_PLT_OPT .EQ. 'FNU')THEN
	    YAXIS='F\d\gn\u(Jy kpc\u2\d)'
	  END IF
!
	  IF(X .EQ. 'DDF2')THEN
	    DO I=1,ND-1
	      YV(I)=YV(I)/ABS(XV(I+1)-XV(I))
	    END DO
	    J=INDEX(YAXIS,')')
	    YAXIS='d'//TRIM(YAXIS)//'dX'
	    WRITE(6,*)'Use CUM option in PGPLOT to get intgegrated flux'
	    CALL CURVE(ND-1,XV,YV)
	  ELSE
	    YAXIS='\gd'//TRIM(YAXIS)
	    DO I=1,ND-1
	      XV(I)=0.5D0*(XV(I+1)+XV(I))
	    END DO
	    WRITE(6,*)'Use SUM option in PGPLOT to get intgegrated flux'
	    CALL CURVE(ND-1,XV,YV)
	  END IF
!
	ELSE IF(X .EQ. 'DF' .OR. X .EQ. 'DDF')THEN
!
	  CALL USR_OPTION(T1,'Lambda',' ','Wavelength in Ang')
	  FREQ=0.299794D+04/T1
          I=GET_INDX_DP(FREQ,NU,NCF)
	  IF(NU(I)-FREQ .GT. FREQ-NU(I+1))I=I+1
!
	  XV(1:ND)=XV_SAV(1:ND)
	  XAXIS=XAXIS_SAV
	  YV(1:ND)=dFR(1:ND,I)
	  IF(Y_PLT_OPT .EQ. 'NU_FNU')THEN
	    T1=FREQ*1.0D-08
	    YV(1:ND)=YV(1:ND)*T1
	    YAXIS='\gnF\dgn\u(ergs\d cm\u-2 \ds kpc\u2\d)'
	  ELSE IF(Y_PLT_OPT .EQ. 'FLAM')THEN
	    T1=FREQ*FREQ*1.0D-01/C_CMS
	    YV(1:ND)=YV(1:ND)*T1
	    YAXIS='F\d\gn\u(ergs\d \ukpc\u^2 \dcm\u-2 \ds\u-1 \d\A\u-1\d)'
	  ELSE IF(Y_PLT_OPT .EQ. 'FNU')THEN
	    YAXIS='F\d\gn\u(Jy kpc\u2\d)'
	  END IF
	  IF(X .EQ. 'DDF')THEN
	    DO I=1,ND-1
	      YV(I)=YV(I)/ABS(XV(I+1)-XV(I))
	    END DO
	    J=INDEX(YAXIS,')')
	    YAXIS='d'//TRIM(YAXIS)//'dX'
	    WRITE(6,*)'Use CUM option in PGPLOT to get intgegrated flux'
	    CALL CURVE(ND-1,XV,YV)
	  ELSE
	    YAXIS='\gd'//TRIM(YAXIS)
	    DO I=1,ND-1
	      XV(I)=0.5D0*(XV(I+1)+XV(I))
	    END DO
	    WRITE(6,*)'Use SUM option in PGPLOT to get intgegrated flux'
	    WRITE(6,*)'Use H curve to plot'
	    CALL CURVE(ND-1,XV,YV)
	  END IF
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
	  CALL USR_HIDDEN(RAD_VEL,'RAD_VEL','0.0D0','Radial velcoity (+ve if away)')
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
	  NCF_MAX=MAX(NCF,ND)
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
! 
! Plot section:
!
	ELSE IF(X(1:2) .EQ. 'GR') THEN
	  CALL GRAMON_PGPLOT(XAXIS,YAXIS,NAME,' ')
	  XAXIS=XAXIS_SAV
!
	ELSE IF(X(1:4) .EQ.'GRNL') THEN
	  CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	  XAXIS=XAXIS_SAV
!
! 
!
	ELSE IF(X(1:2) .EQ. 'LI' .OR. X(1:4) .EQ. 'LIST' .OR. 
	1                             X(1:2) .EQ. 'HE' 
	1          .OR. X(1:4) .EQ. 'HELP')THEN
	  WRITE(6,*)' '
	  WRITE(6,*)' TIT  LX  LY  '
	  WRITE(6,*)' XN:      Set X axis to depth indexV'
	  WRITE(6,*)' XV:      Set X axis to V'
	  WRITE(6,*)' XED:     Set X axis to Log(Ne)'
	  WRITE(6,*)' XROSS:   Set X axis to Log Tau(Ross)'
	  WRITE(6,*)' XLINR:   Set X axis to R/R*'
	  WRITE(6,*)' XLOGV:   Set X axis to Log V'
	  WRITE(6,*)' XLOGR:   Set X axis to Log R/R*'
	  WRITE(6,*)' SP:      Plot spectrum - sums dF(R) over all radii'
	  WRITE(6,*)' '
	  WRITE(6,*)' DF:      Plot dF(R) for a given frequency'
	  WRITE(6,*)' DF2:     Plot dF(R) averaged over a given frequency band'
	  WRITE(6,*)' DDF:     Plots dF/dX (X=default axis) for a given frequency band'
	  WRITE(6,*)' DDF2:    Plots dF/dX (X=default axis) averaged over a given frequency band'
	  WRITE(6,*)' IR:      Plot the contribution, as a funtion of lambda, at a given radius'
	  WRITE(6,*)' RD_OBS:  Read on observed spectrum'
	  WRITE(6,*)' OR:      Plot the radius, as a funtion of lambda, for a given fractional contribution'
	  WRITE(6,*)' OV:      Plot the velocity, as a funtion of lambda, for a given fractional contribution'
	  WRITE(6,*)' ON:      Plot the depth index, as a funtion of lambda, for a given fractional contribution'
	  WRITE(6,*)' '
	  WRITE(6,*)' GR:      Plot data already sent to plot package'
	  WRITE(6,*)' GRNL:    Same as GR but labels are not sent.'
	  WRITE(6,*)' EX:      Exit from program'
	  WRITE(6,*)' '
	  GOTO 1
!
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
