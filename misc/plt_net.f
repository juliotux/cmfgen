!
! Porgram design to plot data from an XzVPRRR file.
!
! NS is read in from MODEL
! The velocity os read in from RVTJ.
!
	PROGRAM PLT_NRR
	USE MOD_USR_OPTION
	USE MOD_USR_HIDDEN
	USE MOD_WR_STRING
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 12-Sep-2016: Extensively modified. Routine is now option driven (as for PLTSPEC and DISPGEN).
!
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: DI(:)
!
	REAL*8, ALLOCATABLE :: TOT_PHOT_RATE(:)
	REAL*8, ALLOCATABLE :: TOT_REC_RATE(:)
	REAL*8, ALLOCATABLE :: COL_ION_RATE(:)
	REAL*8, ALLOCATABLE :: COL_REC_RATE(:)
	REAL*8, ALLOCATABLE :: ADV_REC_RATE(:)
	REAL*8, ALLOCATABLE :: NT_ION_RATE(:)
	REAL*8, ALLOCATABLE :: XRAY_REC_RATE_PS(:)
	REAL*8, ALLOCATABLE :: XRAY_REC_RATE_GS(:)
	REAL*8, ALLOCATABLE :: REC_COEF(:)
	REAL*8, ALLOCATABLE :: NET(:)
!
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
	REAL*8, ALLOCATABLE :: ZV(:)
!
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=80) DIR_NAME
	CHARACTER(LEN=200) STRING
!
        CHARACTER(LEN=80) MAIN_OPT_STR
        CHARACTER(LEN=20) X
        CHARACTER(LEN=10) XOPT
        CHARACTER(LEN=10) XSPEC
        CHARACTER(LEN=20) DEFAULT
        CHARACTER(LEN=80) DESCRIPTION
	CHARACTER(LEN=50) XLABEL,YLABEL
	CHARACTER(LEN=30) UC
	EXTERNAL UC
!
	REAL*8 SUM_VAL
	INTEGER I,J,IBEG
	INTEGER IOS
	INTEGER ND
	INTEGER LEN_DIR
	LOGICAL RADIUS_DONE
	LOGICAL FILE_OPEN
!
! We get the file, and open it here in case we need a directory name.
!
 	FILENAME=' '
	CALL GEN_IN(FILENAME,'File with data to be plotted')
	OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Unable to open the file: ',TRIM(FILENAME)
	  WRITE(6,*)'IOSTAT=',IOS
	  STOP
	END IF
!
! Get default directory.
!
        DIR_NAME=' '            !Valid DIR_NAME if not present.
        LEN_DIR=0
        J=LEN(FILENAME)
        DO WHILE(J .GT. 0)
          IF( FILENAME(J:J) .EQ. ']' .OR.
	1     FILENAME(J:J) .EQ. ':' .OR.
	1     FILENAME(J:J) .EQ. '/'        )THEN
            DIR_NAME=FILENAME(1:J)
            LEN_DIR=J
            J=0
          END IF
          J=J-1
        END DO
!
	IF(LEN_DIR .EQ. 0)THEN
	  STRING='MODEL'
	ELSE
	  STRING=DIR_NAME//'MODEL'
	END IF
	ND=0
	OPEN(UNIT=20,FILE=TRIM(STRING),STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	        READ(STRING,*)ND
	        WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	        EXIT
	      END IF
	    END DO
	  END IF
	  INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
	IF(IOS .NE. 0 .OR. ND .EQ. 0)THEN
	  CALL GEN_IN(ND,'What is the ND value')
	END IF
!
	ALLOCATE(R(ND));     R=0.0D0
	ALLOCATE(V(ND));     V=0.0D0
	ALLOCATE(T(ND));     T=0.0D0
	ALLOCATE(ED(ND));    ED=0.0D0
	ALLOCATE(DI(ND));    DI=0.0D0
!
	ALLOCATE(TOT_PHOT_RATE(ND));      TOT_PHOT_RATE=0.0D0 
	ALLOCATE(TOT_REC_RATE(ND));       TOT_REC_RATE=0.0D0
	ALLOCATE(COL_ION_RATE(ND));       COL_ION_RATE=0.0D0 
	ALLOCATE(COL_REC_RATE(ND));       COL_REC_RATE=0.D0
	ALLOCATE(ADV_REC_RATE(ND));       ADV_REC_RATE=0.0D0
	ALLOCATE(NT_ION_RATE(ND));        NT_ION_RATE=0.0D0
	ALLOCATE(XRAY_REC_RATE_PS(ND));   XRAY_REC_RATE_PS=0.0D0
	ALLOCATE(XRAY_REC_RATE_GS(ND));   XRAY_REC_RATE_GS=0.0D0
	ALLOCATE(REC_COEF(ND));           REC_COEF=0.0D0
	ALLOCATE(NET(ND));                NET=0.0D0
!
	ALLOCATE(XV(ND))
	ALLOCATE(YV(ND))
	ALLOCATE(ZV(ND))
!
	STRING=' '
	DO IBEG=1,ND,10
!
	  RADIUS_DONE=.FALSE.
	  DO WHILE(1 .EQ. 1)
	    IF(STRING .EQ. ' ')THEN
	    ELSE IF(INDEX(STRING,'Radius') .NE. 0)THEN
	      READ(11,*)(R(I),I=IBEG,MIN(IBEG+9,ND))
	      RADIUS_DONE=.TRUE.
	      WRITE(6,*)IBEG,R(IBEG)
	    ELSE IF(INDEX(STRING,'Temperature') .NE. 0)THEN
	      READ(11,*)(T(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Electron Density') .NE. 0)THEN
	       READ(11,*)(ED(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Ion Density') .NE. 0)THEN
	      READ(11,*)(DI(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Photoionization Rates') .NE. 0)THEN
!
! The following assumes there is a blanke line after the photoioization rates.
!
	      DO WHILE(1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(STRING .EQ. ' ')EXIT
	        READ(STRING,*)(YV(I),I=IBEG,MIN(IBEG+9,ND))
	        DO I=IBEG,MIN(IBEG+9,ND)
	          TOT_PHOT_RATE(I)=TOT_PHOT_RATE(I)+YV(I)
	        END DO
	      END DO
	    ELSE IF(INDEX(STRING,'Colisional Ionization Rate') .NE. 0)THEN
	      READ(11,*)(COL_ION_RATE(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Recombination Rates') .NE. 0)THEN
!
! The following assumes there is a blanke line after the recombination rates.
!
	      DO WHILE(1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(STRING .EQ. ' ')EXIT
	        READ(STRING,*)(YV(I),I=IBEG,MIN(IBEG+9,ND))
	        DO I=IBEG,MIN(IBEG+9,ND)
	          TOT_REC_RATE(I)=TOT_REC_RATE(I)+YV(I)
	        END DO
	      END DO
	    ELSE IF(INDEX(STRING,'Colisional Recombination Rate') .NE. 0)THEN
	      READ(11,*)(COL_REC_RATE(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Effective Advection Recombination Rate') .NE. 0)THEN
	      READ(11,*)(ADV_REC_RATE(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Non-Thermal Ionization  Rate') .NE. 0)THEN
	      READ(11,*)(NT_ION_RATE(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Net X-ray recombination rate (to previous ionization state)') .NE. 0)THEN
	      READ(11,*)(XRAY_REC_RATE_PS(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Net X-ray recombination rate to Ar2 g.s.)') .NE. 0)THEN
	      READ(11,*)(XRAY_REC_RATE_GS(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Net Recombination Rate (% of total)') .NE. 0)THEN
	      READ(11,*)(NET(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE IF(INDEX(STRING,'Radiative Recombination Coefficient for explicitly treated levels.') .NE. 0)THEN
	      READ(11,*)(REC_COEF(I),I=IBEG,MIN(IBEG+9,ND))
	    ELSE
	    END IF
	    READ(11,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)EXIT
	    IF(RADIUS_DONE .AND. INDEX(STRING,'Radius') .NE. 0)EXIT
	  END DO
	END DO
	CLOSE(UNIT=11)
!
	IF(V(1) .EQ. 0.0D0)THEN
	  IF(LEN_DIR .EQ. 0)THEN
	    STRING='RVTJ'
	  ELSE
	    STRING=DIR_NAME//'RVTJ'
	  END IF
	  CALL RD_SING_VEC_RVTJ(V,ND,'Velocity',STRING,20,IOS)
	END IF
!
! Set plot defaults.
!
        XV=R
	YV=0.0D0
        XLABEL='R'
        YLABEL=' '
!
! This call resets the .sve algorithm.  Specifically it sets the next
! input answer to be a main option, and all subsequent inputs to be
! sub-options.  
!
 3      CALL SVE_FILE('RESET')
!
        DESCRIPTION=' '                         !As obvious.
        MAIN_OPT_STR='  '
        DEFAULT='GR'
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
! Remove variable changes from main option.
!
        I=INDEX(MAIN_OPT_STR,'(')
        IF(I .EQ. 0)THEN
          X=UC(TRIM(MAIN_OPT_STR))
        ELSE
          X=UC(MAIN_OPT_STR(1:I-1))     !Remove line variables.
        END IF
!
! Remove possile file names etc.
!          
        I=INDEX(X,' ')
        IF(I .EQ. 0)THEN
          X=UC(TRIM(X))
        ELSE
          X=UC(X(1:I-1))        !Remove file names
        END IF
!
! NB: X contains the full option
!     XOPT contains the main part of the option with the species identifier
!     XSPEC is the species identifier (eg. HE2)
!
! _ is presently used to separate the OPTION (e.g. DC) from the species
! (He2). This is the only location it needs ti specified.
! 
!        I=INDEX(X,'_')
!        IF(I .EQ. 0)THEN
!          XSPEC=' '
!          XOPT=X
!        ELSE
!          J=INDEX(X,' ')
!          XOPT=X(1:I-1)
!          XSPEC=X(I+1:J)
!        END IF
!
! Since we are not using species options, we commented out previous few lines.
! However we need to set XOPT, as done below.
!
	XOPT=X
!
! SUM_VAL is used as check whether the rate is available.
!
	SUM_VAL=1.0D+50			!Some non-zero value in case X option
	IF(XOPT .EQ. 'XVEL')THEN
	  XV=V
	  XLABEL='V(km/s)'
	  IF(V(1) .GE. 10000)THEN
	    XV=1.0D-03*V
	    XLABEL='V(Mm/s)'
	  END IF
        ELSE IF(XOPT .EQ. 'XLOGV')THEN
          XV=LOG10(V)
	  XLABEL='Log V(km/s)'
	  IF(V(1) .GE. 10000)THEN
	    XV=XV-3.0D0
	    XLABEL='Log V(Mm/s)'
	  END IF
	ELSE IF(XOPT .EQ. 'XN')THEN
          DO I=1,ND
            XV(I)=I
          END DO
          XLABEL='Index'
        ELSE IF(XOPT .EQ. 'XR')THEN
          XV=R
          XLABEL='R'
        ELSE IF(XOPT .EQ. 'XLOGR')THEN
          XV=LOG10(R)
          XLABEL='Log(R[10\u10\d cm])'
!
	ELSE IF(XOPT .EQ. 'NT')THEN
	  YV=NT_ION_RATE
	  IF(SUM _VAL .NE. 0.0D0)CALL DP_CURVE(ND,XV,YV)
	ELSE IF(XOPT .EQ. 'COL_ION')THEN
	  YV=COL_ION_RATE
	  SUM_VAL=SUM(YV)
	  IF(SUM _VAL .NE. 0.0D0)CALL DP_CURVE(ND,XV,YV)
	ELSE IF(XOPT .EQ. 'COL_REC')THEN
	  YV=COL_REC_RATE
	  SUM_VAL=SUM(YV)
	  IF(SUM _VAL .NE. 0.0D0)CALL DP_CURVE(ND,XV,YV)
	ELSE IF(XOPT .EQ. 'REC_RATE')THEN
	  YV=TOT_REC_RATE
	  SUM_VAL=SUM(YV)
	  IF(SUM _VAL .NE. 0.0D0)CALL DP_CURVE(ND,XV,YV)
!
	ELSE IF(XOPT .EQ. 'PHOT_RATE')THEN
	  YV=TOT_PHOT_RATE
	  SUM_VAL=SUM(YV)
	  IF(SUM _VAL .NE. 0.0D0)CALL DP_CURVE(ND,XV,YV)
	ELSE IF(XOPT .EQ. 'REC_COEF')THEN
	  YV=REC_COEF
	  SUM_VAL=SUM(YV)
	  IF(SUM _VAL .NE. 0.0D0)CALL DP_CURVE(ND,XV,YV)
	ELSE IF(XOPT .EQ. 'NET')THEN
	  YV=NET
	  SUM_VAL=SUM(YV)
	  IF(SUM _VAL .NE. 0.0D0)CALL DP_CURVE(ND,XV,YV)
!
        ELSE IF(XOPT .EQ. 'P' .OR. XOPT .EQ. 'GR')THEN
          CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
        ELSE IF(XOPT .EQ. 'HE' .OR. XOPT .EQ. 'HELP')THEN
	  WRITE(6,'(A1)')' '
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'XR',BLUE_PEN,'Set X axis to R'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'XLOGR',BLUE_PEN,'Set X axis to Log(R)'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'XV',BLUE_PEN,'Set X axis to velocity'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'XLOGV',BLUE_PEN,'Set X axis to Log(velocity)'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'XN',BLUE_PEN,'Set X axis to depth index'
	  WRITE(6,'(A1)')' '
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'COL_ION',BLUE_PEN,'Plot the collisional ionization rate'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'COL_ION',BLUE_PEN,'Plot the collisional recombination rate'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'PHOT_RATE',BLUE_PEN,'Plot the photoionization rate'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'REC_RATE',BLUE_PEN,'Plot the recombination rate'
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'NT',BLUE_PEN,'Plot the non-thermal ionization rate'
	  WRITE(6,'(A1)')' '
	  WRITE(6,'(1X,A1,A,A1,T20,A)')RED_PEN,'EX',BLUE_PEN,'Exit from porgram'
	  WRITE(6,'(A1)')DEF_PEN
        ELSE IF(XOPT .EQ. 'EX')THEN
          STOP
        ELSE
          WRITE(6,*)'Option not recognized'
        END IF
	IF(SUM_VAL .EQ. 0.0D0)THEN
	   WRITE(6,*)'No call to CURVE as data values for requested option are zer'
	END IF
        GOTO 3
!
	END
