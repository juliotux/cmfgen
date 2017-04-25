!
! General purpose line plotting routiine.
!
	SUBROUTINE GRAMON_PGPLOT(XLAB,YLAB,TITL,PASSED_OPT)
	USE NEW_GEN_IN_INTERFACE
	USE MOD_CURVE_DATA
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered:  01-Mar-2016 : Added XN option to XAR option. This allows Y to be plotted against
!                          the index I. By default, a new plot is created [24-Feb-2016].
! Altered:  29-Jun-2015 : Changed RID o check ABS(TAU) which for pp model can -ve.
! Altered:  22-Apr-2015 : Added FILL option to fill the space between two curves that create a polygon.
!                           ANS changed to length 4 (from 3)
! Altered:  17-Feb-2015 : Can now have multi-colored titles.
! Altered:  22-Jan-2015 : Bug fix. SC option for scrolling changed to SCR.
!                         SC is reserved for entering strings by cursor
! Altered:  14-Jan-2014 : Revised LG option
! Altered:  22-Nov-2013 : Added LG option for curve type. This plots the log of the absolute 
!                           value of the data but indicates where the data is -ve.
! Altered:  04-Sep-2013 : Increased MAXPEN (=MAX_PLOTS). Minor cleaning.
! Altered:  31-Aug-2013 : Added long-plot option.
! Altered:  26-Nov-2011 : Curves cycle over pen-colors 2 to 13.
!                         Dashed curve for plots > 13
!                         Marker style ignored when MARK is off (use I for invisible curve).
! Altered:  08-Apr-2006 : In NM option, normalization now correctly handles unequeally 
!                            spaced data.
! Altered:  23-Aug-2005 : Grey pen option installed.
! Altered:  26-Jan-2005 : Change default margins to give more room on borders.
! Altered:  30-May-2003 : Improvements to YAR and XAR options.
! Altered:  11-Feb-2002 : Bug fixed with EW option for reversed data.
! Altered:  04-Apr-2001 : Include option for differen PLT_ST filename.
!                         Code now ouputs name of stored plots when bad plot id entered.
! Altered:  15-Jul-2000 : MOD_CURVE_DATA installed.
!                         Dynamically allocated DAta arrays to allow arbitrary
!                         size plots.
! Altered:  31-May-2000 : Now can handel smaller X and Y ranges, because
!                          exponential format has been included in MON_NUM.
! Altered:  18-Nov-1999 : Maximum number of strings increased to 100.
!                           STR_COL now written to SAVE file.
!                           Only strings inside box are output (unless
!                           PG_LOC is negative.)
! Altered:  05-Mar-1999 : Multiple titles installed. MONBORD_V3 now utilized.
! Altered:  09-Mar-1999 : EXPCHAR_SCALE etc added so that TEXT and plot
!                            symbols have same size relative to the plot,
!                            irrespective of the plotting surface.
! Altered:  30-Jul-1997 : 'W'option installed to allow lines of different
!                            weight. SOme additional work may be needed.
! Altered:  07-Jul-1997 : Revised calls to NEW_GEN_IN_MULT_? installed.
! Altered:  06-May-1997 : Corrections made so that not all curves need
!                            markers.
! Altered:  04-Mar-1997 : Bug fixed in output file name procdure due to 
!                            directory conflicts.
! Finalized 07-Mar-1997 : PGPLOT version.
!                            Based on GRAMON (Mongo)
!
        INTEGER, PARAMETER :: MAXSTR=500
	INTEGER, PARAMETER :: MAXVEC=500
	INTEGER, PARAMETER :: MAXPEN=50            !Should be the same as MAX_PLTS
!
	INTEGER NDEC,GET_INDX_SP
	EXTERNAL SPACING,GET_INDX_SP
	REAL*4 SPACING
!
	CHARACTER(LEN=2) TYPE_CURVE(MAX_PLTS)
	REAL*8 VB_BASE(MAX_PLTS)
!
	LOGICAL DO_ERROR
	CHARACTER*5 LOG_AXIS
!
	REAL*4 XINC,XNUMST,YNUMST
	REAL*4 YINC
	INTEGER IDX,IXTICK
	INTEGER IDY,IYTICK
	REAL*4 XPAR(2),YPAR(2),XT(2),YT(2)
	REAL*4 XMIN,XMAX,YMIN,YMAX
	REAL*4 XPAR_SAV(2)
	REAL*8 YMIN_SAV,YMAX_SAV
!
	REAL*4 YPAR_R_AX(2)
	REAL*4 YNUMST_R_AX
	REAL*4 YINC_R_AX
	LOGICAL NORMAL_R_Y_AXIS
	LOGICAL DO_BORDER
	INTEGER IDY_R_AX,IYTICK_R_AX
	CHARACTER*1 WHICH_Y_AX(MAX_PLTS)
	CHARACTER*80 YLABEL_R_AX
!
	INTEGER IFILL_PLT1(10)
	INTEGER IFILL_PLT2(10)
	INTEGER IFILL_CLR(10)
	INTEGER NFILL
	LOGICAL FILL
!
	INTEGER, PARAMETER :: N_TITLE=10
	CHARACTER*80  XLABEL,YLABEL,TITLE(N_TITLE)
	CHARACTER*(*) XLAB,YLAB,TITL,PASSED_OPT
        CHARACTER*200 FILNAME
        CHARACTER*80 ID_FILNAME
        CHARACTER*80, SAVE :: PLT_ST_FILENAME
	CHARACTER*80 WK_STR,OPTION
	CHARACTER*80 TMP_STR
	CHARACTER*6 TO_TEK
	CHARACTER*2 TO_VT
	CHARACTER*50 PLT_ID,RD_PLT_ID,PLT_ID_SAV
!
! Vector arrays
!
	REAL*4 LINEXST(MAXVEC),LINEYST(MAXVEC)
	REAL*4 LINEXEND(MAXVEC),LINEYEND(MAXVEC)
	LOGICAL VEC,FLAGLINE(MAXVEC),INIT
	LOGICAL QUERYFLAG
	INTEGER VECPEN(MAXVEC)
!
	INTEGER N_LINE_IDS
	LOGICAL OBSERVED_WAVE(5000)
	CHARACTER*10 LINE_ID(5000)
	REAL*4 ID_WAVE(5000)
	REAL*4 ID_WAVE_OFF(5000)
	REAL*4 ID_Y_OFF(5000)
	REAL*4 TAU(5000)
	REAL*4 TAU_CUT
	REAL*4 ID_ORIENT
	REAL*4 ID_SCL
	REAL*4 ID_VEC_BEG
	REAL*4 ID_VEC_END
	REAL*4 ID_EXPCHAR
	REAL*4 CUT_ACCURACY
	INTEGER ID_LOC
	INTEGER ID_LOC_PG
	CHARACTER*10 OMIT_ID(10)
	CHARACTER*10 INC_ID(10)
	INTEGER N_OMIT_ID
	INTEGER N_INC_ID
	LOGICAL AIR_WAVELENGTHS
	LOGICAL WR_ID(5000)
	LOGICAL DRAW_GAUSS_HARD
	EXTERNAL LAM_AIR
	REAL*8 LAM_AIR
	REAL*8 DP_T1
!
! String arrays (not labels or titles)
!
	REAL*4 ORIENTATION(MAXSTR),XSTR(MAXSTR),YSTR(MAXSTR)
	REAL*4 STR_EXP(MAXSTR),LOC_PG(MAXSTR)
	REAL*4 XSTRPOS(MAXSTR),YSTRPOS(MAXSTR)
	CHARACTER*80 STRING(MAXSTR)
	INTEGER ISTR
	INTEGER LOC(MAXSTR)
	INTEGER STR_COL(MAXVEC)
	LOGICAL FLAGSTR(MAXSTR),STR
!
	REAL*4 DXST,DXEND,DYST,DYEND
!
! Pen type, Device ID's
!
! LINE_STYLE   : Type of curve (solid, dsahed etc drawn through data)
! MARKER_STYLE : Type of marker drawn at each datum. If negative, no
!                   curve is drawn between the markers.
!
        INTEGER LINE_STYLE(MAX_PLTS)
        INTEGER LINE_WGT(MAX_PLTS)
	INTEGER MARKER_STYLE(MAX_PLTS)
	LOGICAL HARD,TITONRHS,FIRST,FSTOPEN,DASH,MARK
	LOGICAL FILE_IS_OPEN
	LOGICAL INITIALIZE_ON_EXIT
	LOGICAL, SAVE :: REVERSE_PLOTTING_ORDER=.FALSE.
!
! E, cursor, and continuum parameters.
!
	REAL*4, ALLOCATABLE :: CONT(:)
	REAL*4 EW,CENTROID
	REAL*4 XCUR(2),YCUR(2),SLOPE
	INTEGER PLOT_ID,CURSERR
	INTEGER L_CHAN(2)
	INTEGER, SAVE :: LU_EW=30
	INTEGER, SAVE :: LU_NORM=31
	LOGICAL, SAVE :: FIRST_EW=.TRUE.
	CHARACTER*1 CURSVAL
	LOGICAL CONTINUUM_DEFINED
!
! Functions
!
	LOGICAL END_CURS
	INTEGER PGBEG, PGCURS
	REAL*8 SPEED_OF_LIGHT
	REAL*8 LAM_VAC
!
	REAL*4,    PARAMETER :: RONE=1.0
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: IFOUR=4
!
! Parameters for indicating the size of the plot.
!
	REAL*4 XCM,ASR,TEMPASR,DASR
!
! CENTRAL_LAM must be REAL*8 as LAM_VAC is REAL*8 function.
!
	REAL*8 CENTRAL_LAM
	REAL*8 OLD_CENTRAL_LAM
	REAL*8 C_KMS
	REAL*8 C_VAL
	LOGICAL AIR_LAM
	CHARACTER(LEN=5), SAVE :: VEL_UNIT='km/s'
!
! For XAR and YAR arithmetic options.
!
	CHARACTER*3 XAR_OPERATION,YAR_OPERATION,VAR_OPERATION
	REAL*4 YAR_VAL,XAR_VAL
	INTEGER XAR_PLT,YAR_PLT
	INTEGER VAR_PLT1,VAR_PLT2,VAR_PLT3
!
! Miscellaneous
!
	CHARACTER*4 ANS
	INTEGER LENGTH
	INTEGER LAST_DP
	REAL*4 V1,IT
	REAL*4 EXPCHAR,TICK_FAC,EXPMARK
	REAL*4 EXPCHAR_SCALE,TICK_FAC_SCALE,EXPMARK_SCALE
	REAL*4 XCHAR_SIZE,YCHAR_SIZE
	REAL*4 MARGINX(2),MARGINY(2)
	REAL*4 MEAN,SIGMA
	REAL*4 T1,T2,T3,T4
 	REAL*4 XVAL,YVAL
	REAL*4 XVAL_SAV,YVAL_SAV
	REAL*4, ALLOCATABLE :: TA(:)
!
	INTEGER, SAVE :: PEN_OFFSET
	INTEGER BEG
	INTEGER Q	!Used for pen color
!
	INTEGER, PARAMETER :: T_IN=5
	INTEGER, PARAMETER :: T_OUT=6
!
! Loop variables
	INTEGER I,J,K,L,CNT,IP,OP
	INTEGER IP_ST,IP_END,IP_INC
!
! Color variables
	REAL*4 RED(0:15),BLUE(0:15),GREEN(0:15)
	INTEGER PEN_COL(0:MAXPEN)
	REAL*4 RTEMP,BTEMP,GTEMP
!
! Variables to read in plot in column format
!
	REAL*4 TEMP_VAR(20)
	INTEGER COLUMN(2)
!
! Printer variables.
!
	CHARACTER(LEN=80), SAVE :: PRINTER='pgplot_1.ps/cps'
	CHARACTER(LEN=80), SAVE :: HARD_FILE
	CHARACTER(LEN=20), SAVE :: HARD_TYPE
	CHARACTER(LEN=80) CUR_HARD_FILE
	INTEGER, SAVE :: HARD_CNT
	LOGICAL, SAVE :: FIRST_HARD=.TRUE.
!
	INTEGER PGOPEN,PLT_LINE_WGT
	INTEGER ID
	SAVE ID
	REAL*4 TOTXMM,TOTYMM,SCALEFACY,SCALEFAC
	REAL*4 PRINTX1,PRINTX2,PRINTY1,PRINTY2
!
	LOGICAL LONG_PLOT
	REAL*4 LENGTH_OF_HC_PLOT
	REAL*4 LP_ASR
!
! Variables for options 'WP' (i.e. write plot) and 'RP' (i.e. read plot).
!
	INTEGER IST,IEND,N_REC_SIZE
!
	SAVE XCM,ASR
	SAVE EXPCHAR_SCALE,EXPMARK_SCALE,TICK_FAC_SCALE
	SAVE RED,BLUE,GREEN
	SAVE FSTOPEN,PEN_COL,DASH
	SAVE MARGINX,MARGINY
	SAVE PLT_LINE_WGT
	SAVE LINE_WGT 
	DATA FSTOPEN,DASH/.TRUE.,.FALSE./
	DATA PLT_LINE_WGT/1/
	DATA HARD_CNT/1/
	N_REC_SIZE=1000
!
	IF(NPLTS .EQ. 0)THEN
	  WRITE(T_OUT,*)'Error - No calls made to curve'
	  RETURN
	END IF
	N_LINE_IDS=0
	ID_SCL=1.05D0
	ID_VEC_BEG=1.04D0
	ID_VEC_END=1.01D0
	ID_EXPCHAR=1.0D0
	ID_FILNAME='LINE_ID'
	TAU_CUT=0.1D0
	N_OMIT_ID=0
	N_INC_ID=0
	DO_BORDER=.TRUE.
	DRAW_GAUSS_HARD=.FALSE.
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	LONG_PLOT=.FALSE.
	LENGTH_OF_HC_PLOT=200.0D0       !cm
	VB_BASE=-1000
	FILL=.FALSE.
!
	IF(NPLTS .GT. MAXPEN)THEN
	  WRITE(T_OUT,*)'Error n GRAMON_PLOT -- not enough pen loctions'
	  WRITE(T_OUT,*)'MAXPEN should be set to the same value as MAX_PLTS'
	  RETURN
	END IF
!
! Define character strings for switching between VT and TEK modes.
!
	TO_TEK(1:1)=CHAR(27)
	TO_TEK(2:6)='[?38h'
	TO_VT(1:1)=CHAR(27)
	TO_VT(2:2)=CHAR(3)
!
! Arithmetic options
!
	XAR_OPERATION='*'; XAR_VAL=1.0; XAR_PLT=1
	YAR_OPERATION='*'; YAR_VAL=1.0; YAR_PLT=1
	VAR_OPERATION='*'; VAR_PLT1=1;  VAR_PLT2=1
!
! Set the Aspect ratio of the plot to the default value. No longer ask
! for a new value, since can be set with N option. 
!
	IF(FSTOPEN)THEN
	  ASR=0.0			!Leave as is
	  EXPCHAR_SCALE=1.5   		!Was 1.01 to prevent crashes.
	  EXPMARK_SCALE=1.5   		!Symbol size on plots.
	  TICK_FAC_SCALE=1.5D0
	  PLT_LINE_WGT=1
          PLT_ST_FILENAME=' '
	  LINE_WGT(:)=1
	  PEN_OFFSET=1
	  IFILL_PLT1=0; IFILL_PLT2=0
	END IF
	TITLE(1:N_TITLE)=' '
	CALL GEN_ASCI_OPEN(LU_NORM,'NORM_FACTORS','UNKNOWN','APPEND',' ',IZERO,IOS)
!
	LOG_AXIS=' '
	XLABEL=XLAB
	YLABEL=YLAB
	YLABEL_R_AX=' '
	TITLE(1)=TITL
	STR=.FALSE.
	VEC=.FALSE.
	NORMAL_R_Y_AXIS=.TRUE.
	INITIALIZE_ON_EXIT=.TRUE.
	DO I=1,MAXSTR
	  FLAGSTR(I)=.FALSE.
	  FLAGLINE(I)=.FALSE.
	  STR_EXP(I)=1.0	   !Default (changed using SE option only)
	  STR_COL(I)=1		   !Default (changed using SE option only)
	END DO
	DO_ERROR=.TRUE.
	OPTION=PASSED_OPT
	CALL SET_CASE_UP(OPTION,1,0)
	CENTRAL_LAM=1548.20			!CIV
	OLD_CENTRAL_LAM=0.0D0
	CONTINUUM_DEFINED=.FALSE.
	AIR_WAVELENGTHS=.FALSE.
	PLOT_ID=1
!
	YPAR_R_AX(1)=0.0D0 ; YPAR_R_AX(2)=0.0D0
!
! Determine type of graph.
!
	IF(OPTION(1:4) .EQ. 'MARK')THEN
	  MARK=.TRUE.				!Connect plots
	ELSE IF(OPTION(1:4) .EQ. 'HIST')THEN
	  TYPE_CURVE(1:NPLTS)='H'
	  MARK=.FALSE.
	ELSE IF(OPTION(1:3) .EQ. 'ADJ')THEN  !Histogram, but adjacent points
	  TYPE_CURVE(1:NPLTS)='A'
	  MARK=.FALSE.
	ELSE
	  TYPE_CURVE(1:NPLTS)='L'
	  MARK=.FALSE.
	END IF
	WHICH_Y_AX(1:NPLTS)='L'
!
! Define default line representations (initially not dashed)
!
	DO I=1,MAX_PLTS
	  LINE_STYLE(I)=1
	  MARKER_STYLE(I)=MOD(I,5)
	END DO
	DO I=14,MAX_PLTS
	  LINE_STYLE(I)=2
	END DO
	DASH=.FALSE.
!
! Assign a color index to each pen. Keep previus assignments if they have been
! made. We avoid the white pen as it won't show on a hardcopy.
!
	IF (FSTOPEN) THEN
	  PEN_COL(0)=0
	  DO I=1,MAXPEN
	    PEN_COL(I)=MOD(I-1,14)+1
	  END DO
	END IF
!
! Assign a color index to each vecpen
!
	IF (FSTOPEN) THEN
	  DO I=1,10
	    VECPEN(I)=I+5
	  END DO
	  DO I=11,MAXVEC
	    VECPEN(I)=14
	  END DO
	END IF
!
! Define default type for "Data markers"
!
	DO I=1,NPLTS
	  MARKER_STYLE(I)=MOD(I,12)+1           !1 is not appropriate
	END DO
!
	ANS='P'
	FIRST=.TRUE.
!
! Get absica and ordinate limits.
!
	CALL GET_GRAMON_MIN_MAX(XMIN,XMAX,YMIN,YMAX,TYPE_CURVE,T_OUT) 
!
	XPAR(1)=XMIN
	XPAR(2)=XMAX
	YPAR(1)=YMIN
	YPAR(2)=YMAX
	XPAR_SAV(1)=XPAR(1)		!Indicate value limits evaluated for.
	XPAR_SAV(2)=XPAR(2)
	YMAX_SAV=YMAX; YMIN_SAV=YMIN
!
! Open user set workstation (default is set into the system).
! We only ask for work-station name if first call to GRAMON.
!
	WRITE(T_OUT,4) XMIN,XMAX
4	FORMAT(' Abisca limits  :',1P2E14.4)
	WRITE(T_OUT,11) YMIN,YMAX
11	FORMAT(' Ordinate limits:',1P2E14.4)
	WRITE(T_OUT,'(A)')' '
!
	IF (FSTOPEN) THEN
	  WRITE(T_OUT,*)'NOTE: /XWINDOW supports color'
	  WRITE(T_OUT,*)'    while /XTERM allows auto-switching of the cursor'
	  BEG = PGBEG(IZERO,'?',IONE,IONE)
	  CALL PGQID(ID)
! 
! Initialize viewport.
!
	  CALL PGASK(.FALSE.)                      !TURNS OFF PAGE PROMPT
	  CALL PGVSTD
	  CALL PGENV(XPAR(1),XPAR(2),YPAR(1),YPAR(2),IZERO,IZERO)
	  CALL PGQVP(0.,MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	  IF (MARGINX(1) .LT. (.9)) MARGINX(1)=MARGINX(1)+.05
!
! Increase default border size.
!
	  MARGINX(1)=0.15; MARGINX(2)=0.9
          MARGINY(1)=0.15; MARGINY(2)=0.9
!
! Set preferred defaults for pen colors.
!
	  CALL PGSCR(0,.7,.7,.7)     !set color representations
	  CALL PGSCR(1,0.0,0.0,0.0)
	  CALL PGSCR(2,1.0,0.0,0.0)
	  CALL PGSCR(3,0.0,0.0,1.0)
	  CALL PGSCR(4,0.0,0.6,0.3)
	  CALL PGSCR(5,.5,0.0,.7)
	  CALL PGSCR(13,0.0,1.0,1.0)
	  CALL PGSCR(15,.95,.95,.95)
          DO I=0,15                  !Get these + def color representations.
            CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
!
	  CALL DEFINE_MORE_PENS(MAXPEN)
	  DO I=26,30
	    IFILL_CLR(I-25)=I
	  END DO
!
	END IF
	FSTOPEN=.FALSE.
	HARD=.FALSE.
1000	ANS='P'
!
!	WRITE(T_OUT,*)TO_VT
!
	WRITE(T_OUT,*)' H,P,Z(ZN),A(F),L,CC,CP,N, M,W,D,C,B,'//
	1              ' VC,VF,VE, SC,SF,SE, E'
	CALL NEW_GEN_IN(ANS,'ANS')
	L=LEN_TRIM(ANS)
	CALL SET_CASE_UP(ANS,IONE,IZERO)
        IF(ANS .EQ. 'H')THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'P   - Plot - default'
          WRITE(T_OUT,*)'E   - EXIT from PLOT package'
	  WRITE(T_OUT,*)'Z   - Hardcopy (ZN=Asks for new hard device)'
          WRITE(T_OUT,*)'LP  - Allow a long hard copy postscript plot to be created (with CPS)'
	  WRITE(T_OUT,*)' '
          WRITE(T_OUT,*)'A   - Define Axis Parameters'
          WRITE(T_OUT,*)'2A  - Define labeling of right-hand axis'
          WRITE(T_OUT,*)'F   - Change default axis parameters'
          WRITE(T_OUT,*)'N   - Define aspect ratio, plot margins, character height etc'
	  WRITE(T_OUT,*)'L   - Modify Axis Labels and Titles'
	  WRITE(T_OUT,*)' '
 	  WRITE(T_OUT,*)'D   - Switch dashed lines on/off'
 	  WRITE(T_OUT,*)'DE  - Edit dashed lines one by one'
	  WRITE(T_OUT,*)'W   - Change thickness (weights) of curves'
 	  WRITE(T_OUT,*)'WE  - Edit line weights one by one'
	  WRITE(T_OUT,*)'M   - Switch marking data points on/off'
	  WRITE(T_OUT,*)'C   - Indicate how curves are to be connected (L,H,A,E,I,V,B)'
	  WRITE(T_OUT,*)'B   - Switch error bars on/off'
	  WRITE(T_OUT,*)'CC  - Change Color setting'
	  WRITE(T_OUT,*)'CP  - Change Pen (Color Index)'
	  WRITE(T_OUT,*)'RCP - Reset default color pens'
	  WRITE(T_OUT,*)'GP  - Set default for grey pens'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)'BRD  - Switch border potting on (def) or off'
	  WRITE(T_OUT,*)'FILL - Color in region between 2 curves'
	  WRITE(T_OUT,*)'OFF  - Set offsets when plotting multiple plots'
          WRITE(T_OUT,*)'RPO  - Plot curves in reverse order (switch): does not affect color'
          WRITE(T_OUT,*)'RID  - Read line ID''s'
          WRITE(T_OUT,*)'SID  - Change defaults for writing line ID''s'
	  WRITE(T_OUT,*)' '
!
	  WRITE(T_OUT,*)'LX  - Switch between LINEAR/LOG labeling of X axis'
	  WRITE(T_OUT,*)'LY  - Switch between LINEAR/LOG labeling of Y axis'
	  WRITE(T_OUT,*)'LXY - Switch between LINEAR/LOG baleling of X and Y axes'
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'VC  - Define line vectors using cursor'
	  WRITE(T_OUT,*)'VF  - Define line vectors using file input'
	  WRITE(T_OUT,*)'VE  - Online edit of vectors'
          WRITE(T_OUT,*)'SC  - Define strings using cursor'
          WRITE(T_OUT,*)'SF  - Define strings using file input'
          WRITE(T_OUT,*)'SE  - Online edit of strings'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)'LAM - List wavelengths of common lines'
	  WRITE(T_OUT,*)'VEL - Convert X axis to km/s space'
	  WRITE(T_OUT,*)'XAR - Simple X axis arithmetic'
	  WRITE(T_OUT,*)'YAR - Simple Y axis arithmetic'
	  WRITE(T_OUT,*)'VAR - Simple arithmetic on two plots'
	  WRITE(T_OUT,*)'NM  - Scale average to 1 or to another plot'
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'LOC - Use a cursor to read of (X,Y) coordinates on a plot'
	  WRITE(T_OUT,*)'DC  - Define a straight line continuum for EW'
	  WRITE(T_OUT,*)'EW  - Measure the EW of a single line'
	  WRITE(T_OUT,*)'GF  - Fit a (modfied) gaussian to an absorption or emission line'  
	  WRITE(T_OUT,*)'MGF - Fit multiple gaussian to an absorption or emission complex'  
	  WRITE(T_OUT,*)'EGF - Edit gauss-fit arameters'
	  WRITE(T_OUT,*)'DG  - Draw gauss-fit.'
	  WRITE(T_OUT,*)'WGF - Write gauss-fit parameters to a file'
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'RXY - Read plot from asci file'
	  WRITE(T_OUT,*)'WXY - Write plot to asci file'
	  WRITE(T_OUT,*)'SXY - Write section of data to terminal'
	  WRITE(T_OUT,*)'RP  - Read labeled plots from direct accecs file'
	  WRITE(T_OUT,*)'RPF - Similar to RP but asks for filename'
	  WRITE(T_OUT,*)'WP  - Write labeled plots to direct access file'
	  WRITE(T_OUT,*)'WPF - Similar to WP but asks for filename'
	  READ(T_IN,'(A)')ANS				!can use ANS here.
	  IF(ANS(1:1) .EQ. 'E' .OR. ANS(1:1) .EQ. 'e')GOTO 1000
!
	  WRITE(T_OUT,*)'OLF - Open LOG file'
	  WRITE(T_OUT,*)'OIF - Open input file'
	  WRITE(T_OUT,*)'CLF - Close LOG file'
	  WRITE(T_OUT,*)'CIF - Close inpit file'
	  WRITE(T_OUT,*)'CL  - Clear Graphics Screen'
!
	  WRITE(T_OUT,*)'NOI - Leave data intact on exit (switch)'
	  WRITE(T_OUT,*)'H   - Help'
	  WRITE(T_OUT,*)' '
          GOTO 1000
!
! Exit from Ploting package, saving STRING and VECTOR information.
! If the NOI option has been issued, the plots will still be 
! retained.
!
	ELSE IF(ANS .EQ. 'E')THEN
	  CLOSE(LU_NORM)
	  IF(INITIALIZE_ON_EXIT)THEN
	    DO IP=1,NPLTS
	      IF(ALLOCATED(CD(IP)%XVEC))DEALLOCATE(CD(IP)%XVEC)
	      IF(ALLOCATED(CD(IP)%DATA))DEALLOCATE(CD(IP)%DATA)
	      IF(ALLOCATED(CD(IP)%EMAX))DEALLOCATE(CD(IP)%EMAX)
	      IF(ALLOCATED(CD(IP)%EMIN))DEALLOCATE(CD(IP)%EMIN)
	      ERR(IP)=.FALSE.
	      NPTS(IP)=0
	    END DO
	    NPLTS=0
	  END IF
          IF(STR .OR. VEC)THEN
	    FILNAME='PGOUT'
	    CALL NEW_GEN_IN(FILNAME,'Filname for STRING/VEC storage (no ext)')
	    IF(FILNAME .EQ. ' ')FILNAME='JNK'	!Force output to jnk file as
	    L=LEN_TRIM(FILNAME)			!precaution.
	    IF(FILNAME .NE. ' ')THEN
	      IF(STR)THEN
	        FILNAME=FILNAME(:L)//'.STR'
                L=L+4
                OPEN(UNIT=33,FILE=FILNAME(:L),STATUS='UNKNOWN')
	          DO I=1,MAXSTR
                    IF(FLAGSTR(I))THEN
                      WRITE(33,17)LOC(I),XSTR(I),YSTR(I),ORIENTATION(I),
	1                   STR_EXP(I),STR_COL(I),TRIM(STRING(I))
17	              FORMAT(1X,I1,', ',4(ES12.4,','),I3,',',1X,1H',A,1H')
	           END IF
	         END DO
                CLOSE(UNIT=33)
	      END IF
	      IF(VEC)THEN
	        IF(STR)THEN
	          FILNAME(L-3:L)='.VEC'
	        ELSE
	          FILNAME=FILNAME(:L)//'.VEC'
                  L=L+4
	        END IF
                OPEN(UNIT=33,FILE=FILNAME(:L),STATUS='UNKNOWN')
	          DO I=1,MAXVEC
                    IF(FLAGLINE(I))THEN
	              WRITE(33,18)LINEXST(I),LINEYST(I),
	1                         LINEXEND(I),LINEYEND(I),VECPEN(I)
18	              FORMAT(1X,1P,4E18.8,3X,I3)
	            END IF
	          END DO
                CLOSE(UNIT=33)
	      END IF		!VEC end
            END IF		!Filename check
          END IF		!STR or VEC present
	  RETURN
	ELSE IF(ANS .EQ. 'OLF')THEN
	  FILNAME='LOG_FILE'
	  CALL NEW_GEN_IN(FILNAME,'Log file=')
	  CALL NEW_GEN_IN_OPTS('OPEN_LOG_FILE',FILNAME,20)
          GOTO 1000
	ELSE IF(ANS .EQ. 'CLF')THEN
	  CALL NEW_GEN_IN_OPTS('CLOSE_LOG_FILE','LOG_FILE',20)
          GOTO 1000
	ELSE IF(ANS .EQ. 'OIF')THEN
	  CALL NEW_GEN_IN(FILNAME,'Command input file')
	  CALL NEW_GEN_IN_OPTS('OPEN_FILE_INPUT',FILNAME,21)
          GOTO 1000
	ELSE IF(ANS .EQ. 'CIF')THEN
	  CALL NEW_GEN_IN_OPTS('CLOSE_FILE_INPUT',FILNAME,21)
          GOTO 1000
	END IF
C
	IF(FIRST)THEN
C Define xinc
	  XINC=SPACING(XPAR(1),XPAR(2))
	  YINC=SPACING(YPAR(1),YPAR(2))
	  IXTICK=2
	  IYTICK=2
C Compute Xmin,Xmax,Ymin,Ymax if not input
	  V1=1000
	  XPAR(1)=ABS(XINC)*(AINT(XMIN/ABS(XINC)+V1)-V1)
	  XPAR(2)=ABS(XINC)*(AINT(XMAX/ABS(XINC)-V1)+V1)
	  YPAR(1)=ABS(YINC)*(AINT(YMIN/ABS(YINC)+V1)-V1)
	  YPAR(2)=ABS(YINC)*(AINT(YMAX/ABS(YINC)-V1)+V1)
C
	  XNUMST=XPAR(1)
	  YNUMST=YPAR(1)
C
C Determine number of decimals
	  IDX=NDEC(XINC)
	  IT=NDEC(XPAR(1))
	  IF(IT .GT. IDX)IDX=IT
	  IDY=NDEC(YINC)
	  IT=NDEC(YPAR(1))
	  IF(IT .GT. IDY)IDY=IT
	  FIRST=.FALSE.
	END IF
C
	IF(ANS .EQ. 'A' .OR. ANS .EQ. 'F')THEN
!
! Get absica and ordinate limits.
!
	  CALL GET_GRAMON_MIN_MAX(XMIN,XMAX,YMIN,YMAX,TYPE_CURVE,T_OUT) 
!
	  WRITE(T_OUT,4) XMIN,XMAX
	  CALL NEW_GEN_IN(XPAR,I,ITWO,'XST,XEND')
C
C Look for ordinate limits in X range. This will now generate a plot with 
C roughly the correct scaling even though the oridinates values may be vastly
C different outside the plot window. We only use those plots that will be
C displayed.
C
	    WRITE(6,*)'Getting Y limits'
	    YMIN=1.0E+32
	    YMAX=-1.0E+32
	    DO IP=1,NPLTS
	      IF(TYPE_CURVE(IP) .EQ. 'LG')THEN
	        T1=1.0E+32;T2=-1.0E+32
	        DO I=1,NPTS(IP)
	          IF( (CD(IP)%XVEC(I) .GE. XPAR(1) .AND. CD(IP)%XVEC(I) .LE. XPAR(2)) .OR.
	1           (CD(IP)%XVEC(I) .GE. XPAR(2) .AND. CD(IP)%XVEC(I) .LE. XPAR(1)) )THEN
	             T2=MAX(T2,ABS(CD(IP)%DATA(I)))
	             T1=MIN(T1,ABS(CD(IP)%DATA(I)))
	          END IF
	        END DO
	        YMAX=MAX(LOG10(T2),YMAX)
	        IF(T1 .EQ. 0.0D0)THEN
	          YMIN=-38
	        ELSE
	          YMIN=MIN(YMIN,LOG10(T1))
	        END IF
	      ELSE IF(TYPE_CURVE(IP) .NE. 'I')THEN
	        DO I=1,NPTS(IP)
	          IF( (CD(IP)%XVEC(I) .GE. XPAR(1) .AND. CD(IP)%XVEC(I) .LE. XPAR(2)) .OR.
	1           (CD(IP)%XVEC(I) .GE. XPAR(2) .AND. CD(IP)%XVEC(I) .LE. XPAR(1)) )THEN
	             YMAX=MAX(YMAX,CD(IP)%DATA(I))
	             YMIN=MIN(YMIN,CD(IP)%DATA(I))
	          END IF
	        END DO
	      END IF
	    END DO
	    IF(YMAX .EQ. 0 .AND. YMIN .EQ. 0)THEN
	      WRITE(T_OUT,*)'Poor YMIN & YMAX values'
	      WRITE(T_OUT,*)'YMIN=',YMIN,'  YMAX=',YMAX
	      YMIN=0.0
	      YMAX=1.0
	    END IF
!
	  WRITE(6,*)YMIN,YMAX
	  IF(XPAR(2) .NE. XPAR_SAV(2) .OR. XPAR(1) .NE. XPAR_SAV(1) .OR.
	1    YMAX .NE. YMAX_SAV .OR. YMIN .NE. YMIN_SAV)THEN
	    YINC=SPACING(YMIN,YMAX)
	    IYTICK=2
	    V1=1000
	    YPAR(1)=ABS(YINC)*(AINT(YMIN/ABS(YINC)+V1)-V1)
	    YPAR(2)=ABS(YINC)*(AINT(YMAX/ABS(YINC)-V1)+V1)
	    XPAR_SAV(1)=XPAR(1)		!Indicate value limits evaluated for.
	    XPAR_SAV(2)=XPAR(2)
	    YMIN_SAV=YMIN; YMAX_SAV=YMAX
	  END IF
C
	  WRITE(T_OUT,11)YMIN,YMAX
	  CALL NEW_GEN_IN(YPAR,I,ITWO,'YST,YEND')
C
C Check that limits inserted are not absurd.
C
	  IF(ABS(YPAR(2)-YPAR(1))/MAX(ABS(YPAR(2)),ABS(YPAR(1))) .LT. 1.0E-08)THEN
	    YPAR(1)=0
	    YPAR(2)=1.00
	    WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	  END IF
	  IF(ABS(XPAR(2)-XPAR(1))/MAX(ABS(XPAR(2)),ABS(XPAR(1)))  .LT. 1.0E-08)THEN
	    XPAR(1)=0
	    XPAR(2)=1.00
	    WRITE(T_OUT,*)'Invalid X limits - setting default values'
	  END IF
C
	  XNUMST=XPAR(1)
	  YNUMST=YPAR(1)
C
C Set defaults for XINC and YINC. IXTICK-1 is the number of tick
C marks between numbered ticks.
C
	  XINC=SPACING(XPAR(1),XPAR(2))
	  IXTICK=2
	  YINC=SPACING(YPAR(1),YPAR(2))
	  IYTICK=2
	  WRITE(T_OUT,*)'Indicate spacing between numberd tickmarks'
	  CALL NEW_GEN_IN(XINC,'For X tickmarks')
	  CALL NEW_GEN_IN(YINC,'For Y tickmarks')
C
C Determine number of decimals
C
	  IDX=NDEC(XINC)
	  IT=NDEC(XPAR(1))
	  IF(IT .GT. IDX)IDX=IT
	  IDY=NDEC(YINC)
	  IT=NDEC(YPAR(1))
	  IF(IT .GT. IDY)IDY=IT
C
	  IF(ANS .EQ. 'A')GOTO 1000
C
C This section is for additional axis fiddling. Can change number
C of digits after decimal point, and numbers at which axis labeling
C begins. These parameters are restored to their original values 
C everytime 'A' option is called.
C
	  CALL NEW_GEN_IN(IXTICK,'X minor divisions')
	  IXTICK=MAX(IXTICK,1)
	  CALL NEW_GEN_IN(IYTICK,'Y minor divisions')
	  IYTICK=MAX(IYTICK,1)
	  CALL NEW_GEN_IN(IDX,'# of X Dec digits')
	  CALL NEW_GEN_IN(IDY,'# of Y Dec digits')
	  CALL NEW_GEN_IN(XNUMST,'Number beginning X Axis')
	  CALL NEW_GEN_IN(YNUMST,'Number beginning Y Axis')
C
	  FIRST=.FALSE.
	  GOTO 1000
C
	ELSE IF(ANS .EQ. '2A')THEN
C
C Option to allow a left hand and a right hand axis. The default is a single
C left hand axis. As this is primary for pretty plots, the options are less
C powerful than for the left axis.
C
	  NORMAL_R_Y_AXIS=.TRUE.
	  WRITE(T_OUT,'(A)')' Set Y_AX to R to use with right axis'
	  DO IP=1,NPLTS
	    WRITE(T_OUT,'(I3,'' : '',$)')IP
	    CALL NEW_GEN_IN(WHICH_Y_AX(IP),'Y_AX=')
	    CALL SET_CASE_UP(WHICH_Y_AX(IP),1,1)
	    IF(WHICH_Y_AX(IP) .EQ. 'R')NORMAL_R_Y_AXIS=.FALSE.
	  END DO
C 
	  T1=1.0D+30
	  T2=-1.0D+30
	  DO IP=1,NPLTS
	    IF(WHICH_Y_AX(IP) .EQ. 'R')THEN
	      T3=MINVAL(CD(IP)%DATA)
	      T1=MIN(T1,T3)
	      T3=MAXVAL(CD(IP)%DATA)
	      T2=MAX(T2,T3)
	    END IF
	  END DO
C
	  IF(YPAR_R_AX(1) .EQ. 0 .AND. YPAR_R_AX(2) .EQ. 0)THEN
	    YINC_R_AX=SPACING(T1,T2)
	    IYTICK_R_AX=2
	    V1=1000
	    YPAR_R_AX(1)=ABS(YINC_R_AX)*(AINT(T1/ABS(YINC_R_AX)+V1)-V1)
	    YPAR_R_AX(2)=ABS(YINC_R_AX)*(AINT(T2/ABS(YINC_R_AX)-V1)+V1)
	  END IF
	  WRITE(T_OUT,11)T1,T2
	  CALL NEW_GEN_IN(YPAR_R_AX,I,ITWO,'YST,YEND for RHS axis')
C
C Check that limits inserted are not absurd.
C
	  IF(ABS(YPAR_R_AX(2)-YPAR_R_AX(1))/
	1        MAX(ABS(YPAR_R_AX(2)),ABS(YPAR_R_AX(1)))  .LT. 1.0E-08)THEN
	    YPAR_R_AX(1)=0
	    YPAR_R_AX(2)=1.00
	    WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	  END IF
C
	  YINC_R_AX=SPACING(YPAR_R_AX(1),YPAR_R_AX(2))
	  IYTICK_R_AX=2
	  WRITE(T_OUT,*)'Indicate spacing between numberd tickmarks'
	  CALL NEW_GEN_IN(YINC_R_AX,'For Y tickmarks')
	  CALL NEW_GEN_IN(IYTICK_R_AX,'Y minor divisions')
	  IYTICK=MAX(IYTICK_R_AX,1)
C
	  IDY_R_AX=NDEC(YINC_R_AX)
	  IT=NDEC(YPAR_R_AX(1))
	  IF(IT .GT. IDY_R_AX)IDY_R_AX=IT
	  CALL NEW_GEN_IN(IDY_R_AX,'# of Y Dec digits')
	  YNUMST_R_AX=YPAR_R_AX(1)
	  CALL NEW_GEN_IN(YNUMST_R_AX,'Number beginning Y Axis')
	  CALL NEW_GEN_IN(YLABEL_R_AX,'RHS AXIS LABEL')
C
	  GOTO 1000

	ELSE IF(ANS .EQ. 'L')THEN
	  CALL NEW_GEN_IN(XLABEL,'XLAB')
	  CALL NEW_GEN_IN(YLABEL,'YLAB')
	  DO J=1,N_TITLE
	    CALL NEW_GEN_IN(TITLE(J),'TITLE')
	    IF(TITLE(J) .EQ. ' ')THEN
	      TITLE(J+1:N_TITLE)=' '
	      EXIT
	    END IF
	  END DO
	  CALL NEW_GEN_IN(TITONRHS,'TITONRHS')
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'GL')THEN
	  CALL NEW_GEN_IN(FILNAME,'FILE')
	  CALL GET_TITLES(FILNAME,TITLE,N_TITLE)
	  CALL NEW_GEN_IN(TITONRHS,'TITONRHS')
	  GOTO 1000
!
! Switch to prevent inialization of data curves on exit from routine.
! By default, the plot count is set to zero on exit, and the data is lost on
! a new entry to the plotting package.
!
	ELSE IF(ANS .EQ. 'NOI')THEN
	  INITIALIZE_ON_EXIT=.NOT. INITIALIZE_ON_EXIT
	  IF(INITIALIZE_ON_EXIT)THEN
	    WRITE(T_OUT,*)'Warning: plot data will be initialized on exit.'
	  ELSE
	    WRITE(T_OUT,*)'Plot data will be retained on exit.'
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'CC')THEN
	  CALL CHANGE_COLOR(RED,BLUE,GREEN)
	  GOTO 1000
	ELSE IF(ANS .EQ. 'CP')THEN
	  CALL CHANGE_PEN(PEN_COL,MAXPEN,NPLTS)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SFP')THEN
	  CALL NEW_GEN_IN(PEN_OFFSET,'First colored pen: 0 [black], 1 [red]')
!
! Reset default pen color. Useful after GP option.
!
	ELSE IF(ANS .EQ. 'RCP')THEN
	  CALL PGSCR(0,.7,.7,.7)     !set color representations
	  CALL PGSCR(1,0.0,0.0,0.0)
	  CALL PGSCR(2,1.0,0.0,0.0)
	  CALL PGSCR(3,0.0,0.0,1.0)
	  CALL PGSCR(4,0.0,0.6,0.3)
	  CALL PGSCR(5,.5,0.0,.7)
	  CALL PGSCR(13,0.0,1.0,1.0)
	  CALL PGSCR(15,.95,.95,.95)
          DO I=0,15                  !Get these + def color representations.
            CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
	ELSE IF(ANS .EQ. 'GP')THEN
!
! Set preferred defaults for grey pens.
!
	  CALL PGSCR(0,1.0,1.0,1.0)     !set color representations
	  CALL PGSCR(1,0.0,0.0,0.0)
	  CALL PGSCR(2,0.0,0.0,0.0)
	  CALL PGSCR(3,0.4,0.4,0.4)
	  CALL PGSCR(4,0.8,0.8,0.8)
	  CALL PGSCR(5,0.2,0.2,0.2)
	  CALL PGSCR(6,0.6,0.6,0.6)
          DO I=0,15                  !Get these + def color representations.
            CALL PGQCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
	  GOTO 1000
	ELSE IF(ANS .EQ. 'SCR')THEN
	  T1=100; T2=0.2
	  CALL PGSCRL(T1,T2)
!
!Option to adjust the shape of the box, and the tick marks tec.
!
	ELSE IF(ANS .EQ. 'N')THEN
!
	  CALL NEW_GEN_IN(EXPCHAR_SCALE,'Expand Char.')
	  CALL NEW_GEN_IN(TICK_FAC_SCALE,'Expand TICK')
!
! Define the Aspect ratio of the plot.
!
	  WRITE(T_OUT,*)'Y/X > 0, 0(Device default) , X/Y < 0'
	  CALL NEW_GEN_IN(ASR,'Aspect Ratio')
!
	  IF(LONG_PLOT)THEN
	    WRITE(T_OUT,*)RED_PEN
	    WRITE(T_OUT,*)' As a long plot, 0.05 and 0.95 may be better for left/right margin'
	    WRITE(T_OUT,*)DEF_PEN
	  END IF
700	  CALL NEW_GEN_IN(MARGINX,I,ITWO,'Left and Right Margins (0:1)')
	  IF((MARGINX(1) .GT. 1) .OR. (MARGINX(1) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Left-Margin is incorrect, try again.'
	    GOTO 700
	  END IF
	  IF((MARGINX(2) .GT. 1) .OR. (MARGINX(2).LT. 0)) THEN
	    WRITE(T_OUT,*)'Right-Margin is incorrect, try again.'
	    GOTO 700
	  END IF
!
710	  CALL NEW_GEN_IN(MARGINY,I,ITWO,'Bottom and Top Margins (0:1)')
	  IF((MARGINY(1) .GT. 1.0) .OR. (MARGINY(1) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Bottom-Margin is incorrect, try again.'
	    GOTO 710
	  END IF 
	  IF((MARGINY(2) .GT. 1.0) .OR. (MARGINY(2) .LT. 0)) THEN
	    WRITE(T_OUT,*)'Top-Margin is incorrect, try again.'
	    GOTO 710
	  END IF 
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'W')THEN
	    CALL NEW_GEN_IN(LINE_WGT,I,NPLTS,'Line Weight:: 1,2 etc')
	    IF(LINE_WGT(1) .LT. 0)LINE_WGT(1:NPLTS)=ABS(LINE_WGT(1))
	ELSE IF( ANS .EQ. 'WE')THEN
	  DO IP=1,NPLTS			!Edit individually
	    CALL NEW_GEN_IN(LINE_WGT(IP),'Line weights (1,...,5)')
	  END DO
	ELSE IF( ANS .EQ. 'M')THEN
	  IF(MARK)THEN
            WRITE(T_OUT,*)'Will not mark data points'
	    MARK=.FALSE.
	  ELSE
	    MARK=.TRUE.
            WRITE(T_OUT,*)'Wil now mark data points'
	    CALL NEW_GEN_IN(EXPMARK_SCALE,'Expand Symbol Size')
	    CALL NEW_GEN_IN(MARKER_STYLE,I,NPLTS,'Marker style {+/-}1,...,31)')
	  END IF
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'C')THEN
	  DO IP=1,NPLTS
850	    WRITE(T_OUT,'(I3,'' : '',$)')IP
	    CALL NEW_GEN_IN(TYPE_CURVE(IP),'TC=')
	    CALL SET_CASE_UP(TYPE_CURVE(IP),IONE,IZERO)
	    IF( TYPE_CURVE(IP) .NE. 'L' .AND.              !Normal line
	1       TYPE_CURVE(IP) .NE. 'E' .AND.              !non-monotonic
	1       TYPE_CURVE(IP) .NE. 'EC' .AND.              !non-monotonic
	1       TYPE_CURVE(IP) .NE. 'B' .AND.              !Broken
	1       TYPE_CURVE(IP) .NE. 'I' .AND.              !Invisible
	1       TYPE_CURVE(IP) .NE. 'V' .AND.              !Verticle lines
	1       TYPE_CURVE(IP) .NE. 'VB' .AND.             !Verticle lines
	1       TYPE_CURVE(IP) .NE. 'A' .AND.              !Hist - X vert
	1       TYPE_CURVE(IP) .NE. 'H' .AND.              !Histogram
	1       TYPE_CURVE(IP) .NE. 'LG' )THEN             !Log(ABS)
!
	        WRITE(T_OUT,*)'Invalid connection specifier. Specifiers are:'
	        WRITE(T_OUT,*)'L --- Normal line'
	        WRITE(T_OUT,*)'E --- non-monotonic'
	        WRITE(T_OUT,*)'B --- Broken'
	        WRITE(T_OUT,*)'I --- Invisible'
	        WRITE(T_OUT,*)'V --- Verticle lines'
	        WRITE(T_OUT,*)'A --- Hist - X vert'
	        WRITE(T_OUT,*)'H --- Histogram'
	        WRITE(T_OUT,*)'LG --- Marked ABS(LG)'
	      GOTO 850
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF( ANS .EQ. 'RPO')THEN
	  REVERSE_PLOTTING_ORDER= .NOT. REVERSE_PLOTTING_ORDER
	  IF(REVERSE_PLOTTING_ORDER)THEN
	    WRITE(6,*)'Order of ploting will be reversed'
	  ELSE
	    WRITE(6,*)'Normal order of ploting will be resumed'
	  END IF	    
 	ELSE IF( ANS .EQ. 'D')THEN
	  IF(DASH)THEN
            WRITE(T_OUT,*)'Now all solid line plots '
	    LINE_STYLE(1:NPLTS)=1
	    DASH=.FALSE.
	  ELSE
	    WRITE(T_OUT,*)'Input dashed pens for ',NPLTS
	    DO IP=1,NPLTS,5
	      DO J=1,5
	        IF(IP .LE. NPLTS)LINE_STYLE(IP+J-1)=J
	      END DO
	    END DO
	    CALL NEW_GEN_IN(LINE_STYLE,I,NPLTS,
	1      'Dashed pen types (1,...,5)')
	    DASH=.TRUE.
	  END IF
	  GOTO 1000
!
 	ELSE IF(ANS .EQ. 'BRD')THEN
	  IF(DO_BORDER)THEN
            WRITE(T_OUT,*)'Switching off border'
	    DO_BORDER=.FALSE.
	  ELSE
	    WRITE(T_OUT,*)'Switching on border'
	    DO_BORDER=.TRUE.
	  END IF
	  GOTO 1000
!
! Edit each dashed pen separately.
C
 	ELSE IF( ANS .EQ. 'DE')THEN
	  DO IP=1,NPLTS
	    CALL NEW_GEN_IN(LINE_STYLE(IP),'Dashed pen types (1,...,5)')
	  END DO
	  DASH=.TRUE.
	  GOTO 1000
C
	ELSE IF( ANS .EQ. 'B')THEN
	  IF(DO_ERROR)THEN
            WRITE(T_OUT,*)'Swithcing off error bar plotting'
	    DO_ERROR=.FALSE.
	  ELSE
	    DO_ERROR=.TRUE.
	  END IF
	  GOTO 1000
C
	ELSE IF(ANS .EQ. 'LXY')THEN
	  IF(LOG_AXIS(1:3) .EQ. 'LOG')THEN
            WRITE(T_OUT,*)'Swithcing off ALL LOG axes'
	    LOG_AXIS=' '
	  ELSE
	    LOG_AXIS='LOGXY'
	  END IF
	  GOTO 1000
	ELSE IF( ANS .EQ. 'LX')THEN
	  IF(LOG_AXIS(1:5) .EQ. 'LOGX')THEN
            WRITE(T_OUT,*)'Swithcing off X LOG axes'
	    LOG_AXIS=' '
	  ELSE IF(LOG_AXIS(1:5) .EQ. 'LOGXY')THEN
            WRITE(T_OUT,*)'Swithcing off X LOG axes'
	    LOG_AXIS='LOGY'
	  ELSE IF(LOG_AXIS(1:4) .EQ. 'LOGY')THEN
	    LOG_AXIS='LOGXY'
	  ELSE
	    LOG_AXIS='LOGX'
	  END IF
	  GOTO 1000
	ELSE IF( ANS .EQ. 'LY')THEN
	  IF(LOG_AXIS(1:4) .EQ. 'LOGY')THEN
            WRITE(T_OUT,*)'Swithcing off Y LOG axes'
	    LOG_AXIS=' '
	  ELSE IF(LOG_AXIS(1:5) .EQ. 'LOGXY')THEN
            WRITE(T_OUT,*)'Swithcing off Y LOG axes'
	    LOG_AXIS='LOGX'
	  ELSE IF(LOG_AXIS(1:4) .EQ. 'LOGX')THEN
	    LOG_AXIS='LOGXY'
	  ELSE
	    LOG_AXIS='LOGY'
	  END IF
	  GOTO 1000
!
!
! String handling section.
!
	ELSE IF (ANS .EQ. 'SF')THEN
	  INIT=.FALSE.          !Strings automatically initialized first time.
	  IF(STR)CALL NEW_GEN_IN(INIT,'Initialize string list ?')
	  IF(INIT)THEN
	    DO I=1,MAXSTR
	      FLAGSTR(I)=.FALSE.
	      STR_EXP(I)=1.0	  !Default (changed at editing stage only)
	    END DO
            STR=.FALSE.
	  END IF
!
	  L=INDEX(FILNAME,'.STR')
	  IF(L .NE. 0)FILNAME(L:L+3)='.STR'
	  ISTR=0
923	  CALL NEW_GEN_IN(FILNAME,'Input string file to be read')
	  CALL SET_CASE_UP(FILNAME,1,0)
	  L=LEN_TRIM(FILNAME)
	  IF(L .NE. 0)THEN
	    IF(INDEX(FILNAME,'.') .EQ. 0)FILNAME=FILNAME(:L)//'.STR'
            OPEN(UNIT=33,FILE=TRIM(FILNAME),STATUS='OLD',ERR=923)
	    L=0
            DO ISTR=1,MAXSTR
	      IF( .NOT. FLAGSTR(ISTR))THEN
	        READ(33,*,END=930,ERR=925)LOC(ISTR),XSTR(ISTR),YSTR(ISTR),
	1                   ORIENTATION(ISTR),STR_EXP(ISTR),STR_COL(ISTR),
	1                   STRING(ISTR)
	        FLAGSTR(ISTR)=.TRUE.
	        L=L+1
                GOTO 928
925             WRITE(T_OUT,*)' Error reading string',ISTR
928	        CONTINUE
	      END IF
            END DO
930	    CLOSE(UNIT=33)
	    IF(L .NE. 0)THEN
	      WRITE(T_OUT,932)L,'STRINGS read in from file'
              STR=.TRUE.
	    END IF
          END IF
932	  FORMAT(1X,I3,2X,(A))
!
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SE')THEN
!
! First output strings to user, before doing online editing.
!
	  DO I=1,MAXSTR
	    IF(FLAGSTR(I))THEN
	      WRITE(T_OUT,955)I,TRIM(STRING(I))
955	      FORMAT(1X,I2,3X,A)
	    END IF
	  END DO
!
! Allow online editing of strings. NEW_GEN_IN is used so user can
! see current value. 
!
	  DO WHILE (0 .EQ. 0)
	    ISTR=0			!Default terminates input
	    CALL NEW_GEN_IN(ISTR,'String number - 0 exits')
	    IF(ISTR .LE. 0 .OR. ISTR .GT. MAXSTR)GOTO 1000
	    WRITE(T_OUT,955)ISTR,TRIM(STRING(ISTR))
	    CALL NEW_GEN_IN(FLAGSTR(ISTR),'FLAG')
	    IF(FLAGSTR(ISTR))THEN
	      CALL NEW_GEN_IN(LOC(ISTR),'LOC')
	      CALL NEW_GEN_IN(XSTR(ISTR),'XPOS')
	      CALL NEW_GEN_IN(YSTR(ISTR),'YPOS')
	      CALL NEW_GEN_IN(ORIENTATION(ISTR),'ORI')
	      CALL NEW_GEN_IN(STR_EXP(ISTR),'EXP')
	      CALL NEW_GEN_IN(STR_COL(ISTR),'COL')
	      CALL NEW_GEN_IN(STRING(ISTR),'STRING')
	    END IF
	  END DO
	  GOTO 1000
	  STR=.TRUE.
!
	ELSE IF (ANS .EQ. 'SC')THEN
	  INIT=.FALSE.          !Strings automatically initialized first time.
	  IF(STR)CALL NEW_GEN_IN(INIT,'Initialize string list ?')
	  STR=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXSTR
	      FLAGSTR(I)=.FALSE.
	    END DO
	  END IF
	      WRITE(T_OUT,820)
 820	      FORMAT(1X,'LOC: 1-9, 1=lef bot, 2=cent bot, 3=rit bot')
	      WRITE(T_OUT,821)
 821	      FORMAT(1X,'          4=lef mid, 5=cent mid, 6=rit mid')
	      WRITE(T_OUT,822)
 822	      FORMAT(1X,'          7=lef top, 8=cent top, 9=rit top')
	      WRITE(T_OUT,840)
 840	      FORMAT(1X,'ORI= AN ANGLE FROM 0.0 TO 360')
	  ISTR=1
	  DO WHILE(ISTR .LE. MAXSTR)
	    IF(.NOT. FLAGSTR(ISTR))THEN
	      CURSERR = PGCURS(XSTR(ISTR),YSTR(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
830	      WRITE(T_OUT,810,ADVANCE='NO')
810	      FORMAT(1X,'Desc(LOC,ORI,''STRING''): ')
	      READ(T_IN,'(A)',ERR=830)WK_STR
	      READ(WK_STR,*,ERR=830)LOC(ISTR)
	      IF(LOC(ISTR) .EQ. 0)GOTO 1000
	      READ(WK_STR,*,ERR=830)LOC(ISTR),ORIENTATION(ISTR),STRING(ISTR)
	      FLAGSTR(ISTR)=.TRUE.
	    END IF
	    ISTR=ISTR+1
	  END DO
	  WRITE(T_OUT,*)'Maximum number of strings is',MAXSTR
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'LOC')THEN
	  CALL PGQCS(4,T1,T3)			!T3=Character Height
	  T3=T3*1.5
	  CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	  IF(END_CURS(CURSVAL))GOTO 70
	  T1=XPAR(1)-(XPAR(2)-XPAR(1))/15.0
	  T2=YPAR(2)+(YPAR(2)-YPAR(1))/15.0
	  CALL MON_NUM(T1,T2,XVAL,1,IDX+2)
	  CALL MON_NUM(T1,T2-T3,YVAL,1,IDY+2)
!
	  DO WHILE(1 .EQ. 1)
	    XVAL_SAV=XVAL
	    YVAL_SAV=YVAL
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    CALL PGSCI(0)     ! Background
	    CALL MON_NUM(T1,T2,XVAL_SAV,1,IDX+2)
	    CALL MON_NUM(T1,T2-T3,YVAL_SAV,1,IDY+2)
	    CALL PGSCI(1)     ! Axes
	    IF(END_CURS(CURSVAL))GOTO 70
	    CALL MON_NUM(T1,T2,XVAL,1,IDX+2)
	    CALL MON_NUM(T1,T2-T3,YVAL,1,IDY+2)
	  END DO
70	  CONTINUE
!
	ELSE IF(ANS .EQ. 'SLP')THEN
	  CALL PGQCS(4,T1,T3)			!T3=Character Height
	  T3=T3*1.5
	  T1=XPAR(1)-(XPAR(2)-XPAR(1))/15.0
	  T2=YPAR(2)+(YPAR(2)-YPAR(1))/15.0
!
	  DO WHILE(1 .EQ. 1)
	    CURSERR = PGCURS(XVAL,YVAL,CURSVAL)
	    IF(END_CURS(CURSVAL))EXIT
	    CURSERR = PGCURS(XVAL_SAV,YVAL_SAV,CURSVAL)
	    IF(END_CURS(CURSVAL))EXIT
	    Q=VECPEN(1)
	    CALL PGSCI(IONE)
	    T3=(YVAL_SAV-YVAL)/(XVAL_SAV-XVAL)
	    CALL PGMOVE(XVAL,YVAL)
	    CALL PGDRAW(XVAL_SAV,YVAL_SAV)
	    WRITE(6,*)'Slope of line is:',T3
	  END DO
!
!
!
! Set Line vectors
!
	ELSE IF(ANS .EQ. 'VF')THEN
	  INIT=.FALSE.
	  IF(VEC)CALL NEW_GEN_IN(INIT,'Initialize line drawing list ?')
	  VEC=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXVEC
	      FLAGLINE(I)=.FALSE.
	    END DO
	  END IF
!
	  L=INDEX(FILNAME,'.STR')
	  IF(L .NE. 0)FILNAME(L:L+3)='.VEC'
1100	  CALL NEW_GEN_IN(FILNAME,'Line vector file forinput')
	  CALL SET_CASE_UP(FILNAME,1,0)
	  L=LEN_TRIM(FILNAME)
	  ISTR=1
	  IF(L .NE. 0)THEN
	    IF(INDEX(FILNAME,'.') .EQ. 0)FILNAME=FILNAME(:L)//'.VEC'
            OPEN(UNIT=33,FILE=TRIM(FILNAME),STATUS='OLD',ERR=1100)
	    DO ISTR=1,MAXVEC
	      READ(33,*,END=1110,ERR=1100)LINEXST(ISTR),LINEYST(ISTR),
	1              LINEXEND(ISTR),LINEYEND(ISTR),VECPEN(ISTR)
	      FLAGLINE(ISTR)=.TRUE.
	    END DO
1110	    CONTINUE
	    CLOSE(UNIT=33)
	  END IF
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'VE')THEN
!
! First output vectors to user, before doing online editing.
! NEW_GEN_IN is used so user can see current value. This section
! also used to input vectors by hand.
!
	  DO I=1,MAXVEC
	    IF(FLAGLINE(I))THEN
	      WRITE(T_OUT,'(1X,I2,4(3X,E14.6),I4)')I,
	1               LINEXST(I),LINEXEND(I),
	1               LINEYST(I),LINEYEND(I),VECPEN(I)
	    END IF
	  END DO
!
	  DO WHILE (0 .EQ. 0)
	    ISTR=0			!Default terminates input
	    CALL NEW_GEN_IN(ISTR,'Vector number - 0 exits')
	    IF(ISTR .LE. 0 .OR. ISTR .GT. MAXVEC)GOTO 1150
	    CALL NEW_GEN_IN(FLAGLINE(ISTR),'FLAG')
	    IF(FLAGLINE(ISTR))THEN
	      CALL NEW_GEN_IN(LINEXST(ISTR),'XST')
	      CALL NEW_GEN_IN(LINEYST(ISTR),'YST')
	      CALL NEW_GEN_IN(LINEXEND(ISTR),'XEND')
	      CALL NEW_GEN_IN(LINEYEND(ISTR),'YEND')
	    END IF
	  END DO
 1150	  VEC=.TRUE.
	  QUERYFLAG=.FALSE.
	  CALL NEW_GEN_IN(QUERYFLAG,'Change vector pens?')
	  IF(QUERYFLAG)CALL VECTORPEN(VECPEN,MAXVEC,FLAGLINE)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'VC')THEN
	  INIT=.FALSE.
	  IF(VEC)CALL NEW_GEN_IN(INIT,'Initialize line drawing list?')
	  VEC=.TRUE.
	  IF(INIT)THEN
	    DO I=1,MAXVEC
	      FLAGLINE(I)=.FALSE.
	    END DO
	  END IF
	  ISTR=1
	  DO WHILE (ISTR .LT. MAXVEC)
	    IF( .NOT. FLAGLINE(ISTR) )THEN
	      CURSERR = PGCURS(LINEXST(ISTR),LINEYST(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
	      CURSERR = PGCURS(LINEXEND(ISTR),LINEYEND(ISTR),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
	      IF(LINEXST(ISTR) .EQ. LINEXEND(ISTR) .AND.
	1        LINEYST(ISTR) .EQ. LINEYEND(ISTR) )THEN
	         ISTR=1000
	      ELSE
	        FLAGLINE(ISTR)=.TRUE.
	      END IF
	    END IF
	    ISTR=ISTR+1
	  END DO
	  GOTO 1000
	ELSE IF(ANS .EQ. 'RID')THEN
	  N_LINE_IDS=0
	  J=0
	  CALL NEW_GEN_IN(ID_FILNAME,'File with line IDs')
	  CALL NEW_GEN_IN(TAU_CUT,'Omit lines with central optical depth <')
	  CALL NEW_GEN_IN(AIR_WAVELENGTHS,'Air wavelengths?')
	  CALL SET_CASE_UP(ID_FILNAME,1,0)
          OPEN(UNIT=33,FILE=TRIM(ID_FILNAME),STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    TMP_STR='!'
	    DO WHILE(TMP_STR(1:1) .EQ. '!')
	      READ(33,'(A)')TMP_STR
	    END DO
	    BACKSPACE(33)
	    DO WHILE(J+1 .LE. 5000)
	      READ(33,*,END=1500)LINE_ID(J+1),ID_WAVE(J+1),TAU(J+1),ID_WAVE_OFF(J+1),ID_Y_OFF(J+1)
	      IF(ID_WAVE(J+1) .LT. 0.0D0)THEN
	         ID_WAVE(J+1)=ABS(ID_WAVE(J+1))
	         OBSERVED_WAVE(J+1)=.FALSE.
	      ELSE
	         OBSERVED_WAVE(J+1)=.TRUE.
	      END IF
	     IF(AIR_WAVELENGTHS)THEN
	         DP_T1=ID_WAVE(J+1)
	         ID_WAVE(J+1)=LAM_AIR(DP_T1)
	      END IF
	      IF( (ID_WAVE(J+1)-XPAR(1))*(XPAR(2)-ID_WAVE(J+1)) .GT. 0 .AND.
	1                   ABS(TAU(J+1)) .GT. TAU_CUT)THEN
	         J=J+1
	         N_LINE_IDS=J
	      END IF
	      IF(LINE_ID(J)(2:2) .EQ. 'k')LINE_ID(J)(2:2)='i'
	      IF(LINE_ID(J)(2:2) .EQ. '2')LINE_ID(J)(2:)='II'//LINE_ID(J)(3:)
	      IF(LINE_ID(J)(3:3) .EQ. '2')LINE_ID(J)(3:)='II'//LINE_ID(J)(4:)
	      IF(LINE_ID(J)(2:5) .EQ. 'XSIX')LINE_ID(J)(2:)='XVI'//LINE_ID(J)(6:)
	      IF(LINE_ID(J)(2:5) .EQ. 'XSEV')LINE_ID(J)(2:)='XVII'//LINE_ID(J)(6:)
	      IF(LINE_ID(J)(3:6) .EQ. 'XSIX')LINE_ID(J)(3:)='XVI'//LINE_ID(J)(7:)
	      IF(LINE_ID(J)(3:6) .EQ. 'XSEV')LINE_ID(J)(3:)='XVII'//LINE_ID(J)(7:)
	      IF(LINE_ID(J)(2:4) .EQ. 'SIX')LINE_ID(J)(2:)='VI'//LINE_ID(J)(5:)
	      IF(LINE_ID(J)(2:4) .EQ. 'SEV')LINE_ID(J)(2:)='VII'//LINE_ID(J)(5:)
	      IF(LINE_ID(J)(3:5) .EQ. 'SIX')LINE_ID(J)(3:)='VI'//LINE_ID(J)(6:)
	      IF(LINE_ID(J)(3:5) .EQ. 'SEV')LINE_ID(J)(3:)='VII'//LINE_ID(J)(6:)
	    END DO
	  ELSE
	    WRITE(T_OUT,*)'Unable top open file'
	    GOTO 1000
	  END IF
1500	  CONTINUE
	  CLOSE(UNIT=33)
	  WRITE(T_OUT,*)'Number of lines read in is',N_LINE_IDS
	  CALL NEW_GEN_IN(ID_SCL,'Factor to scale line location')
	  CALL NEW_GEN_IN(ID_EXPCHAR,'Factor to scale size of ID')
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SID')THEN
	  CALL NEW_GEN_IN(ID_SCL,'Factor to scale line location')
	  CALL NEW_GEN_IN(ID_VEC_BEG,'Location to start ID line')
	  CALL NEW_GEN_IN(ID_VEC_END,'Location to end ID line')
	  CALL NEW_GEN_IN(ID_EXPCHAR,'Factor to scale size of ID')
	  N_OMIT_ID=0
!
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)RED_PEN,'Species name are case sensitive and have no spaces'
	  WRITE(T_OUT,*)'Use SiIII etc',DEF_PEN
	  WRITE(T_OUT,*)' '
	  DO I=1,10
	    OMIT_ID(I)='ZZ'
	    CALL NEW_GEN_IN(OMIT_ID(I),'Species ID to OMIT from labels')
	    IF(OMIT_ID(I) .EQ. 'ZZ')EXIT
	    N_OMIT_ID=I
	  END DO
	  N_INC_ID=0
	  DO I=1,10
	    INC_ID(I)='ZZ'
	    CALL NEW_GEN_IN(INC_ID(I),'Species ID to INC in labels')
	    IF(INC_ID(I) .EQ. 'ZZ')EXIT
	    N_INC_ID=I
	  END DO
!
! 
!
! Define a crude straight line continuum about a line, so the EW can be
! computed.
!
	ELSE IF (ANS .EQ. 'DC')THEN
	  IF(NPLTS .EQ. 1)THEN
	    PLOT_ID=1
	  ELSE
	    CALL NEW_GEN_IN(PLOT_ID,'Plot to determine EW for:')
	  END IF
	  CURSERR = PGCURS(XCUR(1),YCUR(1),CURSVAL)
	  CURSERR = PGCURS(XCUR(2),YCUR(2),CURSVAL)
	  SLOPE=(YCUR(2)-YCUR(1))/(XCUR(2)-XCUR(1))
	  IF(ALLOCATED(CONT))DEALLOCATE(CONT)
	  ALLOCATE(CONT(NPTS(PLOT_ID)))
	  DO I=1,NPTS(PLOT_ID)
            CONT(I)=YCUR(1)+SLOPE*(CD(PLOT_ID)%XVEC(I)-XCUR(1))
	  END DO
	  WRITE(T_OUT,*)XCUR(1:2),YCUR(1:2)
	  CALL PGLINE(NPTS(PLOT_ID),CD(PLOT_ID)%XVEC,CONT)
	  CONTINUUM_DEFINED=.TRUE.
	  IP=NPLTS+1
	  CALL NEW_GEN_IN(IP,'Output plot?')
	  IF(IP .NE. 0)THEN
	    I=NPTS(PLOT_ID)
	    TYPE_CURVE(IP)='L'
	    IF(ALLOCATED(CD(IP)%XVEC))THEN
	      DEALLOCATE (CD(IP)%XVEC)
	      DEALLOCATE (CD(IP)%DATA)
	    END IF
            ALLOCATE (CD(IP)%XVEC(I),STAT=IOS)
            IF(IOS .EQ. 0)ALLOCATE (CD(IP)%DATA(I),STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error: unable to allocate new data vectors'
	      WRITE(T_OUT,*)'IOS=',IOS
	      STOP
	    END IF
	    CD(IP)%XVEC(1:I)=CD(PLOT_ID)%XVEC(1:I)
            CD(IP)%DATA(1:I)=CONT(1:I)
            NPTS(IP)=I
            ERR(IP)=.FALSE.
            IF(IP .GT. NPLTS)NPLTS=IP
	  END IF
	  GOTO 1000
C
	ELSE IF(ANS .EQ. 'GF')THEN
	  DRAW_GAUSS_HARD=.TRUE.
	  CALL DO_GAUS_FIT(XPAR(1),XPAR(2))
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'WGF')THEN
	  CALL WR_GAUS_FIT
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'MGF')THEN
	  CALL DO_MULT_GF(XPAR(1),XPAR(2))
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'EG')THEN
	  CALL ED_GAUS_FIT()
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SEW')THEN
	  T1=0.002D0
	  CALL SIMP_EW(CD(1)%DATA(1),CD(1)%XVEC(1),NPTS(1),XPAR(1),XPAR(2),T1)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'DG')THEN
	  CALL DRAW_GAUS()
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'EW')THEN
!
	  IF(CONTINUUM_DEFINED)THEN
!
	  CURSERR = PGCURS(XCUR(1),YCUR(1),CURSVAL)
	  CURSERR = PGCURS(XCUR(2),YCUR(2),CURSVAL)
!
! Find nearest channels to curser positions.
!
	    IP=PLOT_ID
	    L_CHAN(1)=GET_INDX_SP(XCUR(1),CD(IP)%XVEC,NPTS(IP))
	    L_CHAN(2)=GET_INDX_SP(XCUR(2),CD(IP)%XVEC,NPTS(IP))
!
! Draw in line limits
!
	    CALL PGMOVE(CD(IP)%XVEC(L_CHAN(1)),YPAR(1))
	    CALL PGDRAW(CD(IP)%XVEC(L_CHAN(1)),CONT(L_CHAN(1)))
	    CALL PGMOVE(CD(IP)%XVEC(L_CHAN(2)),YPAR(1))
	    CALL PGDRAW(CD(IP)%XVEC(L_CHAN(2)),CONT(L_CHAN(2)))
!
	    EW=0.0
            CENTROID=0.0
	    IF(L_CHAN(1) .LT. L_CHAN(2))THEN
	      T1=0.5D0
	    ELSE
	      K=L_CHAN(1)
	      L_CHAN(1)=L_CHAN(2)
	      L_CHAN(2)=K
	      T1=-0.5D0
	    END IF
	    DO I=L_CHAN(1),L_CHAN(2)-1
	      EW=EW+T1*( (CD(IP)%DATA(I)-CONT(I))/CONT(I) + 
	1               (CD(IP)%DATA(I+1)-CONT(I+1))/CONT(I+1) )*
	1               (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	      CENTROID=CENTROID+T1*( CD(IP)%XVEC(I)*(CD(IP)%DATA(I)-CONT(I))/CONT(I)
	1               + CD(IP)%XVEC(I)*(CD(IP)%DATA(I+1)-CONT(I+1))/CONT(I+1) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	    END DO
	    IF(EW .NE. 0.0)THEN
	      CENTROID=CENTROID/EW
	    ELSE
	      CENTROID=0.0D0
	    END IF
	    WRITE(T_OUT,*)'The equivalent width of the line is =',EW,'X(units)'
	    WRITE(T_OUT,*)'The central postion of the line is   ',CENTROID
	    GOTO 1000
!
	  ELSE
!
	    QUERYFLAG=.FALSE.
	    CALL NEW_GEN_IN(QUERYFLAG,'Compute area?')
	    IF(.NOT. QUERYFLAG)THEN
              WRITE(T_OUT,*)'Continuum not defined'
	      WRITE(T_OUT,*)'Assuming Ic=1 (i.e. normalized data)'
!
! Show the continuum
!
	      T3=1.0
	      CALL PGMOVE(XPAR(1),T3)
	      CALL PGDRAW(XPAR(2),T3)
	    END IF
!
! Get the location of the line.
!
	    DO WHILE(1 .EQ. 1)
	      CURSERR = PGCURS(XCUR(1),YCUR(1),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
	      CURSERR = PGCURS(XCUR(2),YCUR(2),CURSVAL)
	      IF(END_CURS(CURSVAL))GOTO 1000
!
! We do all plot provided the cover the range indicated by the cursors.
! Find nearest channels to curser positions.
!
! Draw in line limits
!
	    IF(QUERYFLAG)THEN
	      T3=0.3*(YPAR(2)-YPAR(1))+YPAR(1)
	    ELSE
	      T3=1.0D0
	    END IF
	    CALL PGMOVE(XCUR(1),YPAR(1))
	    CALL PGDRAW(XCUR(1),T3)
	    CALL PGMOVE(XCUR(2),YPAR(1))
	    CALL PGDRAW(XCUR(2),T3)
!
	    DO IP=1,NPLTS
	      T1=XCUR(1)-CD(IP)%XVEC(1)
	      T2=XCUR(1)-CD(IP)%XVEC(NPTS(IP))
	      T3=XCUR(2)-CD(IP)%XVEC(1)
	      T4=XCUR(2)-CD(IP)%XVEC(NPTS(IP))
	      IF(T1*T2 .LT. 0 .AND. T3*T4 .LT. 0)THEN
	        L_CHAN(1)=GET_INDX_SP(XCUR(1),CD(IP)%XVEC,NPTS(IP))
	        L_CHAN(2)=GET_INDX_SP(XCUR(2),CD(IP)%XVEC,NPTS(IP))
!
	        IF(QUERYFLAG)THEN
!
! Compute area of marked region.
!
	          EW=0.0
                  CENTROID=0.0
	          K=1
	          IF(L_CHAN(1) .LT. L_CHAN(2))THEN
	            T1=0.5D0
	          ELSE
	            K=L_CHAN(1)
	            L_CHAN(1)=L_CHAN(2)
	            L_CHAN(2)=K
	            T1=-0.5D0
	          END IF
	          DO I=L_CHAN(1),L_CHAN(2)-1
	            EW=EW+T1*( CD(IP)%DATA(I) + CD(IP)%DATA(I+1) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	            CENTROID=CENTROID+T1*( CD(IP)%XVEC(I)*CD(IP)%DATA(I)
	1               + CD(IP)%XVEC(I)*CD(IP)%DATA(I+1) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	          END DO
	          IF(EW .NE. 0.0)THEN
	            CENTROID=CENTROID/EW
	          ELSE
	            CENTROID=0.0
	          END IF
	          WRITE(T_OUT,'(A,I3,3X,5X,A,1PE10.3,A,A,1PE14.6)')' Plot ID=',IP,'  AREA=',EW,' X(units);','  Centroid=',CENTROID
	          WRITE(LU_EW,'(A,I3,3X,5X,A,1PE10.3,A,A,1PE14.6)')' Plot ID=',IP,'  AREA=',EW,' X(units);','  Centroid=',CENTROID
	        ELSE
!
! Compute EW
!
	          EW=0.0
                  CENTROID=0.0
	          DO I=L_CHAN(1),L_CHAN(2)-1
	            EW=EW+0.5*( (CD(IP)%DATA(I)-1.0) + (CD(IP)%DATA(I+1)-1.0) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	            CENTROID=CENTROID+0.5*( CD(IP)%XVEC(I)*(CD(IP)%DATA(I)-1.0)
	1               + CD(IP)%XVEC(I)*(CD(IP)%DATA(I+1)-1.0) )*
	1                 (CD(IP)%XVEC(I+1)-CD(IP)%XVEC(I))
	          END DO
	          IF(EW .NE. 0.0)THEN
	            CENTROID=CENTROID/EW
	          ELSE
	            CENTROID=0.0
	          END IF
!
	          IF(FIRST_EW)THEN
	            OPEN(UNIT=LU_EW,FILE='EW_FR_SPEC_PLT',STATUS='UNKNOWN')
	            FIRST_EW=.FALSE.
	            WRITE(LU_EW,'(5A14)')'     Plot ID',' EW(X units)',' X(centroid)',
	1                        '    X(start)','      X(end)'
	          END IF
	          WRITE(T_OUT,'(A,I3,5X,A,1PE10.3,A,A,ES14.6,2ES14.4)')
	1            ' Plot ID=',IP,'  EW=',EW,' X(units);','  Centroid=',CENTROID,
	1                  CD(IP)%XVEC(L_CHAN(1)),CD(IP)%XVEC(L_CHAN(2))
	          IF(IP .EQ. 1)THEN
	            WRITE(LU_EW,'(A)')' '
	            TMP_STR=' '
	            CALL NEW_GEN_IN(TMP_STR,'Comment=')
	            IF(TMP_STR .NE. ' ')THEN
	              WRITE(LU_EW,'(A)')TRIM(TMP_STR)
	              WRITE(LU_EW,'(A)')' '
	              TMP_STR=' '
	              CALL NEW_GEN_IN(TMP_STR,'Comment=')
	              IF(TMP_STR .NE. ' ')WRITE(LU_EW,'(A)')TRIM(TMP_STR)
	            END IF
	          END IF
	          WRITE(LU_EW,'(10X,I4,4ES14.6)')
	1            IP,EW,CENTROID,CD(IP)%XVEC(L_CHAN(1)),CD(IP)%XVEC(L_CHAN(2))
	        END IF
	      END IF
	    END DO
	    END DO		!Multiple plots
	    GOTO 1000
	  END IF
	ELSE IF(ANS .EQ. 'EWG')THEN
	  CALL DO_MANY_EW(TYPE_CURVE,NPLTS,MAX_PLTS)
!
	ELSE IF(ANS .EQ. 'LP')THEN
	  IF(LONG_PLOT)THEN
	    WRITE(T_OUT,*)'Resuming normal hard copy mode'
	    LONG_PLOT=.FALSE.
	  ELSE
	    WRITE(T_OUT,*)'Setting hard-copy to produce long plots'
	    WRITE(T_OUT,*)'This option only work CPS as the device'
	    LONG_PLOT=.TRUE.
	    CALL NEW_GEN_IN(LENGTH_OF_HC_PLOT,'Length of plot surface in cm''s')
	    LP_ASR=2.54D0*8.0D0/LENGTH_OF_HC_PLOT
	    WRITE(T_OUT,*)'Default aspect ratio of pot surface is',LP_ASR
	  END IF
!
!
!
	ELSE IF(ANS .EQ. 'Z' .OR. ANS .EQ. 'ZN')THEN
	  IF (FIRST_HARD .OR. ANS .EQ. 'ZN') THEN
	    WRITE(T_OUT,*)RED_PEN,' '
	    WRITE(T_OUT,*)'Choose a post-script device and file for printing [file/dev] '
	    WRITE(T_OUT,*)'A sensible name format is xxx_1.ps/cps'
	    WRITE(T_OUT,*)'Use a . only for the file extension'
	    WRITE(T_OUT,*)DEF_PEN,' '
	    CALL NEW_GEN_IN(PRINTER,'File and printer: enter ? for list')
	    BEG=PGOPEN(PRINTER)
	    FIRST_HARD=.FALSE.
	    CALL DEFINE_MORE_PENS(MAXPEN)
	  ELSE
	    CALL NEW_GEN_IN(HARD_FILE,'Plot file')
	    PRINTER=TRIM(HARD_FILE)//'/'//HARD_TYPE
	    BEG=PGOPEN(PRINTER)
	    CALL DEFINE_MORE_PENS(MAXPEN)		!Pen definitions are not saved.
	  END IF
!
! Save hard device and set default file name for net plot.
!
	  CALL PGQINF('TYPE',HARD_TYPE,LENGTH)
	  CALL PGQINF('FILE',HARD_FILE,LENGTH)
	  CALL GET_PGI_FILE_NUMBER(HARD_CNT,HARD_FILE)
	  CALL UPDATE_PGI_FILE_NAME(HARD_CNT,HARD_FILE)
	  PRINTER=TRIM(HARD_FILE)//'/'//HARD_TYPE
!
	  HARD=.TRUE.
!
1200	  CALL NEW_GEN_IN(PLT_LINE_WGT,'line weight') 
	  IF (PLT_LINE_WGT .GT. 201) THEN
	    WRITE(T_OUT,*)'value too large'
	    GOTO 1200
	  END IF
	  GOTO 5000
!
	ELSE IF(ANS .EQ. 'CL')THEN
	  IF(.NOT. HARD)THEN			!Terminal
	    CALL PGERAS
	    CALL PGETXT
	    CALL PGPAGE
	    READ(T_IN,'(A)')ANS
	  END IF
	  GOTO 1000
	ELSE IF(ANS .EQ. 'WP' .OR. ANS .EQ. 'WPF')THEN
	  IF(PLT_ST_FILENAME .EQ. ' ')THEN
	    PLT_ST_FILENAME='PLT_ST'
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  ELSE IF(ANS .EQ. 'WPF')THEN
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  END IF
	  OPEN(UNIT=30,FORM='UNFORMATTED',FILE=PLT_ST_FILENAME,STATUS='UNKNOWN'
	1,       IOSTAT=IOS,POSITION='APPEND')
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
	  L=-1
	  IF(NPLTS .NE. 1)CALL NEW_GEN_IN(L,'Plot to output (def=-1=ALL)')
!
	  CNT=0; PLT_ID=' '
	  DO IP=1,NPLTS
	    IF(L .EQ. -1 .OR. L .EQ. IP)THEN
	      CNT=CNT+1
	      CALL NEW_GEN_IN(PLT_ID,'PLT_ID=')
	      WRITE(30)'PLT_ID=',PLT_ID
	      WRITE(30)NPTS(IP)
	      WRITE(30)ERR(IP)
	      DO K=1,(NPTS(IP)+N_REC_SIZE-1)/N_REC_SIZE
	        IST=N_REC_SIZE*(K-1)+1
	        IEND=MIN(IST+N_REC_SIZE-1,NPTS(IP))
	        WRITE(30)(CD(IP)%XVEC(I),I=IST,IEND)
	        WRITE(30)(CD(IP)%DATA(I),I=IST,IEND)
	        IF(ERR(IP))THEN
	          WRITE(30)(CD(IP)%EMAX(I),I=IST,IEND)
	          WRITE(30)(CD(IP)%EMIN(I),I=IST,IEND)
	        END IF
	      END DO
	    END IF
	  END DO
	  WRITE(T_OUT,'(I4,A)')CNT,' plots written to data file'
	  CLOSE(UNIT=30)
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'RP' .OR. ANS .EQ. 'RPF')THEN
	  IF(PLT_ST_FILENAME .EQ. ' ')THEN
	    PLT_ST_FILENAME='PLT_ST'
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  ELSE IF(ANS .EQ. 'RPF')THEN
	    CALL NEW_GEN_IN(PLT_ST_FILENAME,'Name of file with stored plots')
	  END IF
	  OPEN(UNIT=30,FORM='UNFORMATTED',FILE=PLT_ST_FILENAME,STATUS='OLD',
	1          IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
	  PLT_ID=' '
	  PLT_ID_SAV=' '
!
	  DO WHILE(1 .EQ. 1)
	    CALL NEW_GEN_IN(PLT_ID,'PLT_ID=')
	    IF(PLT_ID .EQ. ' ' .OR. PLT_ID .EQ. PLT_ID_SAV)EXIT
	    RD_PLT_ID=' '
	    DO WHILE(RD_PLT_ID .NE. 'PLT_ID='//PLT_ID)
	      READ(30,IOSTAT=IOS)RD_PLT_ID
	      IF(IOS .EQ. -1)THEN
	        WRITE(T_OUT,*)'Unable to identify plot --- IOS',IOS
	        WRITE(T_OUT,*)'Available plots are as follows:'
	        REWIND(UNIT=30)
	        IOS=0
                DO WHILE(IOS .NE. -1)
	          RD_PLT_ID=' '
	          READ(30,IOSTAT=IOS)RD_PLT_ID
	          IF(RD_PLT_ID(1:7) .EQ. 'PLT_ID=')THEN
	             WRITE(T_OUT,'(35X,A)')TRIM(RD_PLT_ID(8:))
	          END IF
	        END DO
	        CLOSE(UNIT=30)
	        GOTO 1000
	      END IF
	    END DO
	    IP=NPLTS+1
	    IF(IP .GT. MAX_PLTS)THEN
	      WRITE(T_OUT,*)'Error in GRAMON_PGPLOT'
	      WRITE(T_OUT,*)'Too many plots have ben stored'
	      GOTO 1000
	    END IF
	    READ(30)NPTS(IP)
	    READ(30)ERR(IP)
	    ALLOCATE (CD(IP)%XVEC(NPTS(IP)))
	    ALLOCATE (CD(IP)%DATA(NPTS(IP)))
	    DO K=1,(NPTS(IP)+N_REC_SIZE-1)/N_REC_SIZE
	      IST=N_REC_SIZE*(K-1)+1
	      IEND=MIN(IST+N_REC_SIZE-1,NPTS(IP))
	      READ(30)(CD(IP)%XVEC(I),I=IST,IEND)
	      READ(30)(CD(IP)%DATA(I),I=IST,IEND)
	      IF(ERR(IP))THEN
	        ALLOCATE (CD(IP)%EMIN(NPTS(IP)))
	        ALLOCATE (CD(IP)%EMAX(NPTS(IP)))
	        READ(30)(CD(IP)%EMAX(I),I=IST,IEND)
	        READ(30)(CD(IP)%EMIN(I),I=IST,IEND)
	      END IF
	    END DO
	    NPLTS=NPLTS+1			!Successfully read.
	    TYPE_CURVE(IP)='L'
	    LINE_WGT(IP)=1
	    IF(IP .NE. 1)THEN
	      LINE_STYLE(IP)=MOD(LINE_STYLE(IP-1),5)+1
              PEN_COL(IP)=PEN_COL(IP-1)+1
            ELSE
	      LINE_STYLE(IP)=1
	      PEN_COL(IP)=2
	    END IF
	    LINE_STYLE(IP)=1
	    IF(DASH)LINE_STYLE(IP)=MOD(IP-1,5)+1
	    PLT_ID_SAV=PLT_ID
	    REWIND(UNIT=30)
	  END DO
!
	  CLOSE(UNIT=30)
	  GOTO 1000
!
!
! Option to read simple plots in column format from file.
!
	ELSE IF(ANS .EQ. 'RXY')THEN
	  FILNAME=' '
	  CALL NEW_GEN_IN(FILNAME,'FILE=')
	  OPEN(UNIT=30,FILE=FILNAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
!
! Use FILNAME as temporary STRING
!
	  CNT=0
	  FILNAME='!'
	  DO WHILE(FILNAME(1:1) .EQ. '!' .OR. FILNAME .EQ. ' ')
	    READ(30,'(A200)',IOSTAT=IOS)FILNAME
	    WRITE(T_OUT,*)TRIM(FILNAME)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error reading file'
	      CLOSE(UNIT=30)
	      GOTO 1000
	    END IF
	    IF(INDEX(FILNAME,'Number of data points:') .NE. 0)THEN
	      I=INDEX(FILNAME,'Number of data points:')
	      READ(FILNAME(I+22:),*)CNT
	      FILNAME='!'
	    END IF
	  END DO
	  BACKSPACE(UNIT=30)
!
	  IP=NPLTS+1
	  NPTS(IP)=CNT
	  CALL NEW_GEN_IN(NPTS(IP),'Number of data points')
	  COLUMN(1)=1; COLUMN(2)=2
	  CALL NEW_GEN_IN(COLUMN,I,ITWO,'Data columns')
	  K=MAX(COLUMN(1),COLUMN(2))
!
	  ALLOCATE (CD(IP)%XVEC(NPTS(IP)),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CD(IP)%DATA(NPTS(IP)),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error allocating data arrays'
	     GOTO 1000
	  END IF
	  DO I=1,NPTS(IP)
	    READ(30,*,IOSTAT=IOS,END=500)(TEMP_VAR(L),L=1,K)
	    CD(IP)%XVEC(I)=TEMP_VAR(COLUMN(1))
	    CD(IP)%DATA(I)=TEMP_VAR(COLUMN(2))
	  END DO
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error in read plot from XY file'
	    WRITE(T_OUT,*)'Data point is',I+1
	    CLOSE(UNIT=30)
	    GOTO 1000
	  END IF
500	  CONTINUE
!
! Data successfully read
!	
	  LAST_DP=LAST_DP+NPTS(IP)
	  NPLTS=NPLTS+1; J=NPLTS
	  TYPE_CURVE(J)='L'
	  LINE_WGT(J)=1
	  IF(J .NE. 1)THEN
	    LINE_STYLE(J)=MOD(LINE_STYLE(J-1),5)+1
            PEN_COL(J)=PEN_COL(J-1)+1
          ELSE
	    LINE_STYLE(J)=1
	    PEN_COL(J)=2
	  END IF
	  CLOSE(UNIT=30)
	  GOTO 1000
!i
	ELSE IF(ANS .EQ. 'WXY')THEN
	  FILNAME=' '
	  CALL NEW_GEN_IN(FILNAME,'FILE=')
	  OPEN(UNIT=30,FILE=FILNAME,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'
	    GOTO 1000
	  END IF
!
	  L=-1
	  CALL NEW_GEN_IN(L,'Plot to output (def=-1=ALL)')
!
	  IF(L .EQ. -1)THEN
	    WRITE(30,*)NPLTS
	    WRITE(30,*)(NPTS(I),I=1,NPLTS)
	    DO J=1,MAXVAL(NPTS)
	      DO IP=1,NPLTS
	        IF(J .LE. NPTS(IP))THEN
	          WRITE(30,'(2X,2ES14.6)',ADVANCE='NO')CD(IP)%XVEC(J),CD(IP)%DATA(J)
	        ELSE
	          WRITE(30,'(2X,2ES14.6)',ADVANCE='NO')0.0D0,0.0D0
	        END IF
	      END DO
	      WRITE(30,'(A)')' '
	    END DO
	    WRITE(T_OUT,*)NPLTS,' plots written to ',TRIM(FILNAME)
	  ELSE
	    WRITE(30,*)NPLTS
	    WRITE(30,*)NPTS(L)
	    IP=L
	    DO J=1,NPTS(IP)
	      WRITE(30,*)CD(IP)%XVEC(J),CD(IP)%DATA(J)
	    END DO
	    WRITE(T_OUT,*)'One plot written to ',TRIM(FILNAME)
	  END IF
	  CLOSE(UNIT=30)
!
! 
!
	ELSE IF(ANS .EQ. 'SXY')THEN
	  DO IP=1,NPLTS
	    DO J=1,NPTS(IP)
	      IF(CD(IP)%XVEC(J) .GE. XPAR(1) .AND. CD(IP)%XVEC(J) .LE. XPAR(2))THEN
	        WRITE(T_OUT,*)IP,J,CD(IP)%XVEC(J),CD(IP)%DATA(J)
	      END IF
	    END DO
	  END DO
! 
!
	ELSE IF(ANS .EQ. 'LAM')THEN
	  WRITE(6,*)RED_PEN
	  WRITE(6,'(A,/)')' All wavelengths are vacuum'//BLUE_PEN
	  CALL WRITE_LINE_LAMBDAS()
	  WRITE(6,*)DEF_PEN
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'CUT')THEN          !Recall ANS in always upper case
	  WRITE(6,'(A)')BLUE_PEN
	  WRITE(6,'(A)')' This option reduces the number of data points in each curve.'
	  WRITE(6,'(A)')' It is useful for creating smaller publication quality plots.'//RED_PEN
	  WRITE(6,'(A)')' This option cannot be undone.'//BLUE_PEN
	  WRITE(6,'(A)')' A negative cut accuracy does nothing'
	  WRITE(6,'(A)')' To compare with original data:'
	  WRITE(6,'(A)')'    Store original data using WP & read after cut with RP, OR'
	  WRITE(6,'(A)')'    do CUT, NOI, and redo plots'
	  WRITE(6,'(A)')DEF_PEN
	  CUT_ACCURACY=0.001D0
	  CALL NEW_GEN_IN(CUT_ACCURACY,'Fractional accuracy to retain in plots')
	  IF(CUT_ACCURACY .GT. 0.0)CALL SHRINK_VECTORS(CUT_ACCURACY)
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'KMS')THEN          !Recall ANS in always upper case
	  VEL_UNIT='km/s'
	  WRITE(6,*)'Using km/s for velocity unit'
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'MMS')THEN
	  VEL_UNIT='Mm/s'
	  WRITE(6,*)'Using Mm/s for velocity unit'
	  GOTO 1000
!
! Convert from Ang to velocity space. Data must have been originally
! in Ang. This option can be done many times, as old Ang scale is
! restored on each call.
!
	ELSE IF (ANS .EQ. 'VEL')THEN
	  C_KMS=1.0D-05*SPEED_OF_LIGHT()
	  C_VAL=C_KMS
	  IF(VEL_UNIT .EQ. 'Mm/s')C_VAL=1.0D-03*C_KMS
	  CALL NEW_GEN_IN(CENTRAL_LAM,'/\(Ang) [-ve: 10^15 Hz]')
	  IF(CENTRAL_LAM .LT. 0)THEN
	    CENTRAL_LAM=1.0D-02*C_VAL/ABS(CENTRAL_LAM)
	  ELSE IF(CENTRAL_LAM .GT. 2000)THEN
	    AIR_LAM=.FALSE.
            CALL NEW_GEN_IN(AIR_LAM,'Wavelength in air?')
	    IF(AIR_LAM)CENTRAL_LAM=LAM_VAC(CENTRAL_LAM)
	  END IF
!
! Undo previous conversion to V space.
!
	  IF(OLD_CENTRAL_LAM .NE. 0)THEN
	    DO IP=1,NPLTS
	      DO J=1,NPTS(IP)
	        CD(IP)%XVEC(J)=OLD_CENTRAL_LAM*
	1                         (1.0+CD(IP)%XVEC(J)/C_VAL)
	      END DO
	    END DO
	  END IF
!
! Puts X-axis in km/s
!
	  IF(CENTRAL_LAM .NE. 0)THEN
	    T1=C_VAL/CENTRAL_LAM
	    DO IP=1,NPLTS
	      DO J=1,NPTS(IP)
	        CD(IP)%XVEC(J)=T1*(CD(IP)%XVEC(J)-CENTRAL_LAM)
	      END DO
	    END DO
	    XLABEL='V(km\d \us\u-1\d)'
	    IF(VEL_UNIT .EQ. 'Mm/s')XLABEL='V(Mm\d \us\u-1\d)'
	  ELSE
	    XLABEL=XLAB
	  END IF
	  OLD_CENTRAL_LAM=CENTRAL_LAM
	  GOTO 1000
! 
!
! Peform simple X-axis arithmetic.
!
	ELSE IF (ANS .EQ. 'XAR')THEN
	  CALL NEW_GEN_IN(XAR_OPERATION,'Operation: *,+,-,/,LG,ALG[=10^x],R[=1/x]')
	  CALL SET_CASE_UP(XAR_OPERATION,IZERO,IZERO)
	  IF(XAR_OPERATION .NE. 'LG' .AND. XAR_OPERATION .NE. 'ALG' .AND. 
	1           XAR_OPERATION .NE. 'XN' .AND. XAR_OPERATION .NE. 'R')THEN
	    CALL NEW_GEN_IN(XAR_VAL,'Value')
	  END IF
	  CALL NEW_GEN_IN(XAR_PLT,'Which plot? (0=ALL,-ve exits)')
	  IP=XAR_PLT
	  IF(IP .LT. 0 .OR. IP .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Invalid plot number'
	    GOTO 1000
	  END IF
	  IP_ST=IP; IF(IP .EQ. 0)IP_ST=1
	  IP_END=IP; IF(IP .EQ. 0)IP_END=NPLTS
	  DO IP=IP_ST,IP_END
	    IF(XAR_OPERATION .EQ. '+')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC+XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. '-')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC-XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. '*')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC*XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. '/')THEN
	      CD(IP)%XVEC=CD(IP)%XVEC/XAR_VAL
	    ELSE IF(XAR_OPERATION .EQ. 'ALG')THEN
	      CD(IP)%XVEC=10.0D0**(CD(IP)%XVEC)
	      IF(XLABEL(1:3) .EQ. 'Log')THEN
	        XLABEL=XLABEL(4:)
	        XLABEL=ADJUSTL(XLABEL)
	      END IF
	    ELSE IF(XAR_OPERATION .EQ. 'LG')THEN
	      T1=TINY(CD(IP)%XVEC(1))
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%XVEC(J) .GT. T1)THEN
	          CD(IP)%XVEC(J)=LOG10(CD(IP)%XVEC(J))
	        ELSE
	          CD(IP)%XVEC(J)=-1000.0
	        END IF
	      END DO
	      IF(XLABEL(1:3) .NE. 'Log')THEN
	        XLABEL='Log '//XLABEL
	      END IF
	    ELSE IF(XAR_OPERATION .EQ. 'R')THEN
	      T1=0.01D0*C_KMS
	      T2=0.01D0*C_KMS
	      CALL NEW_GEN_IN(T1,'Scale factor')
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%XVEC(J) .GT. 0)THEN
	          CD(IP)%XVEC(J)=T1/CD(IP)%XVEC(J)
	        ELSE
	          CD(IP)%XVEC(J)=-1000
	  	END IF
	      END DO
	      IF(XLABEL .EQ. '\gn(10\u15 \dHz)' .AND. T1 .EQ. T2)THEN
	        XLABEL='\gl(\A)'
	      ELSE IF(XLABEL .EQ. '\gl(\A)' .AND. T1 .EQ. T2)THEN
	        XLABEL='\gn(10\u15 \dHz)'
	      END IF
	    ELSE IF(XAR_OPERATION .EQ. 'XN')THEN
	      WRITE(6,*)'Curent plot is',IP
	      K=0
	      CALL NEW_GEN_IN(K,'Output plot: 0 adds new plot?')
	      IF(K .EQ. IP)THEN
	        K=IP
	      ELSE IF(K .EQ. 0)THEN
	        CALL CURVE(NPTS(IP),CD(IP)%XVEC,CD(IP)%DATA)
	        K=NPLTS
	      END IF
	      DO I=1,NPTS(K)
	        CD(K)%XVEC(I)=I
	      END DO
	      TYPE_CURVE(K)=TYPE_CURVE(IP)
	    ELSE
	      WRITE(T_OUT,*)'Invalid operation: try again'
	      GOTO 1000
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'OFF')THEN
	  CALL NEW_GEN_IN(YAR_PLT,'Which plot? (0=ALL,-ve exits)')
	  IP=YAR_PLT
	  IP_ST=IP; IF(IP .EQ. 0)IP_ST=1
	  IP_END=IP; IF(IP .EQ. 0)IP_END=NPLTS
	  IF(IP .LT. 0 .OR. IP .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Invalid plot number'
	    GOTO 1000
	  END IF
	  CALL NEW_GEN_IN(YAR_VAL,'Offset')
	  DO IP=IP_ST,IP_END
	    CD(IP)%DATA=CD(IP)%DATA+YAR_VAL*(IP-IP_ST)
	  END DO
!
! Perform simple Y-axis arithmetic.
!
	ELSE IF (ANS .EQ. 'YAR')THEN
	  CALL NEW_GEN_IN(YAR_OPERATION,'Operation: *,+,-,/,LG,R[=1/Y],ALG[=10^y],ABS')
	  CALL SET_CASE_UP(YAR_OPERATION,IZERO,IZERO)
	  IF(YAR_OPERATION .EQ. 'MX' .OR. YAR_OPERATION .EQ. 'DX')THEN
	    YAR_VAL=0.5D0*(XPAR(1)+XPAR(2))
	    CALL NEW_GEN_IN(YAR_VAL,'X normalization')
	  ELSE IF( YAR_OPERATION .NE. 'LG' .AND. 
	1          YAR_OPERATION .NE. 'ABS' .AND. 
	1          YAR_OPERATION .NE. 'ALG' .AND. 
	1          YAR_OPERATION(1:1) .NE. 'R')THEN
	    CALL NEW_GEN_IN(YAR_VAL,'Value')
	  END IF
	  CALL NEW_GEN_IN(YAR_PLT,'Which plot? (0=ALL,-ve exits)')
	  IP=YAR_PLT
	  IP_ST=IP; IF(IP .EQ. 0)IP_ST=1
	  IP_END=IP; IF(IP .EQ. 0)IP_END=NPLTS
	  IF(IP .LT. 0 .OR. IP .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Invalid plot number'
	    GOTO 1000
	  END IF
	  DO IP=IP_ST,IP_END
	    IF(YAR_OPERATION .EQ. '+')THEN
	      CD(IP)%DATA=CD(IP)%DATA+YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. '-')THEN
	      CD(IP)%DATA=CD(IP)%DATA-YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. '*')THEN
	      CD(IP)%DATA=CD(IP)%DATA*YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. '/')THEN
	      CD(IP)%DATA=CD(IP)%DATA/YAR_VAL
	    ELSE IF(YAR_OPERATION .EQ. 'ABS')THEN
	      CD(IP)%DATA=ABS(CD(IP)%DATA)
	    ELSE IF(YAR_OPERATION .EQ. 'ALG')THEN
	      CD(IP)%DATA=10.0D0**(CD(IP)%DATA)
	      IF(YLABEL(1:3) .EQ. 'Log')THEN
	        YLABEL=YLABEL(4:)
	        YLABEL=ADJUSTL(YLABEL)
	      END IF
	    ELSE IF(YAR_OPERATION .EQ. 'LG')THEN
	      T1=TINY(CD(IP)%DATA(1))
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%DATA(J) .GT. T1)THEN
	          CD(IP)%DATA(J)=LOG10(CD(IP)%DATA(J))
	        ELSE
	          CD(IP)%DATA(J)=-1000.0
	        END IF
	      END DO
	      IF(YLABEL(1:3) .NE. 'Log')THEN
	        YLABEL='Log '//YLABEL
	      END IF
	    ELSE IF(YAR_OPERATION .EQ. 'R')THEN
	      DO J=1,NPTS(IP)
	        IF(ABS(CD(IP)%DATA(J)) .GT. 1.0E-37)THEN
	          CD(IP)%DATA(J)=1.0D0/CD(IP)%DATA(J)
	        ELSE
	          CD(IP)%DATA(J)=1.0E+37
	        END IF
	      END DO
	    ELSE IF(YAR_OPERATION .EQ. 'MX')THEN
	      DO J=1,NPTS(IP)
	       CD(IP)%DATA(J)=CD(IP)%DATA(J)*(CD(IP)%XVEC(J)/YAR_VAL)
	      END DO
	    ELSE IF(YAR_OPERATION .EQ. 'DX')THEN
	      DO J=1,NPTS(IP)
	        IF(CD(IP)%XVEC(J) .NE. 0.0D0)THEN
	          CD(IP)%DATA(J)=CD(IP)%DATA(J)/(CD(IP)%XVEC(J)/YAR_VAL)
	        ELSE
	          CD(IP)%DATA(J)=1.0E+37
	        END IF
	      END DO
	    ELSE
	      WRITE(T_OUT,*)'Invalid operation: try again'
	      GOTO 1000
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'CAX')THEN
	  CALL NEW_GEN_IN(VAR_PLT1,'Plot that is to be altered')
	  CALL NEW_GEN_IN(VAR_OPERATION,'Operation: RX, RY, RXY, RYX, SXY') 
	  CALL SET_CASE_UP(VAR_OPERATION,IZERO,IZERO)
	  VAR_PLT2=VAR_PLT1
	  IF(VAR_OPERATION .NE. 'SXY')CALL NEW_GEN_IN(VAR_PLT2,'Plot with new data')
	  IF(NPTS(VAR_PLT1) .NE. NPTS(VAR_PLT2))THEN
	  ELSE IF(VAR_OPERATION .EQ. 'RX')THEN
	     CD(VAR_PLT1)%XVEC=CD(VAR_PLT2)%XVEC
	  ELSE IF(VAR_OPERATION .EQ. 'RY')THEN
	     CD(VAR_PLT1)%DATA=CD(VAR_PLT2)%DATA
	  ELSE IF(VAR_OPERATION .EQ. 'RXY')THEN
	     CD(VAR_PLT1)%DATA=CD(VAR_PLT2)%XVEC
	  ELSE IF(VAR_OPERATION .EQ. 'RYX')THEN
	     CD(VAR_PLT1)%XVEC=CD(VAR_PLT2)%DATA
	  ELSE IF(VAR_OPERATION .EQ. 'SXY')THEN
	     DO I=1,NPTS(VAR_PLT1)
	       T1=CD(VAR_PLT1)%DATA(I)
	       CD(VAR_PLT1)%DATA(I)=CD(VAR_PLT1)%XVEC(I)
	       CD(VAR_PLT1)%XVEC(I)=T1
	     END DO
	  END IF
!
	ELSE IF (ANS .EQ. 'VAR')THEN
	  CALL NEW_GEN_IN(VAR_PLT1,'Input plot 1?')
	  CALL NEW_GEN_IN(VAR_OPERATION,'Operation: *,+,-,/')
	  CALL SET_CASE_UP(VAR_OPERATION,IZERO,IZERO)
	  CALL NEW_GEN_IN(VAR_PLT2,'Input plot 2?')
	  VAR_PLT3=NPLTS+1
	  CALL NEW_GEN_IN(VAR_PLT3,'Output plot?')
	  TYPE_CURVE(VAR_PLT3)='L'
	  CALL DO_VEC_OP(VAR_PLT1,VAR_PLT2,VAR_PLT3,.TRUE.,VAR_OPERATION)
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'FILL')THEN
	  CALL NEW_GEN_IN(FILL,'Fill enclosed areas?')
	  IF(FILL)THEN
	    I=NFILL+1
	    CALL NEW_GEN_IN(I,'Fill identifier - default is next number')
	    IF(IFILL_PLT1(I) .EQ. 0)IFILL_PLT1(I)=1
	    IF(IFILL_PLT2(I) .EQ. 0)IFILL_PLT2(I)=2
	    CALL NEW_GEN_IN(IFILL_PLT1(I),'Input 1st plot')
	    CALL NEW_GEN_IN(IFILL_PLT2(I),'Input 2nd plot')
	    CALL NEW_GEN_IN(IFILL_CLR(I),'Color for fill region')
	    IF(I .EQ. NFILL+1)NFILL=NFILL+1
	  ELSE
	    NFILL=0
	  END IF
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'CUM')THEN
	  CALL NEW_GEN_IN(IP,'Input plot 1?')
	  OP=NPLTS+1
	  CALL NEW_GEN_IN(OP,'Output plot?')
	  TYPE_CURVE(OP)='L'
	  IF(OP .NE. IP .AND. ALLOCATED(CD(OP)%XVEC))THEN
            DEALLOCATE (CD(OP)%XVEC)
            DEALLOCATE (CD(OP)%DATA)
          END IF
          IF(.NOT. ALLOCATED(CD(OP)%XVEC))THEN
	    ALLOCATE (CD(OP)%XVEC(NPTS(IP)),STAT=IOS)
            IF(IOS .EQ. 0)ALLOCATE (CD(OP)%DATA(NPTS(IP)),STAT=IOS)
	  END IF
	  T1=CD(IP)%DATA(1)
          CD(OP)%DATA(1)=0
          DO I=2,NPTS(IP)
	    T2=CD(IP)%DATA(I)
	    CD(OP)%DATA(I)=CD(OP)%DATA(I-1)+0.5D0*(T1+T2)*
	1              ABS(CD(IP)%XVEC(I)-CD(IP)%XVEC(I-1))
	    T1=T2
	  END DO
	  IF(IP .NE. OP)CD(OP)%XVEC=CD(IP)%XVEC
	  NPTS(OP)=NPTS(IP)
	  ERR(OP)=.FALSE.
	  IF(OP .GT. NPLTS)NPLTS=OP
	  GOTO 1000
!
	ELSE IF (ANS .EQ. 'SUM')THEN
	  CALL NEW_GEN_IN(IP,'Input plot 1?')
	  OP=NPLTS+1
	  CALL NEW_GEN_IN(OP,'Output plot?')
	  TYPE_CURVE(OP)='L'
	  IF(OP .NE. IP .AND. ALLOCATED(CD(OP)%XVEC))THEN
            DEALLOCATE (CD(OP)%XVEC)
            DEALLOCATE (CD(OP)%DATA)
          END IF
          IF(.NOT. ALLOCATED(CD(OP)%XVEC))THEN
	    ALLOCATE (CD(OP)%XVEC(NPTS(IP)),STAT=IOS)
            IF(IOS .EQ. 0)ALLOCATE (CD(OP)%DATA(NPTS(IP)),STAT=IOS)
	  END IF
          CD(OP)%DATA(1)=CD(IP)%DATA(1)
          DO I=2,NPTS(IP)
	    T2=CD(IP)%DATA(I)
	    CD(OP)%DATA(I)=CD(OP)%DATA(I-1)+CD(IP)%DATA(I)
	    T1=T2
	  END DO
	  IF(IP .NE. OP)CD(OP)%XVEC=CD(IP)%XVEC
	  NPTS(OP)=NPTS(IP)
	  ERR(OP)=.FALSE.
	  IF(OP .GT. NPLTS)NPLTS=OP
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'NM')THEN
	  VAR_PLT1=1
	  CALL NEW_GEN_IN(VAR_PLT1,'Reference plot for normalization (0 to norm to 1.0')
	  IF(VAR_PLT1 .LT. 0 .OR. VAR_PLT1 .GT. NPLTS)THEN
	    WRITE(T_OUT,*)'Bad plot number'
	    GOTO 1000
	  END IF
	  IF(VAR_PLT1 .EQ. 0)WRITE(T_OUT,*)'Normalizing plots to 1.0'
	  XT(1)=XPAR(1); CALL NEW_GEN_IN(XT(1),'Beginning of normalization range')
	  XT(2)=XPAR(2)
	  CALL NEW_GEN_IN(XT(2),'End of normalization range')
	  MEAN=0.0D0
	  T2=0.0D0
	  IP=VAR_PLT1
	  IF(VAR_PLT1 .NE. 0)THEN
	    DO J=2,NPTS(VAR_PLT1)-1
	      T1=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	      IF(T1 .GT. 0)THEN
	         MEAN=MEAN+CD(IP)%DATA(J)*ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	         T2=T2+ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	      END IF
	    END DO
	    MEAN=0.5D0*MEAN; T2=0.5D0*T2
	    IF(MEAN .EQ. 0 .OR. T2 .EQ. 0)THEN
	      WRITE(T_OUT,*)'No normalization will be done'
	      WRITE(T_OUT,*)'Bad range of data'
	      GOTO 1000
	    ELSE
	      MEAN=MEAN/T2
	    END IF
	  ELSE
	      MEAN=1.0D0
	  END IF
	  DO IP=1,NPLTS
	    IF(IP .NE. VAR_PLT1 .AND. TYPE_CURVE(IP) .NE. 'I')THEN
	      T1=0.0D0; T2=0.0D0
	      DO J=2,NPTS(IP)-1
	        T3=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	        IF(T3 .GT. 0)THEN
	          T1=T1+CD(IP)%DATA(J)*ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	          T2=T2+ABS(CD(IP)%XVEC(J-1)-CD(IP)%XVEC(J+1))
	        END IF
	      END DO
	      T1=0.5D0*T1; T2=0.5D0*T2
	      IF(T2 .EQ. 0 .OR. T1 .EQ. 0)THEN
	        WRITE(T_OUT,*)IP,' not normalized'
	      ELSE
	        T1=T1/T2
	        T3=MEAN/T1
	        WRITE(T_OUT,*)'Normalization parameter for plot',IP,' is',T3
	        WRITE(LU_NORM,'(A,I3,A,ES12.4,2X,I3,2ES12.4)')'Normalization parameters for plot',IP,' are',
	1                          T3,VAR_PLT1,XT(1),XT(2)
	        DO J=1,NPTS(IP)
	          CD(IP)%DATA(J)=CD(IP)%DATA(J)*T3
	        END DO
	      END IF
	    END IF
	  END DO
	  GOTO 1000
!
	ELSE IF(ANS .EQ. 'SIG')THEN
	  XT(1:2)=XPAR(1:2)
	  CALL NEW_GEN_IN(XT,I,ITWO,'XST,XEND')
	  DO IP=1,NPLTS
	    MEAN=0.0D0; SIGMA=0.0D0; CNT=0
	    DO J=1,NPTS(IP)
	      T1=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	      IF(T1 .GT. 0)THEN
	        CNT=CNT+1
	        MEAN=MEAN+CD(IP)%DATA(J)
	      END IF
	    END DO
	    IF(CNT .GT. 1)THEN
	      MEAN=MEAN/CNT
	      DO J=1,NPTS(IP)
	        T1=(CD(IP)%XVEC(J)-XT(1))*(XT(2)-CD(IP)%XVEC(J))
	        IF(T1 .GT. 0)THEN
	          SIGMA=SIGMA+(CD(IP)%DATA(J)-MEAN)**2
	        END IF
	      END DO
	      SIGMA=SQRT( SIGMA/(CNT-1) )
	      WRITE(T_OUT,*)IP,MEAN,SIGMA,MEAN/SIGMA
	    END IF
	  END DO
	  GOTO 1000

	ELSE IF (ANS .NE. 'P' .AND. ANS .NE. 'R')THEN
	  WRITE(T_OUT,*)' Invalid option - try again'
	  GOTO 1000
	END IF
!
! Actual plotting section. All parameters except the plot size have been set.
!
5000	CONTINUE
!
! Open the appropriate plotting device
!
 6000	IF (HARD) CALL PGBBUF
	IF (HARD .AND. LONG_PLOT)THEN
	  T1=LENGTH_OF_HC_PLOT/2.54D0           !in inches
	  CALL PGPAP(T1,LP_ASR)
	END IF
!
! Initialize plotting parameters to the correct values
!
	CALL PGENV(XPAR(1),XPAR(2),YPAR(1),YPAR(2),0,-1)
	IF(HARD)THEN
	  CALL PGSCR(0,1.0,1.0,1.0)  !change background to white
          DO I=1,15                  !Make sure color are set for hard device.
            CALL PGSCR(I,RED(I),GREEN(I),BLUE(I))
          END DO
	  CALL PGSLW(PLT_LINE_WGT)
	ELSE
	  CALL PGQCR(0,RTEMP,GTEMP,BTEMP)
	  CALL PGSLW(IONE)
	END IF
	CALL PGSVP(MARGINX(1),MARGINX(2),MARGINY(1),MARGINY(2))
	CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
	CALL PGERAS
!
! Set the desired aspect ratio of the plot (useful so that
! aspect ratio is preserved on the screen and hardcopy device)
! 0 implies leave aspect ratio at default value
!
	CALL PGQVP(2,DXST,DXEND,DYST,DYEND)
	DASR=(DYEND-DYST)/(DXEND-DXST)
	TEMPASR=DASR
	IF(ASR .LT. 0)TEMPASR=-1.0/ASR
	IF(ASR .GT. 0)TEMPASR=ASR
!
	TOTXMM=(DXEND-DXST)
	TOTYMM=(DYEND-DYST)
 350	IF(HARD) CALL NEW_GEN_IN(XCM,'Plot size (in cm)')
!
! Commented out to test the box squareing routine
! Gregson Vaux December 1996
!
	IF(HARD .AND. (XCM .NE. 0)) THEN
	  SCALEFAC=XCM*10/TOTXMM
	  SCALEFACY=XCM*10*TEMPASR/TOTYMM
	ELSE IF (DASR .LE. 1.0) THEN
	  SCALEFAC=1
	  SCALEFACY=TEMPASR/DASR
          IF(DASR .LT. TEMPASR)THEN
	    T1=DASR/TEMPASR
	    SCALEFAC=SCALEFAC*T1
	    SCALEFACY=SCALEFACY*T1
	  END IF
	ELSE 
	  SCALEFAC=TEMPASR/DASR
	  SCALEFACY=1
          IF(DASR .GT. TEMPASR)THEN
	    T1=TEMPASR/DASR
	    SCALEFAC=SCALEFAC*T1
	    SCALEFACY=SCALEFACY*T1
	  END IF
	ENDIF
! 
	PRINTX1=MARGINX(1)
	PRINTX2=MARGINX(1)+(MARGINX(2)-MARGINX(1))*SCALEFAC
	PRINTY1=MARGINY(1)
	PRINTY2=MARGINY(1)+(MARGINY(2)-MARGINY(1))*SCALEFACY
        IF(HARD) THEN
	  IF(PRINTX2 .GT. 1.0 .OR. PRINTY2 .GT. 1.0)THEN
	    WRITE(T_OUT,*)' Error - plot too large - try again'
	    GOTO 350
	  ENDIF
	ENDIF
	CALL PGSVP(PRINTX1,PRINTX2,PRINTY1,PRINTY2)
!
! Fiddle with the character size to guarentee that the character size has
! the same relative size to the plot, irrespective of the plot medium.
!
	CALL PGSCH(RONE)
	CALL PGQCS(IFOUR,XCHAR_SIZE,YCHAR_SIZE)
!	T1=ABS(YPAR(2)-YPAR(1))/YCHAR_SIZE/35.0
	T1=(YPAR(2)-YPAR(1))/YCHAR_SIZE/35.0
	EXPCHAR=EXPCHAR_SCALE*T1
	EXPMARK=EXPMARK_SCALE*T1
	TICK_FAC=TICK_FAC_SCALE*T1
	CALL PGSCH(EXPCHAR)
!
! Draw Graphs
!
	IP_ST=1; IP_END=NPLTS; IP_INC=1
	IF(REVERSE_PLOTTING_ORDER)THEN
	  IP_ST=NPLTS; IP_END=1; IP_INC=-1
	END IF
!
! We do the FILL option first so that the curves are drawn on TOP.
!
	IF(FILL)THEN
	  DO I=1,NFILL
	    CALL PGSCI(IFILL_CLR(I))
	    CALL DO_FILL(XPAR,YPAR,IFILL_PLT1(I),IFILL_PLT2(I),.TRUE.)
	  END DO
	END IF
!
	DO IP=IP_ST,IP_END,IP_INC
	  CALL PGSLW(LINE_WGT(IP))
	  CALL PGSLS(LINE_STYLE(IP))
	  Q=PEN_COL(IP+PEN_OFFSET)
	  CALL PGSCI(Q)     ! START WITH COLOR INDEX 2 when PEN_OFFSET is 1 (default)
!
! Checks whether we are installing a right axis.
!
	  IF(WHICH_Y_AX(IP) .EQ. 'R')THEN
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR_R_AX(1),YPAR_R_AX(2))
	  END IF
!
	  IF(TYPE_CURVE(IP) .EQ. 'L' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL PGLINE(NPTS(IP),CD(IP)%XVEC,CD(IP)%DATA)
	  ELSE IF( (TYPE_CURVE(IP) .EQ. 'E' .OR. TYPE_CURVE(IP) .EQ. 'EC') 
	1              .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    IST=1
	    IEND=2
	    T1=CD(IP)%XVEC(NPTS(IEND))-CD(IP)%XVEC(IST)
	    Q=PEN_COL(IP+PEN_OFFSET)-1
	    DO WHILE(IEND .LT. NPTS(IP))
	      DO WHILE(IEND .LT. NPTS(IP))
	         IF( (CD(IP)%XVEC(IEND+1)-CD(IP)%XVEC(IEND))/T1 
	1                                                   .GE. 0)THEN
	           IEND=IEND+1
	         ELSE
	           EXIT
	         END IF
	      END DO
	      IF(TYPE_CURVE(IP) .EQ. 'EC')THEN
	        Q=Q+1;
	        CALL PGSCI(Q)
	      END IF
	      J=IEND-IST+1
	      CALL PGLINE(J,CD(IP)%XVEC(IST),CD(IP)%DATA(IST))
	      IST=IEND+1
	      IEND=IST+1
	      T1=CD(IP)%XVEC(NPTS(IEND))-CD(IP)%XVEC(IST)
	    END DO
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'LG')THEN
	    IST=1
	    IEND=2
	    WRITE(6,*)'Best to use negative marker style'
	    L=ABS(MARKER_STYLE(IP))
	    CALL PGSCH(EXPMARK)
	    IF(ALLOCATED(TA))THEN
	      J=SIZE(TA)
	      IF(J .LT. NPTS(IP))DEALLOCATE(TA)
	    END IF
	    IF(.NOT. ALLOCATED(TA))THEN
	      J=MAXVAL(NPTS(1:NPLTS))
	      ALLOCATE (TA(J))
	    END IF
	    TA(1:NPTS(IP))=0.0D0
	    DO I=1,NPTS(IP)
	      IF(CD(IP)%DATA(I) .NE. 0)THEN
	        TA(I)=LOG10(ABS(CD(IP)%DATA(I)))
	      ELSE
	        TA(I)=-38.0D0
	      END IF
	    END DO
	    DO WHILE(IEND .LT. NPTS(IP))
	      DO WHILE(IEND .LT. NPTS(IP))
	         IF( CD(IP)%DATA(IEND)*CD(IP)%DATA(IEND-1) .GT. 0)THEN
	           IEND=IEND+1
	         ELSE
	           EXIT
	         END IF
	      END DO
	      J=IEND-IST
	      IF(J .NE. 1 .AND. MARKER_STYLE(IP) .GT. 0)CALL PGLINE(J,CD(IP)%XVEC(IST),TA(IST))
	      DO J=IST,IEND-1
	        IF(MARK)CALL PGPT(1,CD(IP)%XVEC(J),TA(J),L)
	        IF(CD(IP)%DATA(J) .LT. 0)CALL PGPT(1,CD(IP)%XVEC(J),TA(J),24)
	      END DO
	      IST=IEND
	      IEND=IST+1
	    END DO
	    CALL PGSCH(EXPCHAR)		!Reset character size
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'H' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL PGBIN(NPTS(IP),CD(IP)%XVEC,CD(IP)%DATA,.TRUE.)
!
! This routine does a histogram plot, but the X axis is assumed
! to represent the vertices of each bin, rather than the central
! position.
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'A' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL HIST_ADJ(CD(IP)%XVEC,CD(IP)%DATA,NPTS(IP))
!
! This routine does a series of verticle lines extending from YMIN to YV(I)
! at each X value.
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'V' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    DO J=1,NPTS(IP)
	      XT(1)=CD(IP)%XVEC(J)
	      XT(2)=CD(IP)%XVEC(J)
	      YT(1)=YPAR(1)
	      YT(2)=CD(IP)%DATA(J)
	      CALL PGLINE(2,XT,YT)
	    END DO
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'VB' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL NEW_GEN_IN(VB_BASE(IP),'Base level')
	    DO J=1,NPTS(IP)
	      XT(1)=CD(IP)%XVEC(J)
	      XT(2)=CD(IP)%XVEC(J)
	      YT(1)=VB_BASE(IP)
	      YT(2)=CD(IP)%DATA(J)
	      CALL PGLINE(2,XT,YT)
	    END DO
!
! This routine does a broken plot. NPTS should be even. Lines are drawn
! from I to I+1 point for I odd only.
!
	  ELSE IF(TYPE_CURVE(IP) .EQ. 'B' .AND. (MARKER_STYLE(IP) .GE. 0 .OR. .NOT. MARK))THEN
	    CALL PGSLS(LINE_STYLE(IP))
	    Q=PEN_COL(IP+1)
	    CALL PGSCI(Q)     ! START WITH COLOR INDEX 2
	    DO J=1,NPTS(IP),2
	      XT(1)=CD(IP)%XVEC(J)
	      XT(2)=CD(IP)%XVEC(J+1)
	      YT(1)=CD(IP)%DATA(J)
	      YT(2)=CD(IP)%DATA(J+1)
	      CALL PGLINE(2,XT,YT)
	    END DO
	  END IF
!
! If this curve was plotted according to the right hand-axis, we need to
! reset the scle to the left axis.
!
	  IF(WHICH_Y_AX(IP) .EQ. 'R')THEN
	    CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
	  END IF
!
	END DO
	CALL PGSLW(PLT_LINE_WGT)
!
! Draw Gaussian fits if needed.
!
	IF(HARD .AND. DRAW_GAUSS_HARD)CALL DRAW_GAUS()
!
! Draw error bars if required.
!
	IF(DO_ERROR)THEN
	  CALL PGSLS(LINE_STYLE(1))
	  DO IP=1,NPLTS
	    IF(ERR(IP) .AND. TYPE_CURVE(IP) .NE. 'I')THEN
	      Q=PEN_COL(IP+1)
	      CALL PGSCI(Q)
	      DO J=1,NPTS(IP)
	        YT(1)=CD(IP)%EMAX(J)
	        YT(2)=CD(IP)%EMIN(J)
	        XT(1)=CD(IP)%XVEC(J)
	        XT(2)=CD(IP)%XVEC(J)
	        CALL PGLINE(2,XT,YT)
	      END DO
	    END IF
	  END DO
	END IF
!
	IF(MARK)THEN
	  CALL PGSCH(EXPMARK)
	  DO IP=1,NPLTS
	    L=ABS(MARKER_STYLE(IP))
!
! Don't draw if invisible curve.
!
	    IF(L .NE. 0 .AND. TYPE_CURVE(IP) .NE. 'I' .AND. 
	1                     TYPE_CURVE(IP) .NE. 'LG')THEN
!	      CALL PGSLS(LINE_STYLE(IP))
	      Q=PEN_COL(IP+1)
	      CALL PGSCI(Q)     ! START WITH COLOR INDEX 2
	      DO J=1,NPTS(IP)
		CALL PGPT(1,CD(IP)%XVEC(J),CD(IP)%DATA(J),L)
	      END DO
	    END IF
	  END DO
	  CALL PGSCH(EXPCHAR)		!Restore to orig status.
	END IF
!
! Draw borders.
!
	I=1
	CALL PGSLS(I)
	Q=PEN_COL(1)
	CALL PGSCI(Q)   ! Draw box in color index 1
!
! *** PGBOX is commented out while MON_BORD is being tested
! September 1996 Gregson Vaux
!
!	CALL PGBOX('ABCNT',0.0,0,'ABCNT',0.0,0)
!
	IF(DO_BORDER)THEN
	  CALL MONBORD_V3(XPAR,XINC,XNUMST,IXTICK,IDX,
	1             YPAR,YINC,YNUMST,IYTICK,IDY,
	1             TICK_FAC,EXPCHAR,
	1             XLABEL,YLABEL,TITLE,N_TITLE,TITONRHS,
	1             LOG_AXIS,OPTION,NORMAL_R_Y_AXIS)
	  IF(.NOT. NORMAL_R_Y_AXIS)THEN
	    CALL DRAW_RIGHT_Y_AXIS(YPAR_R_AX,YINC_R_AX,YNUMST_R_AX,
	1          IYTICK_R_AX,IDY_R_AX,TICK_FAC,
	1          EXPCHAR,YLABEL_R_AX,LOG_AXIS) 
	  END IF
	END IF
!
! Draw strings on graph. The parameters T1 and T2 ensure that only strings
! inside the BOX are plotted. Strings can be ploted outside the box by setting
! LOC_PG negative.
!
	IF(STR)THEN
	  WRITE(6,*)'Calling Justify'
	  CALL JUSTIFY_CONVERT_V2(XSTR,YSTR,LOC,LOC_PG,ORIENTATION,FLAGSTR,
     *    XSTRPOS,YSTRPOS,STRING,MAXSTR)
	  T1=(XSTRPOS(I)-XPAR(1))*(XPAR(2)-XSTRPOS(I))
	  T2=(YSTRPOS(I)-YPAR(1))*(YPAR(2)-YSTRPOS(I))
	  DO I=1,MAXSTR
	    IF(LOC(I) .LT. 0)THEN ; T1=1.0; T2=1.0; END IF
	    IF(FLAGSTR(I))THEN  ! .AND. T1 .GT. 0 .AND. T2 .GT. 0)THEN
	      CALL PGSCI(STR_COL(I))
	      CALL PGSCH(EXPCHAR*STR_EXP(I))
	      WRITE(6,*)'Calling PUT_TEXT'
!	        CALL PGPTXT(XSTRPOS(I),YSTRPOS(I),ORIENTATION(I),LOC_PG(I),STRING(I))
	      CALL PUT_TEXT(XSTRPOS(I),YSTRPOS(I),ORIENTATION(I),LOC_PG(I),STRING(I))
	    END IF
	  END DO
	END IF
!
	IF(N_LINE_IDS .NE. 0)THEN
	  CALL PGSCI(IONE)
	  ID_LOC=4
	  ID_ORIENT=90.0D0
	  ID_LOC_PG=1.0D0
!
	  DO I=1,N_LINE_IDS
	    WR_ID(I)=.TRUE.
	    DO J=1,N_OMIT_ID
	      IF( INDEX(LINE_ID(I),TRIM(OMIT_ID(J))//' ') .NE. 0)THEN
	        WR_ID(I)=.FALSE.
	        EXIT
	      END IF
	    END DO
	    DO J=1,N_INC_ID
	      IF(J .EQ. 1)WR_ID(I)=.FALSE.
	      IF( INDEX(LINE_ID(I),TRIM(INC_ID(J))//' ') .NE. 0 )THEN
	        WR_ID(I)=.TRUE.
	        EXIT
	      END IF
	    END DO
	  END DO
!
	  ID_WAVE_OFF(1:N_LINE_IDS)=ID_WAVE(1:N_LINE_IDS)
	  CALL PGQCS(IFOUR,XCHAR_SIZE,YCHAR_SIZE)
	  WRITE(T_OUT,*)'XCHAR_SIZE=',XCHAR_SIZE
	  DO L=1,8
	    DO I=1,N_LINE_IDS-1
	      IF(WR_ID(I))THEN
	        K=I+1
	        DO WHILE(.NOT. WR_ID(K) .AND. K .LT. N_LINE_IDS)
	          K=K+1
	        END DO
	        T1=(ID_WAVE(I)-XPAR(1))*(XPAR(2)-ID_WAVE(2))
	        IF(T1 .GT. 0)THEN
	          IF(ID_WAVE_OFF(I)+XCHAR_SIZE .GT. ID_WAVE_OFF(K))THEN
	            ID_WAVE_OFF(I)=ID_WAVE_OFF(K)-1.1*XCHAR_SIZE
	          END IF
	        END IF
	      END IF
	    END DO
	  END DO
!
	  DO I=1,N_LINE_IDS
	    T1=ID_SCL*ID_Y_OFF(I)
	    ID_VEC_BEG=0.9*(T1-1.0D0)+1.0D0
	    TMP_STR=' '
	    WRITE(TMP_STR,'(F12.2)')ID_WAVE(I)
	    TMP_STR=TRIM(LINE_ID(I))//'-'//ADJUSTL(TMP_STR)
	    IF(.NOT. OBSERVED_WAVE(I))TMP_STR='*'//TMP_STR
	    CALL JUSTIFY_CONVERT_V2(ID_WAVE_OFF(I),T1,ID_LOC,ID_LOC_PG,ID_ORIENT,.TRUE.,
	1                  XSTRPOS(1),YSTRPOS(1),TMP_STR,IONE)
	    T1=(XSTRPOS(1)-XPAR(1))*(XPAR(2)-XSTRPOS(1))
	    T2=(YSTRPOS(1)-YPAR(1))*(YPAR(2)-YSTRPOS(1))
	    IF(T1 .GT. 0 .AND. T2 .GT. 0)THEN
	      IF(WR_ID(I))THEN
	        CALL PGSCI(ITWO)
	        CALL PGSCH(EXPCHAR*ID_EXPCHAR)
	        CALL PGPTXT(XSTRPOS(1),YSTRPOS(1),ID_ORIENT,ID_LOC_PG,TMP_STR)
	        T1=ID_SCL*ID_Y_OFF(I)
	        CALL PGMOVE(ID_WAVE_OFF(I),ID_VEC_BEG)
	        CALL PGDRAW(ID_WAVE(I),ID_VEC_END)
	      END IF
	    END IF
	  END DO
	END IF
!
! Draw vectors on graphs.
!
	IF(VEC)THEN
	  DO I=1,MAXVEC
	    Q=VECPEN(I)
	    CALL PGSCI(Q)
	    IF(FLAGLINE(I))THEN
	      CALL PGMOVE(LINEXST(I),LINEYST(I))
	      CALL PGDRAW(LINEXEND(I),LINEYEND(I))
	    END IF
	  END DO
	END IF
!
! Centre text output on terminal
!
	CALL PGMOVE((XPAR(2)-XPAR(1))/2.0,(YPAR(2)-YPAR(1))/2.0)
	IF(HARD) THEN
	  CALL PGQINF('FILE',CUR_HARD_FILE,LENGTH)
	  CALL PGEBUF		!send plot to printer file
	  CALL PGCLOS
	  CALL PGSLCT(ID)
	  HARD=.FALSE.
	  IF(LONG_PLOT)CALL MODIFY_PGI_PS(CUR_HARD_FILE,LENGTH_OF_HC_PLOT,LP_ASR)
	END IF
	CALL PGSCR(0,RTEMP,GTEMP,BTEMP)
!
	GO TO 1000
	END
!
! Function to indicate when cursor input is finished. May be system dependent,
! hence function call. Mouse must be used with XWIN.
!
! On VMS 2 and 3rd mouse buttons return control to the screen.
! 1st butoon leaves control in plt window.
!
! Mouse buttons on VMS return A, D, and X.
!
! Thus use 2nd button for input.
!          3rd button to finish cursor input.
!
	LOGICAL FUNCTION END_CURS(CURSVAL)
	CHARACTER*1 CURSVAL
	END_CURS=.TRUE.
!
	IF(CURSVAL .EQ. 'A')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. 'a')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. 'D')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. 'd')THEN
          END_CURS=.FALSE.
	ELSE IF(CURSVAL .EQ. ' ')THEN
          END_CURS=.FALSE.
	END IF
!
	RETURN
	END
