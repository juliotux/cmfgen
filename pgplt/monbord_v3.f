C
C Routine to draw and label axes. Various parameters are used to provide
C  (complete) control over numbering, tick marks, etc.
C
	SUBROUTINE MONBORD_V3(XPAR,XINC,XNUMST,IXTICK,IDX,
	1                  YPAR,YINC,YNUMST,IYTICK,IDY,
	1                  TICK_FAC,EXPCHAR,
	1                  XLABEL,YLABEL,TITLE,N_TITLE,
	1                  TITONRHS,AX_OPT,TOP_AX_OPT,NORMAL_R_AXIS)
C
	IMPLICIT NONE
	REAL*4 FNTICK
!
! Altered 17-Feb-2015 : Can now have multi-colored titles
! Altered 13-Feb-2001 : Titles now correctly treated for flipped axes.
!                       Log axex now handled for flipped axes.
! Altered 31-May-2000 : Y label now positioned correctly for Exponential format.
! Altered 05-Aug-1999 : Multiple title installed. Increased from 2.
!                         As call changed, changed to V3.
!
C Altered 05-Nov-1997: Option not to draw axis included.
C Altered 30-Jul-1997: 10^1 notation along Y axis now written parallel to
C                        X axis.
C Altered  7-Mar-1987; Cleaned. AX_OPT reinstalled.
C
C Altered September 1996 by Gregson Vaux to convert from MONGO to PGPLOT
C         Note: This subroutine is used instead of PGBOX so that the
C               user will have more control over the tick marks and other
C               axis properties
C Altered 12-Apr-1989 - MONGOPAR coomon block puto include statement.
C                       DRAW, RELOCATE, PUTLABEL and GRELOCATE prefixed by MGO.
C Altered 24-Oct-88 - MON_BOX changed to MON_DBOX. Plots box in device
C                     coordinates to avoid roundoff(?) which cause loss
C                     of end axis.
C Altered 21-May-85 - Locating of Y axis label improved and bug fixed.
C Altered 13-May-85 (FNTCK corrected and end axis bars are no longer plotted)
C
	REAL*4 XPAR(2),YPAR(2)		!User coordinates of box
	REAL*4 XNUMST,YNUMST		!Numbers beginning axis labels	
	REAL*4 XINC,YINC		!Spacing between axis coordinates
	REAL*4 TICK_FAC
	REAL*4 EXPCHAR
	REAL*8 LST_X_VAL
C
	INTEGER IXTICK,IYTICK		!Number of ticks between coordinates
	INTEGER IDX,IDY		!Number of digits after decimal point.
C
	LOGICAL TITONRHS		!Title on RHS of plot?
	LOGICAL NORMAL_R_AXIS
C
	INTEGER N_TITLE
	CHARACTER(LEN=*) XLABEL,YLABEL
	CHARACTER(LEN=*) TITLE(N_TITLE)
	CHARACTER(LEN=*) AX_OPT,TOP_AX_OPT
C
C Local paremeters
C
	REAL*4, PARAMETER :: RZERO=0.0
	REAL*4, PARAMETER :: RHALF=0.5
	REAL*4, PARAMETER :: RONE=1.0
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
C
C Local variables.
C
	REAL*4 XBOX(4),YBOX(4)
	REAL*4 XTCK_SIZE,YTCK_SIZE
	REAL*4 XCHAR_SIZE,YCHAR_SIZE
	INTEGER IXST,IYST
C
	REAL*4 XLEN,YLEN			!X and Y range
C
C Miscellaneous variables used locally (i.e. with in the same section of code).
C
	REAL*4 X1,X2,Y1,Y2
	REAL*4 X,Y,DY,T1,ANGLE
	INTEGER I,J,NX,NY
	INTEGER CI_SAV
C
	CHARACTER(LEN=4) STR
	CHARACTER(LEN=12) XNUM,YNUM
	INTEGER XNUM_LEN,YNUM_LEN
	REAL*4 LOG_OFFSET(9)
C
	CHARACTER(LEN=30) UC
	EXTERNAL UC
!
! For printing of titles.
!
	INTEGER PEN_TIT(N_TITLE)
	REAL*4 LEN_TIT(N_TITLE)
	REAL*4 HT_TIT(N_TITLE)
	CHARACTER(LEN=80) LOC_TIT(N_TITLE)
	INTEGER PEN_XLAB
	INTEGER PEN_YLAB
!
	CHARACTER(LEN=80) LOC_XLAB,LOC_YLAB
	LOGICAL MONINSIDE
C
C If TOP_AX_OPT='TOPLAB' a different scale (nonlinear) can be put on the top
C axis, as specifed by SCED, and ED. 
C
	REAL*8 SCED,XED
	INTEGER NXED
	CHARACTER(LEN=30) TOPLABEL
	COMMON /TOPBORD/ SCED(31),XED(31),NXED,TOPLABEL
C
C Two options are available with  logarithmic axis. The default is to
C label the axes with the actual Log values. the alterenate is to
C draw logarithmically space tick marks, and label aixs by actual numbers.
C This is done by specifying AX_OPT=LOGX, LOGY LOGXY. The range in X
C (or Y) must be at least a factor of 10 for this option.
C
        DO I=1,9
	  LOG_OFFSET(I)=LOG10( FLOAT(I) )
	END DO
C
C Estimate appropriate tick sizes (Default is tick same size as 
C character with expand=1)
C
 	CALL PGQTXT(RZERO,RZERO,RZERO,RZERO,'X',XBOX,YBOX)
	XTCK_SIZE=(XBOX(4)-XBOX(2))*TICK_FAC/EXPCHAR
	CALL PGQVP(2,X1,X2,Y1,Y2)
	YTCK_SIZE=XTCK_SIZE*(X2-X1)/(XPAR(2)-XPAR(1))*
	1                  (YPAR(2)-YPAR(1))/(Y2-Y1)
C
C It is possible using a \pn (n=0,1...) on the end of the TITLE to specify the
C title color.
C     
	LOC_TIT(1:N_TITLE)=TITLE(1:N_TITLE)
	PEN_TIT(1:N_TITLE)=0
	PEN_XLAB=0
	PEN_YLAB=0
!
! Get pen color for each title.
!
!	DO J=1,N_TITLE
!	  I=INDEX(TITLE(J),'\p')
!	  IF(I .EQ. 0)I=INDEX(TITLE(J),'\P')
!	  IF(I .NE. 0)THEN
!	   LOC_TIT(J)(I:)=' '
!	   READ(TITLE(J)(I+2:),*)PEN_TIT(J)
!	  END IF
!	END DO
	CALL STRIP_SLASH_P(LOC_TIT,N_TITLE)
!
! Get pen or axis labels
!
	LOC_XLAB=XLABEL
	I=INDEX(LOC_XLAB,'\p')
	IF(I .EQ. 0)I=INDEX(LOC_XLAB,'\P')
	IF(I .NE. 0)THEN
	  LOC_XLAB(I:)=' '
	  READ(XLABEL(I+2:),*)PEN_XLAB
	END IF
!
	LOC_YLAB=YLABEL
	I=INDEX(LOC_YLAB,'\p')
	IF(I .EQ. 0)I=INDEX(LOC_YLAB,'\P')
	IF(I .NE. 0)THEN
	  LOC_YLAB(I:)=' '
	  READ(YLABEL(I+2:),*)PEN_YLAB
	END IF
!
	CALL PGQTXT(RZERO,RZERO,RZERO,RZERO,'XXX',XBOX,YBOX)
	XCHAR_SIZE=(XBOX(4)-XBOX(2))/3
	YCHAR_SIZE=(YBOX(2)-YBOX(4))
C
	XLEN=XPAR(2)-XPAR(1)
	YLEN=YPAR(2)-YPAR(1)
C
C 
C
C Set parameters for drawing tick marks and putting on coordinates depending 
C on the choice of Axes chosen. 
C
	IF(AX_OPT .EQ. 'LOGX')THEN
	  XINC=1.0
	  IF(XPAR(2) .LT. XPAR(1))XINC=-1.0
	  XNUMST= INT( XNUMST+0.00001*SIGN(1.0,XNUMST) )
	  NX=INT((XPAR(2)-XNUMST+0.00001*XINC)/XINC)+1
	  NY=INT(YLEN/YINC*FLOAT(IYTICK)+1.0E-06)
	  IYST=INT( (YNUMST-YPAR(1))*1.0001/ (YINC/IYTICK) )
	ELSE IF(AX_OPT .EQ. 'LOGY')THEN
	  YINC=1.0
	  IF(YPAR(2) .LT. YPAR(1))YINC=-1.0
	  YNUMST= INT( YNUMST+0.00001*SIGN(1.0,YNUMST) )
	  NY=INT((YPAR(2)-YNUMST+0.00001*YINC)/YINC)+1
	  NX=INT(XLEN/XINC*FLOAT(IXTICK)+1.0E-06)
	  IXST=INT( (XNUMST-XPAR(1))*1.0001/ (XINC/IXTICK) )
	ELSE IF(AX_OPT .EQ. 'LOGXY')THEN
	  XINC=1.0
	  YINC=1.0
	  IF(XPAR(2) .LT. XPAR(1))XINC=-1.0
	  IF(YPAR(2) .LT. YPAR(1))YINC=-1.0
C
C Check that XNUMST and YNUMST are multiples of XINC and YINC
C
	  XNUMST= INT( XNUMST+0.00001*SIGN(1.0,XNUMST) )
	  YNUMST= INT( YNUMST+0.00001*SIGN(1.0,YNUMST) )
	  NY=INT((YPAR(2)-YNUMST+0.00001*YINC)/YINC) + 1
	  NX=INT((XPAR(2)-XNUMST+0.00001*XINC)/XINC) + 1
	ELSE
	  NY=INT(YLEN/YINC*FLOAT(IYTICK)+1.0E-06)
	  NX=INT(XLEN/XINC*FLOAT(IXTICK)+1.0E-06)
	  IXST=INT( (XNUMST-XPAR(1))*1.0001/ (XINC/IXTICK) )
	  IYST=INT( (YNUMST-YPAR(1))*1.0001/ (YINC/IYTICK) )
	END IF
C
C If we are switching to LOG axes with full numbers marked we strip Log from 
C the label.
C
	IF(INDEX(AX_OPT,'X') .NE. 0)THEN
	  IF(UC(LOC_XLAB(1:3)) .EQ. 'LOG')THEN
	    IF(UC(LOC_XLAB(1:4)) .EQ. 'LOG(')THEN
	      J=LEN_TRIM(LOC_XLAB)
	      LOC_XLAB=LOC_XLAB(5:J-1)
	    ELSE
	      LOC_XLAB=LOC_XLAB(4:)
	    END IF
	  END IF
	END IF	    
	IF(INDEX(AX_OPT,'Y') .NE. 0)THEN
	  IF(UC(LOC_YLAB(1:3)) .EQ. 'LOG')THEN
	    IF(UC(LOC_YLAB(1:4)) .EQ. 'LOG(')THEN
	      J=LEN_TRIM(LOC_YLAB)
	      LOC_YLAB=LOC_YLAB(5:J-1)
	    ELSE
	      LOC_YLAB=LOC_YLAB(4:)
	    END IF
	  END IF
	END IF
C
C 
C
C Draw a box without labels.
C
! Assume border color is set outside of MONBORD_V3
!
!	CALL PGSCI(IONE)			!Set default color of axes
	CALL PGBOX('BC',RZERO,IZERO,'BC',RZERO,IZERO)
C
C Tick marks and scale
C
C Right edge
C
	IF( NORMAL_R_AXIS .AND.
	1         (AX_OPT .EQ. 'LOGY' .OR. AX_OPT .EQ. 'LOGXY') )THEN
	  DO I=0,NY
	    DO J=1,9
	      Y=YNUMST+(I-1)*YINC+LOG_OFFSET(J)
	      IF( MONINSIDE(Y,YPAR) )THEN
	        CALL PGMOVE(XPAR(2),Y)
	        T1=XTCK_SIZE*( 1.0-0.33*((J+7)/9) )
	        CALL PGDRAW(XPAR(2)-T1,Y)
	      END IF
	    END DO
	  END DO
	ELSE IF(NORMAL_R_AXIS)THEN
	  DO I=0,NY
	    Y=YNUMST+(I-IYST)*YINC/IYTICK
	    IF( MONINSIDE(Y,YPAR) )THEN
	      CALL PGMOVE(XPAR(2),Y)
	      CALL PGDRAW(XPAR(2)-FNTICK(XTCK_SIZE,IYTICK,I-IYST),Y)
	    END IF
	  END DO
	END IF
C
C 
C
C Top axis.
C
C Check if the TOPLABEL option is in effect, before drawing
C a the top scale.
C
	I=LEN(TOP_AX_OPT)
	IF(I .GT. 6)I=6
	IF(AX_OPT .EQ. 'LOGX' .OR. AX_OPT .EQ. 'LOGXY')THEN
	  DO I=0,NX
	    DO J=1,9
	      X=XNUMST+(I-1)*XINC+LOG_OFFSET(J)
	      IF( MONINSIDE(X,XPAR) )THEN
	        CALL PGMOVE(X,YPAR(2))
	        T1=YTCK_SIZE*(1.0-0.33*((J+7)/9))
	        CALL PGDRAW(X,YPAR(2)-T1)
	      END IF
	    END DO
	  END DO
	ELSE IF(TOP_AX_OPT(1:I) .EQ. 'TOPLAB')THEN
          X=(XPAR(1)+XPAR(2))/2.0
	  Y=YPAR(2)+2.0*YTCK_SIZE+YCHAR_SIZE
	  CALL PGPTXT(X,Y,RZERO,RHALF,TOPLABEL)
	  LST_X_VAL=-2.0D0*XPAR(1)
	  DO I=1,NXED
	    X=XED(I)
	    IF( (XPAR(1) .LE. X .AND. XPAR(2) .GE. X) .OR.
	1         (XPAR(1) .GE. X .AND. XPAR(2) .LE. X) )THEN
	      IF(MOD(I,2) .NE. 0)THEN
	        Y=SCED(I)
	        CALL PGMOVE(X,YPAR(2))
	        CALL PGDRAW(X,YPAR(2)-YTCK_SIZE)
	        IF(ABS(X-LST_X_VAL) .GT. ABS(XPAR(2)-XPAR(1))/20.0D0)THEN
	          CALL MON_NUM(X,YPAR(2)+YTCK_SIZE,Y,5,IZERO)               !IDX)
	          LST_X_VAL=X
	        END IF
              ELSE
	        CALL PGMOVE(X,YPAR(2))
	        CALL PGDRAW(X,YPAR(2)-0.67*YTCK_SIZE)
	      END IF
	    END IF
	  END DO
	ELSE
C
C Normal top edge.
C
	  DO I=0,NX
	    X=XNUMST+(I-IXST)*XINC/IXTICK
	    IF( MONINSIDE(X,XPAR) )THEN
	      CALL PGMOVE(X,YPAR(2))
	      CALL PGDRAW(X,YPAR(2)-FNTICK(YTCK_SIZE,IXTICK,I-IXST))
	    END IF
	  END DO
	END IF
C
C 
C Left edge and scale
C
	IF(AX_OPT .EQ. 'LOGY' .OR. AX_OPT .EQ. 'LOGXY')THEN
	  DO I=0,NY
	    DO J=1,9
	      Y=YNUMST+(I-1)*YINC+LOG_OFFSET(J)
	      IF( MONINSIDE(Y,YPAR) )THEN
	        CALL PGMOVE(XPAR(1),Y)
	        T1=XTCK_SIZE*(1.0-0.33*((J+7)/9))
	        CALL PGDRAW(XPAR(1)+T1,Y)
	      END IF
	    END DO
	  END DO
C
C Left Axis labelling.
C
	  X=XPAR(1)-XTCK_SIZE
	  DO I=1,NY
	    Y=YNUMST+(I-1)*YINC-0.4*YCHAR_SIZE
	    IF( (Y-YPAR(1))*(YPAR(2)-Y) .GT. 1.0D-04*ABS(YPAR(2)-YPAR(1)))THEN
	      WRITE(STR,'(I4)')NINT(Y)
	      YNUM='10\u'
	      YNUM_LEN=4
	      DO J=1,4
	        IF(STR(J:J) .NE. ' ')THEN
                  YNUM=YNUM(1:YNUM_LEN)//STR(J:J)
	          YNUM_LEN=YNUM_LEN+1
	        END IF
	      END DO
              YNUM=YNUM(1:YNUM_LEN)//'\d'
	      YNUM_LEN=YNUM_LEN+2
	      CALL PGMOVE(X,Y)
	      ANGLE=0.0
	      CALL PGPTXT(X,Y,ANGLE,RONE,YNUM(1:YNUM_LEN))
	    END IF
	  END DO
	ELSE
	  DO I=0,NY
	    Y=YNUMST+(I-IYST)*YINC/IYTICK
	    IF(MOD(I-IYST,IYTICK) .EQ. 0)THEN
	      IF(MONINSIDE(Y,YPAR))THEN
	        CALL PGMOVE(XPAR(1),Y)
	        CALL PGDRAW(XPAR(1)+XTCK_SIZE,Y)
	      END IF
	      X=XPAR(1)-XTCK_SIZE
	      CALL MON_NUM(X,Y-(YCHAR_SIZE/2.0),Y,4,IDY)
	    ELSE
	      IF(MONINSIDE(Y,YPAR))THEN
	        CALL PGMOVE(XPAR(1),Y)
	        CALL PGDRAW(XPAR(1)+0.67*XTCK_SIZE,Y)
	      END IF
	    END IF
	  END DO
	END IF
C
C 
C
C Bottom Scale.
C
	IF(AX_OPT .EQ. 'LOGX' .OR. AX_OPT .EQ. 'LOGXY')THEN
	  DO I=0,NX
	    DO J=1,9
	      X=XNUMST+(I-1)*XINC+LOG_OFFSET(J)
	      IF( MONINSIDE(X,XPAR) )THEN
	        CALL PGMOVE(X,YPAR(1))
	        T1=YTCK_SIZE*(1.0-0.33*((J+7)/9))
	        CALL PGDRAW(X,YPAR(1)+T1)
	      END IF
	    END DO
	  END DO
C
C Bottom axis labels.
C
	  Y=YPAR(1)-YTCK_SIZE-1.25*YCHAR_SIZE
	  DO I=1,NX
	    X=XNUMST+(I-1)*XINC
	    XNUM='10\u'
	    XNUM_LEN=4
	    WRITE(STR,'(I4)')NINT(X)
	    DO J=1,4
	      IF(STR(J:J) .NE. ' ')THEN
                XNUM=XNUM(1:XNUM_LEN)//STR(J:J)
	        XNUM_LEN=XNUM_LEN+1
	      END IF
	    END DO
            XNUM=XNUM(1:XNUM_LEN)//'\d'
	    XNUM_LEN=XNUM_LEN+2
	    IF( MONINSIDE(X,XPAR) )THEN
	      CALL PGMOVE(X,Y)
	      CALL PGPTXT(X,Y,RZERO,RHALF,XNUM(1:XNUM_LEN))
	    END IF
	  END DO
	ELSE
C
C Normal bottom edge and scale
C
	  DO I=0,NX
	    X=XNUMST+(I-IXST)*XINC/IXTICK
	    IF( MOD(I-IXST,IXTICK) .EQ. 0 )THEN
	      IF( MONINSIDE(X,XPAR) )THEN 
	        CALL PGMOVE(X,YPAR(1))
	        CALL PGDRAW(X,YPAR(1)+YTCK_SIZE)
	      END IF
	      Y=YPAR(1)-YTCK_SIZE
	      CALL MON_NUM(X,Y-(YCHAR_SIZE),X,2,IDX)
	    ELSE
	      IF( MONINSIDE(X,XPAR) )THEN 
	        CALL PGMOVE(X,YPAR(1))
	        CALL PGDRAW(X,YPAR(1)+0.67*YTCK_SIZE)
	      END IF
	    END IF
	  END DO
	END IF
C
C 
C
C Write in titles
C
C Get bounding box for each title.
C	
	DO J=1,N_TITLE
	  CALL PGQTXT(RZERO,RZERO,RZERO,RZERO,TRIM(LOC_TIT(J)),XBOX,YBOX)
	  HT_TIT(J)=(YBOX(3)-YBOX(1))
	  LEN_TIT(J)=(XBOX(4)-XBOX(2))
	END DO
C
	CALL PGQCI(CI_SAV)
	ANGLE=0.0
	IF(TITONRHS)THEN
	  IF(XPAR(1) .GT. XPAR(2))THEN
	    X=XPAR(2)-1*XTCK_SIZE-MINVAL(LEN_TIT)
	  ELSE
	    X=XPAR(2)-1*XTCK_SIZE-MAXVAL(LEN_TIT)
	  END IF
	  Y=YPAR(2)-1.5*YTCK_SIZE
	  DO J=1,N_TITLE
	    IF(LEN_TIT(J) .NE. 0)THEN
	      Y=Y-0.5*HT_TIT(J)
	      IF(PEN_TIT(J) .NE. 0)CALL PGSCI(PEN_TIT(J))
!	      CALL PGPTXT(X,Y,ANGLE,RZERO,LOC_TIT(J))
	      CALL PUT_TEXT(X,Y,ANGLE,RZERO,TITLE(J))
	      IF(J .NE. N_TITLE)Y=Y-0.7*HT_TIT(J+1)
	    END IF
	  END DO
	ELSE
	  X=XPAR(1)+1.2*XTCK_SIZE
	  Y=YPAR(2)-1.5*YTCK_SIZE
	  DO J=1,N_TITLE
	    IF(LEN_TIT(J) .NE. 0)THEN
	      Y=Y-0.5*HT_TIT(J)
!	      IF(PEN_TIT(J) .NE. 0)CALL PGSCI(PEN_TIT(J))
!	      CALL PGPTXT (X,Y,ANGLE,RZERO,LOC_TIT(J))
	      CALL PUT_TEXT (X,Y,ANGLE,RZERO,TITLE(J))
	      IF(J .NE. N_TITLE)Y=Y-0.7*HT_TIT(J+1)
	    END IF
	  END DO
	END IF
	CALL PGSCI(CI_SAV)
C
C X Axis label
C
	X=XPAR(1)+XLEN/2.0
	Y=YPAR(1)-2.0*YTCK_SIZE-(2.0*YCHAR_SIZE)
	CALL PGQCI(CI_SAV)
	IF(PEN_XLAB .NE. 0)CALL PGSCI(PEN_XLAB)
	CALL PGPTXT(X,Y,RZERO,RHALF,LOC_XLAB)
	CALL PGSCI(CI_SAV)
C
C Y Axis label
C
	Y=ABS(YPAR(1))
	IF(Y .LT. ABS(YPAR(2)))Y=ABS(YPAR(2))
	IF(Y .GT. 1.0E+05)THEN
	  DY=5
	ELSE IF(Y .LT. 1.0E-05)THEN
	  DY=5
	ELSE IF(Y .GT. 1.0)THEN
	  DY=AINT(ALOG10(AINT(Y)))+1.0
	ELSE
	  DY=1.0
	END IF
	IF(YPAR(1) .LT. 0 .OR. YPAR(2) .LT. 0)DY=DY+1	!Allow for minus sign
	X=XPAR(1)-(DY+IDY+0.5)*XCHAR_SIZE-2.0*XTCK_SIZE
	Y=YPAR(1)+0.5*YLEN
	ANGLE=90.0
	CALL PGQCI(CI_SAV)
	IF(PEN_YLAB .NE. 0)CALL PGSCI(PEN_YLAB)
	CALL PGPTXT(X,Y,ANGLE,RHALF,LOC_YLAB)
	CALL PGSCI(CI_SAV)
!
	RETURN
	END
