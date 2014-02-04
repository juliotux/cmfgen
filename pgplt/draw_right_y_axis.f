C
C Routine to draw and label the RH Y axis. Various parameters are used 
C to provide (complete) control over numbering, tick marks, etc.
C
	SUBROUTINE DRAW_RIGHT_Y_AXIS(YPAR,YINC,YNUMST,IYTICK,IDY,
	1                  TICK_FAC,EXPCHAR,YLABEL,AX_OPT)
	IMPLICIT NONE
C
	REAL*4 YPAR(2)		!User coordinates of box
	REAL*4 YNUMST		!Numbers beginning axis labels	
	REAL*4 YINC		!Spacing between axis coordinates
	REAL*4 TICK_FAC
	REAL*4 EXPCHAR
C
	INTEGER IYTICK	!Number of ticks between coordinates
	INTEGER IDY		!Number of digits after decimal point.
C
	CHARACTER*(*) YLABEL
	CHARACTER*(*) AX_OPT
	CHARACTER*90 LOC_YLAB
C
C Local parameters
C
	REAL*4, PARAMETER :: RZERO=0.0
	REAL*4, PARAMETER :: RHALF=0.5
	REAL*4, PARAMETER :: RONE=1.0
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER PEN_YLAB
C
C Local variables.
C
	REAL*4 XBOX(4),YBOX(4)
	REAL*4 XTCK_SIZE
	REAL*4 XCHAR_SIZE,YCHAR_SIZE
	REAL*4 YLEN			!X and Y range
C
C Miscellaneous variables used locally (i.e. with in the same section of code).
C
	REAL*4 XPAR(2),YPAR_L_AXIS(2)
C
	REAL*4 X1,X2,Y1,Y2
	REAL*4 X,Y,DY,T1,ANGLE
	INTEGER IYST,I,J,NY
	INTEGER CI_SAV
C
	CHARACTER STR*4,YNUM*12
	INTEGER YNUM_LEN
	REAL*4 LOG_OFFSET(9)
C
	CHARACTER*30 UC
	EXTERNAL UC
	LOGICAL MONINSIDE
C
	CALL PGQWIN(XPAR(1),XPAR(2),YPAR_L_AXIS(1),YPAR_L_AXIS(2))
	CALL PGSWIN(XPAR(1),XPAR(2),YPAR(1),YPAR(2))
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
 	CALL PGQTXT(RZERO,RZERO,RZERO,RZERO,'X',XBOX,YBOX)
	CALL PGQVP(2,X1,X2,Y1,Y2)
C
	CALL PGQTXT(RZERO,RZERO,RZERO,RZERO,'XXX',XBOX,YBOX)
	XCHAR_SIZE=(XBOX(4)-XBOX(2))/3.0D0
	YCHAR_SIZE=(YBOX(2)-YBOX(4))
	YLEN=YPAR(2)-YPAR(1)
C
	LOC_YLAB=YLABEL
	I=INDEX(LOC_YLAB,'\p')
	IF(I .EQ. 0)I=INDEX(LOC_YLAB,'\P')
	IF(I .NE. 0)THEN
	  LOC_YLAB(I:)=' '
	  READ(YLABEL(I+2:),*)PEN_YLAB
	END IF
!
C 
C
C Set parameters for drawing tick marks and putting on coordinates depending 
C on the choice of Axes chosen. 
C
	IF(AX_OPT .EQ. 'LOGY' .OR. AX_OPT .EQ. 'LOGXY')THEN
	  YINC=1.0
	  YNUMST= INT( YNUMST+0.00001*SIGN(1.0,YNUMST) )
	  NY=INT((YPAR(2)-YNUMST+0.00001*YINC)/YINC) + 1
	ELSE
	  NY=INT(YLEN/YINC*FLOAT(IYTICK)+1.0E-06)
	  IYST=INT( (YNUMST-YPAR(1))*1.0001/ (YINC/IYTICK) )
	END IF
C
C If we are switching to LOG axes with full numbers marked we strip Log from 
C the label.
C
	IF(INDEX(AX_OPT,'Y') .NE. 0)THEN
	  IF(UC(LOC_YLAB(1:3)) .EQ. 'LOG')THEN
	    IF(UC(LOC_YLAB(1:4)) .EQ. 'LOG(')THEN
	      J=LEN_TRIM(LOC_YLAB)
	      LOC_YLAB=LOC_YLAB(5:J-1)
	    ELSE
	      LOC_YLAB=LOC_YLAB(4:)
	    END IF
	  END IF
	ELSE
	  LOC_YLAB=LOC_YLAB
	END IF
C
C 
C
	CALL PGSCI(IONE)			
C
C Tick marks and scale: Right edge and scale
C
	IF(AX_OPT .EQ. 'LOGY' .OR. AX_OPT .EQ. 'LOGXY')THEN
	  DO I=0,NY
	    DO J=1,9
	      Y=YNUMST+(I-1)*YINC+LOG_OFFSET(J)
	      IF( MONINSIDE(Y,YPAR) )THEN
	        CALL PGMOVE(XPAR(2),Y)
	        T1=XTCK_SIZE*(1.0-0.33*((J+7)/9))
	        CALL PGDRAW(XPAR(2)-T1,Y)
	      END IF
	    END DO
	  END DO
C
C Right Axis labelling.
C
	  X=XPAR(2)+XTCK_SIZE
	  DO I=1,NY
	    Y=YNUMST+(I-1)*YINC-0.4*YCHAR_SIZE
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
	    CALL PGPTXT(X,Y,ANGLE,RZERO,YNUM(1:YNUM_LEN))
	  END DO
	ELSE
	  DO I=0,NY
	    Y=YNUMST+(I-IYST)*YINC/IYTICK
	    IF(MOD(I-IYST,IYTICK) .EQ. 0)THEN
	      IF(MONINSIDE(Y,YPAR))THEN
	        CALL PGMOVE(XPAR(2),Y)
	        CALL PGDRAW(XPAR(2)-XTCK_SIZE,Y)
	      END IF
	      X=XPAR(2)+XTCK_SIZE
	      CALL MON_NUM(X,Y-(YCHAR_SIZE/2.0),Y,6,IDY)
	    ELSE
	      IF(MONINSIDE(Y,YPAR))THEN
	        CALL PGMOVE(XPAR(2),Y)
	        CALL PGDRAW(XPAR(2)-0.67*XTCK_SIZE,Y)
	      END IF
	    END IF
	  END DO
	END IF
C
C 
C
C Y Axis label
C
	Y=ABS(YPAR(1))
	IF(Y .LT. ABS(YPAR(2)))Y=ABS(YPAR(2))
	IF(Y .GT. 1.0)THEN
	  DY=AINT(ALOG10(AINT(Y)))+1.0
	ELSE
	  DY=1.0
	END IF
	IF(YPAR(1) .LT. 0 .OR. YPAR(2) .LT. 0)DY=DY+1	!Allow for minus sign
	IF(INDEX(AX_OPT,'Y') .NE. 0)THEN
	  DY=2.0
	  X=XPAR(2)+(DY+0.5)*XCHAR_SIZE+2.0*XTCK_SIZE
	  IF(YPAR(1) .LT. 0.1)DY=DY+1
	ELSE
	  X=XPAR(2)+(DY+IDY+0.5)*XCHAR_SIZE+2.0*XTCK_SIZE
	END IF
	Y=YPAR(1)+0.5*YLEN
	ANGLE=270.0
!
	CALL PGQCI(CI_SAV)
	IF(PEN_YLAB .NE. 0)CALL PGSCI(PEN_YLAB)
	CALL PGPTXT(X,Y,ANGLE,RHALF,LOC_YLAB)
	CALL PGSCI(CI_SAV)
C
	CALL PGQWIN(XPAR(1),XPAR(2),YPAR_L_AXIS(1),YPAR_L_AXIS(2))
C
	RETURN
	END
