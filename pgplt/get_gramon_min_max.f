!
! Get aaaaabsica and ordinate maxima and minim for plotting.
!
	SUBROUTINE GET_GRAMON_MIN_MAX(XMIN,XMAX,YMIN,YMAX,TYPE_CURVE,T_OUT)
	USE MOD_CURVE_DATA
!
! Created 06-Sep-2005
!
	REAL*4 XMIN,XMAX,YMIN,YMAX
	INTEGER T_OUT
	CHARACTER*1 TYPE_CURVE(MAX_PLTS)
!
	REAL T1
	INTEGER IP
!
! Look for absica limits
!
	XMIN=CD(1)%XVEC(1)
	XMAX=CD(1)%XVEC(1)
	DO IP=1,NPLTS
	  IF(TYPE_CURVE(IP) .NE. 'I')THEN
	    T1=MINVAl(CD(IP)%XVEC)
	    XMIN=MIN(XMIN,T1)
	    T1=MAXVAL(CD(IP)%XVEC)
	    XMAX=MAX(XMAX,T1)
	  END IF
	END DO
!
! Look for ordinate limits
!
	YMIN=CD(1)%DATA(1)
	YMAX=CD(1)%DATA(1)
	DO IP=1,NPLTS
	  IF(TYPE_CURVE(IP) .NE. 'I')THEN
	    T1=MINVAL(CD(IP)%DATA)
	    YMIN=MIN(YMIN,T1)
	    T1=MAXVAL(CD(IP)%DATA)
	    YMAX=MAX(YMAX,T1)
	  END IF
	END DO
!
! Check range of validity
!
	IF(YMAX .EQ. 0 .AND. YMIN .EQ. 0)THEN
	  YMIN=0.0
	  YMAX=1.0
	  WRITE(T_OUT,*)'Y limits are zero - setting default values'
	ELSE IF(ABS(YMAX-YMIN)/MAX(ABS(YMAX),ABS(YMIN)) .LT. 1.0D-08)THEN
	  YMIN=0.0
	  YMAX=1.0
	  WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	END IF
	IF(XMAX .EQ. 0 .AND. XMIN .EQ. 0)THEN
	  XMIN=0.0
	  XMAX=1.0
	  WRITE(T_OUT,*)'X limits are zero - setting default values'
	ELSE IF(ABS(XMAX-XMIN)/MAX(ABS(XMAX),ABS(XMIN)) .LT. 1.0D-08)THEN
	  XMIN=0.0
	  XMAX=1.0
	  WRITE(T_OUT,*)'Invalid X limits - setting default values'
	END IF
!
	RETURN
	END
