!
! Get absica and ordinate maxima and minima for plotting.
!
	SUBROUTINE GET_GRAMON_MIN_MAX(XMIN,XMAX,YMIN,YMAX,TYPE_CURVE,T_OUT)
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
! Altered 24-Feb-2016: Bug fix - no longer use plot 1 to set initial
!                           values of XMIN and XMAX. 
! Altered 12-Jan-2014: Installed check to see if LG option.
! Created 06-Sep-2005
!
	REAL*4 XMIN,XMAX,YMIN,YMAX
	INTEGER T_OUT
	CHARACTER(LEN=*) TYPE_CURVE(MAX_PLTS)
!
	REAL*4 T1,T2
	INTEGER IP
!
! Look for absica limits
!
	XMIN=0.5D0*HUGE(XMIN)
	XMAX=-XMIN
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
	YMIN=0.5D0*HUGE(YMIN)
	YMAX=-YMIN
	DO IP=1,NPLTS
	  IF(TYPE_CURVE(IP) .EQ. 'LG')THEN
	    T1=MAXVAL(CD(IP)%DATA)
	    T2=-MINVAL(CD(IP)%DATA)
	    YMAX=MAX(T1,T2)
	    YMIN=MINVAL(ABS(CD(IP)%DATA))
	    IF(T1 .EQ. 0.0D0)YMIN=-38.0D0
	  ELSE IF(TYPE_CURVE(IP) .NE. 'I')THEN
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
	  YMIN=0.9*YMIN
	  YMAX=1.1*YMAX
	  WRITE(T_OUT,*)'Invalid Y limits - setting default values'
	END IF
	IF(XMAX .EQ. 0 .AND. XMIN .EQ. 0)THEN
	  XMIN=0.0
	  XMAX=1.0
	  WRITE(T_OUT,*)'X limits are zero - setting default values'
	ELSE IF(ABS(XMAX-XMIN)/MAX(ABS(XMAX),ABS(XMIN)) .LT. 1.0D-08)THEN
	  XMIN=0.9*XMIN
	  XMAX=1.1*XMAX
	  WRITE(T_OUT,*)'Invalid X limits - setting default values'
	END IF
!
	RETURN
	END
