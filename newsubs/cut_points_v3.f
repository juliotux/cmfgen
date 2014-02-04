!===================================================================
!
	SUBROUTINE CUT_POINTS_V3(NU_CUT,CROSS_CUT,NUM_CUT,
	1    NU_SM,CROSS_SM,NUM_SM,CUT_ACCURACY)
!
!===================================================================
!
! This routine cuts points from the given smoothed cross section such that
! linear interpolation can be used to recover intermediate data points to
! a fractional accuracy of CUT_ACCURACY. All local maximum and minimum points 
! are retained.
!
! Altered 20-Apr-2004 : Work arrays passed in call
!                       For use with sm_phot_v3.f
!                       Changed to V3
! Altered  11/Jun/1999 Criterion for cutting points modified.
!                        Required accuracy now passed in call.
! Altered  2/16/96 DLM Change minimum spacing from 5% of total range to
!                      at a minimum of "min_points" in equal steps
!                      of log(Energy)
! Created  9/28/95 DLM
!
      IMPLICIT NONE
!
! On entry *_sm contain the original number of points and data.
! On exit, *_cut contain the new number of data points, and the selected data.
! Both the *_cut and *_sm arrays are assumed to have a maximum length of NUM_SM.
!
      INTEGER NUM_SM
      REAL*8 NU_SM(NUM_SM),CROSS_SM(NUM_SM)
      REAL*8 CUT_ACCURACY
!
      REAL*8 NU_CUT(NUM_SM)
      REAL*8 CROSS_CUT(NUM_SM)
      REAL*8 DERIV(NUM_SM)
      REAL*8 DIST(NUM_SM)
!
      INTEGER I,LOW,HIGH,MID,NUM_CUT,NUM_AREA
      REAL*8 X,Y,M
      REAL*8 COMP_VAL
!
! Find derivatives
!
      DO I=1,NUM_SM-1
        DERIV(I)=(CROSS_SM(I)-CROSS_SM(I+1))/(NU_SM(I)-NU_SM(I+1))
      END DO
!
! Use first (threshold) point
!
      NUM_CUT=1
      NU_CUT(NUM_CUT)=NU_SM(1)
      CROSS_CUT(NUM_CUT)=CROSS_SM(1)
      LOW=1
!
! Some cross sections start above threshold so smoothed cross section will have
! cross_sm=0.  Must find first cross_sm.ne.0 and start calculation there
!
      IF(CROSS_SM(2) .EQ. 0.0D0)THEN
        DO I=3,NUM_SM
          IF(CROSS_SM(I) .NE. 0.0D0)GOTO 400
        ENDDO
        STOP ' ALL CROSS SECTIONS = 0'
 400    LOW=I-1
        NUM_CUT=NUM_CUT+1
        NU_CUT(NUM_CUT)=NU_SM(LOW)
        CROSS_CUT(NUM_CUT)=CROSS_SM(LOW)
        LOW=LOW+1
        NUM_CUT=NUM_CUT+1
        NU_CUT(NUM_CUT)=NU_SM(LOW)
        CROSS_CUT(NUM_CUT)=CROSS_SM(LOW)
      ENDIF
!
! Can now begin point selection in earnest. All maxima are retained.
! Points are retained such that a linear interpolation gives an accuracy of 
! ACCURACY, when the corss-section is > 1.0D-06. Outside this range at 
! least 20 points per decade are retained.
!
! Find next maximum or minimum point (dy/dx=0).
! Inflection points (d2y/d2x=0) are now omitted.
!
 330  IF(LOW .EQ. 1)THEN
        IF((DERIV(1)*DERIV(2)) .LT. 0.0D0)GOTO 300
        DO I=LOW+2,NUM_SM-1
          IF((DERIV(I-1)*DERIV(I)) .LE. 0.0D0)GOTO 300
        END DO
      ELSE
        DO I=LOW+1,NUM_SM-1
          IF((DERIV(I-1)*DERIV(I)) .LE. 0.0D0)GOTO 300
        END DO
      ENDIF
      HIGH=NUM_SM
 300  HIGH=I
!
 320  IF(HIGH-LOW.EQ.1)THEN
        LOW=HIGH
        NUM_CUT=NUM_CUT+1
        NU_CUT(NUM_CUT)=NU_SM(LOW)
        CROSS_CUT(NUM_CUT)=CROSS_SM(LOW)
        IF(LOW.EQ.NUM_SM)GOTO 340
        GOTO 330
      ELSE
!
! Determine the distance from each exact point (i=low,high)
! to the straight line interpolation between LOW and HIGH.
! 
        M=(CROSS_SM(HIGH)-CROSS_SM(LOW))/(NU_SM(HIGH)-NU_SM(LOW))
        DO I=LOW+1,HIGH-1
	  Y=CROSS_SM(LOW)+M*(NU_SM(I)-NU_SM(LOW))
          DIST(I)=ABS(Y-CROSS_SM(I))
        ENDDO
        MID=LOW+1
        DO I=LOW+2,HIGH-1
          IF(DIST(I).GT.DIST(MID))MID=I
        ENDDO
!
! We can omitted the data points IF and ONLY IF every point
! within the interval is within ACCURACY of the interpolating
! line. If not, we choose to include that point which has the
! biggest separation, and start the test procedure all over.
!
! The use of comp_val avoids problems when the cross-sections are zero,
!
        COMP_VAL=MIN(CROSS_SM(LOW),CROSS_SM(HIGH))
	IF(COMP_VAL .EQ. 0.0D0)COMP_VAL=1.0d-06
!
! We ensure that there is at least 20 points per decade of frequency space,
! when the cross-section (which is Mbarns) is less than 1.0D-06.
!
        IF( (DIST(MID)/COMP_VAL) .GT. CUT_ACCURACY)THEN
          HIGH=MID
          GOTO 320		!RESTART TESTING PROCEDURE
        ELSE IF( NU_SM(HIGH)/NU_SM(LOW) .GT. 1.1D0 
     *                         .AND. COMP_VAL .EQ. 1.0d-06)THEN
          HIGH=MID
          GOTO 320		!RESTART TESTING PROCEDURE
        ELSE
          LOW=HIGH
          NUM_CUT=NUM_CUT+1
          NU_CUT(NUM_CUT)=NU_SM(LOW)
          CROSS_CUT(NUM_CUT)=CROSS_SM(LOW)
          IF(LOW.EQ.NUM_SM)GOTO 340
          GOTO 330
        ENDIF
      ENDIF
!
 340  CONTINUE
!
      RETURN
      END
