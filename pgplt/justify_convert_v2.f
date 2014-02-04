!
! Subroutine to convert MONGO justification commands to PGLPOT
! MONGO uses a system with nine different string justifications
! while PGPLOT has only three.
!
      SUBROUTINE JUSTIFY_CONVERT_V2(XSTR,YSTR,LOC,LOC_PG,ORIENTATION,
     *             FLAGSTR,XSTRPOS,YSTRPOS,STRING,MAX_NUM_STR)
!
!Altered 18-Nov-1990 : Changed to V2: MAX_NUM_STR passed in call
!                      Absolute value of LOC used.
!Created August 1996 by Gregson Vaux
!
      INTEGER MAX_NUM_STR
      REAL*4 XSTR(MAX_NUM_STR), YSTR(MAX_NUM_STR)
      REAL*4 XSTRPOS(MAX_NUM_STR),YSTRPOS(MAX_NUM_STR)
      REAL*4 ORIENTATION(MAX_NUM_STR)
      REAL*4 LOC_PG(MAX_NUM_STR)
      LOGICAL FLAGSTR(MAX_NUM_STR)
      INTEGER LOC(MAX_NUM_STR)
      CHARACTER*(*) STRING(MAX_NUM_STR)
!
!Local variables
!
      INTEGER ISTR
      INTEGER ABS_LOC
      REAL*4 XBOX(4),YBOX(4)
!
      DO ISTR=1,MAX_NUM_STR
        IF (FLAGSTR(ISTR)) THEN
	  ABS_LOC=ABS(LOC(ISTR))
          IF ((ABS_LOC==1) .OR. (ABS_LOC==4) .OR. (ABS_LOC==7)) THEN
            LOC_PG(ISTR)=0.0
          ELSE IF ((ABS_LOC==2) .OR. (ABS_LOC==5) .OR. (ABS_LOC==8))
     *      THEN
            LOC_PG(ISTR)=0.5
          ELSE IF ((ABS_LOC==3) .OR. (ABS_LOC==6) .OR. (ABS_LOC==9))
     *      THEN
            LOC_PG(ISTR)=1.0
          END IF
          CALL PGQTXT(XSTR(ISTR),YSTR(ISTR),ORIENTATION(ISTR),
	1               LOC_PG(ISTR),STRING(ISTR),XBOX,YBOX)
          IF (ABS_LOC .LE. 3)THEN
            XSTRPOS(ISTR)=XSTR(ISTR)
            YSTRPOS(ISTR)=YSTR(ISTR)
	  ELSE IF (ABS_LOC .LE. 6) THEN
            XSTRPOS(ISTR)=XSTR(ISTR)-0.5*(XBOX(2)-XBOX(1))
            YSTRPOS(ISTR)=YSTR(ISTR)-0.5*(YBOX(2)-YBOX(1))
          ELSE  IF (ABS_LOC .LE. 9) THEN
            XSTRPOS(ISTR)=XSTR(ISTR)-(XBOX(2)-XBOX(1))
            YSTRPOS(ISTR)=YSTR(ISTR)-(YBOX(2)-YBOX(1))
          END IF
        END IF
      END DO
!
      RETURN 
      END
