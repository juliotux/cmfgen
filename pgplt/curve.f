!
! Program to initialize plotting arrays.
!
	SUBROUTINE CURVE(NUM,X,Y)
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
! Modified 28-Mar-2003 : Check on NUM=0 installed. Minor clean.
! Modified 15-Jul-2000 : Incorporated MOD_CURVE DATA
!                        Dynamically allocated arrays.
! Modified 20-Jul-1994 : Now automatically allows for NUM to be > 2048.
! Previous version was lost - this version created 12-NOV-85
!
	INTEGER NUM
	REAL X(NUM),Y(NUM)
!
	INTEGER J
!
	IF(NUM .EQ. 0)THEN
	  WRITE(LU_ER,*)'Error - no data passed to curve'
	  RETURN
	END IF
	IF(NPLTS+1 .GT. MAX_PLTS)THEN
	  WRITE(LU_ER,*)'Error - too many call to curve'
	  RETURN
	END IF
!
	J=NPLTS+1
	ALLOCATE (CD(J)%XVEC(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%DATA(NUM),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LU_ER,*)'Erorr allocating plt storage in CURVE'
	  WRITE(LU_ER,*)'IOS=',IOS
	  IF(ALLOCATED(CD(J)%XVEC))DEALLOCATE(CD(J)%XVEC)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%DATA)
	END IF
!
	CD(J)%XVEC(1:NUM)=X(1:NUM)
	CD(J)%DATA(1:NUM)=Y(1:NUM)
	NPTS(J)=NUM
	ERR(J)=.FALSE.
	NPLTS=J
!
	RETURN
	END
!
!
!
! Program to initialize plotting arrays - Includes error option.
!
! Created 06-Mar-1990 - Based on curve
!                   
	SUBROUTINE CURVE_AND_ER(NUM,X,Y,SIGMA,OPT)
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
	INTEGER NUM
	REAL X(NUM),Y(NUM),SIGMA(NUM)
	CHARACTER*(*) OPT
!
	REAL T1
	INTEGER I,J
!
	IF(NUM .EQ. 0)THEN
	  WRITE(LU_ER,*)'Error - no data passed to curve'
	  RETURN
	END IF
	IF(NPLTS+1 .GT. MAX_PLTS)THEN
	  WRITE(LU_ER,*)'Error - too many call to curve'
	  RETURN
	END IF
!
	J=NPLTS+1
	ALLOCATE (CD(J)%XVEC(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%DATA(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%EMAX(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%EMIN(NUM),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LU_ER,*)'Erorr allocating plt storage in CURVE'
	  WRITE(LU_ER,*)'IOS=',IOS
	  IF(ALLOCATED(CD(J)%XVEC))DEALLOCATE(CD(J)%XVEC)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%DATA)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%EMAX)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%EMIN)
	END IF
!
	CD(J)%XVEC(1:NUM)=X(1:NUM)
!
	IF(OPT(1:3) .EQ. 'LOG')THEN
	  DO I=1,NUM
	    IF(Y(I) .GT. 0)THEN
	      CD(J)%DATA(I)=LOG10(Y(I))
	      CD(J)%EMAX(I)=LOG10(Y(I)+SIGMA(I))
	      T1=Y(I)-SIGMA(I)
	      IF(T1 .LE. 0)THEN
	        CD(J)%EMIN(I)=LOG10(Y(I))-2.0
	      ELSE
	        CD(J)%EMIN(I)=LOG10(T1)
	      END IF
	    ELSE
	      CD(J)%DATA(I)=-1000.0
	      CD(J)%EMAX(I)=0.0
	      CD(J)%EMIN(I)=0.0
	    END IF
	  END DO
	ELSE
	  CD(J)%XVEC=X
	  CD(J)%DATA=Y
	  CD(J)%EMAX=Y+SIGMA
	  CD(J)%EMIN=Y-SIGMA
	END IF
!
	NPLTS=J
	NPTS(J)=NUM
	ERR(NPLTS)=.TRUE.
!
	RETURN
	END 
!
!
!
! Program to initialize plotting arrays.
!
	SUBROUTINE DP_CURVE(NUM,X,Y)
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
	INTEGER NUM
	REAL*8 X(NUM),Y(NUM)
!
	INTEGER J
!
	IF(NUM .EQ. 0)THEN
	  WRITE(LU_ER,*)'Error - no data passed to curve'
	  RETURN
	END IF
	IF(NPLTS+1 .GT. MAX_PLTS)THEN
	  WRITE(LU_ER,*)'Error - too many call to curve'
	  RETURN
	END IF
!
	J=NPLTS+1
	ALLOCATE (CD(J)%XVEC(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%DATA(NUM),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LU_ER,*)'Erorr allocating plt storage in CURVE'
	  WRITE(LU_ER,*)'IOS=',IOS
	  IF(ALLOCATED(CD(J)%XVEC))DEALLOCATE(CD(J)%XVEC)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%DATA)
	END IF
!
	CD(J)%XVEC(1:NUM)=X(1:NUM)
	CD(J)%DATA(1:NUM)=Y(1:NUM)
	NPTS(J)=NUM
	ERR(J)=.FALSE.
	NPLTS=J
!
	RETURN
	END
!
!
!
! Program to initialize plotting arrays - Includes error option.
!
! Created 06-Mar-1990 - Based on curve
!
	SUBROUTINE DP_CURVE_AND_ER(NUM,X,Y,SIGMA,OPT)
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
	INTEGER NUM
	REAL*8 X(NUM),Y(NUM),SIGMA(NUM)
	CHARACTER*(*) OPT
C
	REAL T1
	INTEGER I,J
!
	IF(NUM .EQ. 0)THEN
	  WRITE(LU_ER,*)'Error - no data passed to curve'
	  RETURN
	END IF
	IF(NPLTS+1 .GT. MAX_PLTS)THEN
	  WRITE(LU_ER,*)'Error - too many call to curve'
	  RETURN
	END IF
!
	J=NPLTS+1
	ALLOCATE (CD(J)%XVEC(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%DATA(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%EMAX(NUM),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(J)%EMIN(NUM),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LU_ER,*)'Erorr allocating plt storage in CURVE'
	  WRITE(LU_ER,*)'IOS=',IOS
	  IF(ALLOCATED(CD(J)%XVEC))DEALLOCATE(CD(J)%XVEC)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%DATA)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%EMAX)
	  IF(ALLOCATED(CD(J)%DATA))DEALLOCATE(CD(J)%EMIN)
	END IF
!
	CD(J)%XVEC(1:NUM)=X(1:NUM)
!
	IF(OPT(1:3) .EQ. 'LOG')THEN
	  DO I=1,NUM
	    IF(Y(I) .GT. 0)THEN
	      CD(J)%DATA(I)=LOG10(Y(I))
	      CD(J)%EMAX(I)=LOG10(Y(I)+SIGMA(I))
	      T1=Y(I)-SIGMA(I)
	      IF(T1 .LE. 0)THEN
	        CD(J)%EMIN(I)=LOG10(Y(I))-2.0
	      ELSE
	        CD(J)%EMIN(I)=LOG10(T1)
	      END IF
	    ELSE
	      CD(J)%DATA(I)=-1000.0
	      CD(J)%EMAX(I)=0.0
	      CD(J)%EMIN(I)=0.0
	    END IF
	  END DO
	ELSE
	  CD(J)%XVEC=X
	  CD(J)%DATA=Y
	  CD(J)%EMAX=Y+SIGMA
	  CD(J)%EMIN=Y-SIGMA
	END IF
!
	NPLTS=J
	NPTS(J)=NUM
	ERR(NPLTS)=.TRUE.
!
	RETURN
	END
