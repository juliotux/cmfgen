!
! Subroutine to write text onto plot. Ths routine allows multiple
! colors in the same string. It currently only works when the strings
! are left justified (i.e., OFFSETs of 1, 4, and 7).
!
! The color of the text is changed by \pn where n=1, to 99.
! Assuming the default definitions, 
!            this is red \p2 this is blue \p3.
! The space before the \p is ignored.
!
! Each string must be lss than 200 characters, and a maximum of 10 color chanes
! is allowed.
!
	SUBROUTINE PUT_TEXT(X,Y,ANGLE,OFFSET,STRING)
	IMPLICIT NONE
!
! Created 02-Feb-2015
!
	REAL*4 X,Y
	REAL*4 ANGLE
	REAL*4 OFFSET
	CHARACTER(LEN=*) STRING
!
	INTEGER, PARAMETER :: NMAX=10
	INTEGER PEN(10)
	CHARACTER(LEN=200) FULL_STR
	CHARACTER(LEN=200) SUB_STR(NMAX)
!
	REAL*4 X1,Y1			!Internal variables controlling location of string (world coord).
	REAL*4 XLEN,YLEN		!Length os string (world coord.)
	REAL*4 XW1,XW2,YW1,YW2		!Location of viewport in world coordinates
	REAL*4 XP1,XP2,YP1,YP2    	!Location of viewport in pixels 
	REAL*4 SCALE_FAC		!Used to adjust string location
	INTEGER I,K
	INTEGER CI_SAV
!
! We use FULL_STR so that we dont change the passed text.
!
	FULL_STR=STRING
	PEN=0; SUB_STR=' '
!
! Check for capiialization of the pen cotrol character.
!
	K=INDEX(FULL_STR,'\P')
	DO WHILE(K .NE. 0)
	  FULL_STR(K+1:K+1)='p'
	  K=INDEX(FULL_STR,'\P')
	END DO
!
! Split text so that each color is a separate string (and stored in SUB_STR).
! Only 10 substrings are allowed.
!	
	DO I=1,NMAX
	  K=INDEX(FULL_STR,'\p')
	  IF(K .EQ. 0)THEN
	    SUB_STR(I)=TRIM(FULL_STR)
	    EXIT
	  ELSE
	    SUB_STR(I)=FULL_STR(1:K-1)
	    FULL_STR=FULL_STR(K+2:)
	    K=INDEX(FULL_STR,' ')
	    IF(K .GE. 3)THEN
	      IF(FULL_STR(2:2) .LT. '0' .AND. FULL_STR(2:0) .GT. '9')K=K-1
	    END IF
	    IF(FULL_STR(1:1) .LT. '0' .AND. FULL_STR(1:1) .GT. '9')K=K-1
	    IF(K .LE. 1)THEN
	      WRITE(6,*)'Invalid pen color in STRING'
	      WRITE(6,*)TRIM(STRING)
	      WRITE(6,*)TRIM(FULL_STR)
	      RETURN
	    END IF
	    READ(FULL_STR(1:K-1),*)PEN(I)
	    FULL_STR=FULL_STR(K:)
	  END IF
	END DO
!
! Save current pen
!
	CALL PGQCI(CI_SAV)
!
! Get size of viewport in world-coordinates, and in pixles. These are use to
! adjust the string location in the Y direction. We need to work in pixel space
! to allow for different viewport aspect ratios.
!
	CALL PGQWIN(XW1,XW2,YW1,YW2)
	CALL PGQVP(3,XP1,XP2,YP1,YP2)
	SCALE_FAC=ABS(YW2-YW1)/ABS(XW2-XW1)*ABS(XP1-XP2)/ABS(YP1-YP2)
!
	X1=X; Y1=Y
	DO I=1,NMAX
	  IF(SUB_STR(I) .EQ. ' ')EXIT
	  IF(PEN(I) .NE. 0)CALL PGSCI(PEN(I))
	  CALL PGPTXT(X1,Y1,ANGLE,OFFSET,SUB_STR(I))
	  K=4; CALL PGLEN(K,TRIM(SUB_STR(I)),XLEN,YLEN)		!Get string length in world coords.
	  X1=X1+COS(ANGLE*3.1415926D0/180.0D0)*XLEN
	  Y1=Y1+SIN(ANGLE*3.1415926D0/180.0D0)*XLEN*SCALE_FAC
	END DO
!
! Restore original pen setting.
!
	CALL PGSCI(CI_SAV)
!
	RETURN
	END
