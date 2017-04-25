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
	SUBROUTINE STRIP_SLASH_P(STRING,N)
	IMPLICIT NONE
!
! Created 02-Feb-2015
!
	INTEGER N
	CHARACTER(LEN=*) STRING(N)
!
	INTEGER, PARAMETER :: NMAX=10
	CHARACTER(LEN=200) SUB_STR(NMAX)
	CHARACTER(LEN=200) FULL_STR
!
	INTEGER I,K,L
!
	DO L=1,N
	  FULL_STR=STRING(L)
	  SUB_STR=' '
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
	        WRITE(6,*)TRIM(STRING(L))
	      END IF
	      FULL_STR=FULL_STR(K:)
	    END IF
	  END DO
	  STRING(L)=SUB_STR(1)
	  DO I=2,NMAX
	    STRING(L)=TRIM(STRING(L))//TRIM(SUB_STR(I))
	  END DO
	END DO
!
	RETURN
	END
