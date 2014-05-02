!
! Routine to compute the depth dependent luminosity arising from:
!    Lines (non-blanketed mode),
!    DIELECTRONIC lines (non-blankted mode)
!    Shock X_RAY emisison
!    Mechanical energy term
!    Deposition of radioactive energy etc
! 
! On entry LINE_LUM should contain r^2 * EMISSIVITY at radius R 
!  (possibly corrected for the LOCAL escape probablity).
!
! We compute the amount of energy emitted between R(I) and R(I+1), and
! store it in LINE_LUM(I). This is computed using the Euler-McLaurin summation
! rule (i.e. the trapazoidal rule with a correction for the first derivatives).
!
	SUBROUTINE LUM_FROM_ETA_V2(LINE_LUM,R,METHOD,ND)
	IMPLICIT NONE
!
! Created 19-Jan-2014 - Based on LUM_FROM_ETA (included METHOD in call).
!
	INTEGER ND
	REAL*8 LINE_LUM(ND)
	REAL*8 R(ND)
	CHARACTER(LEN=*) METHOD
!
	REAL*8 DERIV(ND)
	INTEGER I
!
! We compute the amount of energy emitted between R(I) and R(I+1), and
! store it in LINE(I). This is computed using the Euler-McLaurin summation
! rule (i.e. the trapazoidal rule with a correction for the first derivatives).
!
	IF(METHOD .EQ. ' ')THEN
	  DERIV(1)=(LINE_LUM(1)-LINE_LUM(2))/(R(1)-R(2))
	  DO I=2,ND-1
	    DERIV(I)=(LINE_LUM(I-1)-LINE_LUM(I+1))/(R(I-1)-R(I+1))
	  END DO
	  DERIV(ND)=(LINE_LUM(ND-1)-LINE_LUM(ND))/(R(ND-1)-R(ND))
	ELSE
	  CALL DERIVCHI(DERIV,LINE_LUM,R,ND,METHOD)
	END IF
!
	DO I=1,ND-1
	  LINE_LUM(I)=0.5D0*(R(I)-R(I+1))*( LINE_LUM(I)+LINE_LUM(I+1) 
	1            +(R(I)-R(I+1))*(DERIV(I+1)-DERIV(I))/6.0D0 )
	END DO
	LINE_LUM(ND)=0.0D0
!
	RETURN
	END
