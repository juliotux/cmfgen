C
C Subroutine formulates a string (>= 132 characters) containing the
C transition name, wavelength, continuum flux and line EW.
C
	SUBROUTINE EW_FORMAT(OUTSTR,TRANSITION,LAMBDA,FLUX,EW,SOBOLEV)
	IMPLICIT NONE
C
C Altered 01-Aug-1997 : OUTSTR increased to 132 characters.
C                         Transition name now output at end of string, so
C                         as species with long names will be successfully
C                         identified.
C Altered 24-May-1996 : ERROR_LU, LUER inserted
C                       CNT check inserted.
C
C Created 16-Oct-1989
C
	CHARACTER TRANSITION*(*),OUTSTR*(*)
	REAL*8 LAMBDA,FLUX,EW
	LOGICAL SOBOLEV
C
	INTEGER L,TAB
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
	IF(LEN(OUTSTR) .LT. 132)THEN
	  WRITE(LUER,*)'OUTSTR too small in EW_FORMAT'
	  RETURN
	END IF
C
	TAB=1
	OUTSTR=' '
C
	IF(LAMBDA .LT. 1.0D+05)THEN
	  WRITE(OUTSTR(TAB:),200)LAMBDA
200	  FORMAT(3X,F8.2)
	ELSE
	  WRITE(OUTSTR(TAB:),210)LAMBDA
210	  FORMAT(1PE11.4)
	END IF
C
	TAB=TAB+11
	IF(FLUX .GT. 0.1D0 .AND. FLUX .LT. 1.0D+04)THEN
	  WRITE(OUTSTR(TAB:),300)FLUX
300	  FORMAT(F14.4)
	ELSE
	  WRITE(OUTSTR(TAB:),310)FLUX
310	  FORMAT(1PE14.4)
	END IF
C
	TAB=TAB+15
	IF(ABS(EW) .GT. 0.1D0 .AND. ABS(EW) .LT. 1.0D+04)THEN
	  WRITE(OUTSTR(TAB:),400)EW,SOBOLEV
400	  FORMAT(F12.3,3x,L1)
	ELSE
	  WRITE(OUTSTR(TAB:),410)EW,SOBOLEV
410	  FORMAT(1PE12.3,3X,L1)
	END IF
C
C Allows for a transition name up to 85 characters.
C
	TAB=TAB+20			!4 chracter gap
	L=LEN_TRIM(TRANSITION)
	IF(L .GT. 133-TAB)L=133-TAB
	OUTSTR(TAB:132)=TRIM(TRANSITION(1:L))
C 
	RETURN
	END
