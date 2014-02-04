	REAL*8 FUNCTION INTSIGC(ENR,XION_POT,EMIN,EMAX)
	IMPLICIT NONE
!
	REAL*8 ENR              ! Energy at which we compute the integral of the cross section
	REAL*8 XION_POT         ! (first) ionization potential of species under consideration
	REAL*8 EMIN,EMAX        ! Bounds of integration
	REAL*8 XJ,ONEOVERJ
!
	INTSIGC = 0.0D0
	IF (EMIN.GE.EMAX) THEN
	   WRITE(6,*) 'Emin >= emax in intsigc - we stop ',emin,emax,xion_pot
	   STOP
	ELSE
	   XJ = 0.6D0 * XION_POT
	   ONEOVERJ =  1.0D0/XJ
	   INTSIGC = (DATAN((EMAX-XION_POT)*ONEOVERJ)-DATAN((EMIN-XION_POT)*ONEOVERJ) ) &
	          / DATAN((ENR-XION_POT)*ONEOVERJ/2.0D0)
	ENDIF
	
	RETURN
	END FUNCTION INTSIGC
