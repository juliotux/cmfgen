C
C Routine to compute the vacuum wavelenth from a wavelength in AIR. The
C input and output units are Angstroms.
C
C Formulae is only valid for LAM(AIR) > than 2000Ang.
C
C Results checked against table in Allen: Astrophysical quantities.
C
	FUNCTION LAM_VAC(LAM_AIR)
	IMPLICIT NONE
C
C Created 11-Apr-1997 : Based on LAMVACAIR
C
	REAL*8 LAM_VAC			!Vacuum wavelength in Ang.
	REAL*8 LAM_AIR			!Air wavelength in Ang.
	REAL*8 T1
	INTEGER I
C
	IF(LAM_AIR .LT. 1999.352D0)THEN
	  WRITE(6,*)' ERROR --- in LAM_VAC'
	  WRITE(6,*)' Formulae for conversion of air to vacuum wavelenths',
	1                ' invalid below 2000Ang'
	  STOP
	END IF
C
C Need to iterate, as T1 is unknown vacuum wavelength. Because difference
C between LAM_VAC and LAM_AIR is small, convergence is VERY rapid.
C
	LAM_VAC=LAM_AIR
	DO I=1,3
	  T1=1.0D+08/(LAM_VAC*LAM_VAC)
	  LAM_VAC=LAM_AIR*(  1.0D0+
	1               1.0D-07*( 643.28D0+ 
	1              294981.0D0/(146.0D0-T1)+2554.0D0/(41.0D0-T1) )  )
	END DO
C
	RETURN
	END
