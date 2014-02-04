C
C Function to return the wavelength of a transition which is specified 
C by a frequency (10^15 Hz). The wavelength is returned in Angstroms.
C
C If the wavelenth is less than 2000Ang, it is a vacuum wavelength,
C otherwise it is in air.
C
	FUNCTION LAMVACAIR(FREQ)
	IMPLICIT NONE
	REAL*8 LAMVACAIR,FREQ,T1
C
	LAMVACAIR=2997.92458D0/FREQ
C
	IF(LAMVACAIR .GT. 1999.352D0)THEN
	  T1=1.0D+08/(LAMVACAIR*LAMVACAIR)
	  LAMVACAIR=LAMVACAIR/(  1.0D0+
	1               1.0D-07*( 643.28D0+
	1              294981.0D0/(146.0D0-T1)+2554.0D0/(41.0D0-T1) )  )
	END IF
C
	RETURN
	END
C
	FUNCTION LAM_AIR(LAM)
	IMPLICIT NONE
	REAL*8 LAM_AIR,LAM,T1
C
	LAM_AIR=LAM
C
	IF(LAM_AIR .GT. 1999.352D0)THEN
	  T1=1.0D+08/(LAM_AIR*LAM_AIR)
	  LAM_AIR=LAM_AIR/(  1.0D0+
	1               1.0D-07*( 643.28D0+
	1              294981.0D0/(146.0D0-T1)+2554.0D0/(41.0D0-T1) )  )
	END IF
C
	RETURN
	END
