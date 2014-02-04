	SUBROUTINE ARNAUD_CROSS()
	USE MOD_NON_THERM
	IMPLICIT NONE
!
	INTEGER IT
	INTEGER IKT
	REAL*8 U1
	REAL*8 T1,T2,T3
!
! Fitting formula from Arnaud & Rothenflug 1985:
! http://adsabs.harvard.edu/abs/1985A%26AS...60..425A
! Gives the cross section for direct ionisation.
!
	DO IT=1,NUM_THD
	  IF(THD(IT)%PRES)THEN
	    THD(IT)%CROSS_SEC=0.0D0
	    DO IKT=1,NKT
	      U1 = XKT(IKT) / THD(IT)%ION_POT
	      IF (U1 .GT. 1.0D0) THEN
	        T1=(1.0D0-1.0D0/U1)
	        T2=DLOG(U1)
	        T3 = 1.0D-14 * ( THD(IT)%A_COL*T1 + THD(IT)%B_COL*T1*T1 + &
	              THD(IT)%C_COL*T2 + THD(IT)%D_COL*T2/U1  ) / U1 / THD(IT)%ION_POT**2
	        THD(IT)%CROSS_SEC(IKT)=T3
	      END IF
	    END DO
	  END IF
	END DO
!	
	RETURN
	END
