C
C Auxilary routine to help compute the mean ionic charge (i.e gam)
C of a given atomic species.
C
C Created 10-Apr-1989.
C
	SUBROUTINE UPDATE_GAM(GAM,C2,DC2,ZC2,NC2,ND,CIII_PRES,FIRST)
	IMPLICIT NONE
C
	INTEGER ND,NC2
	REAL*8 GAM(ND),C2(NC2,ND),DC2(ND),ZC2
	LOGICAL CIII_PRES,FIRST
C
C Local variabes.
C
	INTEGER I,J
C
	IF(FIRST)THEN
	  DO J=1,ND
	    GAM(J)=0.0D0
	  END DO
	END IF
C
	DO J=1,ND
	  DO I=1,NC2
	    GAM(J)=GAM(J)+(ZC2-1)*C2(I,J)
	  END DO
	END DO
C
	IF(.NOT. CIII_PRES)THEN
	  DO J=1,ND
	    GAM(J)=GAM(J)+ZC2*DC2(J)
	  END DO
	END IF
C
	FIRST=.FALSE.
	RETURN
	END
