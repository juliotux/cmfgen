C
C Subroutine to compute the opacity due to electron scattering at
C ND depth points.
C
	SUBROUTINE ESOPAC(RKI,ED,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 : Cleaned. (IMPLICT NONE inserted)
C
	INTEGER ND
	REAL*8 RKI(ND),ED(ND)
C
	RKI(1:ND)=6.65D-15*ED(1:ND)
C
	RETURN
	END
