C
C Routine to compute the electron scattering optical depth. The trapazoidal
C rule is used, and at the outer boundary it is assumed the Ne varies as
C r^{-2}. Will need to be altered for an exponential atmosphere.
C
	SUBROUTINE ESTAU(OPT_DEP,R,ED,DTAU,ND)
	IMPLICIT NONE
C
	INTEGER ND,I
	REAL*8 R(ND),ED(ND),DTAU(ND),OPT_DEP(ND)
C
	DO I=1,ND
	  OPT_DEP(I)=0.0D0
	END DO
	CALL ESOPAC(OPT_DEP,ED,ND)
C
	DO I=1,ND-1
	  DTAU(I)=0.5D0*( OPT_DEP(I)+OPT_DEP(I+1) )*( R(I)-R(I+1) )
	END DO
C
	OPT_DEP(1)=OPT_DEP(1)*R(1)
	DO I=2,ND
 	  OPT_DEP(I)=OPT_DEP(I-1)+DTAU(I-1)
	END DO
C
	RETURN
	END




