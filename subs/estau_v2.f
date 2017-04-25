!
! Routine to compute the electron scattering optical depth. The trapazoidal
! rule is used, and at the outer boundary it is assumed the Ne varies as
! r^{-2}. Will need to be altered for an exponential atmosphere.
!
	SUBROUTINE ESTAU_V2(OPT_DEP,R,ED,CLUMP_FAC,DTAU,ND)
	IMPLICIT NONE
!
! Altered 01-Apr-2015: Changed to V2, CLUMP_FAC added to call.
!                        TAU_BND changed to be consistent with output to
!                        MEANOPAC.
!
	INTEGER ND,I
	REAL*8 R(ND),ED(ND),CLUMP_FAC(ND),DTAU(ND),OPT_DEP(ND)
	REAL*8 T1,TAU_BND
!
	DO I=1,ND
	  OPT_DEP(I)=0.0D0
	  DTAU(I)=ED(I)*CLUMP_FAC(I)        !Temporary work variable
	END DO
	CALL ESOPAC(OPT_DEP,DTAU,ND)        !OPT_DEP contain ESEC
!
	T1=LOG(OPT_DEP(1)/OPT_DEP(4))/LOG(R(4)/R(2))
	T1=MAX(T1,2.0D0)
	TAU_BND=OPT_DEP(1)*R(1)/(T1-1.0D0)
!
	DO I=1,ND-1
	  DTAU(I)=0.5D0*( OPT_DEP(I)+OPT_DEP(I+1) )*( R(I)-R(I+1) )
	END DO
!
	OPT_DEP(1)=TAU_BND
	DO I=2,ND
 	  OPT_DEP(I)=OPT_DEP(I-1)+DTAU(I-1)
	END DO
!
	RETURN
	END




