!
! Subroutine to incorporate a free-free line resonance into the
! normal free-free opacity. This is done by modifying the free-free
! gaunt factor for the species under consideration.
!
	SUBROUTINE FF_RES_GAUNT(GFF,FREQ,T,ID,GION,ZION,ND)
	USE PHOT_DATA_MOD
	IMPLICIT NONE
!
	INTEGER ID		!Species identifier
	INTEGER ND		!Number of depth points
	REAL*8 FREQ		!Frequency in units of 10^15 Hz
	REAL*8 GION		!Ion statistical weight
	REAL*8 ZION		!Charge on ion.
	REAL*8 GFF(ND)		!Free-free Gaunt factor
	REAL*8 T(ND)		!Temperature in unitsof 10^4 K
!
	REAL*8 VOIGT
	EXTERNAL VOIGT
!
! Common block with opacity/emissivity constants.
!
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
!
! Local variables
!
	REAL*8 DOP_NU
	REAL*8 A_VOIGT
	REAL*8 V_VOIGT
	REAL*8 PHI_VOIGT 
	REAL*8 CONST
	REAL*8 T1
!
	INTEGER I,K
!
! Check if any free-free resonances for this species.
!
	IF(PD(ID)%NUM_FF .EQ. 0)RETURN
!
! Compute resonance contribution. At present a fixed Doppler width is
! assumed. The Doppler width should be chosen to avoid undersampling of
! the intrinsic line profile with the adopted frequency grid. The intrinsic
! width of the line, set by the autoionization probabilities of the lower and 
! upper levels, is taken into account.
!
	DO I=1,PD(ID)%NUM_FF
	  IF(FREQ .LE. PD(ID)%FF_NU_MAX(I) .AND. FREQ .GE. PD(ID)%FF_NU_MIN(I))THEN
            DOP_NU=PD(ID)%NU_ZERO(I)*PD(ID)%VSM_KMS/2.998D+05
            A_VOIGT=PD(ID)%GAMMA(I)/DOP_NU
            V_VOIGT=(FREQ-PD(ID)%FF_NU_ZERO(I))/DOP_NU
	    CONST=1.489D+15*PD(ID)%FF_GF(I)*(PD(ID)%FF_NU_ZERO(I)**3)/GION/ZION/ZION
            PHI_VOIGT=1.0D-15*VOIGT(A_VOIGT,V_VOIGT)/DOP_NU
            T1=HDKT*PD(ID)%FF_NU_EXCITE(I)
	    WRITE(2,*)FREQ,PD(ID)%FF_NU_ZERO(I)
	    WRITE(2,*)DOP_NU,A_VOIGT,V_VOIGT,CONST,PHI_VOIGT,T1
	    DO K=1,ND
	      GFF(K)=GFF(K)+CONST*PHI_VOIGT*EXP(-T1/T(K))/T(K)
	    END DO
	  END IF
	END DO
!
	RETURN
	END
