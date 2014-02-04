!
! Altered 25-Sep-2011 - Now call GENOPAETA_V10 (changed from V7).
C
C Compute opacity and emissivity. This is a general include file
C provided program uses exactly the same variables. Can be achieved
C by copying declaration statements from CMFGEN. Always advisable
C to use ``IMPLICIT NONE''.
C
C Compute EXP(-hv/kT) and zero CHI, and ETA.
C
	T1=-HDKT*FL
	DO I=1,ND
	  EMHNUKT(I)=EXP(T1/T(I))
	  CHI(I)=0.0D0
	  ETA(I)=0.0D0
	END DO
C
C Compute vector constants for evaluating the level dissolution. These
C constants are the same for all species. These should be stored in a
C common block.             
C
	DO I=1,ND         
	  B_LEV_DIS(I)=( 8.589E+14*(POPION(I)**0.333333D0)/ED(I) )**1.5D0
	  A_LEV_DIS(I)=0.0009*(ED(I)**0.1667)/SQRT(T(I))
	  X_LEV_DIS(I)=(1+A_LEV_DIS(I))**3.15
	END DO
C
C Compute continuum intensity incident from the core assuming a TSTAR
C blackbody.
C
	T1=HDKT*FL/TSTAR
	IF(T1 .GT. 30)THEN
	  IC=TWOHCSQ*(FL**3)*EXP(-HDKT*FL/TSTAR)
	ELSE
	  IC=TWOHCSQ*(FL**3)/(EXP(HDKT*FL/TSTAR)-1.0D0)
	END IF
C
C Compute opacity and emissivity. ESOPAC must be call first since
C CHI is effectively zeroed in that routine.
C
	CALL ESOPAC(ESEC,ED,ND)		!Electron scattering emission factor.
	CALL ESOPAC(CHI,ED,ND)		!Electron scattering opacity.
C
C Free-free and bound-free opacities.
C
	CHI_PAR(:,:)=0.0D0
	ETA_PAR(:,:)=0.0D0
	DO ID=1,NUM_IONS
	  CHI(:)=0.0D0; I=ID
	  ETA(:)=0.0D0
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO J=1,ATM(ID)%N_XzV_PHOT
	      CALL GENOPAETA_V10(I,CHI,ETA,FL,
	1         ATM(ID)%XzV_F,ATM(ID)%XzVLTE_F,ATM(ID)%LOG_XzVLTE_F,ATM(ID)%EDGEXzV_F,
	1          ATM(ID)%GIONXzV_F,ATM(ID)%ZXzV,ATM(ID)%NXzV_F,
	1         ATM(ID+1)%XzV_F,ATM(ID+1)%LOG_XzVLTE_F,
	1         ATM(ID+1)%NXzV_F,
	1         J,ATM(ID)%XzV_ION_LEV_ID(J),ED,T,EMHNUKT,
	1         L_TRUE,ND,LST_DEPTH_ONLY)
	    END DO
	    CHI_PAR(:,ID)=CHI(:)
	    ETA_PAR(:,ID)=ETA(:)
	  END IF
	END DO
!
