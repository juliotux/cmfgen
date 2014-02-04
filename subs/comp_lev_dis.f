C
C Module to store data necssary to compute occupation probabilities
C etc associated with level dissolution.
C
C The vectors should be defined as follows:
C
C  B_LEV_DIS=( 8.589E+14*(POPION**0.333333D0)/ED)**1.5D0
C  A_LEV_DIS=0.0009*(ED**0.1667)/SQRT(T)
C  X_LEV_DIS=(1+A_LEV_DIS)**3.15
C
C NB: B_LEV_DIS is undefined for ED=0, but this should be of no 
C       practical concern, if we specifically handle the case when 
C       level dissoultion is switched off (corresponds to ED=0)
C
C Based on the include file LEV_DIS_BLK.INC. Has dynamic allocation
C and variable indicating whether level dissolution is switched on
C in now explicitly included in the data block.
C
	MODULE MOD_LEV_DIS_BLK
	REAL*8, ALLOCATABLE :: B_LEV_DIS(:)
	REAL*8, ALLOCATABLE :: A_LEV_DIS(:)
	REAL*8, ALLOCATABLE :: X_LEV_DIS(:)
	LOGICAL MOD_DO_LEV_DIS
C
	END MODULE MOD_LEV_DIS_BLK
C
C Subroutine to:
C          (1) Allocate desired memory if not done already.
C          (2) Compute level dissolution constants.
C
	SUBROUTINE COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DIS,ND)
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
C
	INTEGER ND
	REAL*8 ED(ND)
	REAL*8 POPION(ND)
	REAL*8 T(ND)
	LOGICAL DO_LEV_DIS
C
	IF(.NOT. ALLOCATED(B_LEV_DIS))THEN
	  ALLOCATE( B_LEV_DIS(ND) )
	  ALLOCATE( A_LEV_DIS(ND) )
	  ALLOCATE( X_LEV_DIS(ND) )
	END IF
C
	MOD_DO_LEV_DIS=DO_LEV_DIS
	IF(DO_LEV_DIS)THEN
	  B_LEV_DIS=( 8.589D+14*(POPION**0.333333D0)/ED)**1.5D0
	  A_LEV_DIS=0.0009D0*(ED**0.1667D0)/SQRT(T)
	  X_LEV_DIS=(1.0D0+A_LEV_DIS)**3.15D0
	ELSE
	  B_LEV_DIS=1.0D+30		!Large number
	  A_LEV_DIS=0.0D0
	  X_LEV_DIS=1.0D0
	END IF
C
	RETURN
	END
