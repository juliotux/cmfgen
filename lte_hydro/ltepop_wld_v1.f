C
C General routine to:
C
C       Compute the LTE populations of the levels in the FULL atom given
C       ED (electron density), T (electron temperature) and the ion density
C       DI C2. Level dissolution is taken into account.
C              
	SUBROUTINE LTEPOP_WLD_V1(C2LTE,W_C2,EDGEC2,GC2,
	1             ZC2,GION_C2,NC2,DIC2,ED,T,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 : ND_MAX removed (was unused).
C Created 30-May-1995 : Based on LTEPOP.
C             
	INTEGER ND
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperature 10^4K
	REAL*8 DIC2(ND)			!Ion density (Full model atom)
C
	INTEGER NC2
	REAL*8 C2LTE(NC2,ND)
	REAL*8 W_C2(NC2,ND)
	REAL*8 EDGEC2(NC2)
	REAL*8 GC2(NC2)
	REAL*8 GION_C2			!Statistical weight of ion groun state.
	REAL*8 ZC2			!Ion charge
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables.
C
	INTEGER I,K
	REAL*8 X,Y,RGU
C
C Compute the occupation probabilities.
C
	CALL OCCUPATION_PROB(W_C2,EDGEC2,ZC2,NC2,ND)
C
C Compute the LTE populations of the levels in the full atom, taking level
C dissolution into account.
C
	RGU=DLOG(2.07078D-22)
	DO K=1,ND
	  X=HDKT/T(K)
	  Y=ED(K)*DIC2(K)*( T(K)**(-1.5) )/GION_C2
	  IF(Y .GT. 0)THEN
	    DO I=1,NC2
	      C2LTE(I,K)=W_C2(I,K)*GC2(I)*Y*EXP(EDGEC2(I)*X+RGU)
	    END DO
	  ELSE
	    C2LTE(:,K)=0.0D0
	  END IF
	END DO
C
	RETURN
	END
