!
! General routine to:
!
!       Compute the LTE populations of the levels in the FULL atom given
!       ED (electron density), T (electron temperature) and the ion density
!       DI C2. Level dissolution is taken into account.
!              
	SUBROUTINE LTEPOP_WLD_V1(C2LTE,W_C2,EDGEC2,GC2,
	1             ZC2,GION_C2,NC2,DIC2,ED,T,ND)
	IMPLICIT NONE
!
!
! Altered 18-FEb-2010 : Take log of density to extend range of LTE populations.
!                         We use NEW_METHOD to allow easy change to previous
!                         version. Only necessary if something untoward crops up
! Altered 24-May-1996 : ND_MAX removed (was unused).
! Created 30-May-1995 : Based on LTEPOP.
!             
	INTEGER ND
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperature 10^4K
	REAL*8 DIC2(ND)			!Ion density (Full model atom)
!
	INTEGER NC2
	REAL*8 C2LTE(NC2,ND)
	REAL*8 W_C2(NC2,ND)
	REAL*8 EDGEC2(NC2)
	REAL*8 GC2(NC2)
	REAL*8 GION_C2			!Statistical weight of ion groun state.
	REAL*8 ZC2			!Ion charge
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER I,K
	REAL*8 X,Y,RGU
	LOGICAL, PARAMETER :: NEW_METHOD=.TRUE.
!
! Compute the occupation probabilities.
!
	CALL OCCUPATION_PROB(W_C2,EDGEC2,ZC2,NC2,ND)
!
! Compute the LTE populations of the levels in the full atom, taking level
! dissolution into account.
!
	IF(NEW_METHOD)THEN
	  DO K=1,ND
	    X=HDKT/T(K)
	    RGU=2.07078D-22*ED(K)*DIC2(K)*( T(K)**(-1.5D0) )/GION_C2
	    RGU=DLOG(RGU)
	    DO I=1,NC2
	      C2LTE(I,K)=W_C2(I,K)*GC2(I)*EXP(EDGEC2(I)*X+RGU)
	    END DO
	  END DO
	ELSE
	  RGU=DLOG(2.07078D-22)
	  DO K=1,ND
	    X=HDKT/T(K)
	    Y=ED(K)*DIC2(K)*( T(K)**(-1.5D0) )/GION_C2
	    DO I=1,NC2
	      C2LTE(I,K)=W_C2(I,K)*GC2(I)*Y*EXP(EDGEC2(I)*X+RGU)
	    END DO
	  END DO
	END IF
!
	RETURN
	END
