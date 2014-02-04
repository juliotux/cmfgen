!
! Subroutine to compute the LTE populations (at NR depth points)
! given ED (electron density) and DI,GU (density and statistical
! weight of the ground state of the next ionization state 
! respectively).
!
	SUBROUTINE LTEPOP(HNST,ED,DI,G,NUION,T,GU,N,NR)
	IMPLICIT NONE
!
! Altered 18-FEb-2010 : Take log of density to extend range of LTE populations.
!                         We use NEW_METHOD to allow easy change to previous
!                         version. Only necessary if something untoward crops up.
! Altered 05-Dec-1996 : END DO used to terminate DO LOOPS.
! Altered 28-May-1996 : Generica calls for EXP and LOG
! Altered 10-Apr-1989 - Implicit none installed.
! Altered 14-AUG-1984
!
	INTEGER N,NR
	REAL*8 HNST(N,NR),ED(NR),DI(NR),T(NR),G(N),NUION(N),GU
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local varaiables.
!
	INTEGER I,J
	REAL*8 X,Y,RGU
	LOGICAL, PARAMETER :: NEW_METHOD=.TRUE.
!
	IF(NEW_METHOD)THEN
	  DO I=1,NR
	    X=HDKT/T(I)
	    RGU=2.07078D-22*ED(I)*DI(I)*( T(I)**(-1.5D0) )/GU
	    RGU=LOG(RGU)
	    DO J=1,N
	      HNST(J,I)=G(J)*EXP(NUION(J)*X+RGU)
	    END DO
	  END DO
	ELSE
	  RGU=2.07078D-22/GU
	  RGU=LOG(RGU)
	  DO I=1,NR
	    X=HDKT/T(I)
	    Y=ED(I)*DI(I)*( T(I)**(-1.5D0) )
	    DO J=1,N
	      HNST(J,I)=G(J)*Y*EXP(NUION(J)*X+RGU)
	    END DO
	  END DO
	END IF
!
	RETURN
	END
