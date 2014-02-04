C
C Routine to compute the contibution to the electron density by
C a set of ionization species. U and PHI are non-LTE partition
C functions as computed by PAR_FUN. Also returned are dCHARGEdNE
C (so that a Newton-Rapshon technique can be used to solve for Ne
C as equations are non-linear) and the POPULATION of the highest
C ionization species.
C
C Created 29-Aug-1990.
C
	SUBROUTINE EVAL_ED(CHARGE,dCHARGEdNE,U,PHI,Z,ED,POPOXY,
	1                  DION,LOC_CH,dDIONdNe,dLOC_CH_dNe,
	1                  NSPEC,ND,FIRST)
	IMPLICIT NONE
C
	LOGICAL FIRST
	INTEGER ND,NSPEC
	REAL*8 CHARGE(ND),dCHARGEdNE(ND),U(ND,NSPEC),PHI(ND,NSPEC)
	REAL*8 Z(NSPEC),ED(ND),POPOXY(ND)
	REAL*8 DION(ND),LOC_CH(ND),dDIONdNe(ND),dLOC_CH_dNe(ND)
C
	REAL*8 T1
	INTEGER I,J
C
	IF(FIRST)THEN
	  DO I=1,ND
	    CHARGE(I)=0.0D0
	    dCHARGEdNe(I)=0.0D0
	  END DO
	  FIRST=.FALSE.
	END IF
C
	DO I=1,ND
	  DION(I)=0.0D0
	  dDIONdNe(I)=0.0D0
	  LOC_CH(I)=0.0D0
	  dLOC_CH_dNe(I)=0.0D0
	END DO
!
	DO I=1,ND
	  DO J=1,NSPEC-1
	    T1=PHI(I,J)*( U(I,J)/U(I,J+1) )
            DION(I)=T1*ED(I)*( 1.0D0+DION(I) )
	    dDIONdNe(I)=T1*( (NSPEC-1) + ED(I)*dDIONdNe(I) )
	    LOC_CH(I)=T1*ED(I)*( Z(J)+LOC_CH(I) )
	    dLOC_CH_dNe(I)=T1*( Z(J)*(NSPEC-1) +
	1                        ED(I)*dLOC_CH_dNe(I) )
	  END DO
	  DION(I)=DION(I)+1.0D0
	  LOC_CH(I)=LOC_CH(I)+Z(NSPEC)
	END DO
C
	DO I=1,ND
	  T1=POPOXY(I)/DION(I)
	  CHARGE(I)=CHARGE(I)+T1*LOC_CH(I)
	  dCHARGEdNE(I)=dCHARGEdNE(I) + T1*(
	1                  dLOC_CH_dNe(I) - dDIONdNe(I)*(LOC_CH(I)/DION(I)) )
	  DION(I)=T1
	END DO
C
	RETURN
	END
	
