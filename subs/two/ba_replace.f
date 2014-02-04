C
C Routine to automatically replace the ground state equation of a species
C with the ionization/recombination equation.
C
C \eg He2(n=1) with He2/HeIII ioinzation/recombination equation.
C
C Used for stability, as recomination rates mab be very small compared with
C other terms in the equations.
C
	SUBROUTINE BA_REPLACE(BA,STEQ,BAION,STEQION,C2,DC2,EDGEC2,T,
	1                     NC2,EQC2,EQC2ION,NION,NT,NUM_BNDS,ND,
	1                     DESC,C2PRES)
	IMPLICIT NONE
C
C Altered 24-May-1996 : ERROR_LU inserted
C                       Descriptor written to unit 53.
C Altered 13-Jul-1993 : Criterion for inclusion of recombination equation
C                       altered.
C
	INTEGER*4 NC2,EQC2,EQC2ION,NION,NT,NUM_BNDS,ND
	REAL*8 BA(NT,NT,NUM_BNDS,ND),STEQ(NT,ND)
	REAL*8 BAION(NION,NT,NUM_BNDS,ND),STEQION(NION,ND)
	REAL*8 C2(NC2,ND),DC2(ND),EDGEC2(NC2),T(ND)
	LOGICAL C2PRES
	CHARACTER*(*) DESC
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	REAL*8 T1,FAC
	INTEGER*4 I,J,K,L,LIM,ICHRLEN,DIAG
	INTEGER*4 LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
	IF(.NOT. C2PRES)RETURN
	IF(EQC2ION .EQ. 0)RETURN
	LUER=ERROR_LU()
C
C We compare the variation of the ground state equation with the
C that of the next ionization state.
C
C The first term will reflect collisions, or recombinations to the
C lower ionization state. The second term will reflect recombinations
C from the nexit ionization stage.
C
C Formerly we compared the collisional rate with the recombination rate.
C If collisional rate is much larger, we replace one of the equations with
C the ionization balance equation. We assume Omega/gl is unity, and assume
C that both the collison rate and recombination rate are inversely
C proportional to sqrt(T). Recombination coefficient is assumed to be 10^{-12}.
C
C	LIM=0					!Old code.
C	T1=HDKT*(EDGEC2(2)-EDGEC2(1))
C	DO K=1,ND
C	  COL=1.0D+05*C2(1,K)*EXP(T1/T(K))/DC2(K)
C	  IF(COL .GT. 1.0D+05)LIM=K
C	END DO
C
C Fac is the value by which dN1 is to exceed dNION before the
C equation is replaced. The optimal value is unknown.
C
	FAC=1.0D+05
C
	DIAG=(NUM_BNDS+1)/2
	LIM=0
	DO K=1,ND
	  IF(NUM_BNDS .EQ. ND)DIAG=K
	  IF( ABS(BA(EQC2,EQC2,DIAG,K))*C2(1,K) .GT.
	1        FAC*ABS(BA(EQC2,EQC2+NC2,DIAG,K))*DC2(K) )LIM=K
	END DO
	IF(LIM .EQ. 0)RETURN
	K=ICHRLEN(DESC)
	WRITE(LUER,100)DESC(1:K),LIM
	WRITE(53,100)DESC(1:K),LIM
100	FORMAT(1X,A,' g.s. eq. replaced by ionization eq. to d=',I3)
C
C In all cases, we replce the ground state equation.
C
	DO K=1,LIM
 	  T1=0.0D0
	  DO I=1,NC2
	     T1=T1+STEQ(EQC2+I-1,K)
	  END DO
	  WRITE(53,*)T1,STEQION(EQC2ION,K)
	  STEQ(EQC2,K)=STEQION(EQC2ION,K)
	  DO L=1,NUM_BNDS
	    DO J=1,NT
	      BA(EQC2,J,L,K)=BAION(EQC2ION,J,L,K)
	    END DO
	  END DO
	END DO
C
	RETURN
	END
