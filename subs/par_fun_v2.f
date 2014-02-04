C
C Routine to evaluate the population of the highest ionization stage
C assuming fixed departure coefficients. U AND PHI are defined in
C Mihalas (1978), page 110. U and PHI are also evaluated for
C next ionization stage. This will be overwritten if an additional
C ionization stage is present.
C
	SUBROUTINE PAR_FUN_V2(U,PHI,ZPFN,HIGH_POP,
	1             HE2,HE2LTE,W_HE2,DHE2,EDGE,GHE2,GION,ZION,
	1             T,N,ND,NSPEC,NSPEC_MAX,PRES)
	IMPLICIT NONE
C
C Altered 07-MAr-2006 - 2.07D-22 pulled into exponential.
C Altered 25-May-1996 - ERROR_LU inserted.
C Altered 16-Jan-1995 - Occupation probabilities included in calculation of
C                         partition function. CALL altered, but still V2.
C Altered 15-Jan-1995 - HIGH_POP and DHE2 inserted. Allowest population of
C                         the highest ionization state to be set through
C                         successive calls. CALL altered.
C Altered 27-Sep-1990 - HE2 is now converted to departure coefficients
C                       in routine.
C
	INTEGER N,ND,NSPEC,NSPEC_MAX
	REAL*8 HE2(N,ND)
	REAL*8 HE2LTE(N,ND)
	REAL*8 W_HE2(N,ND)
	REAL*8 DHE2(ND)
	REAL*8 EDGE(N)
	REAL*8 GHE2(N)
	REAL*8 GION,ZION
	REAL*8 U(ND,NSPEC_MAX),PHI(ND,NSPEC_MAX),ZPFN(NSPEC_MAX),T(ND)
	REAL*8 HIGH_POP(ND)
	LOGICAL PRES
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
C
	INTEGER I,J
	REAL*8 T1,T2,RGU
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	IF(.NOT. PRES)RETURN
	LUER=ERROR_LU()
C
	NSPEC=NSPEC+1		!Update species number.
	IF(NSPEC+1 .GT. NSPEC_MAX)THEN
	  WRITE(LUER,*)'Error in PAR_FUN ---- NSPEC too small'
	  RETURN
	END IF
C
C Convert HE2 array to departure coefficients.
C
	DO I=1,ND
	  DO J=1,N
	    HE2(J,I)=HE2(J,I)/HE2LTE(J,I)
	  END DO
	  HIGH_POP(I)=DHE2(I)
	  IF(W_HE2(1,I) .NE. 1.0)THEN
	    WRITE(LUER,*)'Error in PAR_FUN --- occupation probability for ',
	1               'ground state must be zero.'
	    WRITE(LUER,*)'N=',N
	    STOP
	  END IF
	END DO
C
	ZPFN(NSPEC)=ZION-1
	ZPFN(NSPEC+1)=ZION
	RGU=2.07078D-22
	RGU=LOG(RGU)
	T1=HDKT*EDGE(1)
	DO I=1,ND
	  U(I,NSPEC)=GHE2(1)
	  U(I,NSPEC+1)=GION
	  PHI(I,NSPEC+1)=0.0D0
	  PHI(I,NSPEC)=HE2(1,I)*EXP( RGU+T1/T(I) )/T(I)/SQRT(T(I))
	  DO J=2,N
	    T2=HDKT*(EDGE(J)-EDGE(1))
	    U(I,NSPEC)=U(I,NSPEC) + GHE2(J)*EXP(T2/T(I))*
	1                W_HE2(J,I)*(HE2(J,I)/HE2(1,I))
	  END DO
	END DO
C
	RETURN
	END
