C
C General routine to:
C
C   1:  Compute the LTE populations of the levels in the super level atom
C       given the POPULATIONS in the FULL level atom.
C   2:  Compute the ion level population for the super-level atom.
C
C Routine is written for any 2 successive ionization stages --- not just
C C2 and CIII.
C
C Notation:
C
C         We use _F to denote populations and variables for the FULL atom,
C            with all terms and levels treated separately.
C	  We use _S to denote populations and variables for the SMALL model
C            atom, with many terms and levels treated as one (i.e using
C            SUPER levels).
C
	SUBROUTINE FULL_TO_SUP(
	1   C2_S,NC2_S,DC2_S,C2_PRES,
	1   C2_F,F_TO_S_MAP_C2,NC2_F,DC2_F,
	1   CIII_S,NCIII_S,CIII_PRES,ND)
	IMPLICIT NONE
C
	INTEGER ND
	REAL*8 DC2_S(ND)		!Ion density (Super levels)
	REAL*8 DC2_F(ND)		!Ion density (Full model atom)
C
	INTEGER NC2_F
	REAL*8 C2_F(NC2_F,ND)
	INTEGER F_TO_S_MAP_C2(NC2_F)
C
	INTEGER NC2_S
	REAL*8 C2_S(NC2_S,ND)
	LOGICAL C2_PRES
C
	INTEGER NCIII_S
	REAL*8 CIII_S(NCIII_S,ND)
	LOGICAL CIII_PRES
C
C Local variables.
C
	INTEGER I,K,L
C
	IF(.NOT. C2_PRES)RETURN
C
	IF(CIII_PRES)THEN
	  DO K=1,ND
	    DC2_S(K)=CIII_S(1,K)
	  END DO
	ELSE
	  DO K=1,ND
	    DC2_S(K)=DC2_F(K)
	  END DO
	END IF
C
	DO K=1,ND
	  DO I=1,NC2_S
	    C2_S(I,K)=0.0D0
	  END DO
	END DO
C
	DO K=1,ND
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    C2_S(L,K)=C2_S(L,K)+C2_F(I,K)
	  END DO
	END DO
C
	RETURN
	END
