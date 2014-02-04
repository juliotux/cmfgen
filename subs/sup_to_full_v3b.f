!
! General routine to
!
!   1:  compute the LTE populations of the levels in the FULL atom given the
!       POPULATIONS in the super level atoms,
!   2:  compute the populations of the levels in the full atom given the
!       POPULATIONS in the supel level atoms,
!   3:  compute the populations of the levels in the super-level atom,
!   4:  compute the occupation probailities, and
!   5:  routine the ion population used to compute the LTE population in the
!       FULL model atom.
!
! for C2.
!
! Routine is written for any 2 successive ionization stages --- not just
! C2 and CIII.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
	SUBROUTINE SUP_TO_FULL_V3b(
	1   C2_F,C2LTE_F,LOG_C2LTE_F,C2LTE_F_ON_S,W_C2_F,EDGEC2_F,
	1   GC2_F,F_TO_S_MAP_C2,INT_SEQ_C2,NC2_F,DIC2_F,GIONC2_F,
	1   C2_S,C2LTE_S,LOG_C2LTE_S,NC2_S,DIC2_S,ZC2,C2_PRES,
	1   EDGECIII_F,GCIII_F,F_TO_S_MAP_CIII,NCIII_F,
	1   CIII_PRES,T,ED,ND)
	IMPLICIT NONE
!
! Altered: 5-Apr-2011: MAX_LOG_LTE_POP parameter introduced.
!                         HNST_F set to zero if HNST_S=0 
!                         Based on SUP_TO_FULL_V3 (original coding early 2011).
!                         Call changed as LOG_C2LTE_S/_F variables introduced.
! Altered 20-Feb-2010 : Changed computation of LTE population to increase 
!                         dynamic range.
! Altered 06-Nov-2009 : No longer gnerates error for truncated interp seq.
! Altered 27-may-1996 : Dynamic arrays installed for GION,EDGE_S, SUM and CNT.
!
! Altered 02-Jan-1996. INT_SEQ_C2 inserted in call.
!                      This variable will us to estimate non constant
!                        departure coeficients in a super level. This
!                        has been done specifically to improve the treatment
!                        of hign n lines in HeII and HeI.
!
! Altered 24-Oct-1995. GIONC2_S deleted from call.
!
	INTEGER ND
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperature 10^4K
	REAL*8 DIC2_S(ND)		!Ion density (Super levels)
	REAL*8 DIC2_F(ND)		!Ion density (Full model atom)
	REAL*8 ZC2			!Ion charge
!
	INTEGER NC2_F
	REAL*8 C2_F(NC2_F,ND)
	REAL*8 C2LTE_F(NC2_F,ND)
	REAL*8 LOG_C2LTE_F(NC2_F,ND)
	REAL*8 C2LTE_F_ON_S(NC2_F,ND)
	REAL*8 W_C2_F(NC2_F,ND)
	REAL*8 EDGEC2_F(NC2_F)
	REAL*8 GC2_F(NC2_F)
	INTEGER F_TO_S_MAP_C2(NC2_F)
	INTEGER INT_SEQ_C2(NC2_F)
	REAL*8 GIONC2_F
!
	INTEGER NC2_S
	REAL*8 C2_S(NC2_S,ND)
	REAL*8 C2LTE_S(NC2_S,ND)
	REAL*8 LOG_C2LTE_S(NC2_S,ND)
!
	LOGICAL C2_PRES,CIII_PRES
!
	INTEGER NCIII_F
	REAL*8 EDGECIII_F(NCIII_F)
	REAL*8 GCIII_F(NCIII_F)
	INTEGER F_TO_S_MAP_CIII(NCIII_F)
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	REAL*8, PARAMETER :: MAX_LOG_LTE_POP=600.0D0
	REAL*8 B,B1,B2,T1
	REAL*8 X,Y,RGU
	INTEGER I,K,L,M,LU_ER
	INTEGER INT_SL
!
! Dynamic arrays.
!
	REAL*8 GION(ND)
	REAL*8 EDGE_S(NC2_S)	
	REAL*8 SCALE_FAC(NC2_S)
	REAL*8 SUM(NC2_S)
	INTEGER CNT(NC2_S)
!
	IF(.NOT. C2_PRES)RETURN
!
! Compute the partition function of the CIII g.s.. NB. If CIII is not present
! the statistical weight of the ion ground state must just be GIONC2_F.
!
	IF(CIII_PRES)THEN
	  DO K=1,ND
	    GION(K)=0.0D0
	    DO I=1,NCIII_F
	      IF(F_TO_S_MAP_CIII(I) .EQ. 1)THEN
	        GION(K)=GION(K)+GCIII_F(I)*EXP( HDKT*(EDGECIII_F(I)-
	1                                         EDGECIII_F(1))/T(K) )
	      END IF
	    END DO
	  END DO
	ELSE
	  DO K=1,ND
	    GION(K)=GIONC2_F
	  END DO
	END IF
!
! Compute the ion density used to compute LTE populations in the full atom.
! This is essentially the ground-state population.
!
	DO K=1,ND
	  DIC2_F(K)=DIC2_S(K)*GIONC2_F/GION(K)
	END DO
!
! Compute the occupation probabilities.
!
	CALL OCCUPATION_PROB(W_C2_F,EDGEC2_F,ZC2,NC2_F,ND)
!
! Since no the effective statistical weight, can now compute the LTE
! populations of the levels in the full atom.
!
	C2LTE_F=0.0D0
	DO K=1,ND
	  X=HDKT/T(K)
	  RGU=2.07078D-22*ED(K)*DIC2_S(K)*( T(K)**(-1.5) )/GION(K)
	  RGU=LOG(RGU)
	  DO I=1,NC2_F
	    LOG_C2LTE_F(I,K)=LOG(W_C2_F(I,K)*GC2_F(I)) + EDGEC2_F(I)*X + RGU
	    IF(LOG_C2LTE_F(I,K) .LT. MAX_LOG_LTE_POP)C2LTE_F(I,K)=EXP(LOG_C2LTE_F(I,K))    !24-Feb-2011
	  END DO
	END DO
!
! Compute the LTE pops in the atom with super levels, after initializing them.
!
	DO K=1,ND
	  DO I=1,NC2_S
	    C2LTE_S(I,K)=0.0D0
	    SCALE_FAC(I)=0.0D0
	  END DO
!
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    IF(SCALE_FAC(L) .EQ. 0.0D0)SCALE_FAC(L)=LOG_C2LTE_F(I,K)
	    C2LTE_S(L,K)=C2LTE_S(L,K)+EXP(LOG_C2LTE_F(I,K)-SCALE_FAC(L))
	  END DO
!
	  DO L=1,NC2_S
	    LOG_C2LTE_S(L,K)=LOG(C2LTE_S(L,K))+SCALE_FAC(L)
            IF(LOG_C2LTE_S(L,K) .LT. MAX_LOG_LTE_POP)THEN                 !Changed 24-Feb-2010
	      C2LTE_S(L,K)=EXP(LOG_C2LTE_S(L,K))
	    ELSE
	      C2LTE_S(L,K)=0.0D0
	    END IF
	  END DO
!
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    IF(C2LTE_S(L,K) .EQ. 0.0D0)C2LTE_F(I,K)=0.0D0
	    C2LTE_F_ON_S(I,K)=EXP(LOG_C2LTE_F(I,K)-LOG_C2LTE_S(L,K))
	  END DO
!    
	END DO
!
! Can now compute populations in full atom.
!
	DO K=1,ND
!
	  DO L=1,NC2_S
	    CNT(L)=0
	    EDGE_S(L)=0.0D0
	  END DO
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    EDGE_S(L)=EDGE_S(L)+EDGEC2_F(I)*C2LTE_F_ON_S(I,K)
	    CNT(L)=CNT(L)+1
	  END DO
!
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    IF(CNT(L) .EQ. 1 .OR. INT_SEQ_C2(I) .EQ. 0)THEN
	      C2_F(I,K)=C2_S(L,K)*C2LTE_F_ON_S(I,K)
	    ELSE
!
! Find the closest level of the same interpolating sequence. We first
! attempt to bracket (in energy) the enegry of the level whose population
! is to be determined.
!
	      INT_SL=0
	      IF(EDGEC2_F(I) .LE. EDGE_S(L))THEN
	        M=I+1
	        DO WHILE(M .LE. NC2_F .AND. INT_SL .EQ. 0)
	          IF(INT_SEQ_C2(M) .EQ. INT_SEQ_C2(I) .AND.
	1             F_TO_S_MAP_C2(I) .NE. F_TO_S_MAP_C2(M))THEN
	            INT_SL=M
	          END IF
	          M=M+1
	        END DO
	      END IF
	      IF(INT_SL .EQ. 0)THEN
	        M=I-1
	        DO WHILE(M .GE. 1 .AND. INT_SL .EQ. 0)
	          IF(INT_SEQ_C2(M) .EQ. INT_SEQ_C2(I) .AND.
	1             F_TO_S_MAP_C2(I) .NE. F_TO_S_MAP_C2(M))THEN
	            INT_SL=M
	          END IF
	          M=M-1
	        END DO
	      END IF
!
! This might occur when we dont use the full atom.
!
	      IF(INT_SL .EQ. 0)THEN
	        C2_F(I,K)=C2_S(L,K)*C2LTE_F_ON_S(I,K)
	      ELSE
	        INT_SL=F_TO_S_MAP_C2(INT_SL)
!
! As this should only be high levels, overflow may not be an issue.
!
	        T1=LOG(EDGE_S(L)/EDGEC2_F(I)) / LOG(EDGE_S(L)/EDGE_S(INT_SL))
	        B1=C2_S(INT_SL,K)/C2LTE_S(INT_SL,K)
	        B2=C2_S(L,K)/C2LTE_S(L,K)
                B=T1*B1 + (1.0D0-T1)*B2
!
! Constrain the interpolation. Hopefully this is not necessary.
!
	        IF(B1 .LE. 1. .AND. B2 .LE. 1 .AND. B .GT. 1)B=1.0
	        IF(B1 .GE. 1. .AND. B2 .GE. 1 .AND. B .LT. 1)B=1.0
	        IF(B .LT. 0)B=B1
	        C2_F(I,K)=C2LTE_F(I,K)*B
	      END IF
	    END IF
	  END DO
!
! The revised departure coefficents now needs adjusting so that the total
! population of the C2 levels in the full atom matches the corresponding
! super level. This correction will generally be small.
!
	  DO L=1,NC2_S
	    SUM(L)=0.0D0
	  END DO
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    SUM(L)=SUM(L)+C2_F(I,K)
	  END DO
	  DO I=1,NC2_F
	    L=F_TO_S_MAP_C2(I)
	    C2_F(I,K)=C2_F(I,K)*(C2_S(L,K)/SUM(L))
	  END DO
!
	END DO
!
	RETURN
	END
