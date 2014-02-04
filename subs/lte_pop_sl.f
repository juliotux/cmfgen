C
C Routine to compute the LTE population for a super level. The derivative with
C respect to dln_/dlnT is also returned.
C
C NB: The LTE populations of the levels in the FULL atom must have been
C       computed previously.
C
C This routine is specifically designed for the handling of super levels.
C That is, we treat the process in a large atom but assume that the populations
C can be described by a smaller set of levels.
C
C Notation:
C
C         We use _F to denote populations and variables for the FULL atom,
C            with all terms and levels treated separately.
C	  We use _S to denote populations and variables for the SMALL model
C            atom, with many terms and levels treated as one (i.e using
C            SUPER levels).
C
	SUBROUTINE LTE_POP_SL(HNST_S,dlnHNST_S_dlnT,N_S,
	1                       HNST_F,EDGE_F,F_TO_S_MAPPING,N_F,SPEC_PRES,
	1                       T,ND)
	IMPLICIT NONE
C
	INTEGER N_S,N_F,ND
C
	REAL*8 HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
C
	REAL*8 HNST_F(N_F,ND)
	REAL*8 EDGE_F(N_F)
	INTEGER F_TO_S_MAPPING(N_F)
	LOGICAL SPEC_PRES
C
	REAL*8 T(ND)
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables.
C
	INTEGER I,L,K
C
	IF(.NOT. SPEC_PRES)RETURN
C
	DO K=1,ND
	  DO L=1,N_S
	    HNST_S(L,K)=0.0D0
	    dlnHNST_S_dlnT(L,K)=0.0D0
	  END DO
	END DO
C
C The LTE population is simply a linear sum over the combined levels.
C
	DO K=1,ND
	  DO I=1,N_F
	    L=F_TO_S_MAPPING(I)
	    HNST_S(L,K)=HNST_S(L,K)+HNST_F(I,K)
	    dlnHNST_S_dlnT(L,K)=dlnHNST_S_dlnT(L,K) +
	1                          HNST_F(I,K)*EDGE_F(I)
	  END DO
	END DO
C
	DO K=1,ND
	  DO L=1,N_S
	    dlnHNST_S_dlnT(L,K)=-1.5D0 -
	1          HDKT*dlnHNST_S_dlnT(L,K)/T(K)/HNST_S(L,K)
	  END DO
	END DO
C
	RETURN
	END
