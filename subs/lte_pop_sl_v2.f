!
! Routine to compute the LTE population for a super level. The derivative with
! respect to dln_/dlnT is also returned.
!
! NB: The LTE populations of the levels in the FULL atom must have been
!       computed previously.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
	SUBROUTINE LTE_POP_SL_V2(HNST_S,LOG_HNST_S,dlnHNST_S_dlnT,N_S,
	1           HNST_F,LOG_HNST_F,HNST_F_ON_S,EDGE_F,F_TO_S_MAPPING,N_F,
	1           SPEC_PRES,T,ND)
	IMPLICIT NONE
!
! Altered: 5-Apr-2011: MAX_LOG_LTE_POP parameter introduced.
!                      HNST_F set to zero if HNST_S=0 
!                      Based on LTE_POP_SL_V1 (original coding early 2011).
!                      Call changed as LOG_HNST_S/F variables introduced.
!
!
	INTEGER N_S,N_F,ND
!
	REAL*8 HNST_S(N_S,ND)
	REAL*8 LOG_HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
!
	REAL*8 HNST_F(N_F,ND)
	REAL*8 HNST_F_ON_S(N_F,ND)
	REAL*8 LOG_HNST_F(N_F,ND)
	REAL*8 EDGE_F(N_F)
	INTEGER F_TO_S_MAPPING(N_F)
	LOGICAL SPEC_PRES
!
	REAL*8 T(ND)
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables.
!
	REAL*8, PARAMETER :: MAX_LOG_LTE_POP=600.0D0
	REAL*8 SCALE_FAC(N_S)
	REAL*8 T1
	INTEGER I,L,K
!
	IF(.NOT. SPEC_PRES)RETURN
!
! The LTE population is simply a linear sum over the combined levels.
! Because of floating overflow, we operate in LOG space. Note that
! HSNT_F/HNST_S should always be well defined, since the F levels are always
! realtively close in energy to the supler level S.
!
	DO K=1,ND
!
	  DO L=1,N_S
	    HNST_S(L,K)=0.0D0
	    dlnHNST_S_dlnT(L,K)=0.0D0
	    SCALE_FAC(L)=0.0D0
	  END DO
!
	  DO I=1,N_F
	    L=F_TO_S_MAPPING(I)
	    IF(SCALE_FAC(L) .EQ. 0.0D0)SCALE_FAC(L)=LOG_HNST_F(I,K)
	    HNST_S(L,K)=HNST_S(L,K)+EXP(LOG_HNST_F(I,K)-SCALE_FAC(L))
	  END DO
!
	  DO L=1,N_S
	    LOG_HNST_S(L,K)=LOG(HNST_S(L,K))+SCALE_FAC(L)
	    IF(LOG_HNST_S(L,K) .LT. MAX_LOG_LTE_POP)THEN
	      HNST_S(L,K)=EXP(LOG_HNST_S(L,K))
	    ELSE
	      HNST_S(L,K)=0.0D0
	    END IF
	  END DO
!
! Make HNST_F consitent zeroing.
!
	  DO I=1,N_F
	    L=F_TO_S_MAPPING(I)
	    IF(HNST_S(L,K) .EQ. 0.0D0)HNST_F(I,K)=0.0D0
	  END DO
!
	  DO I=1,N_F
	    L=F_TO_S_MAPPING(I)
	    T1=EXP(LOG_HNST_F(I,K)-LOG_HNST_S(L,K))
	    HNST_F_ON_S(I,K)=T1
	    dlnHNST_S_dlnT(L,K)=dlnHNST_S_dlnT(L,K) + T1*EDGE_F(I)
	  END DO
!
	  DO L=1,N_S
	    dlnHNST_S_dlnT(L,K) = -1.5D0 - HDKT*dlnHNST_S_dlnT(L,K)/T(K)
	  END DO
	END DO
!
	RETURN
	END
