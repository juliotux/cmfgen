!
! Subroutine to compute the opacity variation due to BOUND-FREE and FREE-FREE
! transitions as a function of the level populations for a general ion.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Routine can handle ionizations to differnt super levels.
!
! Notation:
!
!         We use _F to denote poupulations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote poupulations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
	SUBROUTINE VAR_OP_V10(VCHI,VETA,VCHI_NOT_IMP,VETA_NOT_IMP,
	1             HN_S,HNST_S,LOG_HNST_S,dlnHNST_S_dlnT,N_S,
	1	      HNST_F_ON_S,EDGE_F,N_F,F_TO_S_MAPPING,
	1             DI_S,LOG_DIST_S,dlnDIST_S_dlnT,N_DI,
	1             PHOT_ID,ION_LEV,ED,T,EMHNUKT,IMP_VAR,
	1             NU,Z,ID,IONFF,
	1             EQHN,GS_ION_EQ,NT,ND,LST_DEPTH_ONLY)
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
!
! Altered 23-Oct-2016 - Inckude PHOT_DIS_PARAMETER
! Altered 11-Nov-2011 - Only include FF variation when PHOT_ID=1. When this is true,
!                         we now do it for ALL levels.
! Altered 05-Apr-2011 - Changed to V5.
!                         HNST_F_ON_S (rather than HNST_F) is passed in call.
!                         LOG_HNSR_S passed in call.
!                         Changes done to faciliate modifications allowing lower temperaturs.
!                         Most of editing done 6-Feb-2010
! Altered 21-May-2002 - Bug fix. In first bound-free section, VCHI was always being
!                         accesed, instead of PCHI (which points to VCHI or 
!                         VCHI_NOT_IMP).
! Altered 13-Feb-2002 - VCHI_NOT_IMP, VETA_NOT_IMP inserted.
! Altered 05-May-1998 - Bug fix --- ALPHA_VEC not correctly zeroed when level
!                         dissolution is switched off.
! Altered 15-Dec-1997 - MOD_LEV_DIS_BLK replaces include file. Level
!                         dissolution can be switched off completely.
! Altered 25-Aug-1996 - Bound-free section altered to improve speed.
C                       As major changes, called V4 (13-Dec-1996)
C Altered 28-May-1996 - GFF_VAL no dynamically dimensioned.
C                       PHOT_GEN_BLEND_V2 now called. LST_DEPTH_ONLY option
C                         installed in this call. NB: previos version was
C                         incorrectly returning value at ND=1 when
C                         LST_DEPTH_ONLY was true.
C Altered 25-Nov-1995 - Bug fixed with the variation of the b-f emissivity.
C Altered 17-Oct-1995 - _V3 (Call changed)
C                       Now allows opacity due to ionizations/recobinations to
C                         excited states to be taken automatically into
C                         account.
C Altered 07-Jun-1995 - Bug fix. Computing the variation in the bound-free
C                       opacities and emissivities using DI_F instead of DI.
C                       DI_F left in call, but not used.
C Altered 02-Jun-1995 - LST_DEPTH_ONLY option installed. Allows variation in
C                         ETA and CHI to be computed at the last depth point
C                         only. This is usefule when computing dTdR.
C                         Call was changed, hence call _V2
C Created 16-May-1995 - Based on VAR_GEN_V3
C
	INTEGER ID
	INTEGER N_S
	INTEGER N_F
	INTEGER N_DI		!Number of levels in ion
	INTEGER EQHN
	INTEGER GS_ION_EQ
	INTEGER NT,ND
	LOGICAL IONFF			!Include free-free opacity for level?
	LOGICAL LST_DEPTH_ONLY	!for computing dTdR
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	REAL*8, TARGET :: VCHI(NT,ND)		!VCHI(I,K)=dCHI(K)/dN(I,K)
	REAL*8, TARGET :: VETA(NT,ND)		!VETA(I,K)=dETA(K)/dN(I,K)
	REAL*8, TARGET :: VCHI_NOT_IMP(NT,ND)
	REAL*8, TARGET :: VETA_NOT_IMP(NT,ND)
C
	REAL*8 HN_S(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
	REAL*8 LOG_HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
C
	REAL*8 HNST_F_ON_S(N_F,ND)
	REAL*8 EDGE_F(N_F)
	INTEGER F_TO_S_MAPPING(N_F)
C
C Ion population information.
C
	REAL*8 DI_S(N_DI,ND)
	REAL*8 LOG_DIST_S(N_DI,ND)
	REAL*8 dlnDIST_S_dlnT(N_DI,ND)
C
	LOGICAL IMP_VAR(NT)
C
	INTEGER PHOT_ID		!Photoionization ID
	INTEGER ION_LEV		!target level in ION for ionizations.
C
	REAL*8 NU			!Frequency (10^15 Hz)
	REAL*8 Z			!Charge on ion
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperatuure (10^4 K)
	REAL*8 EMHNUKT(ND)		!exp(-hv/kT)
C
C Local variables.
C
	REAL*8, POINTER :: PCHI(:,:)
	REAL*8, POINTER :: PETA(:,:)
C
	REAL*8 GFF_VAL(ND)		!Used as work vector
	REAL*8 LOG_DI_RAT(ND)		!Used as work vector
	REAL*8 DT_TERM(ND)		!Used as work vector
	REAL*8 HDKT_ON_T(ND)		!Used as work vector
	REAL*8 POP_SUM(ND)
C
	REAL*8 YDIS(ND)			!Constant for computing level dissolution/
	REAL*8 XDIS(ND)			!Constant for computing level dissolution/
	REAL*8 DIS_CONST(N_F)		!Constant appearing in dissolution formula.
	REAL*8 ALPHA_VEC(N_F)		!Photionization cross-section
	REAL*8 VCHI_TMP(N_F,ND)
	REAL*8 SUM_ION; REAL*8 SUM_T1; REAL*8 SUM_T2
	REAL*8 SUM_ION_NOT_IMP; REAL*8 SUM_T1_NOT_IMP; REAL*8 SUM_T2_NOT_IMP
	REAL*8 NEFF,ZION_CUBED,T1,T2
C
	INTEGER ND_LOC
	INTEGER I                     !Used as level index (same) in atom.
	INTEGER L			!index of level in full atom.
	INTEGER K_ST,K		!Used as depth index.
	INTEGER GENLEV		!Level index in VCHI, VETA
	INTEGER EQION			!Ion variable in VCHI,VETA
	INTEGER NO_NON_ZERO_PHOT
C
	REAL*8 TCHI1,TCHI2,TETA1,TETA2,TETA3
	REAL*8 HNUONK,ALPHA
C
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	HNUONK=HDKT*NU
	EQION=GS_ION_EQ+(ION_LEV-1)
!
! If EQION is not important, all free-free and bound-free processes
! are added to the IMP variation.
!
	IF( IMP_VAR(EQION) )THEN
	  PCHI=>VCHI
	  PETA=>VETA
	ELSE
	  PCHI=>VCHI_NOT_IMP
	  PETA=>VETA_NOT_IMP
	END IF
C
C ND_LOC indicates the number of depth points we are going to compute the
C opacity at.
C
C K_ST indicates the depth point to start, and is either 1, or ND.
C
	IF(LST_DEPTH_ONLY)THEN
	  ND_LOC=1
	  K_ST=ND
	ELSE
	  ND_LOC=ND
	  K_ST=1
	END IF
C
C Free-free processes
C
	IF( IONFF .AND. PHOT_ID .EQ. 1)THEN
C
C Compute free-free gaunt factors. Replaces call to GFF in following DO loop.
C
	  IF(LST_DEPTH_ONLY)THEN
	    CALL GFF_VEC(GFF_VAL(ND),NU,T(ND),Z,ND_LOC)
	    POP_SUM(ND)=SUM(DI_S(:,ND))
	  ELSE
	    CALL GFF_VEC(GFF_VAL,NU,T,Z,ND_LOC)
	    POP_SUM=SUM(DI_S,1)
	  END IF
C
	  TCHI1=CHIFF*Z*Z/( NU**3 )
	  TETA1=CHIFF*Z*Z*TWOHCSQ
!
!$OMP PARALLEL DO PRIVATE(ALPHA,TCHI2,TETA2,I,K,L)
	  DO K=K_ST,ND
	    ALPHA=GFF_VAL(K)/SQRT(T(K))
C
	    TCHI2=TCHI1*ALPHA
	    PCHI(NT-1,K)=PCHI(NT-1,K)+POP_SUM(K)*TCHI2*(1.0D0-EMHNUKT(K))
	    PCHI(NT,K)=PCHI(NT,K)+ED(K)*POP_SUM(K)*TCHI2/T(K)*( -0.5D0+(0.5D0-HNUONK/T(K))*EMHNUKT(K) )
C
	    TETA2=TETA1*ALPHA*EMHNUKT(K)
	    PETA(NT-1,K)=PETA(NT-1,K)+TETA2*POP_SUM(K)
	    PETA(NT,K)  =PETA(NT,K)  +TETA2*POP_SUM(K)*ED(K)*(HNUONK/T(K)-0.5D0)/T(K)
!
	    TCHI2=TCHI2*ED(K)*(1.0D0-EMHNUKT(K))
	    TETA2=TETA2*ED(K)
	    DO I=1,N_DI
	      L=EQION+I-1
	      PCHI(L,K)=PCHI(L,K)+TCHI2
	      PETA(L,K)=PETA(L,K)+TETA2
	    END DO
	  END DO
!$OMP END PARALLEL DO
!
	END IF
C 
C
C Now add in BOUND-FREE contributions. We first compute, via a subroutine call,
C a vector containing the cross-section for each level; This will decrease
C the compuation time.
C
C If PHOT_ID=1 and DISSOLUTION is switched ON, ALPHA_VEC contians the 
C threshold cross-section when NU < EDGE.
C
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU,EDGE_F,N_F,PHOT_ID,L_TRUE)
	ELSE
	  CALL SUB_PHOT_GEN(ID,ALPHA_VEC,NU,EDGE_F,N_F,PHOT_ID,L_FALSE)
	END IF
!
	NO_NON_ZERO_PHOT=COUNT(ALPHA_VEC .GT. 0)
	IF(NO_NON_ZERO_PHOT .EQ. 0)RETURN
C
C DIS_CONST is the constant K appearing in the expression for level dissolution.
C A negative value for DIS_CONST implies that the cross-section is zero.
C
C**** NB: If NU < EDGE but there is no dissolution, we MUST set ALPHA_VEC to
C         zero, as in the loops to evaluate VCHI we only check ALPHA_VEC
C         and DIS_CONST.
C
	DIS_CONST(1:N_F)=-1.0D0
	IF(MOD_DO_LEV_DIS .AND. PHOT_ID .EQ. 1)THEN
	  ZION_CUBED=Z*Z*Z
	  DO I=1,N_F
	    IF(NU .LT. EDGE_F(I) .AND. ALPHA_VEC(I) .NE. 0)THEN
	      NEFF=SQRT(3.289395D0*Z*Z/(EDGE_F(I)-NU))
	      IF(NEFF .GT. 2*Z)THEN
	        T1=MIN(1.0D0,16.0D0*NEFF/(1+NEFF)/(1+NEFF)/3.0D0)
	        DIS_CONST(I)=( T1*ZION_CUBED/(NEFF**4) )**1.5D0
	      ELSE
	        ALPHA_VEC(I)=0.0D0
	      END IF
	    END IF
	  END DO
	END IF
C
C Compute dissolution vectors that are independent of level.
C
	IF(MOD_DO_LEV_DIS)THEN
	  DO K=K_ST,ND
	    YDIS(K)=1.091D0*(X_LEV_DIS(K)+4.0D0*(Z-1)*A_LEV_DIS(K))*
	1                 B_LEV_DIS(K)*B_LEV_DIS(K)
	    XDIS(K)=B_LEV_DIS(K)*X_LEV_DIS(K)
	  END DO
	END IF
C
C Compute quantities to save execution time.
C
C NB:  TMP_HNST=HNST(I,K)*(DI(ION_LEV,K)/DIST(ION_LEV,K))*(DIST(1,K)/DI(1,K))
C
C
	DO K=K_ST,ND
	  LOG_DI_RAT(K)=LOG(DI_S(ION_LEV,K)/DI_S(1,K))-LOG_DIST_S(ION_LEV,K)+LOG_DIST_S(1,K)
	  DT_TERM(K)=( 1.5D0 + (dlnDIST_S_dlnT(ION_LEV,K)-dlnDIST_S_dlNT(1,K)) )/T(K)
	  HDKT_ON_T(K)=HDKT/T(K)
	END DO
C
	TETA1=TWOHCSQ*NU*NU*NU
	IF(NO_NON_ZERO_PHOT .LT. 4)THEN
!	IF(NO_NON_ZERO_PHOT .LT. 2*(ND-K_ST+1))THEN
C	IF(NO_NON_ZERO_PHOT .LT. 1000)THEN
	  DO I=1,N_F
	    L=F_TO_S_MAPPING(I)
	    IF(ALPHA_VEC(I) .GT. 0)THEN
	      GENLEV=L+EQHN-1
	      IF( IMP_VAR(GENLEV) )THEN
	        PCHI=>VCHI
	        PETA=>VETA
	      ELSE
	        PCHI=>VCHI_NOT_IMP
	        PETA=>VETA_NOT_IMP
	      END IF
C
C NB: We divide by DI and not DI_F since we want the variation with
C     resepect to DI. Although the LTE_F pops depend directly on DI_F
C     DI is proportional to DI_F. Thus
C
C     d(LTE_F)/dDI = LTE_F/Di_F *(DI_F/DI)=LTE_F/DI
C
	      DO K=K_ST,ND
	        ALPHA=ALPHA_VEC(I)*HNST_F_ON_S(I,K)
	        IF(DIS_CONST(I) .GE. 0.0D0)THEN
	          T1=7.782D0+XDIS(K)*DIS_CONST(I)
	          T2=T1/(T1+YDIS(K)*DIS_CONST(I)*DIS_CONST(I))
	          ALPHA=ALPHA*T2
	          IF(T2 .LT. PHOT_DIS_PARAMETER)ALPHA=0.0D0
	        END IF
	        IF(ALPHA .GT. 0.0D0)THEN
	          PCHI(GENLEV,K)=PCHI(GENLEV,K)+ALPHA
	          TCHI1=ALPHA*EXP(LOG_HNST_S(L,K)+LOG_DI_RAT(K)-HDKT*NU/T(K))
	          TCHI2=DT_TERM(K)+HDKT_ON_T(K)*(EDGE_F(I)-NU)/T(K)
	          PCHI(EQION,K)=PCHI(EQION,K)-TCHI1/DI_S(ION_LEV,K)
	          PCHI(NT-1,K)=PCHI(NT-1,K)-TCHI1/ED(K)
	          PCHI(NT,K)=PCHI(NT,K) + TCHI1*TCHI2 - HN_S(L,K)*ALPHA*
	1           (1.5D0+HDKT_ON_T(K)*EDGE_F(I)+dlnHNST_S_dlnT(L,K))/T(K)
C
C NB. The cross-section ALPHA is contained in TCHI1.
C
	          TETA3=TETA1*TCHI1
	          PETA(EQION,K)=PETA(EQION,K)+TETA3/DI_S(ION_LEV,K)
	          PETA(NT-1,K)=PETA(NT-1,K)+TETA3/ED(K)
	          PETA(NT,K)=PETA(NT,K)-TETA3*TCHI2
	        END IF
	      END DO
	    END IF
	  END DO
C
C 
C
	ELSE
	  VCHI_TMP(:,:)=0.0D0
!
!$OMP PARALLEL DO PRIVATE(SUM_ION, SUM_T1, SUM_T2,
!$OMP1 SUM_ION_NOT_IMP, SUM_T1_NOT_IMP, SUM_T2_NOT_IMP,
!$OMP1 I, L, ALPHA, T1, T2, TCHI1,TCHI2) IF(K_ST .EQ. 1)
	  DO K=K_ST,ND
	    SUM_ION=0.0D0; SUM_T1=0.0D0;   SUM_T2=0.0D0
	    SUM_ION_NOT_IMP=0.0D0; SUM_T1_NOT_IMP=0.0D0;   SUM_T2_NOT_IMP=0.0D0
	    DO I=1,N_F
	      L=F_TO_S_MAPPING(I)
	      IF(ALPHA_VEC(I) .GT. 0.0D0)THEN
	        ALPHA=ALPHA_VEC(I)*HNST_F_ON_S(I,K)
	        IF(DIS_CONST(I) .GE. 0.0D0)THEN
	          T1=7.782D0+XDIS(K)*DIS_CONST(I)
	          T2=T1/(T1+YDIS(K)*DIS_CONST(I)*DIS_CONST(I))
	          ALPHA=ALPHA*T2
	          IF(T2 .LT. PHOT_DIS_PARAMETER)ALPHA=0.0D0
	        END IF
C
	        IF(ALPHA .GT. 0.0D0)THEN
	          VCHI_TMP(I,K)=ALPHA
	          TCHI1=ALPHA*EXP(LOG_HNST_S(L,K)+LOG_DI_RAT(K)-HDKT*NU/T(K))
	          TCHI2=DT_TERM(K)+HDKT_ON_T(K)*(EDGE_F(I)-NU)/T(K)
	          IF(IMP_VAR(L+(EQHN-1)))THEN
	            SUM_ION=SUM_ION+TCHI1
	            SUM_T1=SUM_T1+TCHI1*TCHI2
	            SUM_T2=SUM_T2+HN_S(L,K)*ALPHA*
	1             (1.5D0+HDKT_ON_T(K)*EDGE_F(I)+dlnHNST_S_dlnT(L,K))
	          ELSE
	            SUM_ION_NOT_IMP=SUM_ION_NOT_IMP+TCHI1
	            SUM_T1_NOT_IMP=SUM_T1_NOT_IMP+TCHI1*TCHI2
	            SUM_T2_NOT_IMP=SUM_T2_NOT_IMP+HN_S(L,K)*ALPHA*
	1               (1.5D0+HDKT_ON_T(K)*EDGE_F(I)+dlnHNST_S_dlnT(L,K))
	          END IF
	        END IF
	      END IF
	    END DO
C
	    VCHI(EQION,K)=VCHI(EQION,K)-SUM_ION/DI_S(ION_LEV,K)
	    VCHI(NT-1,K)=VCHI(NT-1,K)-SUM_ION/ED(K)
	    VCHI(NT,K)=VCHI(NT,K)+SUM_T1-SUM_T2/T(K)
!
	    VCHI_NOT_IMP(EQION,K)=VCHI_NOT_IMP(EQION,K)-SUM_ION_NOT_IMP/DI_S(ION_LEV,K)
	    VCHI_NOT_IMP(NT-1,K)=VCHI_NOT_IMP(NT-1,K)-SUM_ION_NOT_IMP/ED(K)
	    VCHI_NOT_IMP(NT,K)=VCHI_NOT_IMP(NT,K)+SUM_T1_NOT_IMP-SUM_T2_NOT_IMP/T(K)
C
	    T1=TETA1
	    VETA(EQION,K)=VETA(EQION,K)+T1*SUM_ION/DI_S(ION_LEV,K)
	    VETA(NT-1,K)=VETA(NT-1,K)+T1*SUM_ION/ED(K)
	    VETA(NT,K)=VETA(NT,K)-T1*SUM_T1
!
	    VETA_NOT_IMP(EQION,K)=VETA_NOT_IMP(EQION,K)+T1*SUM_ION_NOT_IMP/DI_S(ION_LEV,K)
	    VETA_NOT_IMP(NT-1,K)=VETA_NOT_IMP(NT-1,K)+T1*SUM_ION_NOT_IMP/ED(K)
	    VETA_NOT_IMP(NT,K)=VETA_NOT_IMP(NT,K)-T1*SUM_T1_NOT_IMP
	  END DO
!$OMP END PARALLEL DO
!
!	  DO I=1,N_F
!	    L=F_TO_S_MAPPING(I)
!	    GENLEV=L+EQHN-1
!	    IF( IMP_VAR(GENLEV) )THEN
!	      DO K=K_ST,ND
!	        VCHI(GENLEV,K)=VCHI(GENLEV,K)+VCHI_TMP(I,K)
!	      END DO
!	    ELSE
!	      DO K=K_ST,ND
!	        VCHI_NOT_IMP(GENLEV,K)=VCHI_NOT_IMP(GENLEV,K)+VCHI_TMP(I,K)
!	      END DO
!	    END IF
!	  END DO
!
!$OMP PARALLEL DO PRIVATE(GENLEV) IF(K_ST .EQ. 1)
	  DO K=K_ST,ND
	    DO I=1,N_F
	      GENLEV=F_TO_S_MAPPING(I)+EQHN-1
	      IF( IMP_VAR(GENLEV) )THEN
	        VCHI(GENLEV,K)=VCHI(GENLEV,K)+VCHI_TMP(I,K)
	      ELSE
	        VCHI_NOT_IMP(GENLEV,K)=VCHI_NOT_IMP(GENLEV,K)+VCHI_TMP(I,K)
	      END IF
	    END DO
	  END DO
!$OMP END PARALLEL DO
	END IF
!
! Nullify pointer assignments.
!
	NULLIFY (PCHI)
	NULLIFY (PETA)
!
	RETURN
	END
