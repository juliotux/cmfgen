C
C Routine to compute the collisional recobination and ionization rates, and
C the collisional cooling rate for ions with super levels.
C
	SUBROUTINE COLCOOL_SL_V4(CPR,CRR,COOL,CNM,DCNM,
	1             HN_S,HNST_S,dlnHNST_S_dlnT,N_S,
	1             HN_F,HNST_F,A_F,W_F,EDGE_F,G_F,LEVNAME_F,
	1             F_TO_S_MAP,N_F,
	1             ZION,ID,COL_FILE,OMEGA_COL,ED,T,ND)
	IMPLICIT NONE
C
C Altered 24-Dec-1996 : SUB_PHOT replaces PHOT_FUN (superficial).
C Altered 24-May-1995 : N_F_MAX removed using F90
C Altered 01-Feb-1996 : Routine now calls SUBCOL_MULTI_V3 (not _V2). The _V3
C                        routine was introdued when interpolationg between
C                        levels in a given super level.
C                        Changed to _V3 for consistency with SUBCOL.
C Created 07-JUn-1995 :Based on COLGENCOOL
C
	EXTERNAL OMEGA_COL
	EXTERNAL ERROR_LU
	INTEGER ERROR_LU	
C
	INTEGER ID
	INTEGER N_S,N_F,ND
	REAL*8 COOL(ND),CRR(ND),CPR(ND)
	REAL*8 CNM(N_S,N_S),DCNM(N_S,N_S)
	REAL*8 ZION,T(ND),ED(ND)
C
	REAL*8 HN_S(N_S,ND),HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
C
	REAL*8 HN_F(N_F,ND),HNST_F(N_F,ND)	
	REAL*8 W_F(N_F,ND)
	REAL*8 A_F(N_F,N_F),EDGE_F(N_F),G_F(N_F)
	INTEGER F_TO_S_MAP(N_F)
	CHARACTER*(*) COL_FILE,LEVNAME_F(N_F)
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables.
C
	REAL*8 H,TMP_ED
	INTEGER I,J,IONE
	PARAMETER (IONE=1)
C
	REAL*8 OMEGA_F(N_F,N_F)
	REAL*8 dln_OMEGA_F_dlnT(N_F,N_F)
C
	H=6.6261965D-12		!H*1.0E+15  (1.0E+15 due to times frequency)
	TMP_ED=1.0D0
C
	DO I=1,ND			!Which depth
C                            
C Compute collisional cross-sections (and their T derivatives)
C
	  COOL(I)=0.0D0
	  CALL SUBCOL_MULTI_V4(OMEGA_F,dln_OMEGA_F_dlNT,
	1          CNM,DCNM,
	1          HN_S(1,I),HNST_S(1,I),dlnHNST_S_dlnT(1,I),N_S,
	1          HN_F(1,I),HNST_F(1,I),W_F(1,I),EDGE_F,
	1          A_F,G_F,LEVNAME_F,N_F,
	1          ZION,ID,COL_FILE,OMEGA_COL,
	1          F_TO_S_MAP,COOL(I),T(I),TMP_ED,IONE)
C                        
	  CPR(I)=0.0D0
	  CRR(I)=0.0D0
	  COOL(I)=COOL(I)*ED(I)*H
C
	  DO J=1,N_S				!Level
	    CPR(I)=CPR(I)+HN_S(J,I)*ED(I)*CNM(J,J)
	    CRR(I)=CRR(I)+HNST_S(J,I)*ED(I)*CNM(J,J)
	  END DO
C
	END DO
C
	RETURN
	END
