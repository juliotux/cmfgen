C
C Subroutine to compute the value of the statistical equilibrium
C equations and the variation of the statistical equilibrium matrix for
C terms which are radiation field independent.
C
C This routine is specifically designed for the handling of super levels.
C That is, we treat the process in a large atom but assume that the populations
C can be described by a smaller set of levels.
C
C Routine also handles level dissolution.
C
C Notation:
C
C         We use _F to denote populations and variables for the FULL atom,
C            with all terms and levels treated separately.
C	  We use _S to denote populations and variables for the SMALL model
C            atom, with many terms and levels treated as one (i.e using
C            SUPER levels).
C
C
C The COLLISION routine that is called has a special FORM, which is distinct
C from that in STEQGEN_V2.
C
C NB - ZION is the charge on the ion - thus ZHYD=1.0D0
C
C Equation NW     : Radiative equilibrium
C Equation NW-1   : Charge conservation
C Equation EQPOP   : Population conservation (EQPOP-NST+1)
C
C Routine also increments the ionization equilibrium equations. Routine no
C longer works for NUM_BNDS=ND.
C
C At present only collisional ionizations to ground state are considered.
C
C NION is the the first dimension of STEQ[ION]. In general we
C NION would be the total number of ionic species.
C
C EQUAT gives the equation number for the species under consideration.
C
	SUBROUTINE STEQ_MULTI_V3(BA,SE,CNM,DCNM,ED,T,
	1       HN_S,HNST_S,dlnHNST_S_dlnT,N_S,DI_S,
	1       HN_F,HNST_F,W_F,A_F,FEDGE_F,G_F,LEVNAME_F,N_F,
	1       F_TO_S_MAPPING,POP,NEXT_PRES,ZION,
	1       SUB_PHOT,COL_FILE,OMEGA_GEN,
	1       NST,EQPOP,NT,NUM_BNDS,ND,
	1       BAION,STEQION,EQUAT,NION,DST,DEND)
	IMPLICIT NONE
C
C Altered 20-Sep-1999 : TMP_VEC_ED and TMP_VEC_COOL used in call to 
C                                                         SUBCOL_MULTI_V3
C Altered 14-Dec-1996 : SUB_PHOT replaces PHOT_FUN (superficial).
C Altered 15-Jun-1996 : T1 initialized before being passed to SUMBCOL_MULTI_V3.
C Altered 26-May-1996 : N_F_MAX removed. Now use dynamic memoery allocation
C                         for OMEGA_F etc.
C Altered 03-Nov-1995 : Version changed to _V3
C                       HN_F inserted in call to SUBCOL_MULTI_V3 (prev. _V2)
C
C Altered 10-Nov-1995 : Call to CUBCOL_MULTI_V2 updated.
C Altered 27-Oct-1995 : Call altered to handle new SUBCOL routine.
C                       Now _V2.
C Altered 07-Jun-1995 : Bug fix. Wrong values of HNLTE_S etc being
C                        passed to SUBCOL (effectively those at d=1).
C Created 16-May-1995 : Based on STEQGEN_V2
C
	EXTERNAL OMEGA_GEN,SUB_PHOT
C
	INTEGER EQPOP,NST,NT,ND,NUM_BNDS,EQUAT,NION,DST,DEND
	REAL*8 BA(2-NST:NT-NST+1,2-NST:NT-NST+1,NUM_BNDS,ND)
	REAL*8 SE(2-NST:NT-NST+1,ND)
	REAL*8 BAION(NION,2-NST:NT-NST+1,NUM_BNDS,ND),STEQION(NION,ND)
C
C CNM, and DCNM are used as work arrays. DCNM refers to dCNMdT
C
	INTEGER N_S,N_F
	REAL*8 CNM(N_S,N_S),DCNM(N_S,N_S)
C
	REAL*8 T(ND)
	REAL*8 ED(ND)
	REAL*8 DI_S(ND)
C
	REAL*8 HN_S(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
C
	REAL*8 HN_F(N_F,ND)
	REAL*8 HNST_F(N_F,ND)
	REAL*8 W_F(N_F,ND)
	REAL*8 A_F(N_F,N_F)
	REAL*8 FEDGE_F(N_F)
	REAL*8 G_F(N_F)
	CHARACTER*(*) LEVNAME_F(N_F),COL_FILE
	INTEGER F_TO_S_MAPPING(N_F)
	REAL*8 ZION
C
	REAL*8 POP(ND)		!Population of species.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
	LOGICAL NEXT_PRES
C
C Local variables.
C
	INTEGER EQION,IONE
	INTEGER I,J,K,L,M,NW
	REAL*8 T1,T2
	REAL*8 TMP_VEC_ED(1)
	REAL*8 TMP_VEC_COOL(1)
	PARAMETER (IONE=1)
C
	REAL*8 OMEGA_F(N_F,N_F)
	REAL*8 dln_OMEGA_F_dlnT(N_F,N_F)
C
	NW=NT-NST+1
	EQION=N_S+1			!Ion equation : Local dimensions
C
	DO I=DST,DEND			!Which depth
	  M=(NUM_BNDS/2)+1
C
C Compute collisional cross-sections (and their T derivatives)
C We call this routine ND times so the CNM and DCM arrays can be
C smaller (i.e. no ND dimension).
C
C OMEGA_F,dln_OMEGA_dlnT are work arrays only.
C T1 is returned with the toal cooling rate. Not used in this routine.
C We use arrays (even though of length 1) so that some F90 compilers 
C don't give an error message because a scaler is passed a vector.
C
	  TMP_VEC_ED(1)=1.0D0		!Electron density
	  TMP_VEC_COOL(1)=0.0D0		!Initialize cooling rate even 
!                                                        though not used here.
C
	  CALL SUBCOL_MULTI_V3(
	1         OMEGA_F,dln_OMEGA_F_dlnT,
	1         CNM,DCNM,
	1         HN_S(1,I),HNST_S(1,I),dlnHNST_S_dlnT(1,I),N_S,
	1         HN_F(1,I),HNST_F(1,I),W_F(1,I),FEDGE_F,
	1         A_F,G_F,LEVNAME_F,N_F,
	1         ZION,SUB_PHOT,COL_FILE,OMEGA_GEN,
	1         F_TO_S_MAPPING,TMP_VEC_COOL,T(I),TMP_VEC_ED,IONE)
C
	  DO J=1,N_S			!Which S.E. equation
	    DO K=1,N_S			!Which variable
	      IF(K.EQ.J)THEN
	 	T1=0.0D0
	        DO L=1,N_S
		  T1=T1+CNM(J,L)
	        END DO
		BA(J,K,M,I)=BA(J,K,M,I)-T1*ED(I)
	      ELSE
	        BA(J,K,M,I)=BA(J,K,M,I)+ED(I)*CNM(K,J)
	      END IF
	    END DO
	    T1=0.0
	    T2=0.0
C
	    DO L=1,N_S
	      T1=T1+( HN_S(L,I)*CNM(L,J)-HN_S(J,I)*CNM(J,L) )
	      T2=T2+( HN_S(L,I)*DCNM(L,J)-HN_S(J,I)*DCNM(J,L) )
	    END DO
C
	    BA(J,N_S+1,M,I)=BA(J,N_S+1,M,I) +
	1       HNST_S(J,I)*ED(I)/DI_S(I)*CNM(J,J)
	    BA(J,NW-1,M,I)=BA(J,NW-1,M,I) +
	1      T1+CNM(J,J)*(2*HNST_S(J,I)-HN_S(J,I))
	    BA(J,NW,M,I)=BA(J,NW,M,I) +
	1      ED(I)*( T2+(HNST_S(J,I)-HN_S(J,I))*DCNM(J,J)+
	1      CNM(J,J)*HNST_S(J,I)*dlnHNST_S_dlnT(J,I)/T(I) )
	    SE(J,I)=SE(J,I)+(T1+(HNST_S(J,I)-HN_S(J,I))*CNM(J,J))*ED(I)
	  END DO
C
C EQION is the ion equation in local dimensions (Recall that BA and STEQ
C arrays don't begin at 1).  Adding NST-1 to EQION offsets it so
C dimensions begin at 1, and hence is directly comparable to EQPOP.
C
	  IF(EQION+(NST-1) .LT. EQPOP)THEN
	    T1=0.0D0
	    T2=0.0D0
	    DO J=1,N_S
	      T1=T1+(HNST_S(J,I)-HN_S(J,I))*CNM(J,J)
	      T2=T2+HNST_S(J,I)*CNM(J,J)
	      BA(EQION,J,M,I)=BA(EQION,J,M,I)+CNM(J,J)*ED(I)
	      BA(EQION,NW,M,I)=BA(EQION,NW,M,I) -
	1            ED(I)*(  (HNST_S(J,I)-HN_S(J,I))*DCNM(J,J) +
	1            CNM(J,J)*HNST_S(J,I)*dlnHNST_S_dlnT(J,I)/T(I)  )
	    END DO
	    SE(EQION,I)=SE(EQION,I)-T1*ED(I)
	    BA(EQION,N_S+1,M,I)=BA(EQION,N_S+1,M,I) -
	1            T2*ED(I)/DI_S(I)
	    BA(EQION,NW-1,M,I)=BA(EQION,NW-1,M,I) - T1 - T2
	  END IF
C
	  IF(EQUAT .NE. 0)THEN
	    T1=0.0D0
	    T2=0.0D0
	    DO J=1,N_S
	      T1=T1+(HNST_S(J,I)-HN_S(J,I))*CNM(J,J)
	      T2=T2+HNST_S(J,I)*CNM(J,J)
	      BAION(EQUAT,J,M,I)=BAION(EQUAT,J,M,I)-CNM(J,J)*ED(I)
	      BAION(EQUAT,NW,M,I)=BAION(EQUAT,NW,M,I) +
	1            ED(I)*(  (HNST_S(J,I)-HN_S(J,I))*DCNM(J,J) +
	1            CNM(J,J)*HNST_S(J,I)*dlnHNST_S_dlnT(J,I)/T(I)  )
	    END DO
	    STEQION(EQUAT,I)=STEQION(EQUAT,I)+T1*ED(I)
	    BAION(EQUAT,N_S+1,M,I)=BAION(EQUAT,N_S+1,M,I) +
	1            T2*ED(I)/DI_S(I)
	    BAION(EQUAT,NW-1,M,I)=BAION(EQUAT,NW-1,M,I) + T1 + T2
	  END IF
C
C N_S+1  - Pop. conservation
C NW-1   - Charge conservation
C NW     - Radiative Equilibrium
C
	  T1=0.0
	  DO L=1,N_S
	    BA(NW-1,L,M,I)=BA(NW-1,L,M,I)+(ZION-1.0D0)
	    BA(EQPOP-NST+1,L,M,I)=BA(EQPOP-NST+1,L,M,I)+1.0D0
	    T1=T1+HN_S(L,I)
	  END DO
	  SE(EQPOP-NST+1,I)=SE(EQPOP-NST+1,I)+T1
	  SE(NW-1,I)=SE(NW-1,I)+(ZION-1.0D0)*T1
C
C We only include DI in the population and charge conservation equations
C if the higher ionization species is not present. Necessary to do this as
C DI is the ground state of the next species. We also correct the
C conservation equation for POP if the higher ionization species is not
C present.
C
C Note the charge on DI is ZION.
C
	  J=EQPOP-NST+1			!Carbon conservation equation.
	  IF(.NOT. NEXT_PRES)THEN
	    SE(J,I)=SE(J,I)+DI_S(I)-POP(I)
	    BA(J,J,M,I)=BA(J,J,M,I)+1.0D0
	    SE(NW-1,I)=SE(NW-1,I)+DI_S(I)*ZION
	    BA(NW-1,J,M,I)=BA(NW-1,J,M,I)+ZION
	  END IF
C
	END DO
C
	RETURN
	END
