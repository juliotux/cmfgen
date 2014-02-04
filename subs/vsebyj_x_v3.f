C
C Subroutine to increment the variation matrix due to terms which
C depend directly on the intensity J. The Radiative equilibrium equation
C is not altered. This routine is for X-ray ionization only, where
C 2 electrons are ejected.
C
C Routine also increments the ionization equilibrium equations.
C
	SUBROUTINE VSEBYJ_X_V3(BA,WSE_X,
	1             HN_A,HNST_A,dlnHNST_AdlnT,N_A,
	1             HN_B,HNST_B,EDGE_B,N_B,
	1             DI,ED,T,EMHNUKT,NU,ML,RJ,
	1             EQ_A,ION_EQ,SPEC_EQ,NCF,NT,NUM_BNDS,ND,
	1             BAION,EQ_A_BAL,EQ_B_BAL,NION,DST,DEND)
	IMPLICIT NONE
C 
C Altered 27-OCt-1995 : Changed to be compatible with super levels.
C                        Call chnaged --- Now _V3.
C Altered 06-Mar-1995 : Dimensions of WSE_X changed to (N,ND) from (N,NCF).
C                       _V2 append to call.
C Testing 22-Jul-1994 : Minor mods.
C Created 19-Jul-1993 : Based on VSEBYJ_COM and EVALSE_QWVJ.
C
	INTEGER N_A		!Number of levles in ionizations state i
	INTEGER N_B		!Number of levles in ionizations state i+1
	INTEGER EQ_A		!Eqn. for ground state of ion. state i
	INTEGER ION_EQ	!Eqn. # of target species [gs. of (i+2)th ]
	INTEGER SPEC_EQ	!Eqn. # of abundance equation
C
	INTEGER NCF		!Number of freuencies in NU 
	INTEGER ML		!Inidcates current freqency in vector NU
	INTEGER NT		!Total number of levels
	INTEGER ND		!Number of depth points
C
	INTEGER EQ_A_BAL	!Eqn. # for     ith ion. stage in ion. matrix
	INTEGER EQ_B_BAL	!Eqn. # for (i+1)th ion. stage in ion. matrix
        INTEGER NION		!Numer of Eqns. in ionization matrix.
C
C NB --- NION is the total number of ionic species i.e. for
C HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and CV).
C
	INTEGER NUM_BNDS,DST,DEND
	REAL*8 BA(NT,NT,NUM_BNDS,ND)
	REAL*8 BAION(NION,NT,NUM_BNDS,ND)
	REAL*8 WSE_X(N_A,ND)
C
C _A refers to quantities associated with the atom WITH super levels.
C
	REAL*8 HN_A(N_A,ND)		!Pops. of ith ionzation stage
	REAL*8 HNST_A(N_A,ND)		!LTE   "    "  "      "       "
	REAL*8 dlnHNST_AdlnT(N_A,ND)	
C
C _B refers to quantities assiciated with the FULL atom of the next 
C  ionization stage.
C
	REAL*8 HN_B(N_B,ND)		!Pops. of (i+1)th ionization stage
	REAL*8 HNST_B(N_B,ND)		!LTE   "   "    "       "        "
	REAL*8 EDGE_B(N_B)
C
	REAL*8 DI(ND)			!Ion densit for B levels.
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperature in 10^4K.
	REAL*8 EMHNUKT(ND)		!EXP(-hv/kT)
	REAL*8 RJ(ND)			!Mean intensity
	REAL*8 NU(NCF)			!frequency (10^15 Hz)
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables
C
	INTEGER I,J,L,NJ
	REAL*8 T1,T3,BSTIM
	REAL*8 WSE_BY_RJ,DI_FAC,ED_FAC,T_FAC
C           
	T1=TWOHCSQ*( NU(ML)**3 )
	L=(NUM_BNDS/2)+1
	DO I=DST,DEND			!Which depth point.
	  BSTIM=(T1+RJ(I))*EMHNUKT(I)*HNST_B(1,I)/HN_B(1,I)
	  IF(NUM_BNDS .EQ. ND)L=I
	  DO J=1,N_A			!Which equation (for S.E. only)
	    IF(WSE_X(J,I) .NE. 0)THEN
	      NJ=J+EQ_A-1
	      WSE_BY_RJ=WSE_X(J,I)*RJ(I)
	      BA(NJ,NJ,L,I)=BA(NJ,NJ,L,I)-WSE_BY_RJ
C
C NB: In the following, the factor of 2.0 for ED_FAC arrises because 
C the rate is prop. to 
C
C        HNST_A(J,I)*HNST_B(1,I) .
C
C The variation with respect to HN_B(1,I) is intrinsically zero, since 
C HNST_A(J,I)/HN_B(1,I) is independent of HN_B(1,I)
C
	      T3=HNST_A(J,I)*WSE_X(J,I)*BSTIM
	      DI_FAC=T3/DI(I)
	      ED_FAC=2.0D0*T3/ED(I)
	      T_FAC=T3*( HDKT*(NU(ML)-EDGE_B(1)-1.5D0)/T(I) +
	1             dlnHNST_AdlnT(J,I) )/T(I)
C
	      BA(NJ,ION_EQ,L,I)=BA(NJ,ION_EQ,L,I) +DI_FAC
	      BA(NJ,NT-1,L,I)  =BA(NJ,NT-1,L,I)   +ED_FAC
	      BA(NJ,NT,L,I)    =BA(NJ,NT,L,I)     +T_FAC
C
C Include ionizations/recombinations implicitly in the rate equation
C of the target ion (eg He++(gs) for He+ ion/recoms ). The rates are
C not included if the target ion is the final ionization state, as then
C the equation is the density constraint.
C
	      IF(ION_EQ .LT. SPEC_EQ)THEN
	        BA(ION_EQ,NJ,L,I)=BA(ION_EQ,NJ,L,I)+WSE_BY_RJ
	        BA(ION_EQ,ION_EQ,L,I)=BA(ION_EQ,ION_EQ,L,I)-DI_FAC
	        BA(ION_EQ,NT-1,L,I)=BA(ION_EQ,NT-1,L,I)-ED_FAC
	        BA(ION_EQ,NT,L,I)=BA(ION_EQ,NT,L,I)-T_FAC
	      END IF		!ION_EQ .NE. 0
C
C Add in effect of X-rays to ionization/recombination balance equation.
C X-rays ionize from level i to i+2. However, for numerical stability
C we analytically cancel lower phot/recom. rates from rate equations. 
C As a consequence X-ray rates are only included in consecutive ionization
C stages.
C
	      IF(EQ_A_BAL .NE. 0)THEN
	        BAION(EQ_A_BAL,NJ,L,I)=BAION(EQ_A_BAL,NJ,L,I)-WSE_BY_RJ
	        BAION(EQ_A_BAL,ION_EQ,L,I)=BAION(EQ_A_BAL,ION_EQ,L,I)+DI_FAC
	        BAION(EQ_A_BAL,NT-1,L,I)=BAION(EQ_A_BAL,NT-1,L,I)+ED_FAC
	        BAION(EQ_A_BAL,NT,L,I)=BAION(EQ_A_BAL,NT,L,I)+T_FAC
	      END IF		!EQ_A_BAL .NE. 0
	      IF(EQ_B_BAL .NE. 0)THEN
	        BAION(EQ_B_BAL,NJ,L,I)=BAION(EQ_B_BAL,NJ,L,I)-WSE_BY_RJ
	        BAION(EQ_B_BAL,ION_EQ,L,I)=BAION(EQ_B_BAL,ION_EQ,L,I)+DI_FAC
	        BAION(EQ_B_BAL,NT-1,L,I)=BAION(EQ_B_BAL,NT-1,L,I)+ED_FAC
	        BAION(EQ_B_BAL,NT,L,I)=BAION(EQ_B_BAL,NT,L,I)+T_FAC
	      END IF		!EQ_B_BAL .NE. 0
	    END IF		!WSE(J,ML) .NE. 0
	  END DO
	END DO
C
	RETURN
	END
