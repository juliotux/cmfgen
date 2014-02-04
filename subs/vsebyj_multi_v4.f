C
C Subroutine to increment the variation matrix due to terms which
C depend directly on the intensity J. The Radiative equilibrium equation
C is not altered.
C
C Routine also increments the ionization equilibrium equations.
C
	SUBROUTINE VSEBYJ_MULTI_V4(BA,WSE,dWSEdT,
	1             HN,HNST,dlnHNST_dlnT,NLEV,
	1             DI,DIST,dlnDIST_dlnT,N_DI,ION_LEV,
	1             ED,T,JREC,dJRECdT,JPHOT,
	1             FRST_EQ,GS_ION_EQ,SPEC_EQ,
	1             NT,NUM_BNDS,ND,
	1             BAION,EQUAT,NION,DST,DEND)
	IMPLICIT NONE
C
C Altered : 08-Jun-1995 EDGE frequency delted from call.
C                       Change from _V1 to _V2 as call changed.
C Created - May 1995 
C
	INTEGER NLEV		!Numer of levls in HN
        INTEGER N_DI		!Number of levels in target ion
	INTEGER FRST_EQ	!Equation number for species
	INTEGER GS_ION_EQ	!Equation number of g.s target species
	INTEGER ION_LEV	!Super level target in ION
	INTEGER SPEC_EQ	!Equation number of abundance equation
	INTEGER NT		!Total number of levels
	INTEGER ND		!Number of depth points
	INTEGER EQUAT		!Equation number in ioization matrix
        INTEGER NION		!Numer of Eqns. in ionization matrix.
C
C NB --- NION is the total number of ionic species i.e. for
C HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and CV).
C
	INTEGER NUM_BNDS,DST,DEND
C
	REAL*8 BA(NT,NT,NUM_BNDS,ND),BAION(NION,NT,NUM_BNDS,ND)
	REAL*8 WSE(NLEV,ND),dWSEdT(NLEV,ND)
C
C Populations of species undergoing photoionization.
C
	REAL*8 HN(NLEV,ND),HNST(NLEV,ND),dlnHNST_dlnT(NLEV,ND)
C
C Ion populations.
C
	REAL*8 DI(N_DI,ND),DIST(N_DI,ND)
	REAL*8 dlnDIST_dlnT(N_DI,ND)
C
	REAL*8 ED(ND),T(ND)
	REAL*8 JREC(ND)
	REAL*8 dJRECdT(ND)
	REAL*8 JPHOT(ND)
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables
C
	INTEGER J,K,L,NJ,ION_EQ
	REAL*8 T3
	REAL*8 B_RAT
C
C REV_HNST referes to the LTE population  of the level defined with respect
C to the actual destination (target) level.
C
	REAL*8 REV_HNST
	REAL*8 WSE_BY_RJ,DI_FAC,ED_FAC,T_FAC
C
	IF(ION_LEV .EQ. 0)RETURN
C
	DO K=DST,DEND			!Which depth point.
	  L=(NUM_BNDS/2)+1
	  B_RAT=(DI(ION_LEV,K)/DIST(ION_LEV,K))*(DIST(1,K)/DI(1,K))
	  DO J=1,NLEV			!Which equation (for S.E. only)
	    IF(WSE(J,K) .NE. 0)THEN
	      NJ=J+FRST_EQ-1
	      WSE_BY_RJ=WSE(J,K)*JPHOT(K)
	      BA(NJ,NJ,L,K)=BA(NJ,NJ,L,K)-WSE_BY_RJ
C
	      REV_HNST=HNST(J,K)*B_RAT
	      T3=REV_HNST*WSE(J,K)*JREC(K)
	      DI_FAC=T3/DI(ION_LEV,K)
	      ED_FAC=T3/ED(K)
	      T_FAC=T3*( dlnHNST_dlnT(J,K) +
	1             (dlnDIST_dlnT(1,K)-dlnDIST_dlnt(ION_LEV,K)) )/T(K) +
	1             dWSEdT(J,K)*(REV_HNST*JREC(K)-HN(J,K)*JPHOT(K)) + 
	1             REV_HNST*WSE(J,K)*dJRECdT(K)
C
	      ION_EQ=GS_ION_EQ+(ION_LEV-1)
	      BA(NJ,ION_EQ,L,K)=BA(NJ,ION_EQ,L,K) + DI_FAC
	      BA(NJ,NT-1,L,K)  =BA(NJ,NT-1,L,K)   + ED_FAC
	      BA(NJ,NT,L,K)    =BA(NJ,NT,L,K)     + T_FAC
C
C Include ionizations/recombinations implicitly in the rate equation
C of the target ion (eg He++(gs) for He+ ion/recoms ). The rates are
C not included if the target ion is the final ionization state, as then
C the equation is the density constraint.
C
	      IF(ION_EQ .LT. SPEC_EQ)THEN
	        BA(ION_EQ,NJ,L,K)=BA(ION_EQ,NJ,L,K) + WSE_BY_RJ
	        BA(ION_EQ,ION_EQ,L,K)=BA(ION_EQ,ION_EQ,L,K) - DI_FAC
	        BA(ION_EQ,NT-1,L,K)=BA(ION_EQ,NT-1,L,K) - ED_FAC
	        BA(ION_EQ,NT,L,K)=BA(ION_EQ,NT,L,K) - T_FAC 
	      END IF		!ION_EQ .NE. 0
C
C NB: We do not increment the ionization equation for a species by
C     ionizations/recombinations from/to the lower ionization state.
C     This is satisfactory provided there are no transitions between
C     states differing by a charge of 2 --- such as occurs with
C     Auger ionization. In such cases the X-ray ionizations must be 
C     incorporated in a special way. 
C
	      IF(EQUAT .NE. 0)THEN
	        BAION(EQUAT,NJ,L,K)=BAION(EQUAT,NJ,L,K) - WSE_BY_RJ
	        BAION(EQUAT,ION_EQ,L,K)=BAION(EQUAT,ION_EQ,L,K) + DI_FAC
	        BAION(EQUAT,NT-1,L,K)=BAION(EQUAT,NT-1,L,K) + ED_FAC
	        BAION(EQUAT,NT,L,K)=BAION(EQUAT,NT,L,K) + T_FAC
	      END IF		!EQUAT .NE. 0
	    END IF		!WSE(J,K) .NE. 0
	  END DO
	END DO
C
	RETURN
	END
