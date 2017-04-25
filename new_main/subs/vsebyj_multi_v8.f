!
! Subroutine to increment the variation matrix due to terms which
! depend directly on the intensity J. The Radiative equilibrium equation
! is not altered.
!
! Routine also increments the ionization equilibrium equations.
!
	SUBROUTINE VSEBYJ_MULTI_V8(ID,WSE,dWSEdT,
	1             HN,HNST,dlnHNST_dlnT,NLEV,
	1             DI,LOG_DIST,dlnDIST_dlnT,N_DI,ION_LEV,
	1             ED,T,JREC,dJRECdT,JPHOT,FIXED_T,
	1             NUM_BNDS,ND,DST,DEND)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 17-Jun-2016 - Changed to V8, added FIXED_T to call (added to main line: 04-Oct-2016).
! Altered 05-Apr-2011 - Changed to V7.
!                       LOG_DIST (rather than DIST) is passed in call.
!                       Changes done to faciliate modifications allowing lower temperatures.
!                       Major correction don 26-Nov-2010
!
! Altered : 08-Sep-2004 Bug fix. ION_V is now set to ION_EQ, rather than NLEV+ION_LEV.
!                         This is corect because of the way the important variables are
!                         assigned in CREATE_IV_LINKS.
! Altered : 01-Apr-2001 Changed to use STEQ_DATA_MOD
!                       Extensive chnages in call (changed to V6).
! Altered : 08-Jun-1995 EDGE frequency delted from call.
!                       Change from _V1 to _V2 as call changed.
! Created - May 1995 
!
	INTEGER ID		!Number of ionization stage
	INTEGER NLEV		!Numer of levls in HN
        INTEGER N_DI		!Number of levels in target ion
	INTEGER ION_LEV	!Super level target in ION
	INTEGER ND		!Number of depth points
        INTEGER NION		!Numer of Eqns. in ionization matrix.
!
! NB --- NION is the total number of ionic species i.e. for
! HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and CV).
!
	INTEGER NUM_BNDS,DST,DEND
!
	REAL*8 WSE(NLEV,ND),dWSEdT(NLEV,ND)
!
! Populations of species undergoing photoionization.
!
	REAL*8 HN(NLEV,ND),HNST(NLEV,ND),dlnHNST_dlnT(NLEV,ND)
!
! Ion populations.
!
	REAL*8 DI(N_DI,ND)
	REAL*8 LOG_DIST(N_DI,ND)
	REAL*8 dlnDIST_dlnT(N_DI,ND)
!
	REAL*8 ED(ND),T(ND)
	REAL*8 JREC(ND)
	REAL*8 dJRECdT(ND)
	REAL*8 JPHOT(ND)
	LOGICAL FIXED_T
!
! Constants for opacity etc.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables
!
	INTEGER J,K,L
	INTEGER NT
	INTEGER ION_V
	INTEGER ION_EQ
	REAL*8 T3
	REAL*8 B_RAT
	REAL*8 LOG_B_RAT
!
! REV_HNST referes to the LTE population  of the level defined with respect
! to the actual destination (target) level.
!
	REAL*8 REV_HNST
	REAL*8 WSE_BY_RJ,DI_FAC,ED_FAC,T_FAC
!
!
!
	IF(ION_LEV .EQ. 0)RETURN
!
	NT=SE(ID)%N_IV
        ION_EQ=SE(ID)%ION_LEV_TO_EQ_PNT(ION_LEV)
	ION_V=ION_EQ
	L=(NUM_BNDS/2)+1
!
!$OMP PARALLEL DO PRIVATE(J,K,LOG_B_RAT,B_RAT,WSE_BY_RJ,REV_HNST,T3,DI_FAC,ED_FAC,T_FAC)
!
	DO K=DST,DEND			!Which depth point.
	  IF(ION_LEV .NE. 1)THEN
	    LOG_B_RAT=LOG(DI(ION_LEV,K)/DI(1,K))+LOG_DIST(1,K)-LOG_DIST(ION_LEV,K)
	    B_RAT=0.0D0
	    IF(LOG_B_RAT .LT. 780.0D0)B_RAT=EXP(LOG_B_RAT)
	  ELSE
	    B_RAT=1.0D0
	    LOG_B_RAT=0.0D0
	  END IF 

	  DO J=1,NLEV			!Which equation (for S.E. only)
	    IF(WSE(J,K) .NE. 0)THEN
	      WSE_BY_RJ=WSE(J,K)*JPHOT(K)
	      SE(ID)%BA_PAR(J,J,K)=SE(ID)%BA_PAR(J,J,K)-WSE_BY_RJ
!
	      REV_HNST=HNST(J,K)*B_RAT
	      T3=REV_HNST*WSE(J,K)*JREC(K)
	      DI_FAC=T3/DI(ION_LEV,K)
	      ED_FAC=T3/ED(K)
	      SE(ID)%BA_PAR(J,ION_V,K)=SE(ID)%BA_PAR(J,ION_V,K)  + DI_FAC
	      SE(ID)%BA_PAR(J,NT-1,K) =SE(ID)%BA_PAR(J,NT-1,K)   + ED_FAC
!
! Include ionizations/recombinations implicitly in the rate equation
! of the target ion (eg He++(gs) for He+ ion/recoms ). 
!
	      SE(ID)%BA_PAR(ION_EQ,J,K)    =SE(ID)%BA_PAR(ION_EQ,J,K)     + WSE_BY_RJ
	      SE(ID)%BA_PAR(ION_EQ,ION_V,K)=SE(ID)%BA_PAR(ION_EQ,ION_V,K) - DI_FAC
	      SE(ID)%BA_PAR(ION_EQ,NT-1,K) =SE(ID)%BA_PAR(ION_EQ,NT-1,K)  - ED_FAC
!
	      IF(.NOT. FIXED_T)THEN
	        T_FAC=T3*( dlnHNST_dlnT(J,K) +
	1             (dlnDIST_dlnT(1,K)-dlnDIST_dlnt(ION_LEV,K)) )/T(K) +
	1             dWSEdT(J,K)*(REV_HNST*JREC(K)-HN(J,K)*JPHOT(K)) + 
	1             REV_HNST*WSE(J,K)*dJRECdT(K)
	        SE(ID)%BA_PAR(J,NT,K)   =SE(ID)%BA_PAR(J,NT,K)     + T_FAC
	        SE(ID)%BA_PAR(ION_EQ,NT,K)   =SE(ID)%BA_PAR(ION_EQ,NT,K)    - T_FAC 
	      END IF
!
	    END IF		!WSE(J,K) .NE. 0
	  END DO
	END DO
!$OMP END PARALLEL DO
!
	RETURN
	END
