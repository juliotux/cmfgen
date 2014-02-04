C
C Subroutine to increment the statistical equilibrium equations for each 
C depth point given the value of the mean intensity at each depth point.
C
C Subroutine also increments the QFV_R AND QFV_P matrices that describe the 
C variation of the SE quations with respect to RJ.
C
C Routine also increments the ionization equilibrium equations.
C
C Routine is for X-ray ionizations in which 2 electrons are ejected.
C
	SUBROUTINE EVALSE_X_QWVJ_V3(STEQ,QFV_R,QFV_P,WSE_X,
	1                    HN_A,HNST_A,N_A,HN_B,HNST_B,N_B,
	1                    JREC,JPHOT,EQ_A,ION_EQ,SPEC_EQ,NT,ND,
	1                    STEQION,QFVION_R,QFVION_P,EQ_A_BAL,EQ_B_BAL,NION)
	IMPLICIT NONE
C
C Altered 17-Sep-1997 : QFV matix split into QFV_R and QFV_B so that a
C                         constant cross-section can be handelled.
C Altered 25-Feb-1996 : Major bug fix: X-ray rates were being added to
C                         wrong equation (ION_EQ qas okay).
C Altered 06-Mar-1995 : Dimensioniong of WSE changed to (N,ND) from (N,NCF).
C                        _V2 append to name.
C                       Call unchanged.
C Created 19-Jul-1993 : Based on EVALSE_QWVJ
C
	INTEGER N_A		!Number of levles in ionizations state i
	INTEGER N_B		!Number of levles in ionizations state i+1
	INTEGER EQ_A		!Eqn. for ground state of ion. state i
	INTEGER ION_EQ	!Eqn. # of target species [gs. of (i+2)th ]
	INTEGER SPEC_EQ	!Eqn. # of abundance equation
C
	INTEGER NT		!Total number of levels
	INTEGER ND		!Number of depth points
C
	INTEGER EQ_A_BAL	!Eqn. # for     ith ion. stage in ion. matrix
	INTEGER EQ_B_BAL	!Eqn. # for (i+1)th ion. stage in ion. matrix
        INTEGER NION		!Numer of Eqns. in ionization matrix.
C
C NB --- NION is the total number of ionic species i.e. for
C HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and 
C CV [if no CVI]).
C
	REAL*8 STEQ(NT,ND)
	REAL*8 STEQION(NION,ND)
	REAL*8 QFV_R(NT,ND)			!Recombination weight
	REAL*8 QFVION_R(NION,ND)
	REAL*8 QFV_P(NT,ND)			!Photoiozation weight
	REAL*8 QFVION_P(NION,ND)
C
	REAL*8 HN_A(N_A,ND)		!    Pops. of ith ionzation stage
	REAL*8 HNST_A(N_A,ND)		!LTE   "    "  "      "       "
	REAL*8 HN_B(N_B,ND)		!    Pops. of (i+1)th ionization stage
	REAL*8 HNST_B(N_B,ND)		!LTE   "   "    "       "        "
C
C WSE_X is the quadrature weight for X-ray ionization (with 2e ejected) for
C ionization state i [final product is (i+1)].
C
	REAL*8 WSE_X(N_A,ND)
	REAL*8 JREC(ND)			! Int (2h/c2v^3+J)*EXP(-hv/kT)/v dv
	REAL*8 JPHOT(ND)		! Int J/v dv
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables.
C
	INTEGER I,J
	REAL*8 NETR
	REAL*8 SUM_SE
	REAL*8 SUM_VJ_R,SUM_VJ_P
	REAL*8 J_B_ION,B_ION
C
C The net ionization (collisional and radaitive) to the last ionization stage
C must be zero from the sum of the previous equilibrum equations. Hence
C there is no need for a rate equation for the final species - it is 
C preserved for the abundance equation.
C
C Note that the product HNST_A(j, )*HNST_B(1, )/HN_B(1, ) is effectively the
C LTE population  of the state j with respect to the g.s. of the (I+2)th
C ionization stage.
C
	DO J=1,ND
	  SUM_SE=0.0D0
	  SUM_VJ_R=0.0D0
	  SUM_VJ_P=0.0D0
	  B_ION=HNST_B(1,J)/HN_B(1,J)		!1/b
	  J_B_ION=JREC(J)*B_ION
	  DO I=1,N_A
	    NETR=WSE_X(I,J)*( HNST_A(I,J)*J_B_ION-HN_A(I,J)*JPHOT(J) )
	    SUM_SE=SUM_SE+NETR
	    STEQ(I+EQ_A-1,J)=STEQ(I+EQ_A-1,J)+NETR
	    QFV_R(I+EQ_A-1,J)=QFV_R(I+EQ_A-1,J) + WSE_X(I,J)*B_ION*HNST_A(I,J)
	    SUM_VJ_R=SUM_VJ_R + WSE_X(I,J)*B_ION*HNST_A(I,J)
	    QFV_P(I+EQ_A-1,J)=QFV_P(I+EQ_A-1,J) + WSE_X(I,J)*HN_A(I,J)
	    SUM_VJ_P=SUM_VJ_P + WSE_X(I,J)*HN_A(I,J)
	  END DO
C
C Include effects of X-rays on equation of target ion. For K shell 
C ionizations of species  with more than 3 electrons, as considered here,
C 2 electrons are given off. Thus ION_EQ should refer to the g.s. of the
C (i+2)th ionization stage.
C
	  IF(ION_EQ .LT. SPEC_EQ)THEN
	    STEQ(ION_EQ,J)=STEQ(ION_EQ,J)-SUM_SE
	    QFV_R(ION_EQ,J)=QFV_R(ION_EQ,J)-SUM_VJ_R
	    QFV_P(ION_EQ,J)=QFV_P(ION_EQ,J)-SUM_VJ_P
	  END IF
C
C Add in effect of X-rays to ionization/recombination balance equation.
C X-rays ionize from level i to i+2. However, for numerical stability
C we analytically cancel lower phot/recom. rates from rate equations. 
C As a consequence X-ray rates are only included consecutive ionization
C stages.
C
	  IF(EQ_A_BAL .NE. 0)THEN
	    STEQION(EQ_A_BAL,J)=STEQION(EQ_A_BAL,J)+SUM_SE
	    QFVION_R(EQ_A_BAL,J)=QFVION_R(EQ_A_BAL,J)+SUM_VJ_R
	    QFVION_P(EQ_A_BAL,J)=QFVION_P(EQ_A_BAL,J)+SUM_VJ_P
	  END IF
	  IF(EQ_B_BAL .NE. 0)THEN
	    STEQION(EQ_B_BAL,J)=STEQION(EQ_B_BAL,J)+SUM_SE
	    QFVION_R(EQ_B_BAL,J)=QFVION_R(EQ_B_BAL,J)+SUM_VJ_R
	    QFVION_P(EQ_B_BAL,J)=QFVION_P(EQ_B_BAL,J)+SUM_VJ_P
	  END IF
	END DO
C
	RETURN
	END
