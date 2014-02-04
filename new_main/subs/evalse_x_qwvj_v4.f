!
! Subroutine to increment the statistical equilibrium equations for each 
! depth point given the value of the mean intensity at each depth point.
!
! Subroutine also increments the QFV_R AND QFV_P matrices that describe the 
! variation of the SE quations with respect to RJ.
!
! Routine also increments the ionization equilibrium equations.
!
! Routine is for X-ray ionizations in which 2 electrons are ejected.
!
	SUBROUTINE EVALSE_X_QWVJ_V4(ID,WSE_X,
	1                    HN_A,HNST_A,N_A,
	1                    HN_B,HNST_B,N_B,ION_EQ_IN_BA,
	1                    JREC,JPHOT,ND,NION)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 12-Oct-2003 : If JREC is zero, the recombination term is not computed.
!                         This is to avoid floating overflows at low temperatures.
!                         In practice, the X-ray recombination term will be effectively zero
!                           at such temperatures, and hence can be neglected.
! Altered 12-Apr-2001 : Changed to utilize STEQ_DATA_MOD
!                       Changed to V4.
! Altered 17-Sep-1997 : QFV matix split into QFV_R and QFV_B so that a
!                         constant cross-section can be handelled.
! Altered 25-Feb-1996 : Major bug fix: X-ray rates were being added to
!                         wrong equation (ION_EQ qas okay).
! Altered 06-Mar-1995 : Dimensioniong of WSE changed to (N,ND) from (N,NCF).
!                        _V2 append to name.
!                       Call unchanged.
! Created 19-Jul-1993 : Based on EVALSE_QWVJ
!
	INTEGER ID            !Ion identifier
	INTEGER N_A		!Number of levles in ionizations state i
	INTEGER N_B		!Number of levles in ionizations state i+1
	INTEGER ION_EQ_IN_BA	!Eqn. # of target species [gs. of (i+2)th ]
	INTEGER NION		!Total number of IONS
	INTEGER ND		!Number of depth points
!
	REAL*8 HN_A(N_A,ND)		!    Pops. of ith ionzation stage
	REAL*8 HNST_A(N_A,ND)		!LTE   "    "  "      "       "
	REAL*8 HN_B(N_B,ND)		!    Pops. of (i+1)th ionization stage
	REAL*8 HNST_B(N_B,ND)		!LTE   "   "    "       "        "
!
! WSE_X is the quadrature weight for X-ray ionization (with 2e ejected) for
! ionization state i [final product is (i+1)].
!
	REAL*8 WSE_X(N_A,ND)
	REAL*8 JREC(ND)			! Int (2h/c2v^3+J)*EXP(-hv/kT)/v dv
	REAL*8 JPHOT(ND)		! Int J/v dv
!
! Local variables.
!
	INTEGER ION_EQ		!Ion eqation in SE(ID)%BA matrix.
	INTEGER ION_V                 !Location of ion variable in SE(ID)%BA.
	INTEGER I,J
	REAL*8 NETR
	REAL*8 SUM_SE
	REAL*8 SUM_VJ_R,SUM_VJ_P
	REAL*8 J_B_ION,B_ION
!
	ION_EQ=SE(ID)%XRAY_EQ
	ION_V=SE(ID)%LNK_TO_IV(ION_EQ_IN_BA)
	IF(SUM(WSE_X) .EQ. 0D0)RETURN
!
! Note that the product HNST_A(j, )*HNST_B(1, )/HN_B(1, ) is effectively the
! LTE population  of the state j with respect to the g.s. of the (I+2)th
! ionization stage.
!
	DO J=1,ND
	  SUM_SE=0.0D0
	  SUM_VJ_R=0.0D0
	  SUM_VJ_P=0.0D0
	  IF(JREC(J) .NE. 0.0D0)THEN
	    B_ION=HNST_B(1,J)/HN_B(1,J)		!1/b
	    J_B_ION=JREC(J)*B_ION
	  ELSE
	    B_ION=0.0D0
	    J_B_ION=0.0D0
	  END IF
	  DO I=1,N_A
	    NETR=WSE_X(I,J)*( HNST_A(I,J)*J_B_ION-HN_A(I,J)*JPHOT(J) )
	    SUM_SE  =SUM_SE+NETR
	    SUM_VJ_R=SUM_VJ_R + WSE_X(I,J)*B_ION*HNST_A(I,J)
	    SUM_VJ_P=SUM_VJ_P + WSE_X(I,J)*HN_A(I,J)
	    SE(ID)%STEQ(I,J) =SE(ID)%STEQ(I,J)  + NETR
	    SE(ID)%QFV_R(I,J)=SE(ID)%QFV_R(I,J) + WSE_X(I,J)*B_ION*HNST_A(I,J)
	    SE(ID)%QFV_P(I,J)=SE(ID)%QFV_P(I,J) + WSE_X(I,J)*HN_A(I,J)
	  END DO
!
! Include effects of X-rays on equation of target ion. For K shell 
! ionizations of species  with more than 3 electrons, as considered here,
! 2 electrons are given off. Thus ION_EQ should refer to the g.s. of the
! (i+2)th ionization stage.
!
	  SE(ID)%STEQ(ION_EQ,J) =SE(ID)%STEQ(ION_EQ,J) -SUM_SE
	  SE(ID)%QFV_R(ION_EQ,J)=SE(ID)%QFV_R(ION_EQ,J)-SUM_VJ_R
	  SE(ID)%QFV_P(ION_EQ,J)=SE(ID)%QFV_P(ION_EQ,J)-SUM_VJ_P
!
	END DO
!
	RETURN
	END
