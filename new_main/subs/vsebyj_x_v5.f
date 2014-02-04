!
! Subroutine to increment the variation matrix due to terms which
! depend directly on the intensity J. The Radiative equilibrium equation
! is not altered. This routine is for X-ray ionization only, where
! 2 electrons are ejected.
!
! Routine also increments the ionization equilibrium equations.
!
	SUBROUTINE VSEBYJ_X_V5(ID,WSE_X,
	1             HN_A,HNST_A,dlnHNST_AdlnT,N_A,
	1             HN_B,HNST_B,EDGE_B,N_B,DI,ION_EQ_IN_BA,
	1             ED,T,EMHNUKT,NU,ML,RJ,
	1             NCF,ND,NION,DST,DEND)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 12-Apr-2001 : Changed to use STEQ_DATA_MOD.
!                       Changed to V5.
! Altered 27-OCt-1995 : Changed to be compatible with super levels.
!                        Call chnaged --- Now _V3.
! Altered 06-Mar-1995 : Dimensions of WSE_X changed to (N,ND) from (N,NCF).
!                       _V2 append to call.
! Testing 22-Jul-1994 : Minor mods.
! Created 19-Jul-1993 : Based on VSEBYJ_COM and EVALSE_QWVJ.
!
	INTEGER ID		!Ion identification
	INTEGER N_A		!Number of levels in ionizations state i
	INTEGER N_B		!Number of levels in ionizations state i+1
	INTEGER NCF		!Number of freuencies in NU 
	INTEGER ML		!Indicates current freqency in vector NU
	INTEGER ND		!Number of depth points
        INTEGER NION		!Total Number of ions in model.
        INTEGER ION_EQ_IN_BA  !
!
	INTEGER DST,DEND
!
	REAL*8 WSE_X(N_A,ND)		!Quadrature weights (incl. cross. sec.)
!
! _A refers to quantities associated with the atom WITH super levels.
!
	REAL*8 HN_A(N_A,ND)		!Pops. of ith ionzation stage
	REAL*8 HNST_A(N_A,ND)		!LTE   "    "  "      "       "
	REAL*8 dlnHNST_AdlnT(N_A,ND)	
!
! _B refers to quantities assiciated with the FULL atom of the next 
!  ionization stage.
!
	REAL*8 HN_B(N_B,ND)		!Pops. of (i+1)th ionization stage
	REAL*8 HNST_B(N_B,ND)		!LTE   "   "    "       "        "
	REAL*8 EDGE_B(N_B)
!
	REAL*8 DI(ND)			!Ion density for B levels.
	REAL*8 ED(ND)			!Electron density
	REAL*8 T(ND)			!Temperature in 10^4K.
	REAL*8 EMHNUKT(ND)		!EXP(-hv/kT)
	REAL*8 RJ(ND)			!Mean intensity
	REAL*8 NU(NCF)			!frequency (10^15 Hz)
!
! Constants for opacity etc.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables
!
	INTEGER NIV
	INTEGER ION_EQ		!Ion eqation in SE(ID)%BA matrix.
	INTEGER ION_V                 !Location of ion variable in SE(ID)%BA.
	INTEGER I,J
	REAL*8 T1,T3,BSTIM
	REAL*8 WSE_BY_RJ,DI_FAC,ED_FAC,T_FAC
!
	NIV=SE(ID)%N_IV
	ION_EQ=SE(ID)%XRAY_EQ
	ION_V=SE(ID)%LNK_TO_IV(ION_EQ_IN_BA)
!
	T1=TWOHCSQ*( NU(ML)**3 )
	DO I=DST,DEND			!Which depth point.
	  BSTIM=(T1+RJ(I))*EMHNUKT(I)*HNST_B(1,I)/HN_B(1,I)
	  DO J=1,N_A			!Which equation (for S.E. only)
	    IF(WSE_X(J,I) .NE. 0)THEN
	      WSE_BY_RJ=WSE_X(J,I)*RJ(I)
	      SE(ID)%BA_PAR(J,J,I)=SE(ID)%BA_PAR(J,J,I)-WSE_BY_RJ
!
! NB: In the following, the factor of 2.0 for ED_FAC arrises because 
! the rate is prop. to 
!
!        HNST_A(J,I)*HNST_B(1,I) .
!
! The variation with respect to HN_B(1,I) is intrinsically zero, since 
! HNST_A(J,I)/HN_B(1,I) is independent of HN_B(1,I)
!
	      T3=HNST_A(J,I)*WSE_X(J,I)*BSTIM
	      DI_FAC=T3/DI(I)
	      ED_FAC=2.0D0*T3/ED(I)
	      T_FAC=T3*( HDKT*(NU(ML)-EDGE_B(1)-1.5D0)/T(I) +
	1             dlnHNST_AdlnT(J,I) )/T(I)
!
	      SE(ID)%BA_PAR(J,ION_V,I)=SE(ID)%BA_PAR(J,ION_V,I) +DI_FAC
	      SE(ID)%BA_PAR(J,NIV-1,I)=SE(ID)%BA_PAR(J,NIV-1,I) +ED_FAC
	      SE(ID)%BA_PAR(J,NIV,I)  =SE(ID)%BA_PAR(J,NIV,I)   +T_FAC
!
! Include ionizations/recombinations explicitly in the rate equation
! of the target ion (eg He++(gs) for He+ ion/recoms). 
!
	      SE(ID)%BA_PAR(ION_EQ,J,I)    =SE(ID)%BA_PAR(ION_EQ,J,I)    +WSE_BY_RJ
	      SE(ID)%BA_PAR(ION_EQ,ION_V,I)=SE(ID)%BA_PAR(ION_EQ,ION_V,I)-DI_FAC
	      SE(ID)%BA_PAR(ION_EQ,NIV-1,I)=SE(ID)%BA_PAR(ION_EQ,NIV-1,I)-ED_FAC
	      SE(ID)%BA_PAR(ION_EQ,NIV,I)  =SE(ID)%BA_PAR(ION_EQ,NIV,I)  -T_FAC
!
	    END IF		!WSE(J,ML) .NE. 0
	  END DO
	END DO
!
	RETURN
	END
