!
! Subroutine to compute the quadrature weights for the statistical
! equilibrium equations for K shell ionization (by X-rays) of ions 
! with MORE than 3 electrons.
!
! The factor of DEX(-10) in TP1 is due to the definition of PRGEN which is 
! DEX(10) times the photoionization cross section so that CHI*R is constant.
!
	SUBROUTINE QUAD_X_GEN_V5(ZCORE,NUM_ELEC,WSE,WCR,NU_CONT,
	1             HNST_S,N_S,
	1             HNST_F_ON_S,EDGE_F,F_TO_S,N_F,
	1             EDGE_B,N_B,ND)
	IMPLICIT NONE
!
! Altered 05-Apr-2011 - Changed to V5.
!                       HNST_F_ON_S (rather than HNST_F) is passed in call.
!                       HNST_F/HNST_S replaced by HNST_F_ON_S - done to faciliate
!                         modifications allowing lower temperaturs.
!                       Most of editing done early 2011 
! Altered 17-Sep-1997 : Altered so that a constant continuum cross-section
!                         across a band can be handled. WSE now must be
!                         be effectively multiplied by FQW/NU_CONT when it 
!                         is used.
! Altered 30-Jan-1995 : Buf fix. WSE was a factor of 10^{15} too large.
! Altered 27-Oct-1995 : Adapted to allow for the presence of super levels.
!                        Several new variables inserted in call
!                        NOW _V3.
! Altered 06-Mar-1995 : Dimensioning of WSE changed. New WSE for each
!                         frequency.
! Created 20-Jul-1993 : Based on QUADGEN
!
	INTEGER N_S,N_B,N_F,ND
	REAL*8 ZCORE,NUM_ELEC
!
! N_S refers to the model atom with SUPER levels.
!
	REAL*8 WSE(N_S,ND)
	REAL*8 WCR(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
!
! _F refers to populations in the full atom.
!
	REAL*8 HNST_F(N_F,ND)
	REAL*8 HNST_F_ON_S(N_F,ND)
	REAL*8 EDGE_F(N_F)
	INTEGER F_TO_S(N_F)
!
! _B Refers to FULL atom in nest ionization stage.
!
	REAL*8 EDGE_B(N_B)
	REAL*8 NU_CONT
!
	EXTERNAL XCROSS_V2
	REAL*8 XCROSS_V2
!
	INTEGER I_F,I_S,K
	REAL*8 WEIGHT
	REAL*8 FOUR_PI_ON_H
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! NB: WSE represents the statistical weight
!
	WSE(:,:)=0.0D0		!N_S,ND
	WCR(:,:)=0.0D0		!N_S,ND
!
! NB: Because SUM[HNST_F] = HNST_S and the cross section is independent of
! the level under consideration, WSE_S will just be = WEIGHT. [Because
! the cross section switches on a little earlier for the higher levels,
! there is a weak dependance.] Therefor we do not need a dWSE_S_dT.
!
! NB: Constant FOUR_PI_ON_H differs from original constant in QUADGEN
!       and QUAD_X_GEN by factor of 10^{-15} due to the fact the FQW
!       is now in units of Hz, not 10^15 Hz. The factor of 10^{-15}
!       arises from the 1/nu term. [Change made when we removed the
!       frequency index from the WSE array].
!
	FOUR_PI_ON_H=1.8965D+02		!4*PI/H*DEX(-10)/1.0D+15
	WEIGHT=FOUR_PI_ON_H*XCROSS_V2(NU_CONT,ZCORE,NUM_ELEC,
	1                                   IZERO,IZERO,L_FALSE,L_FALSE)
	IF(WEIGHT .EQ. 0)RETURN
!
	DO K=1,ND
	  DO I_F=1,N_F
	    I_S=F_TO_S(I_F)
	    WSE(I_S,K)=WSE(I_S,K) + WEIGHT*HNST_F_ON_S(I_F,K)
	    WCR(I_S,K)=WCR(I_S,K) - (EDGE_F(I_F)+EDGE_B(1))*WEIGHT*HNST_F_ON_S(I_F,K)
	  END DO
	END DO
!
	RETURN
	END
