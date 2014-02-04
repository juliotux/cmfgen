C
C Subroutine to compute the quadrature weights for the statistical
C equilibrium equations for K shell ionization (by X-rays) of ions 
C with MORE than 3 electrons.
C
C The factor of DEX(-10) in TP1 is due to the definition of PRGEN which is 
C DEX(10) times the photoionization cross section so that CHI*R is constant.
C
	SUBROUTINE QUAD_X_GEN_V4(ZCORE,NUM_ELEC,WSE,WCR,NU_CONT,
	1             HNST_S,N_S,
	1             HNST_F,EDGE_F,F_TO_S,N_F,
	1             EDGE_B,N_B,ND)
	IMPLICIT NONE
C
C
C Altered 17-Sep-1997 : Altered so that a constant continuum cross-section
C                         across a band can be handled. WSE now must be
C                         be effectively multiplied by FQW/NU_CONT when it 
C                         is used.
C Altered 30-Jan-1995 : Buf fix. WSE was a factor of 10^{15} too large.
C Altered 27-Oct-1995 : Adapted to allow for the presence of super levels.
C                        Several new variables inserted in call
C                        NOW _V3.
C Altered 06-Mar-1995 : Dimensioning of WSE changed. New WSE for each
C                         frequency.
C Created 20-Jul-1993 : Based on QUADGEN
C
	INTEGER N_S,N_B,N_F,ND
	REAL*8 ZCORE,NUM_ELEC
C
C N_S refers to the model atom with SUPER levels.
C
	REAL*8 WSE(N_S,ND)
	REAL*8 WCR(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
C
C _F refers to populations in the full atom.
C
	REAL*8 HNST_F(N_F,ND)
	REAL*8 EDGE_F(N_F)
	INTEGER F_TO_S(N_F)
C
C _B Refers to FULL atom in nest ionization stage.
C
	REAL*8 EDGE_B(N_B)
	REAL*8 NU_CONT
C
	EXTERNAL XCROSS_V2
	REAL*8 XCROSS_V2
C
	INTEGER I_F,I_S,K
	REAL*8 WEIGHT
	REAL*8 FOUR_PI_ON_H
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
C
C NB: WSE represents the statistical weight
C
	WSE(:,:)=0.0D0		!N_S,ND
	WCR(:,:)=0.0D0		!N_S,ND
C
C NB: Because SUM[HNST_F] = HNST_S and the cross section is independent of
C the level under consideration, WSE_S will just be = WEIGHT. [Because
C the cross section switches on a little earlier for the higher levels,
C there is a weak dependance.] Therefor we do not need a dWSE_S_dT.
C
C NB: Constant FOUR_PI_ON_H differs from original constant in QUADGEN
C       and QUAD_X_GEN by factor of 10^{-15} due to the fact the FQW
C       is now in units of Hz, not 10^15 Hz. The factor of 10^{-15}
C       arises from the 1/nu term. [Change made when we removed the
C       frequency index from the WSE array].
C
	FOUR_PI_ON_H=1.8965D+02		!4*PI/H*DEX(-10)/1.0D+15
	WEIGHT=FOUR_PI_ON_H*XCROSS_V2(NU_CONT,ZCORE,NUM_ELEC,
	1                                   IZERO,IZERO,L_FALSE,L_FALSE)
	IF(WEIGHT .EQ. 0)RETURN
C
	DO K=1,ND
	  DO I_F=1,N_F
	    I_S=F_TO_S(I_F)
	    WSE(I_S,K)=WSE(I_S,K) + WEIGHT*HNST_F(I_F,K)/HNST_S(I_S,K)
	    WCR(I_S,K)=WCR(I_S,K) - (EDGE_F(I_F)+EDGE_B(1))*
	1                          WEIGHT*HNST_F(I_F,K)/HNST_S(I_S,K)
	  END DO
	END DO
C
	RETURN
	END
