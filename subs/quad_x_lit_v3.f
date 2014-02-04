C
C Subroutine to compute the quadrature weights for the statistical
C equilibrium equations due to K shell ionization of Lithium like
C ions. This routine assumes, for example, that ionization if CIV(1s^2nl)
C results in the CV ground state. Should be revised when CV is available,
C since resulting state is CV(1snl).
C
C These cross section get ADDED directly to the normal (ie valence electron)
C quadrature weights for the Lithium like ion.
C
C The factor of DEX(-10) in TP1 is due to the definition of PRGEN which 
C is DEX(10) times the photoionization cross section so that CHI*R 
C is constant.
C
	SUBROUTINE QUAD_X_LIT_V3(ZCORE,NUM_ELEC,WSE,NU_CONT,
	1               N,ND,COMPUTE_NEW_CROSS)
	IMPLICIT NONE
C
C Altered 17-Sep-1997 : Altered so that a constant continuum cross-section
C                         across a band can be handled. WSE now must be
C                         be effectively multiplied by FQW/NU_CONT when it 
C                         is used.
C Altered 26-May-1996 : FOUR_PI_ON_H constant fixed.
C                       ERROR_LU installed.
C Altered 06-Mar-1995 : Dimensionioning of WSE changed. New WSE for each
C                         frequency.
C Created 20-Jul-1993 : Based on QUADGEN
C
	INTEGER N,ND
	REAL*8 ZCORE,NUM_ELEC
	REAL*8 WSE(N,ND)
	REAL*8 NU_CONT
	LOGICAL COMPUTE_NEW_CROSS
C
	INTEGER ERROR_LU
	REAL*8 XCROSS
	EXTERNAL XCROSS,ERROR_LU
C
	INTEGER I,K,LUER
	REAL*8 FOUR_PI_ON_H
	REAL*8 WEIGHT
C
	IF(NUM_ELEC .NE. 3)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in QUAD_X_LIT --- invalid number of electrons'
	  WRITE(lUER,*)'This routine only valid for K shell ionization of',
	1           ' lithium like ions'
	  STOP
	END IF
C
	IF(COMPUTE_NEW_CROSS)THEN
C
C Divide by 10^15 as FQW now in Hz.
C See comments in QUAD_X_GEN_V3.
C
	  FOUR_PI_ON_H=1.8965D+02		!4*PI/H*DEX(-10)/1.0D+15
C
C Now evaluate the quadrature weights. 
C
C The cross section are the same for all levels. 
C
	  WEIGHT=FOUR_PI_ON_H*XCROSS(NU_CONT,ZCORE,NUM_ELEC)
C
C Now set quadrature weights for other levels. We are assuming the
C cross-section, and ionization energy, are independent of the outer
C electron configuration.
C
	  DO K=1,ND
	    DO I=1,N
	      WSE(I,K)=WSE(I,K)+WEIGHT
	    END DO
	  END DO
	END IF
C
	RETURN
	END
