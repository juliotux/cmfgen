C
C This routine solves for the mean intensity as a function of depth using the
C Feautrier Technique. A Schuster or diffusion approaximation is used for the
C lower boundary condition. This routine must be in a loop so that the f
C values are iterated to convergence.
C
	SUBROUTINE JFEAU_IBC(TA,TB,TC,DTAU,R,RJ,Q,F,
	1                     ZETA,THETA,CHI,DBB,IC,HBC_J,HBC_S,
	1                     INBC,THK,DIFF,ND,METHOD)
	IMPLICIT NONE
C
C Altered 24-May-1996 - Call to DP_ZERO deleted; IONE isnatlled in THOMAS call.
C Created 12-JUN-1991 - Based on JFEAUNEW. HBC replaced by HBC_J and HBC_S.
C                       [ IBC - Improved boundary confition. ]
C
C Altered 29-May-1989 - Q now computed in routine.
C Altered 24-Feb-1987 - Q nolonger computed in routine. Method made a string.
C Altered 10-Feb-1987 - The accuracy of the optical depth scale has been
C                       improved by correcting the integral by the first
C                       derivatives. (Was previously done for J but now also
C                       done for U and hence f computation).
C
C Altered 31-Oct-1986 - Schuster boundary condition installed at inner boundary.
C                       Two new variables INBC and INBCNEW now in call. These
C                       are used at the inner boundary (NB IC is intensity
C                       incident on the inner boundary).
C                       Calls to XVECFEAU and TFEAU altered.
C Altered 4-Mar-1986 -  New f Feautrier factors returned in the NEWRK
C                       array. (similarly HBC in HBCNEW).
C Altered 28-FEB-1986 - AQW3 Installed. Integrating bu AQW*(mu)**2.0 d(mu)
C                       gave invalid f values (i.e not 0.33333) at inner
C                       boundary since only trapazoidal weights.
C Created 17-FEB-1986
C
	INTEGER ND,I
	REAL*8 TA(ND),TB(ND),TC(ND),R(ND),ZETA(ND)
	REAL*8 RJ(ND),DTAU(ND),Q(ND),F(ND)
	REAL*8 THETA(ND),CHI(ND)
	REAL*8 DBB,HBC_J,HBC_S,INBC,IC,T1
	CHARACTER*6 METHOD
	LOGICAL DIFF,THK
C
	INTEGER, PARAMETER :: IONE=1
C
	RJ(:)=0.0D0
C
C Compute the Q factors from F.
C
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA, TB are work vectors.
C
C Form "SPHERICAL" optical depth scale.
C
	DO I=1,ND
	  TA(I)=Q(I)*CHI(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
C
C COMPUTE T ( A TRIDIAGONAL MATRIX) AND STORE IT AS THREE VECTORS
C TA,TB AND TC .
C
	T1=HBC_J-HBC_S*THETA(1)
	CALL TFEAU(TA,TB,TC,R,Q,F,THETA,DTAU,T1,INBC,DIFF,ND)
C
C Form the SOURCE vector (section replaces call to XVECFEAU).
C
	RJ(1)=-HBC_S*R(1)*R(1)*ZETA(1)
	DO I=2,ND-1
	  RJ(I)=R(I)*R(I)*ZETA(I)/Q(I)
	END DO
C
C Note well - DBB =dB/dR (and Q(ND)=1.0 by definition)
C
	IF(DIFF)THEN
	  RJ(ND)=R(ND)*R(ND)*DBB/CHI(ND)/3.0D0
	ELSE
	  RJ(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*INBC)
	END IF
C
C Find the solution
C
	CALL THOMAS(TA,TB,TC,RJ,ND,IONE)
C
	RETURN
	END
