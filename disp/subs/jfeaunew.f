C
C This routine solves for the mean intensity as a function of depth using the
C Feautrier Technique. A Schuster or diffusion approaximation is used for the
C lower boundary condition.
C This routine must be in a loop so that the f values are iterated to
C convergence.
C
C Created 17-FEB-1986
C Altered 28-FEB-1986 - AQW3 Installed. Integrating bu AQW*(mu)**2.0 d(mu)
C                       gave invalid f values (i.e not 0.33333) at inner
C                       boundary since only trapazoidal weights.
C Altered 4-Mar-1986 -  New f Feautrier factors returned in the NEWRK
C                       array. (similarly HBC in HBCNEW).
C
C Altered 31-Oct-1986 - Schuster boundary condition installed at inner boundary.
C                       Two new variables INBC and INBCNEW now in call. These
C                       are used at the inner boundary (NB IC is intensity
C                       incident on the inner boundary).
C                       Calls to XVECFEAU and TFEAU altered.
C
C Altered 10-Feb-1987 - The accuracy of the optical depth scale has been
C                       improved by correcting the integral by the first
C                       derivatives. (Was previously done for J but now also
C                       done for U and hence f computation).
C
C Altered 24-Feb-1987 - Q nolonger computed in routine. Method made a string.
C
C Altered 29-May-1989 - Q now computed in routine.
C
	SUBROUTINE JFEAUNEW(TA,TB,TC,DTAU,R,RJ,Q,F,
	1                   ZETA,THETA,CHI,DBB,IC,HBC,
	1                   INBC,THK,DIFF,ND,METHOD)
	IMPLICIT NONE
C
	INTEGER ND,I
	REAL*8 TA(ND),TB(ND),TC(ND),R(ND),ZETA(ND)
	REAL*8 RJ(ND),DTAU(ND),Q(ND),F(ND)
	REAL*8 THETA(ND),CHI(ND)
	REAL*8 DBB,HBC,INBC,IC
	CHARACTER*6 METHOD
	LOGICAL DIFF,THK
C
	CALL DP_ZERO(RJ,ND)
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
	CALL TFEAU(TA,TB,TC,R,Q,F,THETA,DTAU,HBC,INBC,DIFF,ND)
C
C Form the SOURCE vector
C
	CALL XVECFEAU(RJ,R,Q,ZETA,DIFF,DBB,INBC,IC,CHI(ND),ND)
C
C Find the solution
C
	CALL THOMAS(TA,TB,TC,RJ,ND,1)
C
	RETURN
	END
