!
! This routine solves for the mean intensity as a function of depth using the
! Feautrier Technique. A Schuster or diffusion approaximation is used for the
! lower boundary condition. This routine must be in a loop so that the f
! values are iterated to convergence.
!
	SUBROUTINE JFEAU_IBC_V2(TA,TB,TC,DTAU,R,RJ,Q,F,
	1            ZETA,THETA,CHI,DBB,IC,HBC_J,HBC_S,
	1            INBC,THK,INNER_BND_METH,ND,METHOD)
	IMPLICIT NONE
!
! Altered 15-Jan-2015 - Added check for -ve intensities. 
! Altered 07-Jun-2010 - Changed to V2. Changed DIF to INNER_BND_METH.
! Altered 24-May-1996 - Call to DP_ZERO deleted; IONE isnatlled in THOMAS call.
! Created 12-JUN-1991 - Based on JFEAUNEW. HBC replaced by HBC_J and HBC_S.
!                       [ IBC - Improved boundary confition. ]
!
! Altered 29-May-1989 - Q now computed in routine.
! Altered 24-Feb-1987 - Q nolonger computed in routine. Method made a string.
! Altered 10-Feb-1987 - The accuracy of the optical depth scale has been
!                       improved by correcting the integral by the first
!                       derivatives. (Was previously done for J but now also
!                       done for U and hence f computation).
!
! Altered 31-Oct-1986 - Schuster boundary condition installed at inner boundary.
!                       Two new variables INBC and INBCNEW now in call. These
!                       are used at the inner boundary (NB IC is intensity
!                       incident on the inner boundary).
!                       Calls to XVECFEAU and TFEAU altered.
! Altered 4-Mar-1986 -  New f Feautrier factors returned in the NEWRK
!                       array. (similarly HBC in HBCNEW).
! Altered 28-FEB-1986 - AQW3 Installed. Integrating bu AQW*(mu)**2.0 d(mu)
!                       gave invalid f values (i.e not 0.33333) at inner
!                       boundary since only trapazoidal weights.
! Created 17-FEB-1986
!
	INTEGER ND,I
	REAL*8 TA(ND),TB(ND),TC(ND),R(ND),ZETA(ND)
	REAL*8 RJ(ND),DTAU(ND),Q(ND),F(ND)
	REAL*8 THETA(ND),CHI(ND)
	REAL*8 DBB,HBC_J,HBC_S,INBC,IC,T1
	CHARACTER(LEN=6) METHOD
	CHARACTER(LEN=*) INNER_BND_METH
	LOGICAL THK
!
	INTEGER, PARAMETER :: IONE=1
!
	RJ(:)=0.0D0
!
! Compute the Q factors from F.
!
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA, TB are work vectors.
!
! Form "SPHERICAL" optical depth scale.
!
	DO I=1,ND
	  TA(I)=Q(I)*CHI(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
!
! Compute T ( a tridiagonal matrix) and store it as three vectors
! TA,TB and TC .
!
	T1=HBC_J-HBC_S*THETA(1)
	TA(1)=0.0D0
	TC(1)=F(2)*Q(2)*R(2)*R(2)/DTAU(1)
	TB(1)=-R(1)*R(1)*(F(1)*Q(1)/DTAU(1)+T1)
!
	DO I=2,ND-1
	  T1=0.5D0*(DTAU(I-1)+DTAU(I))
	  TA(I)=-F(I-1)*Q(I-1)*R(I-1)*R(I-1)/DTAU(I-1)/T1
	  TC(I)=-F(I+1)*Q(I+1)*R(I+1)*R(I+1)/DTAU(I)/T1
	  TB(I)=R(I)*R(I)*((1.0D0-THETA(I))/Q(I)
	1        +F(I)*Q(I)*(1.0D0/DTAU(I)+1.0D0/DTAU(I-1))/T1)
	END DO
!
! Note that Q(ND)=1.0d0 by definition.
!
        IF(INNER_BND_METH .EQ. 'DIFUSION')THEN
          TA(ND)=-F(ND-1)*Q(ND-1)*R(ND-1)*R(ND-1)/DTAU(ND-1)
          TB(ND)=F(ND)*R(ND)*R(ND)/DTAU(ND-1)
          TC(ND)=0.0D0
        ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
          TA(ND)=-F(ND-1)*Q(ND-1)*R(ND-1)*R(ND-1)/DTAU(ND-1)
          TB(ND)=R(ND)*R(ND)*(F(ND)/DTAU(ND-1)+INBC)
          TC(ND)=0.0D0
	ELSE
          TA(ND)=-F(ND-1)*Q(ND-1)*R(ND-1)*R(ND-1)/DTAU(ND-1)
          TB(ND)=F(ND)*R(ND)*R(ND)/DTAU(ND-1)
          TC(ND)=0.0D0
        END IF
!
! Form the SOURCE vector (section replaces call to XVECFEAU).
!
	RJ(1)=-HBC_S*R(1)*R(1)*ZETA(1)
	DO I=2,ND-1
	  RJ(I)=R(I)*R(I)*ZETA(I)/Q(I)
	END DO
!
! Note well - DBB =dB/dR (and Q(ND)=1.0 by definition)
! Final B.C. is hollow core / zero flux.
!
	IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	  RJ(ND)=R(ND)*R(ND)*DBB/CHI(ND)/3.0D0
	ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	  RJ(ND)=R(ND)*R(ND)*IC*(0.25D0+0.5D0*INBC)
	ELSE
	  RJ(ND)=0.0D0
	END IF
!
! Find the solution
!
	CALL THOMAS(TA,TB,TC,RJ,ND,IONE)
!
	IF(MINVAL(RJ) .LE. 0)THEN
	  WRITE(6,*)'Error in FQCOMP_IBC -- -ve mean intensities'
	  DO I=1,ND
	     WRITE(6,'(I4,6ES14.6)')I,RJ(I),CHI(I),DTAU(I),ZETA(I),THETA(I),F(I),Q(I)
	  END DO
	  STOP
	END IF
!
	RETURN
	END
