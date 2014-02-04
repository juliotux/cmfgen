C
C This routine solves for the mean intensity as a function of depth using the
C Feautrier Technique for a spherical gray atmosphere. A diffusion approaximation
C is used for the lower boundary condition.
C This routine must be in a loop so that the f values are iterated to
C convergence.
C
	SUBROUTINE JGREY(TA,TB,TC,XM,DTAU,R,Z,P,RJ,NEWRJ,NEWRK,Q,F
	1  ,CHI,dCHIdr,AQW,AQW3,LUM,HBC,HBCNEW,NC,ND,NP,METHOD)
	IMPLICIT NONE
C
C Altered 24-May-1996 - Calls to DP_ZERO deleted. IONE inserted in call to
C                          THOMAS.
C Altered  7-Apr-1987 - Bug fixed (not passing CHI to NORDTAU in ray
C                          integration.
C Altered 13-Mar-1987 - Method option installed. Correct units installed for
C                          so that J has the units of B.
C Created 2-DEC-1986 - Based on JFEAU
C
	INTEGER NC,ND,NP,I,NI,LS
	REAL*8 TA(ND),TB(ND),TC(ND),XM(ND),R(ND),Z(ND),CHI(ND),dCHIdr(ND)
	REAL*8 RJ(ND),NEWRJ(ND),NEWRK(ND),DTAU(ND),Q(ND),F(ND)
	REAL*8 AQW(ND,NP),AQW3(ND,NP),P(NP),LUM,HBC,HBCNEW,E1,E2,E3
	REAL*8 DBC,DBB,T1,T2,IC
	CHARACTER*6 METHOD
	LOGICAL DIFF
C
	INTEGER, PARAMETER :: IONE=1
C
	DIFF=.TRUE.
	RJ(:)=0.0D0
	NEWRJ(:)=0.0D0
	NEWRK(:)=0.0D0
	HBCNEW=0.0D0
C
C Form the sphericity factor Q from F
C
	CALL QFROMF(F,Q,R,TA,TB,ND)	!TA work vector
C
C Form "SPHERICAL" optical depth scale.
C
	DO I=1,ND
	  TA(I)=Q(I)*CHI(I)
	END DO
	CALL DERIVCHI(TB,TA,R,ND,METHOD)
	CALL NORDTAU(DTAU,TA,R,R,TB,ND)
C
C Compute the solution vector. Note that the units need to be
C eventually included. Note T1=Lsun/16/(PI*PI)/10**2 (10**2 for 1/R**2).
C
	T1=3.826D+13*LUM/16.0D0/(3.141592654D0)**2.0
	T2=F(1)*Q(1)/HBC
	RJ(1)=T1/HBC/R(1)/R(1)
	DO I=2,ND
	  T2=T2+DTAU(I-1)
	  RJ(I)=T1/R(I)/R(I)/F(I)/Q(I)*T2
	END DO
C
C DBB =3L/16(piR)**2 and is used for the lower boundary diffusion approximation.
C
	DBB=3.0D0*T1/R(ND)/R(ND)
C
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
C
C ENTER LOOP FOR EACH IMPACT PARAMETER P
C
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*SQRT( (R(ND)-P(LS))*(R(ND)+P(LS)) )/R(ND)
	  END IF
C
C Compute Z for this imapct parameter
C
	  IF(NI .GT. 2)THEN
	    CALL ZALONGP(R,Z,P(LS),NI)
 	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdr,NI)
	  END IF
C
C Compute XM vector.
C
	  DO I=1,NI
	    TA(I)=RJ(I)		!Since ZETA=0, THETA=1.0
	  END DO
C
	  IF(NI .GT. 2)THEN
	    CALL XVECD(DTAU,TA,XM,DIFF,DBC,IC,LS,NC,ND,NI)
C
C Compute T ( a tridiagonal matrix) and store it as three vectors
C TA,TB and TC .
C
	    CALL TCOMPD(TA,TB,TC,DTAU,DIFF,LS,NC,ND,NI)
C
C Solve the tridiagonal system of equations.
C
	    CALL THOMAS(TA,TB,TC,XM,NI,IONE)
C
	ELSE IF(NI .EQ. 1)THEN
	  XM(1)=0.0D0
	ELSE IF(NI .EQ. 2)THEN
	  Z(1)=SQRT( (R(1)-P(LS))*(R(1)+P(LS)) )
	  DTAU(1)=0.5D0*Z(1)*(CHI(1)+CHI(2))		!Z(2)=0.0
	  E1=DEXP(-DTAU(1))
	  E2=1.0D0-(1.0D0-E1)/DTAU(1)
	  E3=(1.0D0-E1)/DTAU(1)-E1
	  IF(DTAU(1) .LT. 1.0D-03)THEN
	    E2=DTAU(1)*0.5D0+DTAU(1)*DTAU(1)/6.0D0
	    E3=DTAU(1)*0.5D0-DTAU(1)*DTAU(1)/3.0D0
	  END IF
C
	  XM(2)=TA(2)*E2+TA(1)*E3
          XM(1)=0.5D0*(XM(2)*E1+TA(1)*E2+TA(2)*E3)
	END IF
C
C Update the ANGLE integrations.
C
	  DO I=1,NI
	    NEWRJ(I)=NEWRJ(I)+AQW(I,LS)*XM(I)
	    NEWRK(I)=NEWRK(I)+AQW3(I,LS)*XM(I)
	  END DO
C
	  HBCNEW=HBCNEW+AQW(1,LS)*(XM(1))*Z(1)/R(1)
C
2000	CONTINUE
C
C Compute the new Feautrier factors. These are stored in NEWRK so as not
C to destroy the old factors needed for the linearization.
C
	DO I=1,ND
	  NEWRK(I)=NEWRK(I)/NEWRJ(I)
	END DO
C
C Compute the factor for the outer boundary condition.
C
	HBCNEW=HBCNEW/NEWRJ(1)
C
	RETURN
	END
