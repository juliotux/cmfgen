C
C Compute the flux using NORDULUND weights.
C
	SUBROUTINE NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,
	1    SOURCE,CHI,dCHIdr,HAQW,SOB,
	1    S1,THICK,DIF,DBB,IC,NC,ND,NP,METHOD)
	IMPLICIT NONE
C
C Altered 28-Oct-1996 - COS converted back to ACOS in TOR expression.
C Altered 24-May-1996 - Call to DP_ZERO removed.
C                       DEXP etc replaced by GENRIC calls.
C                       IONE used in call to THOMAS
C
C Altered 07-Jan-1991 - Computation of TA,TB,TC split into two loops
C                       so CRAY compiler recognized that its vectorizable.
C Altered 22-Feb-1987 - Call changed (dCHIdr vector included)
C                       Method is now also a passed parameter.
C Altered 10-Feb-1987 - NI=2 is now handled correctly. Dont have to
C                       consider the case NI=1 since the flux is zero.
C
	INTEGER NC,ND,NP
	REAL*8 TA(ND),TB(ND),TC(ND),XM(ND),DTAU(ND),dCHIdr(ND)
	REAL*8 R(ND),Z(ND),P(NP),SOURCE(ND),CHI(ND),SOB(ND)
	REAL*8 HAQW(ND,NP),S1,DBB,IC
	LOGICAL THICK,DIF
	CHARACTER*6 METHOD
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER  LS,I,NI
	REAL*8 IBOUND,TOR,DBC
	REAL*8 E1,E2,E3,PP
C
C Compute dCHIdr
C
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
C
C Compute flux distribution and luminosity (in L(sun)) of star.
C
	SOB(:)=0.0D0
	DO 8000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	  END IF
	  IF(NI .LT. 2)GOTO 8000
C
	IF(THICK)THEN
	  IF(P(LS) .GT. 0.0D0)THEN
	    TOR=CHI(1)*R(1)*R(1)*(1.570796D0-ACOS(P(LS)/R(1)))/P(LS)
	  ELSE
	    TOR=CHI(1)*R(1)
	  END IF
	  IBOUND=S1*(1.0D0-EXP(-TOR))
	ELSE
	  IBOUND=0.0D0
	END IF
C
C
C Compute Z and the optical depth scale DTAU for this imapct parameter.
C
	  IF(NI .GT. 1)THEN
C	    CALL ZALONGP(R,Z,P(LS),NI)
C	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdr,NI)
	    PP=P(LS)*P(LS)
	    DO I=1,NI
	      Z(I)=SQRT(R(I)*R(I)-PP)
	    END DO
	    DO I=1,NI-1
	      DTAU(I)=0.5D0*(Z(I)-Z(I+1))*(CHI(I)+CHI(I+1)+(Z(I)-Z(I+1))
	1     *(dCHIdR(I+1)*Z(I+1)/R(I+1)-dCHIdR(I)*Z(I)/R(I))/6.0D0)
	    END DO
	  END IF
C
C Compute XM. Compute T ( a tridiagonal matrix) and store it as three vectors
C TA,TB and TC . This code is a combined version of XVECD and TCOMPD.
C
C	    CALL XVECD(DTAU,SOURCE,XM,DIF,DBC,IC,LS,NC,ND,NI)
C	    IF(THICK)XM(1)=-IBOUND
C	    CALL TCOMPD(TA,TB,TC,DTAU,DIF,LS,NC,ND,NI)
C
	  IF(NI .GT. 2)THEN
	    XM(1)=0.0D0
	    IF(THICK)XM(1)=-IBOUND
	    TA(1)=0.0D0
	    TC(1)=1.0D0/DTAU(1)
	    TB(1)=-1.0D0-TC(1)
	    DO I=2,NI-1
	      TC(I)=1.0D0/DTAU(I)
	    END DO
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TB(I)=-0.5D0*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5D0
	    END DO
C
	    IF(LS .LE. NC .AND. DIF)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=TC(NI-1)
	      XM(NI)=DBC
	    ELSE IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-TA(NI)+DTAU(NI-1)/2.0D0
	      XM(NI)=0.5D0*DTAU(NI-1)*SOURCE(NI)
	    ELSE
	      TA(NI)=-TC(NI-1)
	      TB(NI)=1.0D0-TA(NI)
	      XM(NI)=IC
	    END IF
	    TC(NI)=0.0D0
C
C Solve the tridiagonal system of equations.
C
	  CALL THOMAS(TA,TB,TC,XM,NI,IONE)
C
C Compute the total emitted flux through each depth point using
C linear differencing. The following is based on FLUXGRID.
C
C	  CALL FLUXGRID(XM,DTAU,TB,TC,NI)
C
C Compute the midpoint flux
C
	    DO I=1,NI-1
	      TC(I)=(XM(I+1)-XM(I))/DTAU(I)
	    END DO
C
C Interpolate the flux onto the radius grid. Handle the end points by
C extrapolation.
C
	    DO I=2,NI-1
	      TB(I)=(DTAU(I-1)*TC(I)+DTAU(I)*TC(I-1))
	1           /(DTAU(I-1)+DTAU(I))
	    END DO
	    TB(1)=(TC(1)*(2.0D0*DTAU(1)+DTAU(2))-DTAU(1)*TC(2))
	1              /(DTAU(1)+DTAU(2))
	    TB(NI)=(TC(NI-1)*(2.0D0*DTAU(NI-1)+DTAU(NI-2))-
	1               DTAU(NI-1)*TC(NI-2))
	1              /(DTAU(NI-1)+DTAU(NI-2))
C
C End of flux computation on the radius mesh.
C
	  DO I=3,NI
	    SOB(I)=SOB(I)+TB(I)*HAQW(I,LS)
	  END DO
	  SOB(1)=SOB(1)+TB(1)*HAQW(1,LS)
	  SOB(2)=SOB(2)+XM(1)*HAQW(1,LS)
C
	ELSE IF(NI .EQ. 2)THEN
	  E1=EXP(-DTAU(1))
	  E2=1.0D0-(1.0D0-E1)/DTAU(1)
	  E3=(1.0D0-E1)/DTAU(1)-E1
	  IF(DTAU(1) .LT. 1.0D-03)THEN
	    E2=DTAU(1)*0.5D0+DTAU(1)*DTAU(1)/6.0D0
	    E3=DTAU(1)*0.5D0-DTAU(1)*DTAU(1)/3.0D0
	  END IF
C
	  XM(2)=IBOUND*E1+SOURCE(2)*E2+SOURCE(1)*E3
          XM(1)=0.5D0*(IBOUND+XM(2)*E1+SOURCE(1)*E2+SOURCE(2)*E3)
	  TB(1)=XM(1)-0.5D0*IBOUND
	  TB(2)=0.0D0
C
	  SOB(1)=SOB(1)+TB(1)*HAQW(1,LS)
	  SOB(2)=SOB(2)+XM(1)*HAQW(1,LS)
C
	END IF
8000	CONTINUE
C
C Update fluxes by R**2
C
	  SOB(1)=SOB(1)*R(1)*R(1)
	  SOB(2)=SOB(2)*R(1)*R(1)		!U rather than V at d=1.
	  DO I=3,ND
	    SOB(I)=SOB(I)*R(I)*R(I)
	  END DO
C
	RETURN
	END
