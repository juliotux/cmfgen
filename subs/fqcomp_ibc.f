C
C This routine solves for the mean intensity as a function of depth using the
C Feautrier Technique. A Schuster or diffusion approaximation is used for the
C lower boundary condition. This routine may need to be in a loop so that
C the f values are iterated to convergence.
C
	SUBROUTINE FQCOMP_IBC(TA,TB,TC,XM,DTAU,R,Z,P,NEWRJ,NEWRK,
	1            SOURCE,CHI,dCHIdr,AQW,AQW3,DBB,HBC_J,HBC_S,
	1            INBCNEW,IC,THK,DIFF,NC,ND,NP,METHOD)
	IMPLICIT NONE
C
C Altered 26-May-2005 --- R*R-P*P replaved by (R-P)*(R+P)
C Altered 09-Dec-2001 --- PP replaced by P(LS)*P(LS)
C Altered 28-Oct-1997 --- Bug fix: COS corrected back to ACOS in TOR expression.
C Altered 24-May-1996 --- Call to DP_ZERO eliminated.
C                         IONE inserted in CALL to THOMAS
C                         Generic calls for EXP, COS
C Altered 12-Jun-1991 --- HBCNEW replaced by HBC_J and HBC_S. In order
C                         to improve convergence at outer boundary when wind
C                         is very thick. S1 deleted and use SOURCE(1) instead.
C                         Name changed from FQCOMP.
C
C Altered 21-Sep-1990 --- Bug fix: IBOUND was not being zeroed when THK
C                         option was off. Only effected direct computation of
C                         HBCNEW.
C Created 24-FEB-1987
C
	INTEGER NC,ND,NP,I,NI,LS
	REAL*8 TA(ND),TB(ND),TC(ND),XM(ND),R(ND),Z(ND)
	REAL*8 NEWRJ(ND),NEWRK(ND),DTAU(ND),dCHIdR(ND)
	REAL*8 SOURCE(ND),CHI(ND),AQW(ND,NP),AQW3(ND,NP),P(NP)
C
	INTEGER, PARAMETER :: IONE=1
C
	REAL*8 DBB,DBC,IBOUND,TOR,HBC_J,HBC_S,INBCNEW,IC,E1,E2,E3
	LOGICAL DIFF,THK
	CHARACTER*6 METHOD
C
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
C
C Zero all parameters.
C
	NEWRJ(1:ND)=0.0D0
	NEWRK(1:ND)=0.0D0
	HBC_J=0.0D0
	HBC_S=0.0D0
	INBCNEW=0.0D0
C
C Enter loop for each impact parameter P.
C
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*DSQRT( (R(ND)-P(LS))*(R(ND)+P(LS)) )/R(ND)/CHI(ND)
	  END IF
C
	  IF(THK)THEN
	    IF(P(LS) .GT. 0.0D0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796D0-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=SOURCE(1)*(1.0D0-EXP(-TOR))
	  ELSE
	    IBOUND=0.0D0
	  END IF
C
C Compute Z and the optical depth scale DTAU for this imapct parameter.
C
	  IF(NI .GT. 1)THEN
C	    CALL ZALONGP(R,Z,P(LS),NI)
C	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdr,NI)
	    DO I=1,NI
	      Z(I)=DSQRT( (R(I)-P(LS))*(R(I)+P(LS)) )
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
	  IF(NI .GT. 2)THEN
	    XM(1)=-IBOUND
	    TA(1)=0.0D0
	    TC(1)=1.0D0/DTAU(1)
	    TB(1)=-1.0D0-TC(1)
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TC(I)=1.0D0/DTAU(I)
	      TB(I)=-0.5D0*(DTAU(I-1)+DTAU(I))-TA(I)-TC(I)
	      XM(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5D0
	    END DO
C
	    IF(LS .LE. NC .AND. DIFF)THEN
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
	  ELSE IF(NI .EQ. 1)THEN
	    XM(1)=IBOUND
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
	  END IF
C
C Update the FA and FB matrices (see notes).
C
	  DO I=1,NI
	    NEWRJ(I)=NEWRJ(I)+AQW(I,LS)*XM(I)
	    NEWRK(I)=NEWRK(I)+AQW3(I,LS)*XM(I)
	  END DO
C
	  HBC_J=HBC_J+AQW(1,LS)*XM(1)*Z(1)/R(1)
	  HBC_S=HBC_S+AQW(1,LS)*IBOUND*Z(1)/R(1)
	  IF(NI .EQ. ND)THEN
	    INBCNEW=INBCNEW + AQW(ND,LS)*
	1    (XM(ND)-(XM(ND)-XM(ND-1))*Z(ND)/R(ND)/DTAU(ND-1))
	  END IF
C
2000	CONTINUE
C
C Compute the factor for the outer boundary condition.
C
	HBC_J=HBC_J/NEWRJ(1)
	HBC_S=HBC_S/SOURCE(1)
	INBCNEW=INBCNEW/(2.0D0*NEWRJ(ND)-IC)
C
C Compute the new Feautrier factors.
C
	DO I=1,ND
	  NEWRK(I)=NEWRK(I)/NEWRJ(I)
	END DO
C
C Compute the Q factors from F. The Q values are stored in NEWRJ.
C We now compute Q in JFEAUNEW so that we nolonger need to save Q (only f).
C We still compute Q here because
C    (i) Q is passed (as NEWRJ vector) and hence is corrupted by this
C        routine. We still need Q, however, for the varaition
C        routines.
C   (ii) Q is now completely consistent with the updated f value.
C
	CALL QFROMF(NEWRK,NEWRJ,R,TA,TB,ND)	!TA work vector
C
	RETURN
	END
