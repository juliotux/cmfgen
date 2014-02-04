C
C This routine solves for the mean intensity as
C a function of depth. A Schuster or diffusion approaximation
C is used for the lower boundary condition.
C
C Altered 28-Oct-1996 - Bug fix: COS corrected back to ACOS in TOR expression.
C Altered 24-May-1996 - IMPLCIT NONE Installed.
C                       IONE used for calls (eg SIMPTH).
C                       Call to DP_ZERO removed.
C
C Altered  6-Feb-1987 - Method made into a character string; dCHIdr inserted.
C Altered 13-Feb-1987 - (Method option installed to allow the calculation
C                        of an improved optical depth scale.
C Altered  9-Ded-1986 - (Call changed to NEWJSOLD and call to WMAT changed -
C                           Now NEWWMAT. Note that S1 has been deleted from
C                           call to NEWJSOLD).
C Altered  6-NOV-1986 - (NI=1 and NI=2 pts now handled)
C Altered 22-AUG-1982 - (THICK boundary condition added)
C Altered 20-APR-1984 - (THICK boundary condition improved)
C Altered 28-JUL-1982
C
	SUBROUTINE NEWJSOLD(TA,TB,TC,XM,WM,FB,RJ,DTAU,R,Z,P,
	1   ZETA,THETA,CHI,dCHIdR,AQW,
	1   THK,DIFF,DBB,IC,NC,ND,NP,METHOD)
	IMPLICIT NONE
C
	INTEGER NC,ND,NP
	REAL*8 TA(ND),TB(ND),TC(ND),XM(ND),WM(ND,ND)
	REAL*8 FB(ND,ND),RJ(ND),DTAU(ND),R(ND),Z(ND),ZETA(ND)
	REAL*8 THETA(ND),CHI(ND),dCHIdr(ND),AQW(ND,NP),P(NP)
	REAL*8 DBB,IC
C
	LOGICAL DIFF,THK
	CHARACTER*6 METHOD
C
	REAL*8 DBC,IBOUND,TOR
	INTEGER I,LS,NI,KS
	INTEGER, PARAMETER :: IONE=1
C
C Zero RJ vector and set FB matrix equal to a unit matrix of dimension
C (ND*ND) . Compute derivative of the opacity.
C
	FB(:,:)=0.0D0
	DO 40 I=1,ND
	  FB(I,I)=1.0D0
40	CONTINUE
	RJ(:)=0.0D0
	CALL DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
C
C ENTER LOOP FOR EACH IMPACT PARAMETER P
C
	DO 2000 LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)THEN
	    NI=ND
	    DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))
	1   /R(ND)/CHI(ND)
	  END IF
C
	IF(THK)THEN
	  IF(P(LS) .GT. 0)THEN
	    TOR=CHI(1)*R(1)*R(1)*(1.570796D0-ACOS(P(LS)/R(1)))/P(LS)
	  ELSE
	    TOR=CHI(1)*R(1)
	  END IF
	  IBOUND=ZETA(1)*(1.0D0-EXP(-TOR))
	ELSE
	  TOR=0.0D0
	  IBOUND=0.0D0
	END IF
C
C Compute Z and optical depth scale for this imapct parameter.
C
	IF(NI .GT. 1)THEN
	  CALL ZALONGP(R,Z,P(LS),NI)
	  CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	END IF
C
	IF(NI .GT. 2)THEN
C
C Compute T ( A tridiagonal matrix) and store it as three vectors
C TA,TB and TC .
C
	  CALL TCOMPD(TA,TB,TC,DTAU,DIFF,LS,NC,ND,NI)
C
C Compute WM and XM matrices. Thete is no need to pass the thick
C variable to WMAT since the non thick case is equivalent to
C TOR=0.0 .
C
	  CALL NEWWMAT(WM,DTAU,THETA,TOR,LS,NC,NI)
	  CALL XVECD(DTAU,ZETA,XM,DIFF,DBC,IC,LS,NC,ND,NI)
	  XM(1)=-IBOUND				!IBOUND=0 if not thick.
C
C solve the tridiagonal system of equations for the A and B
C matrices. (note after solution WXM contains both A and B)
C
	  CALL THOMAS(TA,TB,TC,WM,NI,NI)
	  CALL SIMPTH(TA,TB,TC,XM,NI,IONE)
C
	ELSE
	  CALL LAST2RAYS(XM,WM,R,Z,P(LS),DTAU,ZETA,THETA,CHI,TOR,NI)
	END IF
C
C Update the FA and FB matrices . (see notes)
C
	  CALL MULTVEC(RJ,RJ,AQW(1,LS),XM,NI)
	  CALL MULT2D(FB,AQW(1,LS),WM,ND,NI,IONE)
2000	CONTINUE
C
C Solve for mean intensity J as a function of depth.
C
	CALL SIMQ(FB,RJ,ND,KS)
C
	RETURN
	END
