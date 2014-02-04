C
C Subroutine to compute the solution for the intensity along the last two
C rays. In this case there is only one or two points, and the normal tridiagonal
C solution does not work. Note that for the thick case, it is assumed that
C Theta << 1. Program will not give a good estimate if the electron scattering 
C opacity dominates when the total optical depth is much greater than unity.
C
	SUBROUTINE LAST2RAYS(XM,WM,R,Z,P,DTAU,ZETA,THETA,CHI,TOR,NI)
	IMPLICIT NONE
C
C Altered 24-May-1996. - DEXP replaced by EXP
C Altered 30-Jan-1989. - Two minor bug fixes:- DTAU and Z are now
C                        passed to the routine - they are no longer
C                        computed internally. Thus DTAU is effected
C                        by the METHOD option in the calling routine.
C                        A [ -THETA(1)*T1 ] term was inadvertantly
C                        omitted from the WM(1,1) term.
C
	INTEGER NI
	REAL*8 XM(NI),WM(NI,NI),THETA(NI),CHI(NI),P
	REAL*8 DTAU(NI),R(NI),Z(NI),ZETA(NI)
	REAL*8 TOR,E1,E2,E3,T1,IBOUND
C
C The normal boundary condition of no incident radiation is obtained when
C TOR is equal to ZERO. IF TOR is not ZERO, the thick boundary condition
C is applied.
C
	IF(TOR .GT. 0.01)THEN
	  T1=(1.0D0-EXP(-TOR))
	ELSE
	  T1=(1.D0-TOR/2.0*(1.0D0-TOR/3.0*(1.0-TOR/4.0)))*TOR
	END IF
	IBOUND=ZETA(1)*T1
C
	IF(NI .EQ. 1)THEN
	  XM(1)=ZETA(1)*T1
	  WM(1,1)=-THETA(1)*T1
	ELSE IF(NI .EQ. 2)THEN
	  E1=EXP(-DTAU(1))
	  E2=1.0D0-(1.0D0-E1)/DTAU(1)
	  E3=(1.0D0-E1)/DTAU(1)-E1
	  IF(DTAU(1) .LT. 1.0D-03)THEN
	    E2=DTAU(1)*0.5+DTAU(1)*DTAU(1)/6.0D0
	    E3=DTAU(1)*0.5-DTAU(1)*DTAU(1)/3.0D0
	  END IF
C
	  XM(2)=IBOUND*E1+ZETA(2)*E2+ZETA(1)*E3
          XM(1)=0.5*(IBOUND+XM(2)*E1+ZETA(1)*E2+ZETA(2)*E3)
C
          WM(2,2)=-THETA(2)*E2
          WM(2,1)=-THETA(1)*E3-THETA(1)*T1*E1
	  WM(1,1)=0.5*( WM(2,1)*E1-THETA(1)*(E2+T1) )
	  WM(1,2)=0.5*(WM(2,2)*E1-THETA(2)*E3)
	END IF
C
	RETURN
	END
