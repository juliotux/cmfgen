C
C Subroutine to compute JBAR and ZNET using the escape probability
C approximation. Six vectors are returned with new values.
C
C    JBAR  = Line mean intensity
C    ZNET  = Net rate in the line
C    VB    = Varition of ZNET with respect to CHIL.
C    VC    = Variation of ZNET with respect to ETAL.
C    BETA  = Sobolev Escape Probability.
C    BETAC = Sobolev Escape Probability for the continuum.
C
C Altered 28-Dec-1987 - (1.0-EXP(-GAM))/GAM is now calculated with
C                       improved accuracy to allow for cancellation when
C                       X is small.
C Altered 4-Jan-87 -  Thick option installed. Note that there are two extra
C                     parameters in the call statement.
C Altered 15-Jan-88 - Two extra vectors are returned in call. These are
C                     BETA, the escape probability, and B_c'u - the
C                     continuum term. These can be used to allow for the
C                     variation of Z with J. SOURCE installed in call rather
C                     than ETA.
C ALtered 20-JAN-88 - S1 removed and METHOD installed in call.
C
C Altered 8-Feb-88 - Bug fix. There was an error in the computaion of VB
C                    which was introduced on 28-Dec-87. Important for
C                    lines with (1-AV/S) .NE. to 1.
C
C Altered 12-May-88 - NV put at 150. Check put on ND
C Altered 26-Aug-91 - When computing VB(I) routine had
C                        +-AV(I)*T2/ETAL(I) )*AQW(I,LS)
C                     (ie two operators). Minus sign is correct.
C
	SUBROUTINE SOBJBAR(SOURCE,CHI,CHIL,ETAL,V,SIGMA,R,P,AQW,
	1                  JBAR,ZNET,VB,VC,BETA,BETAC,
	1                  FL,DIF,DBB,IC,THICK,NLF,NC,NP,ND,METHOD)
	IMPLICIT NONE
	INTEGER NLF,NC,NP,ND,NV
	PARAMETER (NV=150)
	REAL*8 SOURCE(ND),CHI(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP),AQW(ND,NP)
	REAL*8 JBAR(ND),ZNET(ND),VB(ND),VC(ND),BETA(ND),BETAC(ND)
	REAL*8 DBB,IC,FL
	CHARACTER*(*) METHOD
	LOGICAL DIF,THICK
C
	COMMON /SCRATCH/ TA,TB,TC,AV,DTAU,Z,GAM,GAMH,dCHIdR
	REAL*8 TA(NV),TB(NV),TC(NV),AV(NV),DTAU(NV),Z(NV)
	REAL*8 GAM(NV),GAMH(NV),dCHIdR(NV)
	DOUBLE PRECISION EXPONX,d_EXPONX_dX
C
C Local variables.
C
	INTEGER I,LS,NI
	REAL*8 T1,T2,DBC,TOR
C
C Check scratch block okay.
C
	IF(ND .GT. NV)THEN
	  WRITE(6,*)'Error in SOBJBAR - NV is too small'
	  WRITE(6,*)'Maximum number of grid points is 150'
	  STOP
	END IF
C
C Zero arrays which are incremented as we integrate over angle.
C Evaluate the SOBOLEV optical depth (GAMH) without angle factor.
C
	DO I=1,ND
	  JBAR(I)=0.0D0
	  ZNET(I)=0.0D0
	  VB(I)=0.0D0
	  VC(I)=0.0D0
	  BETA(I)=0.0D0
	  BETAC(I)=0.0D0
	  GAMH(I)=CHIL(I)*3.0D-10*R(I)/V(I)/FL    	!C/dex(15)/dex(5)
	END DO
C
	CALL DERIVCHI(dCHIdR,CHI,R,ND,METHOD)
C
C Enter loop to perform integration along each ray.
C
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
C
C Zero AV vector.
C
	  DO I=1,NI
	    AV(I)=0.0
	  END DO
C
	  CALL ZALONGP(R,Z,P(LS),NI)
	  DO I=1,NI
	    GAM(I)=GAMH(I)/(1.0D0+Z(I)*Z(I)/R(I)/R(I)*SIGMA(I))
	  END DO
C
	  IF(DIF .AND. LS .LE. NC)THEN
	   DBC=DBB*DSQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	  END IF
	  CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	  CALL TCOMPD(TA,TB,TC,DTAU,DIF,LS,NC,ND,NI)
	  CALL XVECD(DTAU,SOURCE,AV,DIF,DBC,IC,LS,NC,ND,NI)
C
C SOURCE(1) is the boundary continuum source function.
C
	  IF(THICK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796-DACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    AV(1)=-SOURCE(1)*(1.0D0-DEXP(-TOR))
	  END IF
C
C Solve for the radiation field along ray for this frequency.
C
	  CALL THOMAS(TA,TB,TC,AV,NI,1)
C
C The function EXPONX is given by (1.0-EXP(-X))/X.
C The function d_EXPONX_dX is given by d[ (1.0-EXP(-X))/X ]/dX. These
C functions are called to allow for cancellation when X is small. Note
C GAM(I)/CHIL(I) = d_GAM(I)/d_CHIL(I) .
C
	  DO I=1,NI
	    T1=1.0D0-AV(I)*CHIL(I)/ETAL(I)
	    T2=EXPONX(GAM(I))
	    ZNET(I)=ZNET(I)+T1*T2*AQW(I,LS)
C
	    BETA(I)=BETA(I)+T2*AQW(I,LS)
	    BETAC(I)=BETAC(I)+T2*AV(I)*CHIL(I)*AQW(I,LS)/ETAL(I)
C
	    VB(I)=VB(I)+( T1*d_EXPONX_dX(GAM(I))*GAM(I)/CHIL(I)-
	1     AV(I)*T2/ETAL(I) )*AQW(I,LS)
C
	    VC(I)=VC(I)+AV(I)*CHIL(I)/ETAL(I)/ETAL(I)*T2*AQW(I,LS)
C
	  END DO
	END DO
C
C Compute the mean line intensity
C
	DO I=1,ND
	  JBAR(I)=(1.0D0-ZNET(I))*ETAL(I)/CHIL(I)
	END DO
C
	RETURN
	END
