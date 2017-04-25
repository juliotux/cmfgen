C
C Routine to compute the Equivalent width of a line using the SOBOLEV
C approximation.
C
	SUBROUTINE SOBEW(SOURCE,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                 AQW,HAQW,LINE_FLUX,EW,CONT_FLUX,
	1                 FL,DIF,DBB,IC,THICK,DIE,NC,NP,ND,METHOD)
	IMPLICIT NONE
C
C Altered 26-Oct-2014 : We now write ABS(dmuV/dz) so as to handle non-monotonic velocity laws.
C Altered 19-Nov-1998 : Adapted to allow only 1 or 2 points along a r ray.
C                         (i.e. so can be run with NP=ND+NC).
C                         Changes follow SOBJBAR_SIM.
C                         Computation of TB(1) and TB(NI) [Fluxes] improved.
C Altered 28-OCt-1996 : COS converted back to ACOS in TOR expression.
C Altered 26-May-1996 : Dynamic allocation used for scratch vectors.
C                       Generic calls ued for EXP, SQRT
C                       IONE passed in call to THOMAS.
C
C Altered 13-Nov-1988 : LINEFLUX returned in call. Can be used to illustrate
C                           line origin.
C Created 6-Oct-1988  : Based on SOBJBARANDEW which in turn based on SOBJBAR
C
	INTEGER NC,NP,ND
	REAL*8 SOURCE(ND),CHI(ND),ESEC(ND),CHIL(ND),ETAL(ND)
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP),AQW(ND,NP),HAQW(ND,NP)
	REAL*8 DBB,IC,FL,EW,CONT_FLUX,LINE_FLUX(ND)
	CHARACTER*(*) METHOD
	LOGICAL DIF,THICK,DIE
C
C Use dynamic allocation for required vectors.
C
	REAL*8 TA(ND),TB(ND),TC(ND),AV(ND),DTAU(ND),Z(ND)
	REAL*8 GAM(ND),GAMH(ND),dCHIdR(ND),NOES(ND)
C
	REAL*8 EXPONX
C
C Local variables.
C
	INTEGER, PARAMETER :: IONE=1
C
	INTEGER I,LS,NI
	REAL*8 T1,T2,T3,T4,DBC,TOR
	REAL*8 E1,E2,E3
	REAL*8 IBOUND
C
C Zero arrays which are incremented as we integrate over angle.
C Evaluate the SOBOLEV optical depth without angle factor (GAMH).
C Evaluate the thermal opacity.
C
	DO I=1,ND
	  LINE_FLUX(I)=0.0D0
	  NOES(I)=CHI(I)-ESEC(I)
	  GAMH(I)=CHIL(I)*3.0D-10*R(I)/V(I)/FL    	!C/dex(15)/dex(5)
	END DO
	EW=0.0D0
	CONT_FLUX=0.0D0
C
	IF(DIE)THEN
	  DO I=1,ND
	    GAMH(I)=0.0D0
	  END DO
	END IF
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
	  AV(1:NI)=0.0
C
C SOURCE(1) is the boundary continuum source function. NB: If THICK it is
C           not obvious what you mean by the FLUX going to the observer.
C           Is it Iplus or Iplus-Ibound? Hopefully the outer boundary is
C           at sufficiently large radii that it makes no difference.
C
	  IF(THICK)THEN
	    IF(P(LS) .GT. 0)THEN
	      TOR=CHI(1)*R(1)*R(1)*(1.570796-ACOS(P(LS)/R(1)))/P(LS)
	    ELSE
	      TOR=CHI(1)*R(1)
	    END IF
	    IBOUND=SOURCE(1)*(1.0D0-EXP(-TOR))
	  ELSE 
	    IBOUND=0
	  END IF
C
	  IF(NI .EQ. 1)THEN
	    AV(1)=IBOUND
	    TB(1)=0.0D0			!Flux at grid point.
	  ELSE IF(NI .EQ. 2)THEN
	    CALL ZALONGP(R,Z,P(LS),NI)
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    E1=EXP(-DTAU(1))
	    E2=1.0D0-(1.0D0-E1)/DTAU(1)
	    E3=(1.0D0-E1)/DTAU(1)-E1
	    IF(DTAU(1) .LT. 1.0D-03)THEN
	      E2=DTAU(1)*0.5+DTAU(1)*DTAU(1)/6.0D0
	      E3=DTAU(1)*0.5-DTAU(1)*DTAU(1)/3.0D0
	    END IF
	    AV(2)=IBOUND*E1+SOURCE(2)*E2+SOURCE(1)*E3
            AV(1)=0.5*(IBOUND+AV(2)*E1+SOURCE(1)*E2+SOURCE(2)*E3)
            TB(2)=0.0D0
	    TB(1)=AV(1)-IBOUND
	  ELSE
C
	    CALL ZALONGP(R,Z,P(LS),NI)
C
C NB: GAM is the angle dependent SOBOLEV optical depth.
C
	    DO I=1,NI
	      GAM(I)=GAMH(I)/ABS(1.0D0+Z(I)*Z(I)/R(I)/R(I)*SIGMA(I))
	    END DO
C
	    IF(DIF .AND. LS .LE. NC)THEN
	      DBC=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	    END IF
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	    CALL TCOMPD(TA,TB,TC,DTAU,DIF,LS,NC,ND,NI)
	    CALL XVECD(DTAU,SOURCE,AV,DIF,DBC,IC,LS,NC,ND,NI)
	    AV(1)=-IBOUND
C
C Solve for the radiation field along ray for this frequency.
C
	    CALL THOMAS(TA,TB,TC,AV,NI,IONE)
C
C Compute the midpoint flux in order that we can compute the
C line EW.  TB and TC are not required further.
C
	    DO I=1,NI-1
	      TC(I)=(AV(I+1)-AV(I))/DTAU(I)
	    END DO
C
C Interpolate the flux onto the radius grid. We handle the end points by
C using the bondary conditions.
C
	    DO I=2,NI-1
	      TB(I)=(DTAU(I-1)*TC(I)+DTAU(I)*TC(I-1))
	1           /(DTAU(I-1)+DTAU(I))
	    END DO
	    TB(1)=AV(1)-IBOUND
	    TB(NI)=0.0D0		     !Symmetry
	    IF(LS .LE. NC .AND. DIF)THEN
	      TB(NI)=DBC
	    ELSE IF(LS .LE. NC)THEN
	      TB(NI)=IC-AV(NI)
	    END IF
	  END IF			!NI > 2
C
C Update the continuum flux.
C
	  CONT_FLUX=CONT_FLUX+TB(1)*HAQW(1,LS)
C
C Compute the thermal optical depth scale. Can use DTAU since
C this is not required anymore.
C
	  DTAU(1)=NOES(1)*R(1)
	  DO I=2,NI
	    DTAU(I)=DTAU(I-1)+0.5D0*(NOES(I-1)+NOES(I))*(Z(I-1)-Z(I))
	  END DO
C
C The function EXPONX is given by (1.0-EXP(-X))/X.
C
	  DO I=1,NI
	    T2=EXPONX(GAM(I))
C
C NB: T3 is the optical depth from the oberver to the point I on the NEAR side.
C        of the atmosphere.
C
	    IF(DTAU(I) .GT. 50)THEN
	      T3=0.0D0
	    ELSE
	      T3=EXP(-DTAU(I))
	    END IF
C
C NB: T4 is the optical depth from the oberver to the point I on the READ side
C        of the atmosphere.
C
	    T4=DTAU(NI)+DTAU(NI)-DTAU(I)
	    IF(LS .LE. NC .OR. T4 .GT. 50)THEN
	      T4=0.0D0
	    ELSE
	      T4=EXP(-T4)
	    END IF
C
	    LINE_FLUX(I)=LINE_FLUX(I)
	1         +0.5D0*T2*AQW(I,LS)
	1         *(  (T3+T4)-( (AV(I)+TB(I))*T3+(AV(I)-TB(I))*T4 )
	1         *CHIL(I)/ETAL(I)  )
C
	  END DO
C
	END DO			!LS
C
C Evaluate the continuum intensity (Jy (d=1kpc)) and compute line
C emission function such that integral gives the equivalent width
C in Angstroms.
C
	CONT_FLUX=CONT_FLUX*R(1)*R(1)*13.19868
	T1=13.19868*(2.997924D-12/FL/FL)/CONT_FLUX
	DO I=1,ND
	  LINE_FLUX(I)=T1*LINE_FLUX(I)*ETAL(I)*( R(I)**3 )
	END DO
C
C Compute the line flux, and then evaluate the continuum intensity
C (Jy (d=1kpc)) and the line equivalent width (Angstroms).
C
	EW=0.0D0
	T2=(LINE_FLUX(1)-LINE_FLUX(2))/LOG(R(1)/R(2))
	DO I=1,ND-2
	  T1=T2
	  T2=(LINE_FLUX(I)-LINE_FLUX(I+2))/LOG(R(I)/R(I+2))
	  EW=EW+LOG(R(I)/R(I+1))*( LINE_FLUX(I)+LINE_FLUX(I+1)
	1      +LOG(R(I)/R(I+1))*(T2-T1)/6.0 )
	END DO
	T1=T2
	T2=(LINE_FLUX(ND-1)-LINE_FLUX(ND))/LOG(R(ND-1)/R(ND))
	EW=EW+LOG(R(ND-1)/R(ND))*( LINE_FLUX(ND-1)+LINE_FLUX(ND)
	1            +LOG(R(ND-1)/R(ND))*(T2-T1)/6.0 )
C
	EW=EW*0.5D0
C
	RETURN
	END
