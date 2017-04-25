!
! Routine to compute the equivalent width of a line and the 
! radiative force of a line using the SOBOLEV approximation.
! The radiative force is expressed as a multiple of the
! force on a free electron.
!
	SUBROUTINE SOBEW_GRAD_V2(SOURCE,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                 FORCE_MULT,RLUM,
	1                 AQW,HAQW,LINE_FLUX,EW,CONT_FLUX,
	1                 FL,INNER_BND_METH,DBB,IC,THICK,DIE,NC,NP,ND,METHOD)
	IMPLICIT NONE
!
! Altered 26-Oct-2014 : We now write ABS(dmuV/dz) so as to handle non-monotonic velocity laws.
! Altered 07-Jun-2010: Changed to THOMAS_RH to handle later stages of SN.
!                        Replaced DIF variable by INNER_BND_METH.
! Created 24-Oct-2001: Based on SOBEW
!
	INTEGER NC,NP,ND
	REAL*8 SOURCE(ND),CHI(ND),ESEC(ND),CHIL(ND),ETAL(ND)
	REAL*8 FORCE_MULT(ND)
	REAL*8 RLUM 
	REAL*8 V(ND),SIGMA(ND),R(ND),P(NP),AQW(ND,NP),HAQW(ND,NP)
	REAL*8 DBB,IC,FL,EW,CONT_FLUX,LINE_FLUX(ND)
	CHARACTER*(*) METHOD
	CHARACTER*(*) INNER_BND_METH
	LOGICAL THICK,DIE
!
! Use dynamic allocation for required vectors.
!
	REAL*8 TA(ND),TB(ND),TC(ND),AV(ND),DTAU(ND),Z(ND)
	REAL*8 GAM(ND),GAMH(ND),dCHIdR(ND),NOES(ND)
	REAL*8 GLINE(ND)
!
	REAL*8 EXPONX
!
! Local variables.
!
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER I,LS,NI
	REAL*8 T1,T2,T3,T4,TOR
	REAL*8 E1,E2,E3
	REAL*8 IBOUND
	REAL*8 JFLUX
!
	IF( INNER_BND_METH .NE. 'DIFFUSION' .AND.
	1   INNER_BND_METH .NE. 'SCHUSTER' .AND.
	1   INNER_BND_METH .NE. 'HOLLOW' .AND.
	1   INNER_BND_METH .NE. 'ZERO_FLUX' )THEN
	  WRITE(6,*)'Error in SOBEW_GRAD_V2'
	  WRITE(6,*)'Unrecognized INNER_BND_METH: ',TRIM(INNER_BND_METH)
	  STOP
	END IF
!
! Zero arrays which are incremented as we integrate over angle.
! Evaluate the SOBOLEV optical depth without angle factor (GAMH).
! Evaluate the thermal opacity.
!
	DO I=1,ND
	  LINE_FLUX(I)=0.0D0
	  NOES(I)=CHI(I)-ESEC(I)
	  GAMH(I)=CHIL(I)*3.0D-10*R(I)/V(I)/FL    	!C/dex(15)/dex(5)
	END DO
	EW=0.0D0
	CONT_FLUX=0.0D0
	GLINE(1:ND)=0.0D0
	JFLUX=0.0D0
!
	IF(DIE)THEN
	  DO I=1,ND
	    GAMH(I)=0.0D0
	  END DO
	END IF
!
	CALL DERIVCHI(dCHIdR,CHI,R,ND,METHOD)
!
! Enter loop to perform integration along each ray.
!
	DO LS=1,NP
	  NI=ND-(LS-NC-1)
	  IF(LS .LE. NC)NI=ND
!
! Zero AV vector.
!
	  AV(1:NI)=0.0
!
! SOURCE(1) is the boundary continuum source function. NB: If THICK it is
!           not obvious what you mean by the FLUX going to the observer.
!           Is it Iplus or Iplus-Ibound? Hopefully the outer boundary is
!           at sufficiently large radii that it makes no difference.
!
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
!
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
!
	    CALL ZALONGP(R,Z,P(LS),NI)
!
! NB: GAM is the angle dependent SOBOLEV optical depth.
!
	    DO I=1,NI
	      GAM(I)=GAMH(I)/ABS(1.0D0+Z(I)*Z(I)/R(I)/R(I)*SIGMA(I))
	    END DO
!
! Compute optical depth increments.
!
	    CALL NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
!
! Compute tridiagonal components.
! NB: TB(I)[tri matrix] = -[TB(I)+TA(I)+TC(I)]
! This section was modified from TCOMPD.
!
	    TA(1)=0.0D0
	    TB(1)=1.0D0
	    AV(1)=-IBOUND
	    DO I=1,NI-1
	      TC(I)=1.0D0/DTAU(I)
	    END DO
	    DO I=2,NI-1
	      TA(I)=TC(I-1)
	      TB(I)=0.5D0*(DTAU(I-1)+DTAU(I))
	      AV(I)=-SOURCE(I)*(DTAU(I-1)+DTAU(I))*0.5D0
	    END DO
!
	    IF(LS .GT. NC)THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-DTAU(NI-1)/2.0D0
	      AV(NI)=0.5D0*DTAU(NI-1)*SOURCE(NI)
	    ELSE IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=0.0D0
	      AV(NI)=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	    ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-1.0D0
	      AV(NI)=IC
	    ELSE
	      TA(NI)=-TC(NI-1)
	      TB(NI)=-0.5D0*DTAU(NI-1)
	      AV(NI)=0.5D0*DTAU(NI-1)*SOURCE(NI)
	    END IF
	    TC(NI)=0.0D0
!
!
! Solve for the radiation field along ray for this frequency.
!
	    CALL THOMAS_RH(TA,TB,TC,AV,NI,IONE)
!
! Compute the midpoint flux in order that we can compute the
! line EW.  TB and TC are not required further.
!
	    DO I=1,NI-1
	      TC(I)=(AV(I+1)-AV(I))/DTAU(I)
	    END DO
!
! Interpolate the flux onto the radius grid. We handle the end points by
! using the bondary conditions.
!
	    DO I=2,NI-1
	      TB(I)=(DTAU(I-1)*TC(I)+DTAU(I)*TC(I-1))/(DTAU(I-1)+DTAU(I))
	    END DO
	    TB(1)=AV(1)-IBOUND
	    IF(LS .LE. NC)THEN
	      IF(INNER_BND_METH .EQ. 'DIFFUSION')THEN
	        TB(NI)=DBB*SQRT(R(ND)*R(ND)-P(LS)*P(LS))/R(ND)/CHI(ND)
	       ELSE IF(INNER_BND_METH .EQ. 'SCHUSTER')THEN
	        TB(NI)=IC-AV(NI)
	      ELSE
	        TB(NI)=0.0D0            !Hollow core and zero_flux options.
	      END IF
	    ELSE
	      TB(NI)=0.0D0		!Symmetry
	    END IF
!
	  END IF			!NI > 2
!
! Update the continuum flux.
!
	  CONT_FLUX=CONT_FLUX+TB(1)*HAQW(1,LS)
	  JFLUX=JFLUX+AV(1)*AQW(1,LS)
!
! Compute the thermal optical depth scale. Can use DTAU since
! this is not required anymore.
!
	  DTAU(1)=NOES(1)*R(1)
	  DO I=2,NI
	    DTAU(I)=DTAU(I-1)+0.5D0*(NOES(I-1)+NOES(I))*(Z(I-1)-Z(I))
	  END DO
!
! The function EXPONX is given by (1.0-EXP(-X))/X.
!
	  DO I=1,NI
	    T2=EXPONX(GAM(I))
!
! The radiative force due to a single line is
!
!      (4pi chil)/(rohc) Int v(mu) Beta(mu) mu dmu
!
! The following just evaluates the integral part of this expression.
!
	    GLINE(I)=GLINE(I)+TB(I)*T2*HAQW(I,LS)
!
! NB: T3 is the optical depth from the oberver to the point I on the NEAR side.
!        of the atmosphere.
!
	    IF(DTAU(I) .GT. 50)THEN
	      T3=0.0D0
	    ELSE
	      T3=EXP(-DTAU(I))
	    END IF
!
! NB: T4 is the optical depth from the oberver to the point I on the READ side
!        of the atmosphere.
!
	    T4=DTAU(NI)+DTAU(NI)-DTAU(I)
	    IF(LS .LE. NC .OR. T4 .GT. 50)THEN
	      T4=0.0D0
	    ELSE
	      T4=EXP(-T4)
	    END IF
!
	    LINE_FLUX(I)=LINE_FLUX(I)
	1         +0.5D0*T2*AQW(I,LS)
	1         *(  (T3+T4)-( (AV(I)+TB(I))*T3+(AV(I)-TB(I))*T4 )
	1         *CHIL(I)/ETAL(I)  )
!
	  END DO
!
	END DO			!LS
!
! Evaluate the continuum intensity (Jy (d=1kpc)) and compute line
! emission function such that integral gives the equivalent width
! in Angstroms.
!
	CONT_FLUX=CONT_FLUX*R(1)*R(1)*13.19868
	T1=13.19868*(2.997924D-12/FL/FL)/CONT_FLUX
	DO I=1,ND
	  LINE_FLUX(I)=T1*LINE_FLUX(I)*ETAL(I)*( R(I)**3 )
	END DO
!
! Compute the line flux, and then evaluate the continuum intensity
! (Jy (d=1kpc)) and the line equivalent width (Angstroms).
!
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
!
	EW=EW*0.5D0
!
	T1=16.0D+20*(3.14159)**2/RLUM/3.845D+33
	DO I=1,ND
	  FORCE_MULT(I)=FORCE_MULT(I)+GLINE(I)*T1*R(I)*R(I)*CHIL(I)/ESEC(I)
	END DO
!
	RETURN
	END
