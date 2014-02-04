C
C Routine to adjust the Radiative equilibrium equation for adiabatic cooling.
C This routine:
C
C  (1) Increments STEQ if adiabatic cooling is to be allowed for
C  (2) Incrments the variation matrix [BA] if it is being computed, and
C       if adiabatic cooling is being included.
C  (3) Evaluates the adiabatic cooling terms, splitting them into 2 terms ---
C       the dTdR term, and the velocity term. This is done for diagnostic
C       purposes.
C
 	SUBROUTINE EVAL_ADIABATIC(BA,STEQ,AD_CR_V,AD_CR_DT,
	1                           COMPUTE_BA,INCL_ADIABATIC,
	1                           POPION,ED,T,R,VEL,SIGMA,DTDR,DIFFW,
	1                           WORK,DIAG_INDX,NUM_BNDS,NT,ND)
 	IMPLICIT NONE
C
C Altered  15-Mar-2001 : Minor bug fiex in BA matrix at d=1.
C Created 256-Jul-1994
C
	INTEGER DIAG_INDX,NUM_BNDS,NT,ND
	REAL*8 BA(NT,NT,NUM_BNDS,ND)	
	REAL*8 STEQ(NT,ND)
	REAL*8 AD_CR_V(ND)
	REAL*8 AD_CR_DT(ND)
	REAL*8 POPION(ND)
	REAL*8 ED(ND)
	REAL*8 T(ND)
	REAL*8 R(ND)
	REAL*8 VEL(ND)
	REAL*8 SIGMA(ND)
	REAL*8 WORK(ND)
	REAL*8 DIFFW(NT)		!d(dT/dR)/d?
	REAL*8 DTDR			!dT/dR (Program units)
	LOGICAL COMPUTE_BA,INCL_ADIABATIC
C
	REAL*8 BOLTZMANN_CONSTANT,FUN_PI
	EXTERNAL BOLTZMANN_CONSTANT,FUN_PI
C
	INTEGER I,J,L,GET_DIAG
	REAL*8 SCALE,T1,PI
C
	GET_DIAG(I)=(NUM_BNDS/ND)*(I-DIAG_INDX)+DIAG_INDX
C
C Compute dlnT/dlnR using a simple difference formulae, which is
C second order accurate for equally space data. dlnT/dlnR is dimensionless.
C
	DO I=2,ND-1
	  AD_CR_DT(I)=LOG(T(I-1)/T(I+1))/LOG(R(I-1)/R(I+1))
	END DO
	AD_CR_DT(1)=LOG(T(1)/T(2))/LOG(R(1)/R(2))
	AD_CR_DT(ND)=R(ND)*DTDR/T(ND)
C
C For historical resions STEQ contains Int[chi.J - eta]dv. Rather than multiply
C this term everywhere by 4pi, we divide the adiabatic cooling rate by that
C term. Not also that the R units are unimportant, since R.Chi is dimensionless.
C The 10^9 arises since T is in units of 10^4K, and V in units of 10^5 km/s.
C
	IF(INCL_ADIABATIC)THEN
	  PI=FUN_PI()
	  SCALE=1.0D+09*BOLTZMANN_CONSTANT()/4.0D0/PI
	  DO I=1,ND
 	    WORK(I)=(POPION(I)+ED(I))*SCALE*VEL(I)*
	1                (1.5D0*AD_CR_DT(I)+SIGMA(I)+3.0D0)/R(I)
 	  END DO
	  DO I=1,ND
	    STEQ(NT,I)=STEQ(NT,I)-WORK(I)*T(I)
	  END DO
	  IF(COMPUTE_BA)THEN
	    L=GET_DIAG(1)
	    BA(NT,NT,L,1)=BA(NT,NT,L,1)-WORK(1)
	    BA(NT,NT-1,L,1)=BA(NT,NT-1,L,1)-WORK(1)*T(1)/(POPION(1)+ED(1))
	    IF(NUM_BNDS .GE. 3)THEN
	      BA(NT,NT,L,1)=BA(NT,NT,L,1)-
	1                         1.5D0*(POPION(1)+ED(1))*SCALE*VEL(1)/
	1                             R(1)/LOG(R(1)/R(2))
	      BA(NT,NT,L+1,1)=BA(NT,NT,L+1,1)+
	1                         1.5D0*(POPION(1)+ED(1))*SCALE*VEL(1)*
	1                             T(1)/R(1)/T(2)/LOG(R(1)/R(2))
  	    END IF
	    DO I=2,ND-1
	      L=GET_DIAG(I)
	      BA(NT,NT,L,I)=BA(NT,NT,L,I)-WORK(I)
	      BA(NT,NT-1,L,I)=BA(NT,NT-1,L,I)-WORK(I)*T(I)/(POPION(I)+ED(I))
	      IF(NUM_BNDS .GE. 3)THEN
	        BA(NT,NT,L-1,I)=BA(NT,NT,L-1,I)-
	1                           1.5D0*(POPION(I)+ED(I))*SCALE*VEL(I)*
	1                             T(I)/R(I)/T(I-1)/LOG(R(I-1)/R(I+1))
	        BA(NT,NT,L+1,I)=BA(NT,NT,L+1,I)+
	1                           1.5D0*(POPION(I)+ED(I))*SCALE*VEL(I)*
	1                             T(I)/R(I)/T(I+1)/LOG(R(I-1)/R(I+1))
	      END IF
	    END DO
	    T1=1.5D0*POPION(1)*SCALE*VEL(1)/R(1)
	    L=GET_DIAG(ND)
	    DO J=1,NT
	      BA(NT,J,L,ND)=BA(NT,J,L,ND)-T1*DIFFW(J)
	    END DO
	    BA(NT,NT,L,ND)=BA(NT,NT,L,ND)-WORK(ND)
	    BA(NT,NT-1,L,ND)=BA(NT,NT-1,L,ND)-WORK(ND)*
	1                                       T(ND)/(POPION(ND)+ED(ND))
	  END IF
	END IF
C
C
C Now compute the adiabatic cooling rate (in ergs/cm^3/sec) for diagnostic
C purposes. The rate is output to the COOLGEN file.
C
C The factor of 0.1 arises from T . V / R [ ( 10^4 . 10^5)/ (10^10) ]
C
C We split the adiabatic terms into 2 parts: The velocity term, and the
C dTdR term. This split is useful for diagnostic purposes.
C
	SCALE=0.1D0*BOLTZMANN_CONSTANT()
	DO I=1,ND
	  AD_CR_V(I)= SCALE*(POPION(I)+ED(I))*
	1                VEL(I)*T(I)/R(I)*(3.0D0+SIGMA(I))
	  AD_CR_DT(I)=SCALE*(POPION(I)+ED(I))*VEL(I)*T(I)/R(I)*
	1                (1.5D0*AD_CR_DT(I))
	END DO
C
	RETURN
	END
