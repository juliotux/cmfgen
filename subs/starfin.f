C
C Routine to compute radius points to be used in the comoving frame
C integration. The radius points are chosen to be equally spaced in
C LOG(Tau) where Tau is assumed to be dominated by free-free
C processes and is consequently proportional to the integral of the
C density squared. The velocity is computed from the velocity law
C given by Mihalas with RN(=Radius at which V=0.5*Vinf), ALPHA
C and BETA as parameters which describe the law. ALPHA and BETA 
C are computed from using V=0.5*VINF at RN and V(RP)=VRP. The
C subroutine also returns the Velocity and Sigma(=R/V*dlnV/dlnR-1.)
C This version allows a two component velocity law.
C
	SUBROUTINE STARFIN(R,V,SIGMA,RMAX,RP1,RN1,VRP1,VINF1,
	1             EPS1,GAM1,RP2,RN2,VRP2,VINF2,EPS2,GAM2,ND,TA,TB,TC)
	IMPLICIT NONE
C
C Altered 26-May-1996 : IMPLICIT NONE installed.
C                       Generic calls for LOG and EXP
C Altered 10-MAR-1986 : Tau at boundary computed differently so that
C                        57 will have the same depth points as 30.
C Created 10-FEB-1983
C
	INTEGER ND
	REAL*8 R(ND),V(ND),SIGMA(ND)
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 RMAX,RP1,RN1,VRP1,VINF1
	REAL*8 EPS1,GAM1,RP2,RN2,VRP2,VINF2,EPS2,GAM2
C
C Local varibles.
C
	REAL*8 VF1,VF2
	REAL*8 X,Y,Z
	REAL*8 T1,T2,T3,T4
	REAL*8 DLNR,DLT
	REAL*8 ALPHA1,BETA1
C
	INTEGER I,MND
C
	X=RP2/RP1		!To correct for old velocity law (Change made
	Y=VRP2			! in main program
	Z=RN2/RP2
C
	VF1=0.5D0*(1.0D0+EPS1+GAM1)*VINF1
	VF2=0.5D0*(1.0D0+EPS2+GAM2)*VINF2
C
C Compute ALPHA1 and BETA1
C
	T4=VF1/VRP1-1.0D0
	T1=0.69314718
	T2=10000
	T3=RN1/RP1-1.0
	DO I=1,20
	  BETA1=LOG((T4-EPS1*EXP(T1*T3))/GAM1)/T3
	  ALPHA1=-LOG(0.5+(GAM1/2.0D0-0.5D0-GAM1*EXP(-BETA1))/EPS1)
	  IF(ABS(BETA1-T2)/BETA1+ABS(ALPHA1-T1)/ALPHA1 .LT. 0.0001)
	1 GOTO 100
	  T1=ALPHA1
	  T2=BETA1
	END DO
100	CONTINUE
C
C Compute opacity scale with a radius scale which is equally
C incrmented in LOG(R)
C
	MND=ND-2
	T1=LOG(RMAX)
	T2=LOG(RP1)
	DLNR=(T1-T2)/(ND-1)
	DO I=1,ND
	  TA(I)=EXP(T1-(I-1)*DLNR)	  	  !Radius
	  T3=RN1/TA(I)-1.0D0
	  V(I)=VF1/(1.0+EPS1*EXP(ALPHA1*T3)+GAM1*EXP(BETA1*T3))
	  V(I)=V(I)*(1.0D0+X*(1.0D0-0.999D0*(RP1/TA(I))**Y)**Z)
	  TB(I)=(1E+10/(TA(I)*TA(I)*V(I)))**2	  !Opacity
	END DO
C
C Compute optical depth scale.
C
C	T1=TB(1)*(TA(1)-TA(2))
	T1=TB(1)*TA(1)/3.0D0
	TC(1)=LOG(T1)
	DO I=2,ND
	  T1=T1+(TB(I)+TB(I-1))*(TA(I-1)-TA(I))
	  TC(I)=LOG(T1)
	END DO
	DLT=(TC(ND)-TC(1))/(MND-1)
	DO I=2,MND-1
	  TB(I)=TC(1)+(I-1)*DLT
	  END DO
	TB(1)=TC(1)
	TB(MND)=TC(ND)
C
C COMPUTE NEW RADIUS VALUES
C
	CALL LININT(TB,R,MND,TC,TA,ND)
C
C COMPUTE V AND SIGMA FOR THE NEW RADIUS VALUES.
C
	TA(1)=R(1)
	TA(2)=R(1)-(R(1)-R(2))/50.0
	DO I=2,MND-1
	  TA(I+1)=R(I)
	END DO
	TA(ND)=R(MND)
	TA(ND-1)=TA(ND)+(R(MND-1)-R(MND))/50.0D0
	DO I=1,ND
	  R(I)=TA(I)
	END DO
	DO I=1,ND
	  T1=EPS1*EXP(ALPHA1*(RN1/R(I)-1.0D0))
	  T2=GAM1*EXP(BETA1*(RN1/R(I)-1.0D0))
	  TA(I)=VF1/(1.0D0+T1+T2)
	  TB(I)=1.0D0+X*(1.0D0-0.999D0*(RP1/R(I))**Y)**Z
	  V(I)=TA(I)*TB(I)
	  SIGMA(I)=(VF1*RN1*(ALPHA1*T1+BETA1*T2)*TB(I)/R(I)/R(I)
	1 /(1.0D0+T1+T2)**2.0+
	1 TA(I)*X*Y*Z*0.999D0/RP1*(1.0D0-0.999D0*(RP1/R(I))**Y)**(Z-1)
	1 *(RP1/R(I))**(Y+1))*R(I)/V(I)-1.0D0
	END DO
C
	RETURN
	END
