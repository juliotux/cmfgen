C
C Routine to compute new radius grid on a given ray. Grid is equally spaced in 
C LOG(TAU) where TAU is evaluated using the variable CHI (usually the opacity
C or density). At the boudaries, a fine grid spacing can be adopted.
C
	SUBROUTINE NEW_R_SCALE_V2(R_NEW,ND,R,NI,CHI,P,RAMP_IN,RAMP_OUT)
	IMPLICIT NONE
C
C Altered 10-Oct-2004 --- Changed SQRT to (R-P)*(R+P)
C Altered 21-AUg-1997 --- Called V2, RAMP_IN and RAMP_OUT variables installed.
C                           Done to fix numerical problems with CMF_FORM_SOL.
C                           To evaluate TAU at outer boundary, CHI is now 
C                             assumed to vary as a powr law (rather than 1/r^2).
C                             Based on NEW_R_SCALE_V1.
C
	INTEGER ND,NI
	REAL*8 R_NEW(ND)	!New radius grid
	REAL*8 R(NI)		!Old radius grid
	REAL*8 CHI(NI)		!Opacity
	REAL*8 P		!Impact parameter (scaler)
	LOGICAL RAMP_IN		!Finer spacing at inner boundary?
	LOGICAL RAMP_OUT	!Finer spacing at outer boundary?
C
	REAL*8 Z(NI)
	REAL*8 TAU(NI)
C
	REAL*8 TAU_NEW(ND)
	REAL*8 Z_NEW(ND)
C
C Variables to compute TAU at outer boundary.
C
	REAL*8 ALPHA,DELR
	REAL*8 R_MAX,R_BIG,R_SMALL
	REAL*8 Z_BIG,Z_SMALL
	REAL*8 CHI_BIG,CHI_SMALL
	INTEGER IEND
C
C Variables to compute new TAU scale.
C
	INTEGER I,IST,ND_DIV,lU_ER
	REAL*8 DTAU
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C                             
C Compute path length variable Z.
C
	DO I=1,NI
	  Z(I)=SQRT( (R(I)-P)*(R(I)+P) )
	END DO
C
C Evaluate Tau at outer boundary, assuming that CHI scales as a power law,
C and at least as fast as 1/r^2. This is crude, but does not need to be
C very accurate.
C
	IEND=MIN(4,NI)
	ALPHA=LOG(CHI(IEND)/CHI(1))/LOG(R(1)/R(IEND))
	IF(ALPHA .GT. 2)ALPHA=2
C
C Gives nominally a few % accuracy for a power law.
C
	R_MAX=R(1)*100.0D0**(1.0D0/(ALPHA-1.0D0))
	DELR=EXP( LOG(R_MAX/R(1))/50.0D0 )
C
	TAU(1)=0.0D0
	R_BIG=R(1)
	Z_BIG=Z(1)
	DO WHILE(R_BIG .LT. R_MAX)
	  R_SMALL=R_BIG
	  R_BIG=R_SMALL*DELR
	  CHI_SMALL=CHI_BIG
	  CHI_BIG=CHI(1)*(R(1)/R_BIG)**ALPHA
	  Z_SMALL=Z_BIG
	  Z_BIG=SQRT( (R_BIG-P)*(R_BIG+P) )
	  TAU(1)=TAU(1)+0.5D0*(Z_BIG-Z_SMALL)*(CHI_BIG+CHI_SMALL)
	END DO
C
C Now compute the optical depth.
C
	DO I=2,NI
	  TAU(I)=TAU(I-1)+(CHI(I)+CHI(I-1))*(Z(I-1)-Z(I))*0.5D0
	END DO
C
C Compute new tau scale equally spaced in LOG TAU. We can ramp the step size
C at either boundary.
C
	ND_DIV=ND-1			!-1 as ND-1 intervals
	IF(RAMP_IN)ND_DIV=ND_DIV-3	!Insert 3 addit. points at inner bound.
	IF(RAMP_OUT)ND_DIV=ND_DIV-3  	!Insert 3 addit. points at outer bound.
	IF(ND_DIV .LT. 3)THEN
	  LU_ER=ERROR_LU()
	  WRITE(LU_ER,*)'Error in NEW_R_SCALE_V2'
	  WRITE(LU_ER,*)'ND is too small'
	  STOP
	END IF
C
	TAU(1:NI)=LOG( TAU(1:NI) )
	DTAU=(TAU(NI)-TAU(1))/ND_DIV
	TAU_NEW(1)=TAU(1)
	IF(RAMP_OUT)THEN
	  TAU_NEW(2)=TAU(1)+DTAU/64.0D0
	  TAU_NEW(3)=TAU(1)+DTAU/16.0D0
	  TAU_NEW(4)=TAU(1)+DTAU/4.0D0
	  IST=5
	ELSE
	  IST=2
	END IF
C
C IF RAM_IN is true, this does too many evaluations but the values are
C reset anyway.
C
	DO I=IST,ND-1
	  TAU_NEW(I)=TAU_NEW(1)+(I-IST+1)*DTAU
	END DO
C
	TAU_NEW(ND)=TAU(NI)
	IF(RAMP_IN)THEN
	  TAU_NEW(ND-1)=TAU(NI)-DTAU/64.0D0
	  TAU_NEW(ND-2)=TAU(NI)-DTAU/16.0D0
	  TAU_NEW(ND-3)=TAU(NI)-DTAU/4.0D0
	END IF
C
C Compute new Z scale.
C
	CALL LININT(TAU_NEW,Z_NEW,ND,TAU,Z,NI)
C
C Compute new R scale.
C
	R_NEW(1)=R(1)
	R_NEW(ND)=R(NI)
	DO I=2,ND-1
	  R_NEW(I)=SQRT(P*P+Z_NEW(I)*Z_NEW(I))
	END DO
C
	RETURN
	END
