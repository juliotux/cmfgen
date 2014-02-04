C
C Routine to compute radius points to be used in the comoving frame
C integration. The radius points are chosen to be equally spaced in
C LOG(Tau) where Tau is assumed to be dominated by free-free
C processes and is consequently proportional to the integral of the
C density squared.
C
C This version allows a 2nd acceleration zone in the outer wind,
C specified by a second Beta like velocity law.
C
C The Teminal velcoity is VINF1+VEXT
C
	SUBROUTINE STARPCYG_TST(R,V,SIGMA,RMAX,RP,
	1                 SCLHT,VCORE,VPHOT,VINF1,BETA1,EPS1,
	1                 VINF2,BETA2,EPS2,ND,TA,TB,TC,RDINR,LU)
	IMPLICIT NONE
C
C Altered 07-Jul-1997 - We check RDINR file for '!Format date'
C Altered 05-Jun-1996 - T1 in second loop defining radius grid in vector TA
C                         was not set.
C Altered 26-May-1996 - Generic call for LOG installed.
C                       ERROR_LU installed.
C Altered 10-May-1995 - Bug in SIGMA for velocity law 3 at inner boundary.
C Altered 10-Apr-1989 - Reading of R values changed. Program
C                       allows for ION, V etc on the R line.
C                       LU for RDINR incorporated into call.
C Altered 05-Mar-1987 - R valus can be read in from file)
C Created 26-Feb-1987 - Based on STARNEW)
C
	INTEGER ND,LU
	REAL*8 R(ND),V(ND),SIGMA(ND),TA(ND),TB(ND),TC(ND)
C
	REAL*8 RMAX,RP
	REAL*8 SCLHT
	REAL*8 VCORE
	REAL*8 VPHOT
	REAL*8 VINF1
	REAL*8 BETA1
	REAL*8 EPS1
	REAL*8 VINF2
	REAL*8 BETA2
	REAL*8 EPS2
C
	REAL*8 RP1,RP2,VEXT
	REAL*8 V_RAT
	REAL*8 RPHOT
	REAL*8 AV,dAdlnR
	REAL*8 BV,dBdlnR
	INTEGER I,J,LOOP,MND,NUMSCL,NOLD,NDOLD
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	REAL*8 T1,DLNR,DLT
	LOGICAL RDINR
	CHARACTER*80 STRING
C
	MND=ND-2
	SCLHT=RP*SCLHT*EPS1
	VEXT=VINF2-VINF1
	IF(VEXT .LT. 1.0D-06*VINF1)VEXT=0.0D0
	RP1=RP*EPS1
	RP2=RP*EPS2
C
	IF(BETA1 .LT. 1.0D0 .AND. EPS1 .GE. 1.0D0)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V2 --- Invalid EPS1'
          WRITE(LUER,*)'EPS1 should be approximately 0.999 for BETA1<1'
	  STOP
	END IF
	IF(BETA2 .LT. 1.0D0 .AND. EPS2 .GE. 1.0D0)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V2 --- Invalid EPS2'
          WRITE(LUER,*)'EPS2 should be approximately 0.999 for BETA2<1'
	  STOP
	END IF
	IF(EPS1 .LT. 0.5)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V2 --- Invalid EPS1/BETA1'
          WRITE(LUER,*)'EPS1 should lie between 0.5 and 1.0'
	  STOP
	END IF
	IF(EPS2 .LT. 0.5)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V2 --- Invalid EPS2'
          WRITE(LUER,*)'EPS2 should lie between 0.5 and 1.0'
	  STOP
	END IF

C
	IF(EPS2 .EQ. 1.0D0)THEN
	  V_RAT=VPHOT/VCORE - 1.0D0
	ELSE
	  V_RAT= ( VPHOT+ VEXT*(1-EPS2)**BETA2 )/VCORE -1.0D0
	END IF
C
	IF(RDINR)THEN
	  OPEN(UNIT=LU,STATUS='OLD',FILE='RDINR')
C
C Check whether the file has a record containing 'Format date'. Its presence
C effects the way we read the file.
C
	  I=0
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
	    I=I+1
	    READ(LU,'(A)')STRING
	  END DO
	  IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU)
C
	  READ(LU,*)TA(1),TA(1),NOLD,NDOLD
C Check relative values.
	  IF(ND .NE. NDOLD)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error-NDOLD and ND are not equal in RDINR'
	    WRITE(LUER,*)'NDOLD=',NDOLD,' ND=',ND
	    STOP
	  END IF
C TA is used for everything but R which is all we want.
	  DO I=1,ND
	    READ(LU,*)R(I),TA(I),TA(I),TA(I)
	    READ(LU,*)(TA(J),J=1,NOLD)
	  END DO
	  R(1)=RMAX
C
C Compute Velocity and SIGMA
C
	  DO I=1,ND
	    IF(R(I) .EQ. RP1)THEN
	      IF(BETA1 .EQ. 1.0D0)THEN
	        TA(I)=VPHOT
	        dAdlnR=(VINF1-VPHOT)*RP1/R(I)
	      ELSE
	        TA(I)=VPHOT
	        dAdlnR=0.0D0
	      END IF
	      TB(I)=1.0D0+V_RAT
	      dBdlnR=V_RAT*R(I)/SCLHT
	    ELSE IF(R(I) .LT. RP1)THEN
	      TA(I)=VPHOT
	      TB(I)=1.0D0+V_RAT*EXP( (RP1-R(I))/SCLHT )
	      dAdlnR=0.0D0
	      dBdlnR=V_RAT*EXP( (RP1-R(I))/SCLHT )*R(I)/SCLHT
	    ELSE
	      TA(I)=VPHOT+(VINF1-VPHOT)*(1.0D0-RP1/R(I))**BETA1
	      TB(I)=1.0D0+V_RAT*EXP((RP1-R(I))/SCLHT)
	      dAdlnR=(RP1/R(I))*BETA1*(VINF1-VPHOT)*(1.0D0-RP1/R(I))**(BETA1-1)
	      dBdlNR=V_RAT*EXP((RP1-R(I))/SCLHT)*(R(I)/SCLHT)
	    END IF
	    IF(R(I) .EQ. RP2)THEN
	      IF(BETA2 .EQ. 1.0D0)THEN
	        dAdlnR=dAdlnR+VEXT*RP2/R(I)
	      END IF
	    ELSE IF(R(I) .LT. RP2)THEN
	    ELSE
	      TA(I)=TA(I)+VEXT*(1.0D0-RP2/R(I))**BETA2
	      dAdlnR=dAdlnR+(RP2/R(I))*BETA2*VEXT*(1.0D0-RP2/R(I))**(BETA2-1)
	    END IF
	    V(I)=TA(I)/TB(I)
	    SIGMA(I)=dAdlnR/TA(I)+dBdlnR/TB(I)-1.0D0
	  END DO
	  R(ND)=RP
	  CLOSE(UNIT=LU)
	  RETURN
	END IF
C
C Compute opacity scale with a radius scale which is equally
C incremented in LOG(R). Because of the exponetial with scale
C height SCLHT, extra points are inserted in the zone where it is
C important. Note that RPHOT as used to denote the radius beyond which
C the exponential density component becomes insignificant.
C
	IF(V_RAT .GT. 1)THEN
          NUMSCL=LOG(V_RAT)+1.0
	  RPHOT=RP+SCLHT*NUMSCL
	  T1=LOG(RMAX)
	  DLNR=LOG(RMAX/RPHOT)/(ND-1-NUMSCL)
	  DO I=1,ND-1-NUMSCL
	    TA(I)=EXP(T1-(I-1)*DLNR)	  	  !Radius
	  END DO
	  DO I=ND-NUMSCL,ND-1
	    TA(I)=RPHOT-(I-ND+NUMSCL)*SCLHT
	  END DO
	ELSE
	  T1=LOG(RMAX)
	  DLNR=LOG(RMAX/RP)/(ND-1)
	  DO I=1,ND-1
	    TA(I)=EXP(T1-(I-1)*DLNR)	  	  !Radius
	  END DO
	END IF
!
! Set the improtant boundary values.
!
	TA(1)=RMAX
	TA(ND)=RP
C
C Now compute the velocity, and density structure.
C
	DO I=1,ND
	  IF(TA(I) .EQ. RP1)THEN
	    AV=VPHOT
	  ELSE IF(TA(I) .LT. RP1)THEN
	    AV=VPHOT
	    BV=1.0D0+V_RAT*EXP( (RP1-TA(I))/SCLHT )
	  ELSE
	    AV=VPHOT+(VINF1-VPHOT)*(1.0D0-RP1/TA(I))**BETA1
	    BV=1.0D0+V_RAT*EXP((RP1-TA(I))/SCLHT)
	  END IF
	  IF(TA(I) .GT. RP2)THEN
	    AV=AV+VEXT*(1.0D0-RP2/TA(I))**BETA2
	  END IF
	  V(I)=AV/BV
	  TB(I)=(1E+10/(TA(I)*TA(I)*V(I)))**2	  !Opacity
	END DO
C
C Do this procedure twice, as the exponential depth scale wont
C be well resolved on the first pass.
C
C
	DO LOOP=1,2
C
C Compute optical depth scale.
C
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
C Compute V and SIGMA for the new radius values.
C
	  TA(1)=R(1)
	  TA(2)=R(1)-(R(1)-R(2))/20.0
	  DO I=2,MND-1
	    TA(I+1)=R(I)
	  END DO
	  TA(ND)=R(MND)
	  TA(ND-1)=TA(ND)+(R(MND-1)-R(MND))/20.0D0
C
	  DO I=1,ND
	    R(I)=TA(I)
	    IF(R(I) .EQ. RP1)THEN
	      IF(BETA1 .EQ. 1.0D0)THEN
	        TA(I)=VPHOT
	        dAdlnR=(VINF1-VPHOT)*RP1/R(I)
	      ELSE
	        TA(I)=VPHOT
	        dAdlnR=0.0D0
	      END IF
	      dBdlnR=V_RAT*R(I)/SCLHT
	      TB(I)=1.0D0+V_RAT
	    ELSE IF(R(I) .LT. RP1)THEN
	      TA(I)=VPHOT
	      TB(I)=1.0D0+V_RAT*EXP( (RP1-R(I))/SCLHT )
	      dAdlnR=0.0D0
	      dBdlnR=V_RAT*EXP( (RP1-R(I))/SCLHT )*R(I)/SCLHT
	    ELSE
	      TA(I)=VPHOT+(VINF1-VPHOT)*(1.0D0-RP1/R(I))**BETA1
	      TB(I)=1.0D0+V_RAT*EXP((RP1-R(I))/SCLHT)
	      dAdlnR=(RP1/R(I))*BETA1*(VINF1-VPHOT)*(1.0D0-RP1/R(I))**(BETA1-1)
	      dBdlNR=V_RAT*EXP((RP1-R(I))/SCLHT)*(R(I)/SCLHT)
	    END IF
	    IF(R(I) .EQ. RP2)THEN
	      IF(BETA2 .EQ. 1.0D0)THEN
	        dAdlnR=dAdlnR+VEXT*RP2/R(I)
	      END IF
	    ELSE IF(R(I) .LT. RP2)THEN
	    ELSE
	      TA(I)=TA(I)+VEXT*(1.0D0-RP2/R(I))**BETA2
	      dAdlnR=dAdlnR+(RP2/R(I))*BETA2*VEXT*(1.0D0-RP2/R(I))**(BETA2-1)
	    END IF
	    V(I)=TA(I)/TB(I)
	    SIGMA(I)=dAdlnR/TA(I)+dBdlnR/TB(I)-1.0D0
	  END DO
	  R(ND)=RP
C
	  IF(LOOP .EQ. 2)RETURN
C
	  DO I=1,ND
	    TA(I)=R(I)
	    TB(I)=(1E+10/(TA(I)*TA(I)*V(I)))**2.0	    !Opacity
	  END DO
	END DO
C
	END
