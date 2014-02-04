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
	SUBROUTINE STARPCYG_V4(R,V,SIGMA,RMAX,RP,
	1                 SCLHT,VCORE,VPHOT,VINF1,BETA1,EPS1,
	1                 VINF2,BETA2,EPS2,REF_RADIUS,
	1                 NBND_INS,CONS_FOR_TAU_SCL,EXP_FOR_TAU_SCL,
	1                 ND,TA,TB,TC,RDINR,LU)
	IMPLICIT NONE
C
C Altered 13-Jul-2001 - REF_RADIUS inserted to allow the core radius to be changed
C                         while allowing the same V(r) law to be kept.
C                         Changed to version 4.
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
        REAL*8 REF_RADIUS 
!
	INTEGER NBND_INS
	REAL*8 CONS_FOR_TAU_SCL
	REAL*8 EXP_FOR_TAU_SCL
!
        REAL*8 RP1,RP2,RP3,VEXT
	REAL*8 S1,S2
	REAL*8 V_RAT
	REAL*8 RPHOT
	INTEGER I,J,LOOP,MND,NUMSCL,NOLD,NDOLD
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	REAL*8 T1,DLNR,DLT
	LOGICAL RDINR
	CHARACTER*80 STRING
C
	MND=ND-2*NBND_INS
	VEXT=VINF2-VINF1
	IF(VEXT .LT. 1.0D-06*VINF1)VEXT=0.0D0
	SCLHT=REF_RADIUS*SCLHT
	RP1=REF_RADIUS*EPS1
	RP2=REF_RADIUS*EPS2
        RP3=REF_RADIUS
!
! Check whether the passed parameters are valid.
!
	IF(BETA1 .LT. 1.0D0 .AND. RP1 .GE. RP)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V3 --- Invalid EPS1'
          WRITE(LUER,*)'EPS1 should be approximately 0.999 for BETA1<1'
	  STOP
	END IF
	IF(BETA2 .LT. 1.0D0 .AND. RP2 .GE. RP)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V3 --- Invalid EPS2'
          WRITE(LUER,*)'EPS2 should be approximately 0.999 for BETA2<1'
	  STOP
	END IF
	IF(NBND_INS .LT. 1 .OR. NBND_INS .GT. 3)THEN
	  LUER=ERROR_LU()
          WRITE(LUER,*)'Error in STARPCYG_V3 --- Invalid NBND_INS'
          WRITE(LUER,*)'NBND_INS should be 1, 2, or 3'
	END IF

!
! Determine V_RAT using the fact that V(r)=VCORE at REF_RADIUS.
!
        V_RAT=VPHOT
	IF(EPS1 .LT. 1.0D0)THEN
          V_RAT=V_RAT+(VINF1-VPHOT)*(1.0D0-EPS1)**BETA1
	END IF
	IF(EPS2 .LT. 1.0D0)THEN
          V_RAT=V_RAT+VEXT*(1-EPS2)**BETA2
	END IF
	V_RAT=V_RAT/VCORE-1.0D0
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
	    TA(I)=VPHOT
	    S1=0.0D0; S2=0.0D0
	    IF(RP1/R(I) .LT. 1.0D0)THEN
	       TA(I)=TA(I)+(VINF1-VPHOT)*(1.0D0-RP1/R(I))**BETA1
	       S1=RP1*BETA1*(VINF1-VPHOT)*(1.0D0-RP1/R(I))**(BETA1-1)
	    END IF
	    IF(RP2/R(I) .LT. 1.0D0)THEN
	       TA(I)=TA(I)+ VEXT*(1.0D0-RP2/R(I))**BETA2
	       S2=RP2*BETA2*VEXT*(1.0D0-RP2/R(I))**(BETA2-1)    
	    END IF
	    TB(I)=1.0D0+V_RAT*EXP((RP3-R(I))/SCLHT)
	    V(I)=TA(I)/TB(I)
	    SIGMA(I)=(S1+S2)/R(I)/TA(I)
	1          +V_RAT*EXP((RP3-R(I))/SCLHT)*(R(I)/SCLHT/TB(I))-1.0D0
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
C
C Now compute the velocity, and density structure.
C We neglect the BETA component for R =< RPI. This is
C okay for BETA >= 1, since the derivative of tbe BETA
C velocity component at RP1 is zero.
C
	TA(1)=RMAX
	TA(ND)=RP
	DO I=1,ND
	  V(I)=VPHOT
	  IF(RP1/TA(I) .LT. 1.0D0)THEN
	     V(I)=V(I)+(VINF1-VPHOT)*(1.0D0-RP1/TA(I))**BETA1
	  END IF
	  IF(RP2/TA(I) .LT. 1.0D0)THEN
	     V(I)=V(I)+ VEXT*(1.0D0-RP2/TA(I))**BETA2
	  END IF
	  V(I)=V(I)/( 1.0D0+V_RAT*EXP((RP3-TA(I))/SCLHT) )
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
!
! Modify optical depth scale by velocity law. This
! is to allow better depth spacing.
!
	  DO I=1,ND
	    TC(I)=TC(I)-LOG(CONS_FOR_TAU_SCL+V(I)**EXP_FOR_TAU_SCL)
	  END DO
!
! Now compute the TAU scale equally spaced on log of the modified
! TAU scale.
!
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
	  DO I=2,MND-1
	    TA(I+NBND_INS)=R(I)
	  END DO
	  TA(1)=R(1)
	  TA(ND)=R(MND)
	  IF(NBND_INS .EQ. 1)THEN
	    TA(2)=R(1)-(R(1)-R(2))/20.0
	    TA(ND-1)=TA(ND)+(R(MND-1)-R(MND))/20.0D0
	  ELSE IF(NBND_INS .EQ. 2)THEN
	    TA(2)=R(1)-(R(1)-R(2))/10.0D0
	    TA(ND-1)=TA(ND)+(R(MND-1)-R(MND))/10.0D0
	    TA(3)=R(1)-(R(1)-R(2))/3.0D0
	    TA(ND-2)=TA(ND)+(R(MND-1)-R(MND))/3.0D0
	  ELSE IF(NBND_INS .EQ. 3)THEN
	    TA(2)=R(1)-(R(1)-R(2))/20.0D0
	    TA(ND-1)=TA(ND)+(R(MND-1)-R(MND))/20.0D0
	    TA(3)=R(1)-(R(1)-R(2))/8.0D0
	    TA(ND-2)=TA(ND)+(R(MND-1)-R(MND))/8.0D0
	    TA(4)=R(1)-(R(1)-R(2))/3.0D0
	    TA(ND-3)=TA(ND)+(R(MND-1)-R(MND))/3.0D0
	  END IF
C
	  DO I=1,ND
	    R(I)=TA(I)
	    TA(I)=VPHOT
	    S1=0.0D0; S2=0.0D0
	    IF(RP1/R(I) .LT. 1.0D0)THEN
	       TA(I)=TA(I)+(VINF1-VPHOT)*(1.0D0-RP1/R(I))**BETA1
	       S1=RP1*BETA1*(VINF1-VPHOT)*(1.0D0-RP1/R(I))**(BETA1-1)
	    END IF
	    IF(RP2/R(I) .LT. 1.0D0)THEN
	       TA(I)=TA(I)+ VEXT*(1.0D0-RP2/R(I))**BETA2
	       S2=RP2*BETA2*VEXT*(1.0D0-RP2/R(I))**(BETA2-1)    
	    END IF
	    TB(I)=1.0D0+V_RAT*EXP((RP3-R(I))/SCLHT)
	    V(I)=TA(I)/TB(I)
	    SIGMA(I)=(S1+S2)/R(I)/TA(I)
	1          +V_RAT*EXP((RP3-R(I))/SCLHT)*(R(I)/SCLHT/TB(I))-1.0D0
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
