C
C Routine to compute radius points to be used in the comoving frame
C integration. The radius points are chosen to be equally spaced in
C LOG(Tau) where Tau is assumed to be dominated by free-free
C processes and is consequently proportional to the integral of the
C density squared.
C
	SUBROUTINE STARPCYG(R,V,SIGMA,RMAX,RP,GAMMA,SCLHT,
	1                 VCORE,VINF,VPHOT,ND,TA,TB,TC,RDINR,LU)
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
	INTEGER ND,LU,I,J,LOOP,MND,NUMSCL,NOLD,NDOLD
	REAL*8 R(ND),V(ND),SIGMA(ND),TA(ND),TB(ND),TC(ND)
	REAL*8 RMAX,RP,GAMMA,SCLHT,VCORE,VINF,VPHOT,BETA,RPHOT
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	REAL*8 T1,DLNR,DLT
	LOGICAL RDINR
	CHARACTER*80 STRING
C
	MND=ND-2
	SCLHT=RP*SCLHT
	BETA=(VPHOT/VCORE-1.0D0)
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
C Compute Velocity and SIGMA
	  DO I=1,ND-1
	    TA(I)=VPHOT+(VINF-VPHOT)*(1.0D0-RP/R(I))**GAMMA
	    TB(I)=1.0D0+BETA*EXP((RP-R(I))/SCLHT)
	    V(I)=TA(I)/TB(I)
	    SIGMA(I)=GAMMA*(VINF-VPHOT)*(1.0D0-RP/R(I))**(GAMMA-1)
	1                            *(RP/R(I))/TA(I)
	1          +BETA*EXP((RP-R(I))/SCLHT)*(R(I)/SCLHT/TB(I))-1.0D0
	  END DO
	  R(ND)=RP
	  V(ND)=VCORE
	  SIGMA(ND)=GAMMA*(VINF/VPHOT-1.0D0) +
	1                BETA*RP/SCLHT/(1.0D0+BETA)-1.0D0
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
	IF(BETA .GT. 1)THEN
          NUMSCL=LOG(BETA)+1.0
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
C
	DO I=1,ND-1
	  V(I)=VPHOT+(VINF-VPHOT)*(1.0D0-RP/TA(I))**GAMMA
	  V(I)=V(I)/(1.0D0+BETA*EXP((RP-TA(I))/SCLHT))
	  TB(I)=(1E+10/(TA(I)*TA(I)*V(I)))**2	  !Opacity
	END DO
C
C Set the improtant boundary values.
	TA(1)=RMAX
	TA(ND)=RP
	V(ND)=VCORE
	TB(ND)=(1E+10/(RP*RP*VCORE))**2  	!Opacity
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
C COMPUTE V AND SIGMA FOR THE NEW RADIUS VALUES.
C
	  TA(1)=R(1)
	  TA(2)=R(1)-(R(1)-R(2))/20.0
	  DO I=2,MND-1
	    TA(I+1)=R(I)
	  END DO
	  TA(ND)=R(MND)
	  TA(ND-1)=TA(ND)+(R(MND-1)-R(MND))/20.0D0
C
	  DO I=1,ND-1
	    R(I)=TA(I)
	    TA(I)=VPHOT+(VINF-VPHOT)*(1.0D0-RP/R(I))**GAMMA
	    TB(I)=1.0D0+BETA*EXP((RP-R(I))/SCLHT)
	    V(I)=TA(I)/TB(I)
	    SIGMA(I)=GAMMA*(VINF-VPHOT)*(1.0D0-RP/R(I))**(GAMMA-1)
	1                            *(RP/R(I))/TA(I)
	1          +BETA*EXP((RP-R(I))/SCLHT)*(R(I)/SCLHT/TB(I))-1.0D0
	  END DO
	  R(ND)=RP
	  V(ND)=VCORE
	  SIGMA(ND)=GAMMA*(VINF/VPHOT-1.0D0) +
	1                BETA*RP/SCLHT/(1.0D0+BETA)-1.0D0
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
