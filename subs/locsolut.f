C
	SUBROUTINE LOCSOLUT(POPS,STEQ,TA,LU,NT,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 : DOUBLE PRECISION declaration removed.
C
	INTEGER ND,NT,LU
	REAL*8 STEQ(NT,ND),POPS(NT,ND),TA(ND)
C
	INTEGER I,J,K,M,MIN,MAX
	REAL*8 T1,T2,T3,SCALE
C
C Read in solution matrix (STEQ) , and scaling vector (TA).
C
	OPEN(UNIT=LU,FILE='SOLUTION',STATUS='OLD')
	  M=(ND-1)/10 + 1
	  DO J=1,M
	    MIN=(J-1)*10+1
	    MAX=J*10
	    IF(MAX .GT. ND)MAX=ND
	    READ(LU,*)(TA(I),I=MIN,MAX)
	    DO K=1,NT
	      READ(LU,*)(STEQ(K,I),I=MIN,MAX)
	    END DO
	  END DO
	CLOSE(UNIT=LU)
C
C In case one, all solution values are scaled. In case two only
C the local values are scaled.
C
	DO I=1,ND
	  IF(TA(I) .EQ. 100)THEN			!Compute local scaling.
	    SCALE=1.0D0
	    T1=0.95D0
	    T2=-10.0D0
	    DO J=1,NT-1
	      IF(STEQ(J,I) .GT. T1)T1=STEQ(J,I)		!NOTE + MEANS DECREASE
	      IF(STEQ(J,I) .LT. T2)T2=STEQ(J,I)		!NOTE - MEANS INCREASE
	    END DO
	    IF(T1 .GT. 0.95D0)SCALE=0.95D0/T1
	    T2=-10.0D0/T2
	    IF(T2 .LT. SCALE)SCALE=T2
C
C Check if T is a variable.
C
	    T3=0.2D0		!20%
	    IF(T3 .LT. ABS(STEQ(NT,I)))T3=ABS(STEQ(NT,I))
	    IF(T3 .GT. 1.0D-25)T3=0.2D0/T3
	    IF(T3 .LT. SCALE)SCALE=T3
C
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0D0-STEQ(J,I)*SCALE)
	    END DO
	  ELSE
	    SCALE=TA(I)
C
C Update the population levels (and the temperature) .
C
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0D0-STEQ(J,I)*SCALE)
	    END DO
	  END IF
	END DO
C
	RETURN
	END
