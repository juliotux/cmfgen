	FUNCTION EXPN(N,X)
	IMPLICIT NONE
	REAL*8 EXPN,EXPI,X,Z
	EXTERNAL EXPI
	INTEGER I,J,K,N,IFAIL
C
	IFAIL=0
	IF(N .EQ. 1 .AND. X .LT. 10)THEN
	  EXPN=EXP(X)*EXPI(X)
	  RETURN
	ELSE IF(N .EQ. 1)THEN
	  EXPN=(X*X + 4.03640*X + 1.15198)/(X*X + 5.03637*X + 4.1916)/X
	  RETURN
	END IF
C	  
	IF(X .LT. 5)THEN
	  EXPN=EXP(X)*EXPI(X)
	  DO I=2,N
	    EXPN=(1.0D0-X*EXPN)/(I-1)
	  END DO
	ELSE IF(N .GT. 0)THEN
C
C We use recurrence relation in forward direction. 
C
	  EXPN=EXP(X)*EXPI(X)
	  DO I=2,N
	    EXPN=(1.0D0-X*EXPN)/(I-1)
	  END DO
	ELSE
C
C We use recurrence relation in reverse direction.
C
	  J=N+10
	  EXPN=(1.0/(X+J))
	  DO K=J-1,N,-1
	    EXPN=(1.0D0-K*EXPN)/X
	  END DO
	END IF
c
	RETURN
	END
