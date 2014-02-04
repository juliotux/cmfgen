C
C Routine to compute Qk,d and Qk,d+1/2 for the coomoving-frame
C transfer equation.
C
C Altered 27-Apr-1989 - Cleaned - Implicit none installed.
C
	SUBROUTINE QKIM(Q,QH,GAM,GAMH,TCHI,PF,ML,NI,NLF)
	IMPLICIT NONE
	INTEGER ML,NI,NLF
C
	REAL*8 Q(NI),QH(NI)
	REAL*8 GAM(NI),GAMH(NI),TCHI(NI),PF(NLF)
C
C Local varaibles.
C
	INTEGER I
	REAL*8 DNU
C
	IF(ML .NE. 1)THEN
	  DNU=PF(ML-1)-PF(ML)
	  DO 10 I=1,NI-1
	    QH(I)=GAMH(I)*2.0D0/((TCHI(I)+TCHI(I+1))*DNU)
	    Q(I)=GAM(I)/(TCHI(I)*DNU)
10	  CONTINUE
	  QH(NI)=0.0D0
	  Q(NI)=GAM(NI)/(TCHI(NI)*DNU)
	ELSE
C
C Assume first frquency-thus DNU is infinite.
C
	  DO 100 I=1,NI
	    Q(I)=0.0D0
	    QH(I)=0.0D0
100	  CONTINUE
	END IF
C
C
	RETURN
	END
