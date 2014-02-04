C
C Routine to compute the matrices which depend on the integration
C at the previous frequency. It is assumed that ta,tb and tc
C (the tridiagonal matrix to be inverted) have already been
C modified by a call to THOMAS2D.
C
	SUBROUTINE UPDATEU(TX,TVX,KI,U,VB,VC,TA,TB,TC,WHAT,NI,NM)
C
C ALTERED 28-MAY-1996 : IMPLICIT NONE installed.
C                       Cleaned.
C Altered 27-JUL-1982 (-U was replaced everywhere by U)
C
	IMPLICIT NONE
C
	INTEGER NI,NM
	REAL*8 TX(NI,NI,NM),TVX(NI-1,NI,NM)
	REAL*8 KI(NI,3,NM),U(NI),VB(NI),VC(NI)
	REAL*8 TA(NI),TB(NI),TC(NI)
	LOGICAL WHAT
C
	INTEGER I,J,K
C
	DO 500 K=1,NM
C
C Check to see if the matrices need to be modified by the
C calculations at the preceeding frequency.
C
	  IF(WHAT)THEN
	    DO J=1,NI
	      TX(1,J,K)=U(1)*TX(1,J,K)
	      DO I=2,NI-1
	        TX(I,J,K)=U(I)*TX(I,J,K)-VB(I)*TVX(I-1,J,K)-VC(I)*TVX(I,J,K)
	      END DO
	      TX(NI,J,K)=U(NI)*TX(NI,J,K)-VB(NI)*TVX(NI-1,J,K)
	     END DO
	  ELSE
C
	    DO I=1,NI
	      DO J=1,NI
	        TX(J,I,K)=0.0
	      END DO
	    END DO
	    DO I=1,NI
	      DO J=1,NI-1
	        TVX(J,I,K)=0.0
	      END DO
	    END DO
C
	  END IF
C
	  TX(1,1,K)=KI(1,2,K)+TX(1,1,K)
	  TX(1,2,K)=KI(1,3,K)+TX(1,2,K)
	  DO I=2,NI-1
	    TX(I,I-1,K)=KI(I,1,K)+TX(I,I-1,K)
	    TX(I,I,K)=KI(I,2,K)+TX(I,I,K)
	    TX(I,I+1,K)=KI(I,3,K)+TX(I,I+1,K)
	  END DO
	  TX(NI,NI-1,K)=KI(NI,1,K)+TX(NI,NI-1,K)
	  TX(NI,NI,K)=KI(NI,2,K)+TX(NI,NI,K)
C
C Solve the simultaneous equations.
C
	  CALL SIMPTH(TA,TB,TC,TX(1,1,K),NI,NI)
500	CONTINUE
C
	RETURN
	END
