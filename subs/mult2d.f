	SUBROUTINE MULT2D(FB,TB,W,ND,NI,NM)
	IMPLICIT NONE
C
C Altered 05-Dec-1996 : END DO used to terminate DO LOOPs.
C ALtered 24-May-1996 : Implict none installed.
C
	INTEGER ND,NI,NM
	REAL*8 FB(ND,ND,NM),TB(NI),W(NI,NI,NM)
C
	INTEGER I,J,K
C
	DO K=1,NM
	  DO J=1,NI
	    DO I=1,NI
	      FB(I,J,K)=FB(I,J,K)+W(I,J,K)*TB(I)
	    END DO
	  END DO
	END DO
C
	RETURN
	END
