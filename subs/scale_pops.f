C
C Routine to scale the populations so that the species conservation 
C equation is satisfied. Here, SUM, which  must be supplied, is the 
C current species population obtained by summing over all sates.
C
C Created 11-Apr-1989
C
	SUBROUTINE SCALE_POPS(HYD,DHYD,POPHYD,SUM,NHYD,ND)
	IMPLICIT NONE
C
	INTEGER NHYD,ND
	REAL*8 HYD(NHYD,ND),DHYD(ND),POPHYD(ND),SUM(ND)
C
	INTEGER I,J
	REAL*8 T1
C
	IF(POPHYD(ND) .EQ. 0)RETURN
C
	DO J=1,ND
	  T1=POPHYD(J)/SUM(J)
	  DO I=1,NHYD
	    HYD(I,J)=HYD(I,J)*T1
	  END DO
	  DHYD(J)=DHYD(J)*T1
	END DO
C
	RETURN
	END
