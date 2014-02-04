C
C Routine to update the matrices which describe the V flux
C equation at each frequency. Subroutine assumes there is no
C T variation.
C note that:-
C		1)VT CONTAINS PHI(v,T)
C
	SUBROUTINE UPVNOT(TVX,TX,RKB,RKC,VT,GB,H,NI,NM)
	IMPLICIT NONE
C
C Altered 28-May-1996 : IMPLICIT NONE installed.
C Altered 13-SEP-1982 : (variable NM)
C Created 25-JUL-1982
C
	INTEGER NI,NM
	REAL*8 TVX(NI-1,NI,NM)
	REAL*8 TX(NI,NI,NM),GB(NI),H(NI),RKB(NI),RKC(NI)
	REAL*8 VT(NI)
C
	INTEGER I
C
	CALL MARUPV(TVX,TX,GB,H,NI,NM)
C
	DO I=1,NI-1
	 TVX(I,I,1)=TVX(I,I,1)+RKB(I)*VT(I)
	 TVX(I,I+1,1)=TVX(I,I+1,1)+RKC(I)*VT(I+1)
	END DO
	 IF(NM .EQ. 4)THEN
	  DO I=1,NI-1
	   TVX(I,I,3)=TVX(I,I,3)+RKB(I)
	   TVX(I,I+1,3)=TVX(I,I+1,3)+RKC(I)
	  END DO
	 END IF
C
	RETURN
	END
