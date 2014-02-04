C
C Subroutine to compute the angular quadrature weights for the
C computation of the mean intensity. The weights are based on
C the technique illustrated by Nordulund.
C
	SUBROUTINE NORDANGQW(QW,R,P,TA,TB,TC,NC,ND,NP,MOMWEIGHT)
	IMPLICIT NONE
C
C Altered 06-Sep-2002 - Check before taking SQRT installed (for INTEL compiler).
C Altered 24-May-1996 - IMPLICIT NONE installed.
C Created 25-Nov-1986
C
	INTEGER NC,ND,NP
	REAL*8 QW(ND,NP),R(ND),TA(NP),TB(NP),TC(NP),P(NP)
C
C Local variables
C
	INTEGER I,J,NW
	REAL*8 T1
C
	DO 10 I=1,NP
	  TA(I)=P(I)*P(I)
10	CONTINUE
C
C Quadrature weights are to be stored in QW(I,J) where I=1,ND
C signifies which radius and J=1,NW signfies which ray.
C
	DO 300 I=1,ND
	  NW=NC+ND-I+1
	  IF(NW .GT. NP)NW=NP
	  T1=R(I)*R(I)
	  DO 100 J=1,NW
	    TB(J)=0.0D0
	    IF(R(I) .NE. P(J))TB(J)=SQRT(T1-TA(J))/R(I)
100	  CONTINUE
	  CALL MOMWEIGHT(TB,TC,NW)
	  DO 200 J=1,NW
	    QW(I,J)=TC(J)
200	  CONTINUE
300	CONTINUE
C
	RETURN
	END
