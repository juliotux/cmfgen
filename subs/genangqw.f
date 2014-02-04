C
C Subroutine to compute the angular quadrature weights for the
C computation of the angular moments of the radiation field
C assuming spherical geometry.
C
	SUBROUTINE GENANGQW(QW,R,P,TA,TB,TC,NC,ND,NP,MOMWEIGHT,AT_HALF)
	IMPLICIT NONE
C
C 23-Dec-2004 : Only print out warning message when AT_HALF is FALSE.
C                 When AT_HALF is true, we always do extrapolation.
C 04-Sep-2002 : Check before doing SQRT installed (for INTEL compiler).
C 26-May-1996 : Call to DP_ZERO removed
C               ERROR_LU installed.
C               MOMWEIGHT declared external.
C 18-Jul-1991 - Commented out IMPLICIT NONE To overcome CRAY compiler bug.
C               Alternatively it could have been fixed by declaring
C               MOMWEIGHT to be EXTERNAL, and INTEGER.
C 11-Jan-1989 - Renamed to GENANGQW - Option AT_HALF installed.
C 25-Nov-1986 - Created
C
	LOGICAL AT_HALF
	INTEGER NC,ND,NP
	REAL*8 QW(ND,NP),R(ND),TA(NP),TB(NP),TC(NP),P(NP)
C
	INTEGER ERROR_LU,LUER
	EXTERNAL MOMWEIGHT,ERROR_LU
C
	INTEGER I,J,NW,MAXND,COUNTER
	REAL*8 T1,T2
C
C Initialization section
C
	QW(:,:)=0.0D0
	DO 10 I=1,NP
	  TA(I)=P(I)*P(I)
10	CONTINUE
	COUNTER=0
C
C***********************************************************************
C If AT_HALF is TRUE the angular quadrature weight are computed at the
C midpoints of the mesh.
C
C Quadrature weights are to be stored in QW(I,J) where I=1,ND-1
C signifies the radius midpoint (=I+0.5) and J=1,NW signifies which ray.
C This section of the routine is primarily used to compute the H and N
C quadrature weights for the routines where H and N are discretized
C at the mid points of the mesh.
C***********************************************************************
C***********************************************************************
C If AT_HALF is FALSE the angular quadrature weights are computed on
C the mesh.
C
C Quadrature weights are to be stored in QW(I,J) where I=1,ND
C signifies which radius and J=1,NW signifies which ray.
C
C***********************************************************************
C
	IF(AT_HALF)THEN
	  MAXND=ND-1
	ELSE
	  MAXND=ND
	END IF
C
	DO 300 I=1,MAXND
	  NW=NC+MAXND-I+1
	  IF(NW .GT. NP)NW=NP
	  IF(AT_HALF)THEN
	    T2=0.5D0*( R(I)+R(I+1) )
	    T1=T2*T2
	  ELSE
	    T2=R(I)
	    T1=T2*T2
	  END IF
	  DO 100 J=1,NW
	    TB(J)=0.0D0
	    IF(R(I) .NE. P(J))TB(J)=SQRT(T1-TA(J))/T2
100	  CONTINUE
	  CALL MOMWEIGHT(TB,TC,NW)
	  IF(TB(NW) .NE. 0.0D0)COUNTER=COUNTER+1
	  DO 200 J=1,NW
	    QW(I,J)=TC(J)
200	  CONTINUE
300	CONTINUE
C
	IF(COUNTER .NE. 0 .AND. .NOT. AT_HALF)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,500)COUNTER
500	  FORMAT(1X,'Warning -',I3,' extrapolations to zero in GENANGQW')
	END IF
C
	RETURN
	END
