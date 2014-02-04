C
C Altered 28-May-1996 : DOUBLE PRECISION declarations removed.
C                       ERROR_LU installed.
C
	FUNCTION RAT(a,n)
	IMPLICIT NONE
	REAL*8 RAT
	INTEGER A,N,I
C
	IF(A .NE. 0)THEN
	  RAT=1.0D0
	  DO I=0,N-1
	    RAT=RAT*(a+I)
	  END DO
	ELSE
	  RAT=1.0D0
	END IF
C
	RETURN
	END



	FUNCTION HYPE(a1,b1,c,z)
	IMPLICIT NONE
	REAL*8 HYPE,RAT,Z
	INTEGER a1,b1,c,a,b,I
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	IF( B1 .GT. 0 .OR. A1 .GT. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Invalid variable in Hypergeometric'
	  STOP
	END IF
C
	IF(B1 .gt. A1)THEN
	  B=A1
	  A=B1
	ELSE
	  A=A1
	  B=B1
	END IF
C
	IF(A .NE. 0)THEN
	  HYPE=0
	  DO I=ABS(A),1,-1
	    HYPE=( HYPE+RAT(a,I)*RAT(b,I)/RAT(C,I) )* Z/FLOAT(I)
	  END DO
	  HYPE=HYPE+1.0D0
	ELSE
	  HYPE=1.0D0
	END IF
C
	RETURN
	END


	FUNCTION GAMRAT(A,B)
	IMPLICIT NONE
	REAL*8 GAMRAT
	INTEGER A,B,I
C
	GAMRAT=1.0D0
	DO I=B+1,A
	  GAMRAT=GAMRAT*I
	END DO
C
	RETURN
	END

	FUNCTION FAC(L)
	IMPLICIT NONE
	INTEGER I,L
	REAL*8 FAC
C
	FAC=1.0D0
	DO I=1,L
	  FAC=FAC*I
	END DO
C
	RETURN
	END

	FUNCTION RADSQ(N1,L1,NP1,LP1)
	IMPLICIT NONE
	REAL*8 RADSQ,GAMRAT,FAC,HYPE
	REAL*8 T1,T2,Z,MAX,MIN
	INTEGER N,L,NP
	INTEGER N1,L1,NP1,LP1
C
	IF(L1 .GT. LP1)THEN
	  L=L1
	  N=N1
!	  LP=LP1
	  NP=NP1
	ELSE
	  L=LP1
	  N=NP1
!	  LP=L1
	  NP=N1
	END IF	  
C
	MAX=N
	MIN=NP
	IF(MAX .LT. NP)THEN
	  MAX=NP
	  MIN=N
	END IF
C
	RADSQ=0.25D0/FAC(2*L-1)
	IF( MOD(NP-L,2) .NE. 0)RADSQ=-RADSQ
C
	T1=(4.0D0*MIN/MAX)**(L+1)
	1        *( (N-NP)/MAX )**(N+NP-2*L-2)
	1        /( (N+NP)/MAX )**(N+NP)
C
	Z=-4.0D0*N*NP/(N-NP)/(N-NP)
	T2=FLOAT(N-NP)/FLOAT(N+NP)
	T2=T2*T2
C
	RADSQ=RADSQ*SQRT( GAMRAT(N+L,N-L-1)*GAMRAT(NP+L-1,NP-L) )
	1      *T1
	1      *( HYPE(L-N+1,L-NP,2*L,Z) - 
	1      T2*HYPE(L-N-1,L-NP,2*L,Z) )
C
	RADSQ=RADSQ*RADSQ
C
	RETURN
	END