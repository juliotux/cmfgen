C
C Subroutine to compute the quadrature weights for the J integration.
C A trapazoidal rule is used.
C
	SUBROUTINE JTRPWGT(X,W,N)
	IMPLICIT NONE
C
C Aletered 24-May-1996 - Call to DP_ZERO removed
C                        EROR_LU etc installed.
C Created  17-May-1989 - Based on HWEIGHT
C
	INTEGER N,I
	REAL*8 X(N),W(N),H,T1,SUM
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	W(:)=0.0D0
C
	DO I=1,N-1
	  H=0.5D0*(X(I)-X(I+1))
	  W(I)=W(I)+H
	  W(I+1)=W(I+1)+H
	END DO
C
C Assumes that j'=0
C
	IF(X(N) .NE. 0)THEN
C
C Integral from X(N) to 0
C
	  T1=X(N)/( X(N-1)*X(N-1)-X(N)*X(N) )
	  W(N-1)=W(N-1)-2.0D0*T1*X(N)*X(N)/3.0D0
	  W(N)=W(N)+T1*(X(N-1)*X(N-1)-X(N)*X(N)/3.0D0)
C
	END IF
C
C Ensure that the weights have the correct normalization (but dont
C perform the normalization). Two checks are done to insure that
C the correct answer is given for a linear varaition. Because of the
C assumption that du/dmu(mu=0)=0, we have to fiddle with the last check.
C
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - weights require normalization in JWEIGHT'
	  WRITE(LUER,*)'Scale=',SUM
	END IF
C
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	IF(X(N) .EQ. 0)THEN
	  T1=0.5D0
	ELSE
	  T1=X(N)/( X(N-1)*X(N-1)-X(N)*X(N) )
	  T1=T1*( X(N)*(X(N-1)*X(N-1)-X(N)*X(N)/3.0D0)
	1        -2.0D0*X(N)*X(N)*X(N-1)/3.0D0 )
	  T1=0.5D0*(1.0D0-X(N)*X(N))+T1
	END IF
	SUM=T1/SUM
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)THEN
	   LUER=ERROR_LU()
	   WRITE(LUER,*)' Warning - weights require normalization in JWEIGHT'
	END IF
C
	RETURN
	END
