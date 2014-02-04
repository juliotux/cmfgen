C
C Subroutine to compute the quadrature weights for a modified cubic rule.
C Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
C The routine makes specific assumptions concerning the behaviour of the
C flux at the boundary, and would need to be moified to obtain  other
C quantities other than 1st moment of the intensity. The program assumes
C that the first point corresponds to mu=1.0 .
C
C Note that these weights should not be normalized in the usual fashion. 
C Physically, we dont expect V (the flux) to be constant with respect to mu.
C For small mu, we expect a that V is proportional to mu.
C
C D1, and R2 are work vectors.
C
	SUBROUTINE HTRPWGT(X,W,N)
	IMPLICIT NONE
C
C Altered 24-May-1996 - Call to DP_ZERO removed.
C                       ERROR_LU, LUER installed.
C Altered 25-Feb-1987 - Normalization implemented.
C Altered 14-Jan-87 - Program nolonger assumes that X(N)=0, since can
C                     assume X(N+1)=0 with W(N+1)=0.0 with out loss of
C                     accuracy.
C Created 25-Nov-1986 (Based on KWEIGHT)
C
	INTEGER N,I
	REAL*8 X(N),W(N),T1,T2,SUM
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	W(:)=0.0D0
	DO I=1,N-1
	  T1=0.5D0*(X(I)+X(I+1))
	  T2=(X(I)*X(I)+X(I)*X(I+1)+X(I+1)*X(I+1))/3.0D0
	  W(I)=W(I)+T2-X(I+1)*T1
	  W(I+1)=W(I+1)-T2+T1*X(I)
	END DO
C
C Assumes that V(mu=0)=0 and mu.V'(mu)=0 at mu=0. 
C
	IF(X(N) .NE. 0.0D0)THEN
C
C Since V(mu=0) is zero, we dont actually need a ray with mu=0 since the weight
C is automatically zero.
C
C Integral from X(N-1) to X(N)
C
	  W(N)=W(N)+X(N)*X(N)/3.0D0
C
	END IF
C
C Ensure that the weights have the correct normalization (but dont
C perform the normalization). Two checks are done to insure that
C the correct answer is given for a linear varaition. Because of the
C assumption that V(mu=0)=0, we have to fiddle with the last check.
C
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	SUM=1.0D0/SUM/3.0D0
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - weights require normalization in HWEIGHT'
	  WRITE(LUER,*)'Scale=',SUM
	END IF
C
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	IF(X(N) .EQ. 0)THEN
	  T1=0.5D0
	ELSE
	  T1=0.5D0*(1.0D0-X(N)*X(N))+(X(N)**2)/3.0D0
	END IF
	SUM=T1/SUM
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights require normalization in HWEIGHT'
	END IF
C
	RETURN
	END
