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
	SUBROUTINE HWEIGHT(X,W,N)
	IMPLICIT NONE
C
C Altered 24-May-1996 - Call to DP_ZERO deleted.
C                       ERROR_LU,LUER installed.
C Altered 25-Feb-1987 - Normalization implemented.
C Altered 14-Jan-87 - Program nolonger assumes that X(N)=0, since can
C                     assume X(N+1)=0 with W(N+1)=0.0 with out loss of
C                     accuracy.
C Created 25-Nov-1986 (Based on KWEIGHT)
C
	INTEGER N,I
	REAL*8 X(N),W(N),H,HN,RF,RE,SUM
C
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
	W(:)=0.0D0
C
	H=X(1)-X(2)
	HN=X(2)-X(3)
	RF=X(1)-X(3)
	W(1)=0.5D0*H*(X(1)-H/6.0D0)
	1       -H*H/RF/12.0D0*((2.0D0+HN/H)*X(1)-X(2))
	W(2)=0.5D0*H*(X(2)+H/6.0D0)
	1       +H*H/RF/12.0D0*(2.0D0+H/HN+HN/H)*X(1)
	W(3)=-H*H/RF/12.0D0*(X(2)+H/HN*X(1))
C
	DO I=3,N-1
	  H=X(I-1)-X(I)
	  RF=X(I-2)-X(I)
	  RE=X(I-1)-X(I+1)
	  W(I-2)=W(I-2)-X(I-1)*H*H/12.0D0/RF
	  W(I-1)=W(I-1)+0.5D0*H*(X(I-1)-H/6.0D0)
	1              +X(I)*H*H/12.0D0/RE
	  W(I)=W(I)+0.5D0*H*(X(I)+H/6.0D0)
	1              +X(I-1)*H*H/12.0D0/RF
	  W(I+1)=W(I+1)-X(I)*H*H/12.0D0/RE
	END DO
C
C Assumes that V(mu=0)=0 and mu.dV(mu)=0 at mu=0. 
C
	IF(X(N) .EQ. 0.0D0)THEN
	  H=X(N-1)
	  W(N-1)=W(N-1)+H*H/3.0D0
	ELSE
C
C Since V(mu=0) is zero, we dont actually need a ray with mu=0 since the weight
C is automatically zero.
C
C Integral from X(N-1) to X(N)
C
	  H=X(N-1)-X(N)
	  RF=X(N-2)-X(N)
	  RE=X(N-1)				!As X(N+1)=0
	  W(N-2)=W(N-2)-X(N-1)*H*H/12.0D0/RF
	  W(N-1)=W(N-1)+0.5D0*H*(X(N-1)-H/6.0D0)
	1              +X(N)*H*H/12.0D0/RE
	  W(N)=W(N)+0.5D0*H*(X(N)+H/6.0D0)
	1              +X(I-1)*H*H/12.0D0/RF
C
C Integral from X(N) to X(N+1)
C
	  H=X(N)
	  W(N)=W(N)+H*H/3.0D0
	  WRITE(lUER,*)'Warning - Extrapolation to zero required in HWEIGHT'
	END IF
C
C Ensure that the weights have the correct normalization (but dont
C perform the normalization). 
C
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	SUM=1.0D0/SUM/3.0D0
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)
	1  WRITE(LUER,*)' Warning - weights require normalization in HWEIGHT'
C
C
	RETURN
	END
