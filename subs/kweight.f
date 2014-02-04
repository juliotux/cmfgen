C
C Subroutine to compute the quadrature weights for a modified cubic rule.
C Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
C The routine makes specific assumptions concerning the behaviour of the
C mean intensity at the boundary, and would need to be moified to obtain
C other quantities other than 2nd moment of the intensity. The program assumes
C that the first point corresponds to mu=1.0 .
C
	SUBROUTINE KWEIGHT(X,W,N)
	IMPLICIT NONE
C
C Altered  24-May-1996 Call to DP_ZERO removed
C                      ERROR_LU etc installed.
C Modified 25-Feb-1987 (Extrapolation to mu=0 included for so that
C                      the angle integration at the outer boundary is
C                      handled correctly. Program assumes that U is of the
C                      for U= a + b*mu*mu .
C Created 25-Nov-1986 (Based on NORDWEIGHT)
C
	INTEGER N,I
	REAL*8 X(N),W(N),H,HN,RE,RF,SUM
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	W(:)=0.0D0
C
	H=X(1)-X(2)
	HN=X(2)-X(3)
	RF=X(1)-X(3)
C
	W(1)=0.5D0*H*X(1)*(X(1)-H/3.0D0)
	1        -H*H/RF/12.0D0*((2.0D0+HN/H)*(1.0D0-H*H/10.0D0)
	1        -X(2)*X(2)+H*H/10.0D0)
	W(2)=0.5D0*H*X(2)*(X(2)+H/3.0D0)
	1        +H*H/RF/12.0D0
	1        *(2.0D0+H/HN+HN/H)*(1.0D0-H*H/10.0D0)
	W(3)=H*H/RF/12.0D0*(H*H/10.0D0-X(2)*X(2)
	1        -H/HN*(1.0D0-H*H/10.0D0))
C
	DO I=3,N-1
	  H=X(I-1)-X(I)
	  RF=X(I-2)-X(I)
	  RE=X(I-1)-X(I+1)
	  W(I-2)=W(I-2)-(X(I-1)*X(I-1)-H*H/10.0D0)*H*H/12.0D0/RF
	  W(I-1)=W(I-1)+0.5D0*H*X(I-1)*(X(I-1)-H/3.0D0)
	1          +(X(I)*X(I)-H*H/10.0D0)*H*H/12.0D0/RE	
	  W(I)=W(I)+0.5D0*H*X(I)*(X(I)+H/3.0D0)
	1          +(X(I-1)*X(I-1)-H*H/10.0D0)*H*H/12.0D0/RF	
	  W(I+1)=W(I+1)-(X(I)*X(I)-H*H/10.0D0)*H*H/12.0D0/RE
	END DO
C
C Assume that dU/dmu=0 at mu=0, and that U depend on MU*MU. Could
C also have this assumption in previous equation as well.
C
	IF(X(N) .EQ. 0)THEN
	  H=X(N-1)-X(N)
	  W(N-1)=W(N-1)+0.5D0*H*(X(N-1)*X(N-1)-X(N-1)*H/3.0D0)
	1        -H/6.0D0*(X(N-1)*X(N-1)-H*H/10.0D0)
	  W(N)=W(N)+H/6.0D0*(X(N-1)*X(N-1)-H*H/10.0D0)
	ELSE
	  H=15.0D0*(X(N-1)-X(N))*(X(N-1)+X(N))
	  W(N-1)=W(N-1)+2.0D0*(X(N-1)**5)/H
	  W(N)=W(N)+(X(N-1)**3)*(3.0D0*X(N-1)*X(N-1)-5.0D0*X(N)*X(N))/H
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - Extrapolation to zero required in KWEIGHT'
	END IF
C
C Ensure that the weights have the correct normalization (shouldnt be 
C necessary.
C
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)
	END DO
	SUM=1.0D0/3.0D0/SUM
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights required normalization in KWEIGHT'
	END IF
	DO I=1,N
	  W(I)=W(I)*SUM
	END DO
C
	RETURN
	END
