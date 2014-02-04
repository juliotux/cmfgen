C
C Subroutine to compute the quadrature weights for a modified cubic rule.
C Reference - A Nordlund (Methods in Radiative Transfer, W Kallofen).
C The routine makes specific assumptions concerning the behaviour of the
C mean intensity at the boundary, and would need to be moified to obtain
C other quantities other than 2nd moment of the intensity. The program assumes
C that the first point corresponds to mu=1.0 .
C
C The first derivatives are estimated by 
C                               dv(i)= (V(I-1)-V(I))/(MU(I-1)-MU(I))
C
	SUBROUTINE NWEIGHT(X,W,N)
	IMPLICIT NONE
C
C Altered 25-May-1996 - CALL to DP_ZERO deleted
C                       ERROR_LU inserted.     
C Created 11-JAN-1988 - Based on KWEIGHT.
C
	INTEGER N,I
	REAL*8 X(N),W(N),H,R,RE,RF,T1,SUM
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	W(:)=0.0D0
C
C Firstly evaluate weights for integration from I=1 to I=2.
C R, RF and T1 are used in calculating the corection due
C the derivative of V at d=1.
C
	H=X(1)-X(2)
	R=X(2)-X(3)
C
	RF=( (X(1)**3)-0.3D0*X(1)*H*H ) * H*(R+2*H)/(R+H)/12.0D0
	T1=( (X(1)**3)-0.3D0*X(1)*H*H ) * (H**3)/R/(R+H)/12.0D0
	RE=H*H*( (X(2)**3)-0.3D0*H*H*X(2) )/(X(1)-X(3))/12.0D0
C
	W(1)=W(1)+0.5D0*H*(  (X(1)**3) -0.5D0*H*(X(1)**2)
	1         +(H**3)/60.0D0  )+RE-RF
	W(2)=W(2)+0.5D0*H*(  (X(2)**3) +0.5D0*H*(X(2)**2)
	1         -(H**3)/60.0D0  )+RF+T1
	W(3)=W(3)-RE-T1
C
	DO I=3,N-1
	  H=X(I-1)-X(I)
	  RF=H*H*(X(I-1)**3-0.3D0*H*H*X(I-1))/(X(I-2)-X(I))/12.0D0
	  RE=H*H*(X(I)**3-0.3D0*H*H*X(I))/(X(I-1)-X(I+1))/12.0D0
	  W(I-2)=W(I-2)-RF
	  W(I-1)=W(I-1)+0.5D0*H*(  (X(I-1)**3) -0.5D0*H*(X(I-1)**2)
	1         +(H**3)/60.0D0  )+RE
	  W(I)=W(I)+0.5D0*H*(  (X(I)**3) +0.5D0*H*(X(I)**2)
	1         -(H**3)/60.0D0  )+RF
	  W(I+1)=W(I+1)-RE
	END DO
C
C Assume that V(mu)=0 at mu=0, and that it depends linearly on mu.
C When X(N) .ne. 0 we assume that v'(n)=v(n)/mu*
C
	IF(X(N) .EQ. 0)THEN
	  W(N-1)=W(N-1)+0.2D0*(X(N-1)**4)
	ELSE
	  H=X(N-1)-X(N)
	  RF=H*H*(X(N-1)**3-0.3D0*H*H*X(N-1))/(X(N-2)-X(N))/12.0D0
	  RE=H*H*(X(N)**3-0.3D0*H*H*X(N))/X(N)/12.0D0
	  W(N-2)=W(N-2)-RF
	  W(N-1)=W(N-1)+0.5D0*H*(  (X(N-1)**3) -0.5D0*H*(X(N-1)**2)
	1         +(H**3)/60.0D0  )
	  W(N)=W(N)+0.5D0*H*(  (X(N)**3) +0.5D0*H*(X(N)**2)
	1         -(H**3)/60.0D0  )+RF+RE
C
C The above was for the integral from n-1 to n. Now need to extend the
C integral to MU=0.
C
	  W(N)=W(N)+0.2D0*(X(N)**4)
	END IF
C
C Ensure that the weights have the correct normalization (shouldnt be 
C necessary. We assume v=mu for the normalization.
C
	SUM=0.0D0
	DO I=1,N
	  SUM=SUM+W(I)*X(I)
	END DO
	SUM=0.2D0/SUM
	IF(ABS(SUM-1.0D0) .GT. 1.0D-12)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)' Warning - weights required normalization in NWEIGHT'
	  WRITE(LUER,*)'SUM=',SUM
	END IF
	DO I=1,N
	  W(I)=W(I)*SUM
	END DO
C
	RETURN
	END
