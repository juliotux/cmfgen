C
C Subroutine to reduce the size of the extended two dimensional dCHI matrix
C used in the computation of dJ. May also be used for dETA. 
C The integer array GRID gives the positions of the old grid points in
C the new grid. The use of GRID allows for different numbers of grid points
C inserted between adjacet pixels. INDX is used to indicate which depths
C the interpolated opacity is dependant on.
C
	SUBROUTINE REGRID_dCHI(F2DA,CHI,ND,GRID,
	1                          F2DAEXT,CHIEXT,NDEXT,COEF,INDX)
	IMPLICIT NONE
C
C Altered 26-May-1996 - DOUBLE PRECISION removed.
C                       Call to DP_ZERO removed.
C
C Altered  4-May-1988 - Bug fix - was not incrememting F2DA(I,K+3) correctly.
C Altered 26-Apr-1988 - Bug fix. F2DA now zeroed.
C Created  5-Apr-1988
C
	INTEGER ND,NDEXT,INDX(NDEXT),GRID(ND)
	REAL*8 F2DA(ND,ND),CHI(ND),T1
	REAL*8 F2DAEXT(NDEXT,NDEXT)
	REAL*8 CHIEXT(NDEXT),COEF(0:3,NDEXT)
C
	INTEGER I,J,K
C
	F2DA(:,:)=0.0D0
	DO J=1,NDEXT
	  K=INDX(J)
	  DO I=1,ND
	    T1=F2DAEXT(GRID(I),J)*CHIEXT(J)
	    F2DA(I,K)=F2DA(I,K)+T1*COEF(0,J)/CHI(K)
	    F2DA(I,K+1)=F2DA(I,K+1)+T1*COEF(1,J)/CHI(K+1)
	    F2DA(I,K+2)=F2DA(I,K+2)+T1*COEF(2,J)/CHI(K+2)
	    F2DA(I,K+3)=F2DA(I,K+3)+T1*COEF(3,J)/CHI(K+3)
	  END DO
	END DO
C
	RETURN
	END
