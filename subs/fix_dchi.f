C
C Subroutine to reduce the size of the extended two dimensional dCHI 
C matrix used in the computation of dJ. It also used for dETA. The
C integer array GRID gives the positions of the old grid points in the
C new grid. The use of GRID allows for different numbers of grid points
C inserted between adjacent pixels. INDX is used to indicate which depths
C the interpolated opacity is dependant on.
C
C The interpolation may be in the dependent variable (INTERP_TYPE='LIN') or 
C in LOG of the dependent varaiable (INTERP_TYPE='LOG'). The  order of the
C interpolation is assumed to be 3 or less (i.e. CHIEXT at any point depends
C on CHIat 4 points [at most]).
C
C NB J is computed on a grid of size NDEXT. Opacities and polulations are
C known on an original grid of size ND.
C
	SUBROUTINE FIX_dCHI(F2DA,RHS_dHdCHI,CHI,ETA,ND,
	1                   F2DAEXT,RHS_dHdCHIEXT,CHIEXT,ETAEXT,NDEXT,
	1                   COEF,INDX,INTERP_TYPE)
	IMPLICIT NONE
C
C Created  15-May-1997
C
	INTEGER ND,NDEXT
	INTEGER INDX(NDEXT)
C
	REAL*8 F2DA(NDEXT,ND,2)
	REAL*8 RHS_dHdCHI(NDEXT-1,ND)
	REAL*8 CHI(ND),ETA(ND)
C
	REAL*8 F2DAEXT(NDEXT,NDEXT,2)
	REAL*8 RHS_dHdCHIEXT(NDEXT-1,NDEXT)
	REAL*8 CHIEXT(NDEXT),ETAEXT(NDEXT)
	REAL*8 COEF(0:3,NDEXT)
	CHARACTER*(*) INTERP_TYPE
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
	INTEGER I,J,K,L
	INTEGER, PARAMETER :: LRANGE=3
	REAL*8 T1
C
	F2DA(:,:,:)=0.0D0			!NDEXT,ND,2
	RHS_dHdCHI(:,:)=0.0D0			!NDEXT-1,ND
C
	IF(INTERP_TYPE(1:3) .EQ. 'LOG')THEN
C
C NB F2DAEXT(I,K) is only non-zero for K values near I.
C
	  DO L=-LRANGE,LRANGE
	    DO I=MAX(1,1-L),MIN(NDEXT,NDEXT-L)
	      J=I+L
	      K=INDX(J)
	      T1=F2DAEXT(I,J,1)*CHIEXT(J)
	      F2DA(I,K,1)=F2DA(I,K,1)+T1*COEF(0,J)/CHI(K)
	      F2DA(I,K+1,1)=F2DA(I,K+1,1)+T1*COEF(1,J)/CHI(K+1)
	      F2DA(I,K+2,1)=F2DA(I,K+2,1)+T1*COEF(2,J)/CHI(K+2)
	      F2DA(I,K+3,1)=F2DA(I,K+3,1)+T1*COEF(3,J)/CHI(K+3)
	    END DO
	  END DO
C
C Dependence on ETA is diagonal, thus no need to loop over all elements.
C
	  DO I=1,NDEXT
	    J=I
	    K=INDX(J)
	    T1=F2DAEXT(I,J,2)*ETAEXT(J)
	    F2DA(I,K,2)=F2DA(I,K,2)+T1*COEF(0,J)/ETA(K)
	    F2DA(I,K+1,2)=F2DA(I,K+1,2)+T1*COEF(1,J)/ETA(K+1)
	    F2DA(I,K+2,2)=F2DA(I,K+2,2)+T1*COEF(2,J)/ETA(K+2)
	    F2DA(I,K+3,2)=F2DA(I,K+3,2)+T1*COEF(3,J)/ETA(K+3)
	  END DO
C
	  DO L=-LRANGE,LRANGE
	    DO I=MAX(1,1-L),MIN(NDEXT-1,NDEXT-1-L)
	      J=I+L
	      K=INDX(J)
	      T1=RHS_dHdCHIEXT(I,J)*CHIEXT(J)
	      RHS_dHdCHI(I,K)=RHS_dHdCHI(I,K)+T1*COEF(0,J)/CHI(K)
	      RHS_dHdCHI(I,K+1)=RHS_dHdCHI(I,K+1)+T1*COEF(1,J)/CHI(K+1)
	      RHS_dHdCHI(I,K+2)=RHS_dHdCHI(I,K+2)+T1*COEF(2,J)/CHI(K+2)
	      RHS_dHdCHI(I,K+3)=RHS_dHdCHI(I,K+3)+T1*COEF(3,J)/CHI(K+3)
	    END DO
	  END DO
C
	ELSE IF(INTERP_TYPE(1:3) .EQ. 'LIN')THEN
C
	  DO L=-LRANGE,LRANGE
	    DO I=MAX(1,1-L),MIN(NDEXT,NDEXT-L)
	      J=I+L
	      K=INDX(J)
	      T1=F2DAEXT(I,J,1)
	      F2DA(I,K,1)=F2DA(I,K,1)+T1*COEF(0,J)
	      F2DA(I,K+1,1)=F2DA(I,K+1,1)+T1*COEF(1,J)
	      F2DA(I,K+2,1)=F2DA(I,K+2,1)+T1*COEF(2,J)
	      F2DA(I,K+3,1)=F2DA(I,K+3,1)+T1*COEF(3,J)
	    END DO
	  END DO
C
C Dependence on ETA is diagonal, thus no need to loop over all elements.
C
	  DO I=1,NDEXT
	    J=I
	    K=INDX(I)
	    T1=F2DAEXT(I,J,2)
	    F2DA(I,K,2)=F2DA(I,K,2)+T1*COEF(0,J)
	    F2DA(I,K+1,2)=F2DA(I,K+1,2)+T1*COEF(1,J)
	    F2DA(I,K+2,2)=F2DA(I,K+2,2)+T1*COEF(2,J)
	    F2DA(I,K+3,2)=F2DA(I,K+3,2)+T1*COEF(3,J)
	  END DO
C
	  DO L=-LRANGE,LRANGE
	    DO I=MAX(1,1-L),MIN(NDEXT-1,NDEXT-1-L)
	      J=I+L
	      K=INDX(J)
	      T1=RHS_dHdCHIEXT(I,J)
	      RHS_dHdCHI(I,K)=RHS_dHdCHI(I,K)+T1*COEF(0,J)
	      RHS_dHdCHI(I,K+1)=RHS_dHdCHI(I,K+1)+T1*COEF(1,J)
	      RHS_dHdCHI(I,K+2)=RHS_dHdCHI(I,K+2)+T1*COEF(2,J)
	      RHS_dHdCHI(I,K+3)=RHS_dHdCHI(I,K+3)+T1*COEF(3,J)
	    END DO
	  END DO
C
	ELSE
	  I=ERROR_LU()
	  WRITE(I,*)'Error in FIX_dCHI --- invalid interpolation type'
	  STOP
	END IF
C
	RETURN
	END
