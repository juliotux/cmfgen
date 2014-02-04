C
C Routine to interpolate V, T and SIGMA onto a new radius
C
	SUBROUTINE EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NX,
	1                     V,T,SIGMA,ND)
C
	IMPLICIT NONE
	INTEGER NX,ND,INDX(NX)
	REAL*8 VEXT(NX),SIGMAEXT(NX),TEXT(NX),COEF(0:3,NX)
	REAL*8 V(ND),SIGMA(ND),T(ND)
C
	INTEGER I,J
C
C Do intepolation in the log plane - need to add 1 to Sigma as it can be
C negative.
C
	DO I=1,NX
	  VEXT(I)=0.0D0
	  SIGMAEXT(I)=0.0D0
	  TEXT(I)=0.0D0
	  DO J=0,3
	    VEXT(I)=VEXT(I)+COEF(J,I)*LOG( V(J+INDX(I)) )
	    TEXT(I)=TEXT(I)+COEF(J,I)*LOG( T(J+INDX(I)) )
	    SIGMAEXT(I)=SIGMAEXT(I)+COEF(J,I)*
	1                  LOG( SIGMA(J+INDX(I))+1.0D0)
	  END DO
	  VEXT(I)=EXP(VEXT(I))
	  SIGMAEXT(I)=EXP(SIGMAEXT(I))-1.0D0
	  TEXT(I)=EXP(TEXT(I))
	END DO
C
	RETURN	
	END
