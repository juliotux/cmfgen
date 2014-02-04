C
C The function EXPONX is given by (1.0-EXP(-X))/X.
C This functions is called to allow for cancellation when X is small.
C
	FUNCTION  EXPONX(X)
	IMPLICIT NONE
C
C ALtered 24-May-1996 :  DOUBLE PRECISION replaced by REAL*8
C Altered 26-May-1988 : Exponential no longer computed if X is greater than
C                        40. Necessary to overcome a bug with Dec software.

	REAL*8 EXPONX,X
C
	IF( ABS(X) .LT. 1.0D-03 )THEN
	  EXPONX=1.0D0-X*(  0.5D0-X/6.0D0*( 1.0D0-X/4.0D0 )  )
	ELSE IF(X .LT. 40)THEN
	  EXPONX=( 1.0D0-EXP(-X) )/X
	ELSE
	  EXPONX=1.0D0/X
	END IF
C
	RETURN
	END


C
C The function d_EXPONX_dX is given by d[ (1.0-EXP(-X))/X ]/dX. This
C function is called to allow for cancellation when X is small.
C
C Altered 26-May-88 - Exponential no longer computed if X is greater than
C                     40. Necessary to overcome a bug with Dec software.
C
	FUNCTION  d_EXPONX_dX(X)
	IMPLICIT NONE
	DOUBLE PRECISION d_EXPONX_dX,X,Y
C
	IF( ABS(X) .LT. 1.0D-03 )THEN
	  d_EXPONX_dX=-0.5D0+X*( 1.0D0-X*(0.375D0-0.1D0*X) )/3.0D0
	ELSE
	  IF(X .LT. 40)THEN
	    Y=EXP(-X)
	  ELSE
	    Y=0.0D0
	  END IF
	  d_EXPONX_dX=( Y-(1.0D0-Y)/X )/X
	END IF
C
	RETURN
	END
