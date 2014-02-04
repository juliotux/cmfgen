C
C Subroutine to interpolate an array onto a new grid. The grid vector must be 
C either a monotonically decreasing or increasing function. A modified cubic
C polynomial is used to do the interpolation. Instead of using
C the excact cubic estiamtes for the first derivative at the two nodes,
C we use revised estimates which insure that the interpolating function
C is mononotonic in the interpolating interval.
C
C The techniques is somewhat similar to that suggested by Nordulund.
C
C Disadvantages: The interpolating weights can only be defined when the
C                function is known. In principal could use these modified
C                first derivatives to compute an accurate integration
C                formulae. However, the integration weights cannot be defined
C                independently of the function values, as desired in many
C                situations.
C
C Ref: Steffen. M, 1990, A/&A, 239, 443-450
C
	SUBROUTINE MON_INTERP(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 : ERROR_LU installed
C Created 01-Apr-1992 : Code may need recoding for optimal speed, and for
C                         vectorization.
C
	INTEGER NQ,LIN_END,NX,NV,ND
	REAL*8 QZ(NQ,LIN_END),QZR(NX)
	REAL*8 VARRAY(NV,LIN_END),R(ND)
C
	REAL*8 ONE
	PARAMETER (ONE=1.0D0)
	INTEGER I,J,M
	REAL*8 T1
	REAL*8 HI,HIM1,HIP1
	REAL*8 SI,SIM1,SIP1
	REAL*8 A,B,C,D,DYI,DYIP1,SGN
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C The array R may be either monotonically increasing, or decreasing.
C
	SGN=SIGN(ONE,R(ND)-R(1))
	IF( (SGN*QZR(1) .LT. SGN*R(1)) .OR.
	1   (SGN*QZR(NX) .GT. SGN*R(ND)) )THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in MON_INTERP - values outside range'
	  WRITE(LUER,'(1X,2(A,ES20.12),A,I4)')'  R(1)=',R(1),'    R(ND)= ',R(ND),'ND=',ND
	  WRITE(LUER,'(1X,2(A,ES20.12),A,I4)')'QZR(1)=',QZR(1),'  QZR(NX)=',QZR(NX),'NX=',NX
	  STOP
	END IF
	I=1
C
C M is the Index in new interpolated array
C
	DO M=1,NX
500	  IF( SGN*QZR(M) .LE. SGN*R(I+1))THEN
	    IF(I .EQ. 1)THEN
	      HI=R(2)-R(1)
              HIP1=R(3)-R(2)
              DO J=1,LIN_END
                SI=(VARRAY(2,J)-VARRAY(1,J))/HI
                SIP1=(VARRAY(3,J)-VARRAY(2,J))/HIP1
                DYI=SI +(SI-SIP1)*HI/(HI+HIP1)
                DYIP1=(SI*HIP1+SIP1*HI)/(HI+HIP1)
	        DYI=( SIGN(ONE,SI)+SIGN(ONE,DYI) )*
	1            MIN(ABS(SI),0.5D0*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,SIP1) )*
	1            MIN(ABS(SI),ABS(SIP1),0.5D0*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0D0*SI)/HI/HI
	        B=(3.0D0*SI-2.0D0*DYI-DYIP1)/HI
	        C=DYI
	        D=VARRAY(I,J)
                QZ(M,J)=((A*T1+B)*T1+C)*T1+D
	      END DO
	    ELSE IF(I .EQ. ND-1)THEN
	      HI=R(ND)-R(ND-1)
              HIM1=R(ND-1)-R(ND-2)
              DO J=1,LIN_END
                SIM1=(VARRAY(ND-1,J)-VARRAY(ND-2,J))/HIM1
                SI=(VARRAY(ND,J)-VARRAY(ND-1,J))/HI
                DYI=(SIM1*HI+SI*HIM1)/(HIM1+HI)
                DYIP1=SI+(SI-SIM1)*HI/(HIM1+HI)
	        DYI=( SIGN(ONE,SIM1)+SIGN(ONE,SI) )*
	1            MIN(ABS(SIM1),ABS(SI),0.5D0*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,DYIP1) )*
	1            MIN(ABS(SI),0.5D0*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0D0*SI)/HI/HI
	        B=(3.0D0*SI-2.0D0*DYI-DYIP1)/HI
	        C=DYI
	        D=VARRAY(I,J)
                QZ(M,J)=((A*T1+B)*T1+C)*T1+D
	      END DO
	    ELSE
	      HI=R(I+1)-R(I)
              HIM1=R(I)-R(I-1)
              HIP1=R(I+2)-R(I+1)
              DO J=1,LIN_END
                SIM1=(VARRAY(I,J)-VARRAY(I-1,J))/HIM1
                SI=(VARRAY(I+1,J)-VARRAY(I,J))/HI
                SIP1=(VARRAY(I+2,J)-VARRAY(I+1,J))/HIP1
                DYI=(SIM1*HI+SI*HIM1)/(HIM1+HI)
                DYIP1=(SI*HIP1+SIP1*HI)/(HI+HIP1)
	        DYI=( SIGN(ONE,SIM1)+SIGN(ONE,SI) )*
	1            MIN(ABS(SIM1),ABS(SI),0.5D0*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,SIP1) )*
	1            MIN(ABS(SI),ABS(SIP1),0.5D0*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0D0*SI)/HI/HI
	        B=(3.0D0*SI-2.0D0*DYI-DYIP1)/HI
	        C=DYI
	        D=VARRAY(I,J)
                QZ(M,J)=((A*T1+B)*T1+C)*T1+D
	      END DO
	    END IF
	  ELSE
	    I=I+1
	    GOTO 500
	  END IF
C
	END DO
C
	RETURN
	END
