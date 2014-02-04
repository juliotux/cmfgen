!
! Subroutine to perform a simle arithmetic opperation on two arrays.
! The operation is perfomed only on that section contained in both arrays.
! Arrays may be in any numerical order.
!
	  SUBROUTINE DO_VEC_OP(IN1,IN2,OUT,LIN_INT,OPERATION)
	  USE MOD_CURVE_DATA
	  IMPLICIT NONE
!
! Altered 15-May-2002: Bug fixed for unequally spaced data.
!
	  INTEGER IN1
	  INTEGER IN2
	  INTEGER OUT
	  LOGICAL LIN_INT
	  CHARACTER(LEN=*) OPERATION
!         
	  EXTERNAL SP_EQUAL
	  LOGICAL SP_EQUAL
!
! Local variables
!
	  REAL*4, ALLOCATABLE :: XV(:)
	  REAL*4, ALLOCATABLE :: YV(:)
	  REAL*4, ALLOCATABLE :: ZV(:)
!
	  REAL*4 T1
	  LOGICAL UNEQUAL
	  INTEGER I,J,L
	  INTEGER SIGN
	  INTEGER N1,N2
	  INTEGER N
	  INTEGER IL,IU
	  INTEGER LOW_LIM,UP_LIM
!
	  INTEGER, PARAMETER :: IONE=1
	  INTEGER, PARAMETER :: T_OUT=6
!
! Check validity of maps:
!
	  IF(NPTS(IN1) .EQ. 0 .OR. NPTS(IN2) .EQ. 0)THEN
	    WRITE(T_OUT,*)'Invalid Plot ID''s in DO_VEC_OP'
	    WRITE(T_OUT,*)'Plot ID=',IN1,'  NPTS=',NPTS(IN1)
	    WRITE(T_OUT,*)'Plot ID=',IN2,'  NPTS=',NPTS(IN2)
	    RETURN
	  END IF
!
	  T1=1.0D-07
	  I=1
	  UNEQUAL=.FALSE.
	  IF(NPTS(IN2) .NE. NPTS(IN1))UNEQUAL=.TRUE.
	  DO WHILE(.NOT. UNEQUAL .AND. I .LE. NPTS(IN2))
	    IF( SP_EQUAL(CD(IN2)%XVEC(I),CD(IN1)%XVEC(I),T1) )THEN
	      I=I+1
	    ELSE
	      UNEQUAL=.TRUE.
	    END IF
	  END DO
!
	  IF(ALLOCATED(YV))DEALLOCATE(YV)
	  ALLOCATE (YV(NPTS(IN1)),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error in DO_VEC_OP --- unable to allocate YV'
	    WRITE(T_OUT,*)'IOS=',IOS
	    RETURN
	  END IF
!
! We now use a conventional sign convention. SGN is positive if
! X(2) > X(1) etc. If we multiply XVEC by SIGN, we reverse its
! order.
!
	  SIGN=1
	  IF(CD(IN1)%XVEC(1) .GT. CD(IN1)%XVEC(NPTS(IN1)))SIGN=-1
	  N1=NPTS(IN1)
	  N2=NPTS(IN2)
!
	  IF(.NOT. UNEQUAL)THEN
	    LOW_LIM=1
	    UP_LIM=NPTS(IN1)
	    YV(:)=CD(IN2)%DATA(:)
	  ELSE IF(UNEQUAL .AND. LIN_INT)THEN
!
	    IF(ALLOCATED(XV))DEALLOCATE(XV)
	    IF(ALLOCATED(ZV))DEALLOCATE(ZV)
	    ALLOCATE (XV(MAX(N1,N2)),STAT=IOS)
	    IF(IOS .EQ. 0)ALLOCATE (ZV(MAX(N1,N2)),STAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error in DO_VEC_OP --- unable to allocate XV,ZV'
	      WRITE(T_OUT,*)'IOS=',IOS
	      RETURN
	    END IF
!
! Check to see whether we need to flip the array.
!
	    T1= (CD(IN2)%XVEC(N2)-CD(IN2)%XVEC(1))*(CD(IN1)%XVEC(N1)-CD(IN1)%XVEC(1))
	    IF(T1 .LT. 0)THEN
	      DO I=1,N2
	        XV(I)=CD(IN2)%XVEC(N2-I+1)
	        ZV(I)=CD(IN2)%DATA(N2-I+1)
	      END DO
	    ELSE
	      DO I=1,N2
	        XV(I)=CD(IN2)%XVEC(I)
	        ZV(I)=CD(IN2)%DATA(I)
	      END DO
	    END IF
!
	    L=1
	    LOW_LIM=1
	    UP_LIM=NPTS(IN1)
	    DO I=1,NPTS(IN1)
	      IF(SIGN*CD(IN1)%XVEC(I) .LT. SIGN*XV(1))THEN
	        YV(I)=0.0
	        LOW_LIM=I+1
	      ELSE IF(SIGN*CD(IN1)%XVEC(I) .GT. SIGN*XV(N2))THEN
	        YV(I)=0.0
	        UP_LIM=I-1
	        EXIT
	      ELSE 
	        DO WHILE (SIGN*CD(IN1)%XVEC(I) .GT. SIGN*XV(L+1))
	          L=L+1           
	        END DO
	        T1=(CD(IN1)%XVEC(I)-XV(L+1))/(XV(L)-XV(L+1))
	        YV(I)=(1.0D0-T1)*ZV(L+1)+T1*ZV(L)
	      END IF
	    END DO
	  ELSE IF(UNEQUAL)THEN
!
! Check to see whether we need to flip the array.
!
	    T1= (CD(IN2)%XVEC(N2)-CD(IN2)%XVEC(1))*(CD(IN1)%XVEC(N1)-CD(IN1)%XVEC(1))
	    IF(T1 .LT. 0)THEN
	      DO I=1,N2
	        XV(I)=CD(IN2)%XVEC(N2-I+1)
	        ZV(I)=CD(IN2)%DATA(N2-I+1)
	      END DO
	    ELSE
	      DO I=1,N2
	        XV(I)=CD(IN2)%XVEC(I)
	        ZV(I)=CD(IN2)%DATA(I)
	      END DO
	    END IF
!
C
C We will use monotonic cubic interpolation. We first verify the range.
C I & J are temporary variables for the callt o MON_INTERP. I denotes the 
C first element. Initially J denotes the last element, then the numer of
C elements that can be interpolated.
C
	    I=1
	    DO WHILE(SIGN*CD(IN1)%XVEC(I) .LT. SIGN*XV(1))
	      I=I+1
	    END DO
	    J=NPTS(IN1)
	    DO WHILE(SIGN*CD(IN1)%XVEC(J) .GE. SIGN*XV(N2))
	      J=J-1
	    END DO
	    J=J-I+1
C
	    YV(1:NPTS(IN1))=0.0D0
	    CALL SP_MON_INTERP(YV(I),J,IONE,CD(IN1)%XVEC(I),J,ZV,N2,XV,N2)
	    LOW_LIM=I
	    UP_LIM=I+J-1
	  END IF
!
	  IL=lOW_LIM
	  IU=UP_LIM
	  IF(OPERATION .EQ. '*')THEN
	    YV(IL:IU)=CD(IN1)%DATA(IL:IU)*YV(IL:IU)
	  ELSE IF(OPERATION .EQ. '+')THEN
	    YV(IL:IU)=CD(IN1)%DATA(IL:IU)+YV(IL:IU)
	  ELSE IF(OPERATION .EQ. '-')THEN
	    YV(IL:IU)=CD(IN1)%DATA(IL:IU)-YV(IL:IU)
	  ELSE IF(OPERATION .EQ. '/')THEN
!
! The test for overflow needs improving.
!
	    T1=1.0E-12*HUGE(YV(1))	! T1=1.0E-02*HUGE(YV(1))
	    DO J=IL,IU
	      IF( YV(J) .EQ. 0.0D0)THEN
	        YV(J)=T1
	      ELSE IF( ABS(YV(J)) .GT. ABS(CD(IN1)%DATA(J)) )THEN
	        YV(J)=CD(IN1)%DATA(J)/YV(J)
	      ELSE IF( LOG10(ABS(CD(IN1)%DATA(J)))-LOG10(ABS(YV(J))) .LT. LOG(T1) )THEN
	        YV(J)=CD(IN1)%DATA(J)/YV(J)
	      ELSE
	        YV(J)=T1
	      END IF
	    END DO
	  ELSE
	    WRITE(T_OUT,*)'Invalid operation in DO_VEC_OP'
	    RETURN
	  END IF
!
	IF(.NOT. ALLOCATED(XV))ALLOCATE (XV(1:IU))
	XV(IL:IU)=CD(IN1)%XVEC(IL:IU)
	IF(ALLOCATED(CD(OUT)%XVEC))THEN
	  DEALLOCATE (CD(OUT)%XVEC)
	  DEALLOCATE (CD(OUT)%DATA)
	END IF
	N=IU-IL+1
	ALLOCATE (CD(OUT)%XVEC(N),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CD(OUT)%DATA(N),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error in DO_VEC_OP --- unable to allocate',
	1                ' plot storage vectors'
	  WRITE(T_OUT,*)'IOS=',IOS
	  RETURN
	END IF
!
	CD(OUT)%XVEC(1:N)=XV(IL:IU)
	CD(OUT)%DATA(1:N)=YV(IL:IU)
	NPTS(OUT)=N
	ERR(OUT)=.FALSE.
	IF(OUT .GT. NPLTS)NPLTS=OUT
	IF(ALLOCATED(XV))DEALLOCATE(XV)
	IF(ALLOCATED(YV))DEALLOCATE(YV)
	IF(ALLOCATED(ZV))DEALLOCATE(ZV)
!
	RETURN
	END 
!
!
C
C Logical function to determine whether two values are equal to within
C 100Z % . If one of the arguments are zero, EQUAL is set false unless
C both are equal to zero in which case it is set true. Neither X, Y or
C Z are altered.
C
	FUNCTION SP_EQUAL(X,Y,Z)
	IMPLICIT NONE
C
C Altered 24-May-1996 - File now contains DP version only (i.e. not SP_EQUAL)
C Altered 19-Jul-1991 - Rearrangement of LOGICAL descriptor for CRAY.
C Altered 14-Apr-1989 - Now divide by the larger (absolute) of X and Y.
C                       This routine should never give a floating point
C                       overflow.
C Altered  4-NOV-86 (Bug for X or Y=0 fixed).
C
	LOGICAL SP_EQUAL
	REAL*4 X,Y,Z
C
	SP_EQUAL=.FALSE.
	IF(X .EQ. 0.0D0 .AND. Y .EQ. 0.0D0)THEN
	  SP_EQUAL=.TRUE.
	ELSE IF( ABS(Y) .GT. ABS(X) )THEN
	  IF( (1.0D0-X/Y) .LE. Z )SP_EQUAL=.TRUE.
	ELSE
	  IF( (1.0D0-Y/X) .LE. Z )SP_EQUAL=.TRUE.
	END IF
C
	RETURN
	END
C
C!
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
	SUBROUTINE SP_MON_INTERP(QZ,NQ,LIN_END,QZR,NX,VARRAY,NV,R,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 : ERROR_LU installed
C Created 01-Apr-1992 : Code may need recoding for optimal speed, and for
C                         vectorization.
C
	INTEGER NQ,LIN_END,NX,NV,ND
	REAL*4 QZ(NQ,LIN_END),QZR(NX)
	REAL*4 VARRAY(NV,LIN_END),R(ND)
C
	REAL*4 ONE
	PARAMETER (ONE=1.0D0)
	INTEGER I,J,M
	REAL*4 T1
	REAL*4 HI,HIM1,HIP1
	REAL*4 SI,SIM1,SIP1
	REAL*4 A,B,C,D,DYI,DYIP1,SGN
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
	1            MIN(ABS(SI),0.5*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,SIP1) )*
	1            MIN(ABS(SI),ABS(SIP1),0.5*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0*SI)/HI/HI
	        B=(3.0*SI-2.0*DYI-DYIP1)/HI
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
	1            MIN(ABS(SIM1),ABS(SI),0.5*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,DYIP1) )*
	1            MIN(ABS(SI),0.5*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0*SI)/HI/HI
	        B=(3.0*SI-2.0*DYI-DYIP1)/HI
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
	1            MIN(ABS(SIM1),ABS(SI),0.5*ABS(DYI))
	        DYIP1=( SIGN(ONE,SI)+SIGN(ONE,SIP1) )*
	1            MIN(ABS(SI),ABS(SIP1),0.5*ABS(DYIP1))
	        T1=(QZR(M)-R(I))
                A=(DYI+DYIP1-2.0*SI)/HI/HI
	        B=(3.0*SI-2.0*DYI-DYIP1)/HI
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
