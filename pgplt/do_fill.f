!
! Subroutine to perform a simle arithmetic opperation on two arrays.
! The operation is perfomed only on that section contained in both arrays.
! Arrays may be in any numerical order.
!
	  SUBROUTINE DO_FILL(XPAR,YPAR,IN1,IN2,LIN_INT)
	  USE MOD_CURVE_DATA
	  IMPLICIT NONE
!
! Altered 15-May-2002: Bug fixed for unequally spaced data.
!
	  REAL*4  XPAR(2),YPAR(2)
	  INTEGER IN1,IN2
	  LOGICAL LIN_INT
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
	  REAL*4, ALLOCATABLE :: XP(:)
	  REAL*4, ALLOCATABLE :: YP(:)
	  REAL*4, ALLOCATABLE :: ZP(:)
!
	  REAL*4 T1
	  LOGICAL SWITCH
	  LOGICAL UNEQUAL
	  INTEGER I,J,K,L
	  INTEGER SIGN
	  INTEGER N1,N2
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
	  N1=NPTS(IN1)
	  N2=NPTS(IN2)
	  IF(N2 .GT. N1)THEN
	    I=IN1; IN1=IN2; IN2=I
	    I=N1;  N1=N2;   N2=N1
	  END IF
!
	  ALLOCATE (XV(N1),YV(N1),ZV(N1))
	  ALLOCATE (XP(2*N1),YP(2*N1),ZP(2*N1))
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
	  WRITE(6,*)CD(IN1)%XVEC(LOW_LIM),XPAR(1),XPAR(2)
	  WRITE(6,*)CD(IN1)%XVEC(UP_LIM),XPAR(1),XPAR(2)
	  DO WHILE(CD(IN1)%XVEC(LOW_LIM) .LT. XPAR(1))
	    LOW_LIM=LOW_LIM+1
	  END DO
	  DO WHILE(CD(IN1)%XVEC(UP_LIM) .GT. XPAR(2))
	    UP_LIM=UP_LIM-1
	  END DO
	  WRITE(6,*)CD(IN1)%XVEC(LOW_LIM),XPAR(1),XPAR(2)
	  WRITE(6,*)CD(IN1)%XVEC(UP_LIM),XPAR(1),XPAR(2)
!
	  IL=lOW_LIM
	  IU=UP_LIM
	  K=0
	  T1=1.0D-02
	  DO I=IL,IU
	    SWITCH=.FALSE.
	    IF(K .GT. 2)THEN
	      IF( (YP(K)-ZP(K))*(YP(K-1)-ZP(K-1)) .LT. 0)SWITCH=.TRUE.
	    END IF
	    IF( SP_EQUAL(CD(IN1)%DATA(I),YV(I),T1) .OR. SWITCH .OR. I .EQ. IU)THEN
	      K=K+1
	      XP(K)=CD(IN1)%XVEC(I)
	      YP(K)=CD(IN1)%DATA(I)
	      ZP(K)=YP(K)
	      IF(K .GT. 2)THEN
	         DO J=1,K
	           XP(2*K-J+1)=XP(J)
	           YP(2*K-J+1)=ZP(J)
	         END DO
	         K=2*K+1
	         XP(K)=XP(1)
	         YP(K)=ZP(1)
	         CALL PGPOLY(K,XP,YP)
	      END IF
	      K=0
	    ELSE
	      K=K+1
	      XP(K)=CD(IN1)%XVEC(I)
	      YP(K)=CD(IN1)%DATA(I)
	      ZP(K)=YV(I)
	    END IF
	  END DO
!
	DEALLOCATE (XV,YV,ZV)
	DEALLOCATE (XP,YP,ZP)
!
	RETURN
	END
