
C
C Subroutine to compute the increment in optical depth along a ray with impact 
C parameter P. The Euler-Mclaurian summation formula is used to provide an 
C approximate correction to the trapazoidal rule using the first derivatives
C supplied in the call.
C
C NOTE: 
C    [1] If the routine is call with dCHIdR=0, the trapazoidal integration
C          rule is recovered.
C                   
	SUBROUTINE  NORDTAU(DTAU,CHI,Z,R,dCHIdR,NI)
	IMPLICIT NONE
C
C Altered 25-NOV-1986 (Based on TAU)
C
	INTEGER NI,I
	REAL*8 DTAU(NI),CHI(NI),Z(NI),R(NI),dCHIdR(NI)
C
	DO I=1,NI-1
	  DTAU(I)=0.5D0*(Z(I)-Z(I+1))*(CHI(I)+CHI(I+1)+(Z(I)-Z(I+1))
	1   *( dCHIdR(I+1)*Z(I+1)/R(I+1)-dCHIdR(I)*Z(I)/R(I) )/6.0D0)
	END DO
C
	RETURN
	END
C 
C
C Compute the first derivative of the opacity with radius. For use with NORDTAU
C although they could also be used to define a cubic interpolation function. 
C
C The first derivatives are obtained via several methods, as indicated by the
C variable METHOD.
C
C METHOD =  LOGLOG: Derivative in LOG-LOG plane is obtained using the secant
C                      connected to the 2 adjacent points (from Nordulund).
C           LOGLIN: As for LOGLOG, but dependent variable (CHI) is in LOG Plane
C                      while R is in the linear plane.
C           LINEAR: As for LOGLOG, but both the dependent Variable (CHI) and 
C                      the independent variable R are in the linear plane.
C
C           LOGMON: Derivatives at each node are derived by fitting a quadratic
C                      polynomial to it and its adjacent neighbors in the
C                      LOG-LOG plane. These derivatives are converted to the
C                      linear-linear plane, and then adjusted so the the cubic
c                      defined by them is monotonic. When used with NORDTAU it
C                      guarantees that the optical depth increments are always
C                      positive (if CHI is positive).
C	    LINMON: Derivatives chosen so that fitting cure is monotonic.
C                   The fitting is performed in the Linear-Linear plane.
C
C NOTE:
C    [1] Originally, it was assumed that R(1)-R(2) << R(2)-R(3) and 
C R(ND-1)-R(ND) << R(ND-1)-R(ND-2) when computing the derivative at the two 
C positions closet to each boundary. Program now determines whether this is 
C the case.
C
	SUBROUTINE DERIVCHI(dCHIdr,CHI,R,ND,METHOD)
	IMPLICIT NONE
C
C Altered 02-Mar-1999 : In LOGMON the derivatives are now adjusted in the
C                         linear-linear plane to ensure monotocity of the
C                         fitting cubic. 
C Altered 15-Aug-1996 : LINLOG and LINMOM installed in an effort to provided
C                        an improved integration rule that would not produce
C                        negative optical depth increments.
C ALtered 24-May-1996 ; ERROR_LU installed.
C ALtered 25-JAn-88 - Boundary condition checked.
C Altered 20-Feb-1987 (Method option installed)
C
	INTEGER ND,I
	REAL*8 CHI(ND),dCHIdR(ND),R(ND),LIM
	CHARACTER*6 METHOD
C                        
	REAL*8 H(ND),SLOPE(ND)
	REAL*8, PARAMETER :: ONE=1.0D0
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	LIM=3.0D0
C
	IF(METHOD .EQ. 'LOGLOG')THEN
	  DO I=2,ND-1
	    dCHIdR(I)=LOG(CHI(I-1)/CHI(I+1)) /
	1               LOG(R(I-1)/R(I+1)) *CHI(I)/R(I)
	  END DO
	  dCHIdR(1)=LOG(CHI(1)/CHI(2)) / LOG(R(1)/R(2))
	  IF( (R(1)-R(3))/(R(1)-R(2)) .GT. LIM )
	1                  	dCHIdR(2)=dCHIdR(1)*CHI(2)/R(2)
	  dCHIdR(1)=dCHIdR(1)*CHI(1)/R(1)
	  dCHIdR(ND)=LOG(CHI(ND-1)/CHI(ND)) / LOG(R(ND-1)/R(ND))
	  IF( (R(ND-2)-R(ND))/(R(ND-1)-R(ND)) .GT. LIM )
	1		   	dCHIdR(ND-1)=dCHIdR(ND)*CHI(ND-1)/R(ND-1)
	  dCHIdR(ND)=dCHIdR(ND)*CHI(ND)/R(ND)
!
!
	ELSE IF(METHOD .EQ. 'LOGLIN')THEN
	  DO I=2,ND-1
	    dCHIdR(I)=(LOG(CHI(I-1))-LOG(CHI(I+1)))
	1             /(R(I-1)-R(I+1))*CHI(I)
	  END DO
	  dCHIdR(1)=(LOG(CHI(1))-LOG(CHI(2)))/(R(1)-R(2))
	  IF( (R(1)-R(3))/(R(1)-R(2)) .GT. LIM )
	1			dCHIdR(2)=dCHIdR(1)*CHI(2)
	  dCHIdR(1)=dCHIdR(1)*CHI(1)
	  dCHIdR(ND)=(LOG(CHI(ND-1))-LOG(CHI(ND)))/(R(ND-1)-R(ND))
	  IF( (R(ND-2)-R(ND))/(R(ND-1)-R(ND)) .GT. LIM )
	1			dCHIdR(ND-1)=dCHIdR(ND)*CHI(ND-1)
	  dCHIdR(ND)=dCHIdR(ND)*CHI(ND)
!
!
	ELSE IF(METHOD .EQ. 'LINEAR')THEN
	  DO I=2,ND-1
	    dCHIdR(I)=(CHI(I-1)-CHI(I+1))/(R(I-1)-R(I+1))
	  END DO
	  dCHIdR(1)=(CHI(1)-CHI(2))/(R(1)-R(2))
	  IF( (R(1)-R(3))/(R(1)-R(2)) .GT. LIM )dCHIdR(2)=dCHIdR(1)
	  dCHIdR(ND)=(CHI(ND-1)-CHI(ND))/(R(ND-1)-R(ND))
	  IF( (R(ND-2)-R(ND))/(R(ND-1)-R(ND)) .GT. LIM )
	1			dCHIdR(ND-1)=dCHIdR(ND)
!
!
	ELSE IF(METHOD .EQ. 'LINMON')THEN
C
C In this method a cubic is fitted between adjacent points. The first 
C derivative at each node is estimated from the parabola passing through that 
C node and the adjacent points. The derivatives are then adjusted to ensure 
C that the curve in every interval is monotonic. 
C (After Steffen, 1990, A&A, 239, 443-450).
C
C Compute the interval between grid points, and the slope in each interval
C They are computed as vectors as the values are used in 2 intervals.
C It also allows vectorization.
C
	  DO I=1,ND-1
	    H(I)=R(I+1)-R(I)
	    SLOPE(I)=(CHI(I+1)-CHI(I))/H(I)
	  END DO
C
          dCHIdR(1)=SLOPE(1) +(SLOPE(1)-SLOPE(2))*H(1)/(H(1)+H(2))
	  dCHIdR(1)=( SIGN(ONE,SLOPE(1))+SIGN(ONE,dCHIdR(1)) )*
	1            MIN(ABS(SLOPE(1)),0.5D0*ABS(dCHIdR(1)))
	  DO I=2,ND-1
            dCHIdR(I)=(SLOPE(I-1)*H(I)+SLOPE(I)*H(I-1))/(H(I-1)+H(I))
	    dCHIdR(I)=( SIGN(ONE,SLOPE(I-1))+SIGN(ONE,SLOPE(I)) )*
	1            MIN(ABS(SLOPE(I-1)),ABS(SLOPE(I)),0.5D0*ABS(dCHIdR(I)))
	  END DO
	  dCHIdR(ND)=SLOPE(ND-1) +
	1             (SLOPE(ND-1)-SLOPE(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
	  dCHIdR(ND)=( SIGN(ONE,SLOPE(ND-1))+SIGN(ONE,dCHIdR(ND)) )*
	1            MIN(ABS(SLOPE(ND-1)),0.5D0*ABS(dCHIdR(ND)))
!
!
	ELSE IF(METHOD .EQ. 'LOGMON')THEN
!
! In this method a cubic is fitted between adjacent points (in the log-log
! plane). The first derivative at each node is estimated from the parabola 
! passing through that node and the adjacent points. The derivatives are then 
! converted to the linear-linear plane, and adjusted to ensure that the curve 
! in every interval within this plane is monotonic. This ensures that 
! integrals computed using these derivatives remain positive.
! (After Steffen, 1990, A&A, 239, 443-450).
!
! Altered 10-Feb-1999: Derivatives now adjusted in Linear-Linear plane.
!
! Compute the interval between grid points, and the slope in each interval
! They are computed as vectors as the values are used in 2 intervals.
! It also allows vectorization.
!
	  DO I=1,ND-1
	    H(I)=LOG(R(I+1)/R(I))
	    SLOPE(I)=LOG(CHI(I+1)/CHI(I))/H(I)
	  END DO
!
          dCHIdR(1)=SLOPE(1) +(SLOPE(1)-SLOPE(2))*H(1)/(H(1)+H(2))
	  DO I=2,ND-1
            dCHIdR(I)=(SLOPE(I-1)*H(I)+SLOPE(I)*H(I-1))/(H(I-1)+H(I))
	  END DO
	  dCHIdR(ND)=SLOPE(ND-1) +
	1              (SLOPE(ND-1)-SLOPE(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
!
! Change from derivative in the LOG-LOG plane to the LINEAR-LINEAR plane.
!
	  DO I=1,ND
	    dCHIdR(I)=CHI(I)*dCHIdR(I)/R(I)
	  END DO
!
! Compute slopes in linear-linear plane. These are used to used to ensure
! cubic is monotonic in each interval.
!
	  DO I=1,ND-1
	    SLOPE(I)=(CHI(I+1)-CHI(I))/(R(I+1)-R(I))
	  END DO
!
! Now adjust the derivatives so that curve is monotonic in each interval.
!
	  dCHIdR(1)=( SIGN(ONE,SLOPE(1))+SIGN(ONE,dCHIdR(1)) )*
	1            MIN(ABS(SLOPE(1)),0.5D0*ABS(dCHIdR(1)))
	  DO I=2,ND-1
	    dCHIdR(I)=( SIGN(ONE,SLOPE(I-1))+SIGN(ONE,SLOPE(I)) )*
	1            MIN(ABS(SLOPE(I-1)),ABS(SLOPE(I)),0.5D0*ABS(dCHIdR(I)))
	  END DO
	  dCHIdR(ND)=( SIGN(ONE,SLOPE(ND-1))+SIGN(ONE,dCHIdR(ND)) )*
	1            MIN(ABS(SLOPE(ND-1)),0.5D0*ABS(dCHIdR(ND)))
C	
	ELSE IF(METHOD(1:4) .EQ. 'ZERO')THEN
	  DO I=1,ND
	    DCHIDR(I)=0.0D0
	  END DO
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in DCHIDR - invalid method'
	  STOP
	END IF
C
	RETURN
	END
C
C 
C
	SUBROUTINE d_DERIVCHI_dCHI(dCHIdr,CHI,R,ND,METHOD)
	USE MOD_TRAP_DERIVATIVES
	IMPLICIT NONE
C
C Altered 02-Mar-1999 : Module MOD_TRAP_DERIVATIVES replaces COMMON
C                         BLOCK TRAPDERIVATIVES.
C                         In LOGMON the derivatives are now adjusted in the
C                         linear-linear plane to ensure monotocity of the
C                         fitting cubic. 
C ALtered 24-May-1996 ; ERROR_LU installed.
C
	INTEGER ND,I
	REAL*8 CHI(ND),dCHIdR(ND),R(ND),LIM
	CHARACTER*6 METHOD
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	REAL*8 SLOPE(ND),H(ND)
	REAL*8 LIN_SLOPE(ND),LIN_H(ND)
	REAL*8 T1,ORIG_dCHIdR
	REAL*8, PARAMETER :: ONE=1.0D0
C
	LIM=3.0D0
	LUER=ERROR_LU()
	IF(ND_TRAP .LT. ND)THEN
	  IF(ALLOCATED(A))THEN
	    DEALLOCATE(A); DEALLOCATE(B); DEALLOCATE(C)
	  END IF
	  ALLOCATE(A(ND))
	  ALLOCATE(B(ND))
	  ALLOCATE(C(ND))
	  ND_TRAP=ND
	END IF
C
	IF(METHOD .EQ. 'LOGLOG')THEN
	  DO I=2,ND-1
	    A(I)=CHI(I)/CHI(I-1)/R(I)/( LOG(R(I-1))-LOG(R(I+1)) )
	    B(I)=dCHIdR(I)/CHI(I)
	    C(I)=-CHI(I)/CHI(I+1)/R(I)/( LOG(R(I-1))-LOG(R(I+1)) )
	  END DO
	  A(1)=0.0D0
	  B(1)=dCHIdR(1)/CHI(1)+1.0D0/R(1)/( LOG(R(1))-LOG(R(2)) )
	  C(1)=-CHI(1)/CHI(2)/R(1)/( LOG(R(1))-LOG(R(2)) )
	  IF( (R(1)-R(3))/(R(1)-R(2)) .GT. LIM )THEN
	    A(2)=CHI(2)/CHI(1)/R(2)/( LOG(R(1))-LOG(R(2)) )
	    B(2)=dCHIdR(2)/CHI(2)-1.0D0/R(2)/( LOG(R(1))-LOG(R(2)) )
	    C(2)=0.0D0
	  END IF
	  IF( (R(ND-2)-R(ND))/(R(ND-1)-R(ND)) .GT. LIM )THEN
	    A(ND-1)=0.0D0
	    B(ND-1)=dCHIdR(ND-1)/CHI(ND-1)+1.0D0/R(ND-1)/
	1                   ( LOG(R(ND-1))-LOG(R(ND)) )
	    C(ND-1)=-CHI(ND-1)/CHI(ND)/R(ND-1)/
	1                   ( LOG(R(ND-1))-LOG(R(ND)) )
	  END IF
	  A(ND)=CHI(ND)/CHI(ND-1)/R(ND)/
	1                   ( LOG(R(ND-1))-LOG(R(ND)) )
	  B(ND)=dCHIdR(ND)/CHI(ND) - 1.0D0/R(ND)/
	1                   ( LOG(R(ND-1))-LOG(R(ND)) )
	  C(ND)=0.0D0
!
!
	ELSE IF(METHOD .EQ. 'LOGLIN')THEN
	  DO I=2,ND-1
	    A(I)=CHI(I)/CHI(I-1)/(R(I-1)-R(I+1))
	    B(I)=dCHIdR(I)/CHI(I)
	    C(I)=-CHI(I)/CHI(I+1)/(R(I-1)-R(I+1))
	  END DO
	  A(1)=0.0D0
	  B(1)=dCHIdR(1)/CHI(1)+1.0D0/(R(1)-R(2))
	  C(1)=-CHI(1)/CHI(2)/(R(1)-R(2))
	  IF( (R(1)-R(3))/(R(1)-R(2)) .GT. LIM )THEN
	    A(2)=CHI(2)/CHI(1)/(R(1)-R(2))
	    B(2)=dCHIdR(2)/CHI(2)-1.0D0/(R(1)-R(2))
	    C(2)=0.0D0
	  END IF
	  IF( (R(ND-2)-R(ND))/(R(ND-1)-R(ND)) .GT. LIM )THEN
	    A(ND-1)=0.0D0
	    B(ND-1)=dCHIdR(ND-1)/CHI(ND-1)+1.0D0/(R(ND-1)-R(ND))
	    C(ND-1)=-CHI(ND-1)/CHI(ND)/(R(ND-1)-R(ND))
	  END IF
	  A(ND)=CHI(ND)/CHI(ND-1)/(R(ND-1)-R(ND))
	  B(ND)=dCHIdR(ND)/CHI(ND) - 1.0D0/(R(ND-1)-R(ND))
	  C(ND)=0.0D0
!
!
	ELSE IF(METHOD .EQ. 'LINEAR')THEN
	  DO I=2,ND-1
	    A(I)=1.0D0/(R(I-1)-R(I+1))
	    B(I)=0.0D0
	    C(I)=-1.0D0/(R(I-1)-R(I+1))
	  END DO
	  A(1)=0.0D0
	  B(1)=1.0D0/(R(1)-R(2))
	  C(1)=-B(1)
	  IF( (R(1)-R(3))/(R(1)-R(2)) .GT. LIM )THEN
	    A(2)=B(1)
	    B(2)=C(1)
	    C(2)=0.0D0
	  END IF
	  IF( (R(ND-2)-R(ND))/(R(ND-1)-R(ND)) .GT. LIM )THEN
	    A(ND-1)=0.0D0
	    B(ND-1)=1.0D0/(R(ND-1)-R(ND))
	    C(ND-1)=-B(ND-1)
	  END IF
	  A(ND)=1.0D0/(R(ND-1)-R(ND))
	  B(ND)=-A(ND)
	  C(ND)=0.0D0
!
!
	ELSE IF(METHOD .EQ. 'LOGMON')THEN
!
! In this method a cubic is fitted between adjacent points (in the log-log
! plane). The first derivative at each node is estimated from the parabola 
! passing through that node and the adjacent points. The derivatives are then 
! converted to the linear-linear plane, and adjusted to ensure that the curve 
! in every interval within this plane is monotonic. This ensures that 
! integrals computed using these derivatives remain positive.
! (After Steffen, 1990, A&A, 239, 443-450).
!
! Altered 10-Feb-1999: Derivatives now adjusted in Linear-Linear plane.
!
! Compute the interval between grid points, and the slope in each interval
! They are computed as vectors as the values are used in 2 intervals.
! It also allows vectorization.
!
	  DO I=1,ND-1
	    H(I)=LOG(R(I+1)/R(I))
	    SLOPE(I)=LOG(CHI(I+1)/CHI(I))/H(I)
	    LIN_H(I)=R(I+1)-R(I)
	    LIN_SLOPE(I)=(CHI(I+1)-CHI(I))/LIN_H(I)
	  END DO
	  A(1:ND)=0.0D0
	  B(1:ND)=0.0D0
	  C(1:ND)=0.0D0
!
! At the boundary, the interval close to the boundary is generally much
! smaller than the subsequent interval. Thus at the boundary the slope will
! generally be dominated by the boundary slope, and hence we will only
! consider it in computing the derivatives. This means that that we don't
! have to consider the effect (for example) of CHI(3) on dCHIdr(1) which
! would necessitate an extra vector.
!
! We first compute the derivatives assuming no adjustments have been made
! to ensure monotocity.
!
	  C(1)=ONE/H(1)/CHI(2)
	  B(1)=-ONE/H(1)/CHI(1)
	  DO I=2,ND-1
	    T1=1.0D0/(H(I-1)+H(I))
	    A(I)=-T1*H(I)/H(I-1)/CHI(I-1)
	    C(I)=T1*H(I-1)/H(I)/CHI(I+1)
	    B(I)=T1*(H(I)/H(I-1)-H(I-1)/H(I))/CHI(I)
	  END DO
	  B(ND)=1.0D0/H(ND-1)/CHI(ND)
	  A(ND)=-1.0D0/H(ND-1)/CHI(ND-1)
!
! Convert for LOG-LOG to LINEAR-LINEAR plane (i.e. We multiply dCHIdR[i]
! by CHI[i]/R[i] )
!
	  DO I=1,ND
	    A(I)=A(I) * (CHI(I)/R(I))
	    B(I)=dCHIDR(I)/CHI(I) + B(I) * (CHI(I)/R(I))
	    C(I)=C(I) * (CHI(I)/R(I))
	  END DO
!
! Now modify the A,B and C if we have modified dCHIdR to ensure monotocity.
!
          ORIG_dCHIdR=SLOPE(1) +(SLOPE(1)-SLOPE(2))*H(1)/(H(1)+H(2))
	  IF(dCHIdR(1) .EQ. 0.0D0)THEN
	    B(1)=0.0D0
	    C(1)=0.0D0
	  ELSE IF(ABS(SLOPE(1)) .LT. 0.5D0*ABS(ORIG_dCHIdR))THEN
	    B(1)=-2.0D0/LIN_H(1)
	    C(1)=2.0D0/LIN_H(1)
	  END IF
C
	  DO I=2,ND-1
            ORIG_dCHIdR=
	1       (SLOPE(I-1)*H(I)+SLOPE(I)*H(I-1))/(H(I-1)+H(I))
	1      * CHI(I)/R(I)
	    IF(dCHIdR(I) .EQ. 0.0D0)THEN
	      A(I)=0.0D0
	      B(I)=0.0D0
	      C(I)=0.0D0
	    ELSE IF( MIN(ABS(LIN_SLOPE(I-1)),ABS(LIN_ SLOPE(I)))
	1                             .LT. 0.5D0*ABS(ORIG_dCHIdR) )THEN
	      IF( ABS(LIN_SLOPE(I-1)) .LT. ABS(LIN_SLOPE(I)) )THEN
	        A(I)=-2.0D0/LIN_H(I-1)
	        B(I)=2.0D0/LIN_H(I-1)
	        C(I)=0.0D0
	      ELSE
	        A(I)=0.0D0
	        B(I)=-2.0D0/LIN_H(I)
	        C(I)=2.0D0/LIN_H(I)
 	      END IF
	    END IF
	  END DO
!
	  ORIG_dCHIdR=SLOPE(ND-1) +
	1              (SLOPE(ND-1)-SLOPE(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
	  IF(dCHIdR(ND) .EQ. 0.0D0)THEN
	    A(ND)=0.0D0
	    B(ND)=0.0D0
	  ELSE IF(ABS(SLOPE(ND-1)) .LT. 0.5D0*ABS(ORIG_dCHIdR))THEN
	    A(ND)=-2.0D0/LIN_H(ND-1)
	    B(ND)=2.0D0/LIN_H(ND-1)
	  END IF
!
!
	ELSE IF(METHOD .EQ. 'LINMON')THEN
C
C In this method a cubic is fitted between adjacent points. The first 
C derivative at each node is estimated from the parabola passing through that 
C node and the adjacent points. The derivatives are then adjusted to ensure 
C that the curve in every interval is monotonic. 
C (After Steffen, 1990, A&A, 239, 443-450).
C
C Compute the interval between grid points, and the slope in each interval
C They are computed as vectors as the values are used in 2 intervals.
C It also allows vectorization.
C
	  DO I=1,ND-1
	    H(I)=R(I+1)-R(I)
	    SLOPE(I)=(CHI(I+1)-CHI(I))/H(I)
	  END DO
	  A(1:ND)=0.0D0
	  B(1:ND)=0.0D0
	  C(1:ND)=0.0D0
C
C At the boundary, the interval close to the boundary is generally much
C smaller than the subsequent interval. Thus at the boundary the slope will
C generally be dominated by the boundary slope, and hence we will only
C consider it in computing the derivatives. This means that that we do not
C have to consider the effect (for example) of CHI(3) on dCHIdr(1) which
C would necessitate an extra vector.
C
          ORIG_dCHIdR=SLOPE(1) +(SLOPE(1)-SLOPE(2))*H(1)/(H(1)+H(2))
	  T1=SIGN(ONE,SLOPE(1))+SIGN(ONE,ORIG_dCHIdR)
	  IF( ABS(SLOPE(1)) .LT. 0.5D0*ABS(ORIG_dCHIdR) )THEN
	    C(1)=2.0D0/H(1)
	    B(1)=-C(1)
	  ELSE IF(T1 .NE. 0.0D0)THEN
	    C(1)=ONE/H(1)
	    B(1)=-C(1)
	  END IF
C
	  DO I=2,ND-1
            ORIG_dCHIdR=(SLOPE(I-1)*H(I)+SLOPE(I)*H(I-1))/(H(I-1)+H(I))
	    T1=SIGN(ONE,SLOPE(I-1))+SIGN(ONE,SLOPE(I))
	    IF( ABS(SLOPE(I-1)) .LT. 0.5D0*ABS(ORIG_dCHIdR) .OR.
	1	ABS(SLOPE(I)) .LT. 0.5D0*ABS(ORIG_dCHIdR) )THEN
	      IF( ABS(SLOPE(I-1)) .LT. ABS(SLOPE(I)) )THEN
	        A(I)=-T1/H(I-1)*SIGN(ONE,SLOPE(I-1))
	        B(I)=T1/H(I-1)*SIGN(ONE,SLOPE(I-1))
	      ELSE
	        B(I)=-T1/H(I)*SIGN(ONE,SLOPE(I))
	        C(I)=T1/H(I)*SIGN(ONE,SLOPE(I))
 	      END IF
	    ELSE
	      T1=T1*SIGN(0.5D0,ORIG_dCHIdR)/(H(I-1)+H(I))
	      A(I)=-T1*H(I)/H(I-1)
	      C(I)=T1*H(I-1)/H(I)
	      B(I)=-A(I)-C(I)
	    END IF
	  END DO
C
	  ORIG_dCHIdR=SLOPE(ND-1) +
	1               (SLOPE(ND-1)-SLOPE(ND-2))*H(ND-1)/(H(ND-2)+H(ND-1))
	  T1=SIGN(ONE,SLOPE(ND-1))+SIGN(ONE,ORIG_dCHIdR)
	  IF( ABS(SLOPE(ND-1)) .LT. 0.5D0*ABS(ORIG_dCHIdR) )THEN
	    B(ND)=2.0D0/H(ND-1)
	    A(ND)=-B(ND)
	  ELSE IF( T1 .NE. 0)THEN
	    B(ND)=1.0D0/H(ND-1)
	    A(ND)=-B(ND)     
	  END IF
C
	ELSE IF(METHOD(1:4) .EQ. 'ZERO')THEN
	  DO I=1,ND
	   A(I)=0.0D0
	   B(I)=0.0D0
	   C(I)=0.0D0
	  END DO
	ELSE
	  WRITE(LUER,*)'Error in d_DERRIVCHI_dCHI - invalid method'
	  STOP
	END IF
C
	RETURN
	END
