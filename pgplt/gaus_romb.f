!
! Subroutine to return the integral (in ANSWER) of a Gaussian line function
! of the form:
!
!           HEIGHT*EXP(- [ABS(X/SCALE)]**EXPONENT )
!
! TOLERANCE is the absolute tolerance desired, such that the absolute error
!
!	 < PARAMS(1)*PARAMS(2)*TOL
!
! Routine uses Romberg integratio, with up to 8 levels of refinement.
! Routine should provide accuracies of better than 1.0E-04 (and much
! more generally) for
!
! EXPONENT > 0.5
!
	SUBROUTINE GAUS_ROMB(ANSWER,HEIGHT,SCALE,EXPONENT,TOLERANCE)
	IMPLICIT NONE
!
! Altered 02-Apr-2008 : Use convention that EWs are +ve for absorption lines
! Created 28-Sep-2007
!
	REAL*8 ANSWER
	REAL*8 HEIGHT
	REAL*8 SCALE
	REAL*8 EXPONENT
	REAL*8 TOLERANCE
!
! R(I,1) stores the ith trapazoidal integration.
!
	INTEGER, PARAMETER :: M=8
	REAL*8 R(M,M)
!
	REAL*8 H		!Step size
	REAL*8 HMAX		!Maximumstep size adopted
	REAL*8 RANGE		!Range of integration
!
	INTEGER J,K,L		!Loop indices
	INTEGER N		!Number of steps
	LOGICAL TESTING
!
	EXPONENT=ABS(EXPONENT)
	IF(EXPONENT .LT. 0.5)THEN
	  WRITE(6,*)'Error in GAUS_ROMB'
	  WRITE(6,*)'Unable to integrate the modified Gaussian if'//
	1                      ' its exponent is < 0.5'
	  ANSWER=1000
	  RETURN
	END IF
!
	TESTING=.TRUE.
	RANGE=NINT(23.0**(1.0D0/EXPONENT))+1		!Set at function=1.0E-10
	HMAX=1  					!RANGE/10.0D0				!Maximum step size
	IF(TESTING)THEN
	  WRITE(6,*)'RANGE=',RANGE
	  WRITE(6,*)'HMAX=',HMAX
	END IF
!
	DO K=1,M
	  H=HMAX/2**(K-1)
	  R(K,1)=0.5
	  N=NINT(RANGE/H)
	  IF(TESTING)THEN
	    WRITE(6,'(A,I8)')' Refinement number:',K
	    WRITE(6,'(A,I8)')' Number of steps is:',N
	    WRITE(6,'(A,ES16.6)')' Step size is:',H
	  END IF
	  DO L=2,N
	    R(K,1)=R(K,1)+EXP(-((L-1)*H)**EXPONENT )
	  END DO
	  R(K,1)=R(K,1)+0.5D0*EXP(-(N*H)**EXPONENT )
	  R(K,1)=2.0D0*R(K,1)*H                          !*1.0D0/SQRT(ACOS(-1.0D0))
	  IF(TESTING)THEN
	    WRITE(6,'(A,ES14.8)')'Current unscaled trapazoidal answer:',R(K,1)
	  END IF
	  DO L=2,K
	    J=K-L+1
	    R(J,L)=R(J+1,L-1)+(R(J+1,L-1)-R(J,L-1))/(4**(L-1)-1)
	  END DO
	  IF(K .NE. 1)THEN
	    IF(TESTING)THEN
	      WRITE(6,'(A,ES14.6)')'Current accuracy is',ABS(R(1,K)-R(1,K-1))
	    END IF
	    IF( ABS(R(1,K)-R(1,K-1)) .LT. TOLERANCE)EXIT
	  END IF	
	END DO
!
! We now use the convention that EWs are +ve for absorption lines.
!
	K=MIN(K,M)
	ANSWER=-R(1,K)*ABS(SCALE)*HEIGHT
	IF(TESTING)THEN
	  WRITE(6,'(A,ES14.6)')'Final unscaled integeral is:',R(1,K)
	  WRITE(6,'(A,ES14.6)')'          Final integral is:',ANSWER
	END IF
!
	RETURN
	END
