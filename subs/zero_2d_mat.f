!
! Simple do loop to zero a 2D array, with paralleization over the outer loop.
!
	SUBROUTINE ZERO_2D_MAT(A,N,M)
	IMPLICIT NONE
!
! Created 28-Oct-2102
!
	INTEGER N,M
	REAL*8 A(N,M)
	INTEGER I,J
!
!$OMP PARALLEL DO
	DO J=1,M
	  DO I=1,N
	    A(I,J)=0.0D0
	  END DO
	END DO
!
	RETURN
	END
