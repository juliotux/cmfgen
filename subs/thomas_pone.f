C
C Subroutine to solve a tridiagonal system of N1 (.GE. 3) simultaneous
C equations which are tridiagonal in nature for N2 R.H.Sides.
C
C The equations are assumed to have the form.
C
C    A(i).X(i-1) + B(i).X(i) + C(i).X(i+1) = D(i)
C
C By definition A(1)=C(N1)=0. The solutions are returned in D.
C
C This routine should not be used for solving the transfer equation when
C the optical depth steps are less than approximately 10^{-5}.
C
	SUBROUTINE THOMAS_PONE(A,B,C,D,N1,N2)
	IMPLICIT NONE
C
C Altered 08-Feb-2004 : Forced constants to be double precission.
C Altered 29-Sep-1997 : Scaler loop code installed to improve speed on a
C                         an Alphastation.
C Altered 29-May-1996 : Loops reversed in forward elimination and the 
C                         backward substitution to allow CRAY vectorization.
C Altered 21-Feb-1995 :  Cleaned
C
	INTEGER N1,N2
	REAL*8 A(N1),B(N1),C(N1),D(N1,N2)
C
C Local variables
C
	INTEGER I,J
C
C Change the following statement to TRUE if running on a VECTOR machine.
C
	LOGICAL, PARAMETER :: VECTOR_MACHINE=.FALSE.
C
C Compute quantities that will be used repeatedly if the same tridiagonal
C system is used for many R.H. Sides.
C
	B(1)=1.0D0/B(1)
	C(1)=-C(1)*B(1)
	A(1)=-A(1)*B(1)
	B(2)=1.0D0/(B(2)+A(2)*C(1))
	C(2)=-(C(2)+A(2)*A(1))*B(2)
	DO I=3,N1
	  B(I)=1.0D0/(B(I)+A(I)*C(I-1))
	  C(I)=-C(I)*B(I)
	END DO
C
C Entry for THOMAS algorithm when B AND C have already been modified.
C
	ENTRY SIMPTH_PONE(A,B,C,D,N1,N2)
C
	IF(N2 .EQ. 1)THEN
C
C Forward elimination.
C
	  D(1,1)=D(1,1)*B(1)
	  DO I=2,N1
	    D(I,1)=(D(I,1)-A(I)*D(I-1,1))*B(I)
	  END DO
C
C Perform the back substitution. NB D(N1,1)=D(N1,1) is first step.
C
	  DO I=N1-1,1,-1
	    D(I,1)=D(I,1)+C(I)*D(I+1,1)
	  END DO
	  D(1,1)=D(1,1)+A(1)*D(3,1)
	ELSE IF(VECTOR_MACHINE)THEN
C
C This section is for a VECTOR machine. Loop over inner index is outer
C loop to remove a dependancy. This in inefficient on scaler machines
C as array is not accessed sequentially.
C
C Forward elimination.
C
	  DO J=1,N2
	    D(1,J)=D(1,J)*B(1)
	  END DO
	  DO I=2,N1
	    DO J=1,N2
	      D(I,J)=(D(I,J)-A(I)*D(I-1,J))*B(I)
	    END DO
	  END DO
C
C Perform the back substitution. NB D(N1,J)=D(N1,J) is first step.
C
	  DO I=N1-1,1,-1
	     DO J=1,N2
	       D(I,J)=D(I,J)+C(I)*D(I+1,J)
	    END DO
	  END DO
	ELSE
C
C This section of code is for a scaler machine where the dependence
C on a previous computation does not matter. More efficient than previous
C code as array is accessed in correct manner.
C
C Forward elimination.
C
	  DO J=1,N2
	    D(1,J)=D(1,J)*B(1)
	    DO I=2,N1
	      D(I,J)=(D(I,J)-A(I)*D(I-1,J))*B(I)
	    END DO
	  END DO
C
C Perform the back substitution. NB D(N1,J)=D(N1,J) is first step.
C
	  DO J=1,N2
	    DO I=N1-1,1,-1
	       D(I,J)=D(I,J)+C(I)*D(I+1,J)
	    END DO
	  END DO
	END IF
C
	RETURN
	END
