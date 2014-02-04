C
C Subroutine to solve a tridiagonal system of N1 simultaneous
C Equations which are tridiagonal in nature for N2 R.H.sides.
C
C The equations to be solved are assumed to have the form
C
C    A(i).X(i-1) - [H(i)+A(i)+C(i)].X(i) + C(i).X(i+1) = D(i)
C
	SUBROUTINE THOMAS_RH(A,H,C,D,N1,N2)
	IMPLICIT NONE
C
C Altered 29-Sep-1997 : Scaler loop code installed to improve speed on a
C                         an Alphastation.
C Altered 21-Jun-1996 : Bug Fix: H now set to DIV to allow SIMPTH entry
C                         (Thus DIV does not need to be saved between calls.)
C Altered 29-May-1996 : Loops reversed in forward elimination and the 
C                         backward substitution to allow CRAY vectorization.
C Altered 17-May-1996 : DIV now dimension by N1 (using F90 advantages).
C
	INTEGER N1,N2
	REAL*8 A(N1),H(N1),C(N1),D(N1,N2)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
	INTEGER I,J
	REAL*8 DIV(N1)
C
C Change the following statement to TRUE if running on a VECTOR machine.
C
	LOGICAL, PARAMETER :: VECTOR_MACHINE=.FALSE.
C
C Compute quantities that will be used repeatedly if the same tridiagonal
C system is used for many R.H. Sides.
C
	C(N1)=0				!As used.
	DIV(1)=1.0/(C(1)+H(1))
	C(1)=C(1)*DIV(1)
	H(1)=H(1)*DIV(1)
	DO I=2,N1
	  DIV(I)=1.0D0/(A(I)*H(I-1)+H(I)+C(I))
	  H(I)=(A(I)*H(I-1)+H(I))*DIV(I)
	  C(I)=C(I)*DIV(I)
	END DO
	H(1:N1)=DIV(1:N1)
C
C Entry for Thomas algorithm when H,C have been previously modified.
C
	ENTRY SIMPTH_RH(A,H,C,D,N1,N2)
C
C NB: Generally N2 is 1 or ND, hence the check.
C
	IF(N2 .EQ. 1)THEN
C
C Forward elimination.
C
	  D(1,1)=D(1,1)*H(1)
	  DO I=2,N1
	    D(I,1)=(D(I,1)+A(I)*D(I-1,1))*H(I)
	  END DO
C
C Backward substitution.
C
	  D(N1,1)=-D(N1,1)
	  DO I=N1-1,1,-1
	    D(I,1)=C(I)*D(I+1,1)-D(I,1)
	  END DO
C
	ELSE IF(VECTOR_MACHINE)THEN
C
C This section is for a VECTOR machine. Loop over inner index is outer
C loop to remove a dependency. This in inefficient on scaler machines
C as array is not accessed sequentially.
C
C Forward elimination.
C
	  DO J=1,N2
	    D(1,J)=D(1,J)*H(1)
	  END DO
	  DO I=2,N1
	    DO J=1,N2
	      D(I,J)=(D(I,J)+A(I)*D(I-1,J))*H(I)
	    END DO
	  END DO
C
C Backward substitution.
C
	  DO J=1,N2
	    D(N1,J)=-D(N1,J)
	  END DO
	  DO I=N1-1,1,-1
	    DO J=1,N2
	      D(I,J)=C(I)*D(I+1,J)-D(I,J)
	    END DO
	  END DO
	ELSE 
C
C This section of code is for a scaler machine where the dependence
C on a previous computation does not matter. More efficient than previous
C code as array is accessed in correct manner.
C
	  DO J=1,N2
	    D(1,J)=D(1,J)*H(1)
	    DO I=2,N1
	      D(I,J)=(D(I,J)+A(I)*D(I-1,J))*H(I)
	    END DO
	  END DO
C
C Backward substitution.
C
	  DO J=1,N2
	    D(N1,J)=-D(N1,J)
	    DO I=N1-1,1,-1
	      D(I,J)=C(I)*D(I+1,J)-D(I,J)
	    END DO
	  END DO
	END IF
C
	RETURN
	END
