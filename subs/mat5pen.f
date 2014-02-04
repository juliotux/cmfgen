C
C Evaluates the following matrix equation.
C
C                        A = A - B. C + D. (F . C - E)
C
C where:
C               A = real A[N,NS] - Retuned with solution
C               B = real B[N,N]  - Unchanged
C               C = real C[N,NS] - Unchanged
C               D = real D[N,N]  - Unchanged
C               E = real E[N,NS] - Unchanged
C               F = real F[N,N] -  Unchanged
C               WRKMAT = real A[N,NS] - Working matrix (corrupted)
C               VEC = real VEC[N] - Working vector (corrupted)
C               ENONZERO - If FALSE, E is assumed to be zero, and E is not
C                             accessed.
C               DNONZERO - If FALSE, D is assumed to be zero, and D, F and E
C                             are not accessed.
C
C For use in the solution of blockbanded matrices.
C
	SUBROUTINE MAT5PEN(A,B,C,D,E,F,WRKMAT,VEC,N,NR,ENOTZERO,DNOTZERO)
	IMPLICIT NONE
C
C Altered 24-May-1996 : ERROR_LU inserted.
C Created 12-Feb-1989
C
	LOGICAL ENOTZERO,DNOTZERO
	INTEGER N,NR
	REAL*8 A(N,NR),C(N,NR),E(N,NR),WRKMAT(N,NR)
	REAL*8 B(N,N),D(N,N),F(N,N),VEC(N)
C
	INTEGER I,L
C
	CHARACTER*1 NO_TRANS
	REAL*8 DP_ONE,DP_ZERO,DP_NEG_ONE
	PARAMETER (DP_ZERO=0.0D0)
	PARAMETER (DP_ONE=1.0D0)
	PARAMETER (DP_NEG_ONE=-1.0D0)
	PARAMETER (NO_TRANS='N')
C
	IF(DNOTZERO)THEN
C
C Computes WRKMAT = F. C
C
	  CALL DGEMM(NO_TRANS,NO_TRANS,N,NR,N,DP_ONE,F,N,C,N,DP_ZERO,
	1            WRKMAT,N)
C
C Computes WRKMAT=WRKMAT - E  [ (F.C-E) ]
C
	  IF(ENOTZERO)THEN
	    DO L=1,NR
	      DO I=1,N
	       WRKMAT(I,L)=WRKMAT(I,L)-E(I,L)
	      END DO
	    END DO
	  END IF
C
C Computes A = A + D . WRKMAT [ A+ D(F.C-E) ]
C
	  CALL DGEMM(NO_TRANS,NO_TRANS,N,NR,N,DP_ONE,D,N,WRKMAT,N,DP_ONE,
	1            A,N)
	END IF
C
C
C Computes A= A-B. C
C
	CALL DGEMM(NO_TRANS,NO_TRANS,N,NR,N,DP_NEG_ONE,B,N,C,N,DP_ONE,
	1            A,N)
C
	RETURN
	END
