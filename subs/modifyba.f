C
C Altered 04-Mar-88 - Routine to eliminate matrices below diagonal matrices,
C                     and to put diaginal matrix in U form with diagonal
C equal to unity. This then makes the back substitution trivial. The number
C of matrices eliminated below the diagonal is controlled by the variable
C MATELIM. When MATELIM=0, the final sloution should be identical to the
C block iterative technige of NAGGSIT. MATELIM=1 cooresponds to a tridiagonal
C approximation.
C
C
	SUBROUTINE MODIFYBA(BA,STEQ,C,CMORE,N,ND,MATELIM,IFLAG)
	IMPLICIT NONE
C
C Altered 24-MAy-1996 : DOUBLE PRECISION declaration removed
C                       ERROR_LU installed.
C
	INTEGER N,ND,MATELIM
	LOGICAL IFLAG
	REAL*8 BA(N,N,ND,ND),STEQ(N,ND),C(N),CMORE(N,5)
C
	INTEGER I,J,K,L,M,JM,MAT
	REAL*8 TOL,BIG,RBIG,T1
	PARAMETER (TOL=1.0D-30)
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
	IF(MATELIM .GT. 5)THEN
	  WRITE(LUER,*)'MATELIM is too large (should be < 6)'
	  IFLAG=.FALSE.
	END IF
C
C Note that C and CMORE are working arrays.
C
	IFLAG=.TRUE.
C
C Enter main loop to start the elimination.
C
	DO K=1,ND					!Depth
	  DO I=1,N					!Variable
C
C Find maximum coefficient in column I
C
	    BIG=0.0D0
	    JM=I
	    DO J=I,N					!Equation
	      IF(ABS(BIG).LT.ABS(BA(J,I,K,K)))THEN
		BIG=BA(J,I,K,K)
		JM=J
	      END IF
	    END DO
C
C Check that pivot is > TOL.
 	    IF (ABS(BIG) .LE. TOL)THEN
	      IFLAG=.FALSE.
	      RETURN
	    END IF
C
	    BA(JM,I,K,K)=BA(I,I,K,K)
	    RBIG=1.0D0/BIG
	    BA(I,I,K,K)=1.0D0
C
C Store column elements
C
	    DO J=I+1,N
	      C(J)=BA(J,I,K,K)
	    END DO
	    DO MAT=1,MATELIM
	      IF(K+MAT .LE. ND)THEN
                DO J=1,N
	          CMORE(J,MAT)=BA(J,I,K,K+MAT)
	        END DO
	      END IF
	    END DO
C
C Perform elimination on the right hand side. First step is to interchange
C the relvant rows.
C
	    T1=STEQ(JM,K)*RBIG
	    STEQ(JM,K)=STEQ(I,K)
	    STEQ(I,K)=T1
	    DO M=I+1,N
	      STEQ(M,K)=STEQ(M,K)-C(M)*STEQ(I,K)
	    END DO
C
	    DO MAT=1,MATELIM
	      IF(K+MAT .LE. ND)THEN
	        DO M=1,N
	          STEQ(M,K+MAT)=STEQ(M,K+MAT)-CMORE(M,MAT)*STEQ(I,K)
	        END DO
	      END IF
	    END DO
C
C Perform elimination for each column
C
	    DO L=1,ND
	      IF(L .EQ. K)THEN
	        DO J=I+1,N
	          T1=BA(JM,J,K,K)*RBIG
	          BA(JM,J,K,K)=BA(I,J,K,K)
	          BA(I,J,K,K)=T1
	          DO M=I+1,N
	            BA(M,J,K,K)=BA(M,J,K,K)-C(M)*BA(I,J,K,K)
	          END DO
	          BA(J,I,K,K)=0.0D0			!Need to set zero.
	        END DO
	      ELSE
	        DO J=1,N
	          T1=BA(JM,J,L,K)*RBIG
	          BA(JM,J,L,K)=BA(I,J,L,K)
	          BA(I,J,L,K)=T1
	          DO M=I+1,N
	            BA(M,J,L,K)=BA(M,J,L,K)-C(M)*BA(I,J,L,K)
	          END DO
	        END DO
	      END IF
C
C Eliminate matrix immediately below diagonal N*N Matrix.
C
	      DO MAT=1,MATELIM
	        IF(K+MAT .LE. ND)THEN
	          DO J=1,N
	            DO M=1,N
	              BA(M,J,L,K+MAT)=BA(M,J,L,K+MAT)-
	1                                CMORE(M,MAT)*BA(I,J,L,K)
	            END DO
	          END DO
	        END IF
	      END DO	!MAT
C
	    END DO	!L
	  END DO	!I
	END DO		!K
C
	RETURN
	END
