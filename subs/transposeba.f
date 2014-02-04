C
C Routine to transpose the (N*ND)**2 Matrix so that the simultaneous
C equations can be solved LU decomposition.
C
C BA and NEWBA may refer to the same memory space, as can TBA and VJ.
C TBA must be of size N*N or more, and VJ must be at least of size N*N*ND.
C
	SUBROUTINE TRANSPOSEBA(BA,NEWBA,TBA,VJ,N,ND)
C
	IMPLICIT NONE
	INTEGER N,ND
	INTEGER I,J,K,L,KP,LP
C
	REAL*8 BA(N,N,ND,ND)
	REAL*8 NEWBA(N,ND,N,ND)
	REAL*8 TBA(N,N)
	REAL*8 VJ(N,ND,N)
C
C Interchange last two indices in BA array. At the same time we reorder
C the eqations so that the deepest part of the atmosphere occurs first
C in the matrix.
C
	DO L=1,ND-1
	  DO K=1,ND-L
	    LP=ND-K+1
	    KP=ND-L+1
	    DO J=1,N
	      DO I=1,N
	        TBA(I,J)=BA(I,J,K,L)
	        BA(I,J,K,L)=BA(I,J,KP,LP)
	        BA(I,J,KP,LP)=TBA(I,J)
	      END DO
	    END DO
	 END DO
	END DO
C
C Now intechange second and third indices. In addition we now reverse
C the order of the unknowns at a given depth.
C 
C
	DO K=1,ND
C
	  DO L=1,ND	
	    DO J=1,N
	      DO I=1,N
	        VJ(N-I+1,L,N-J+1)=BA(I,J,L,K)
	      END DO
	    END DO
	  END DO
C
	  DO L=1,ND	
	    DO J=1,N
	      DO I=1,N
	        NEWBA(I,L,J,K)=VJ(I,L,J)
	      END DO
	    END DO
	  END DO
C
	END DO
C
	RETURN
	END
