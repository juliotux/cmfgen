	SUBROUTINE BLKBAND(BA,STEQ,WRKMAT,VEC,SOL_TYPE,FLAG,
	1                      DIAG_IND,N,NUM_BNDS,ND)
	IMPLICIT NONE
C
C Altered 24-May-1996 - LUER inserted
C                       String continuations across lines removed.
C                       NMAX removed using F90
C
C Altered  7-Jul-1994 - NAG ROUTINES REPLACED BY LINPACK ROUTINES.
C                       (Editing initiall took place on 14-APR-1994 , but
C                        7-July was final debugging seesion for DIAG and
C                        TRIDIAG sections).
C
C Altered 29-May-1989 - FLAG wasn't being set to TRUE for successful DIAGONAL
C                       solution.
C Created 16-Feb-1989 - Based on BABLKBAND - Tested successfully.
C
C
C Subroutine to solve a "Block-Banded" system of simultaneous equations
C                 BA. X. = STEQ
C The solution X is returned in STEQ. BA is corrupted.
C
C The "BANDED" matrix can either be 'Diagonal', 'Tridiagonal', or
C 'Pentadiagonal'.
C
C Routine is currently designed to operate on the full
C variation matrix (BA in CMFGEN etc) which is of dimension N.N.ND.ND,
C or on BA(n,n,NUM_BNDS,nd) where NUM_BNDS refers to the number of bands.
C Easily modified if the banded components have been stored seaprately
C as "BLOCK VECTORS".
C
C Let A[k], B[k], C[k], D[k], and E[k] be the sub-matrices (dimension N*N) of
C the large Block-Pentadiagonal matrix. Then
C                               A[k]=BA( , ,DI-2,K)
C                               B[k]=BA( , ,DI-1,K)
C                               C[k]=BA( , ,DI,K)
C                               D[k]=BA( , ,DI+1,K)
C                               E[k]=BA( , ,DI+2,K)
C
C where DI=K if ND=NUM_BNDS, else DI=DIAG_IND.
C
C Likewise we have              L[k] =STEQ( , K)
C
C Elimination defined by:
C                  phi[1]= C^{-1) . L[1]
C                  alpha[1]= C^{-1) . D[1]
C                  gamma[1]= C^{-1) . E[1]
C                  SOL[1] = phi[1] - alpha[1].SOL[2] - gamma[1].SOL[3]
C Then
C                  m[k]=C[k]-B[k].alpha[k-1] - A[k].( gamma[k-2]
C                                            - alpha[k-2].alpha[k-1] )
C                  phi[k]=m[k]^{-1}(L[k] - B[k].phi[k-1]- A[k].( phi[k-2]
C                                            - alpha[k-2].phi(k-1) )
C                  alpha[k]=m[k]^{-1}(D[k]-B[k].gamma[k-1] +
C                                            + A[k].alpha[k-2].gamma[k-1]
C                  gamma[k]=m[k]^{-1} . E[k]
C with
C                  SOL[k]= phi[k] - alpha[k]. SOL[k+1] - gamma[k].SOL[k+2]
C
C Note that for k=2, we take A[k]==0.
C           for k=ND-1, we hace E[ND-1]==0 and hence gamma[ND-1]==0
C           for k=ND, we have D[ND]=E[ND]=0 and hence alpha[ND]=gamma[ND]==0.
C Thus
C                  m[ND]=C[ND]-B[ND].alpha[ND-1] - A[ND].( gamma[ND-2]
C                                            - alpha[ND-2].alpha[ND-1] )
C                  phi[ND]=m[ND]^{-1}(L[ND] - B[ND].phi[ND-1]
c                                            - A[ND].( phi[ND-2]
C                                            - alpha[ND-2].phi(ND-1) )
C and              SOL(ND)=phi[ND]
C
C The back substitution is simly obtained from the exprssion for SOL[k]
C since SOL[ND] is now known.
C
C The tridiagonal case can be obtained from these expressinons by setting
C A[k]=E[k]=0, and A[1]=D[ND]=0. When the tridiagonal option is indicated,
C these matrices are ASSUMED to be zero.
C 
C**************************************************************************
C SUBROUTINES CALLED
C ******************
C
C          DGEMV(TRANS,M,N,ALPHA,A,IDA,X,INCX,BETA,Y,INCY)
C
C BLAS2 routine.
C
C Compute Y <--- alpha A.X + beta Y
C
C Where
C
C  	A = Double preciosion A(IDA,N)
C       M = Number of rows in A
C	N = Number of columns in A
C       IDA = First dimension of A [ > max(1,m) ]
C
C	ALPHA, BETA are REAL
C
C       X=REAL   X(1+(N-1}INCX) if TRANS='N'
C       X=REAL   X(1+{M-1}INCY) if TRANS='T'
C	INCX = Stride for X (>0, generally 1)
C
C       Y=REAL   Y(1+(M-1}INCX) if TRANS='N'
C       Y=REAL   Y(1+(N-1}INCX) if TRANS='T'
C	INCY = Stride for Y (>0, generally 1)
C 
C ****************************************************************************
C
C          DGETRF(M,N,A,LDA,IPRIV,INFO)
C
C LAPACK routine.
C
C Perform the LU decomposition of the general matrix A where:
C
C		M   = Number of rows in A
C		N   = Number of columns in A
C       	IDA = First dimension of A [ > max(1,m) ]
C  		A   = Double precision A(IDA,N)
C
C	      IPIV  = Integer work vector with pivot
C                   = IPIV(min[m,n])
C
C             Integer INFO = 0 (successful exit)
C                          = -i (i th argument has illeagal value)
C                          = i (Pivot A(i,i) is exactly zero)
C
C LU decomposition is stored in A.
C
C******************
C
C                   DGETRS(TRANS,N,NRHS,A,LDA,IPIV,B,LDB,INFO)
C
C LAPACK routine.
C
C Solves for the solution of A^{-1} . B . DGETRF must have
C previously been called, and A and IPIV must be passed to DGETRS unchanged.
C The solution is returned in B.
C
C		N      = Order of A
C		NRHS   = Number of RHS
C       	IDA    = First dimension of A [ > max(1,m) ]
C  		A      = Double precision A(IDA,N)
C
C	      IPIV  = Integer work vector with pivot
C                   = IPIV(min[m,n])
C
C  		B      = Double precision B(IDB,NRHS)
C
C             Integer INFO = 0 (successful exit)
C                          = -i (i th argument has illeagal value)
C 
C*******************
C
C	    DGEEQU(M,N,A,LDA,ROW_SF,COL_SF,ROW_CND,COL_CND,MAX_VAL,IFAIL)
C
C LAPACK routine
C
C Compute the Equibriation factors for the M by N matrix A whose first
C dimension is LDA. These factors help improve the condition of the matrix.
C
C
C The equilibrated matrix is found from
C
C         E(I,J)=ROW_SF(I)*A(I,J)*COL_SF(J)
C
C To solve  a set of simultaneous equtions multiply RHS by ROW_SF before
C solution, and by COL_SF after solution.
C
C
C*******************
C                  MAT5PEN(A,B,C,D,E,F,WRKMAT,VEC,N,NS,ENONZERO,DNONZERO)
C
C Evaluates A = A - B. C + D. (F . C - E)
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
C
C******************************************************************************
C
	INTEGER N,ND,NUM_BNDS,DIAG_IND
	LOGICAL FLAG
	CHARACTER*(*) SOL_TYPE
	REAL*8 BA(N,N,NUM_BNDS,ND),STEQ(N,ND)
	REAL*8 WRKMAT(N,N),VEC(N)
C
	REAL*8 ROW_SF(N),COL_SF(N)
	REAL*8 ROW_CND,COL_CND,MAX_VAL
C
C NB: We cheat and use VEC as an intger array in the called subroutines. It
C     doesn't matter as its only a work vector.
C
	INTEGER I,J,K,IFAIL,NSNG,KM2,DIM2
	INTEGER GET_DIAG,DI,DIKM1,DIKM2
	PARAMETER (NSNG=1)
	CHARACTER*10 DESC
	LOGICAL ANONZERO,ENONZERO,TRI
C
	CHARACTER*1 NO_TRANS
	REAL*8 DP_ONE,DP_NEG_ONE,DP_ZERO
	INTEGER INT_ONE
	PARAMETER (DP_NEG_ONE=-1.0D0)
	PARAMETER (DP_ONE=1.0D0)
	PARAMETER (INT_ONE=1)
	PARAMETER (NO_TRANS='N')
C
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
C
C This function computes the index L on BA( , ,?,K) corresponding
C to the local depth variable (i.e that at K). It is equivalent
C to IF (NUM_BNDS .EQ. ND)THEN L=K ELSE L=DIAG END IF
C
	GET_DIAG(K)=(NUM_BNDS/ND)*(K-DIAG_IND)+DIAG_IND
C
	LUER=ERROR_LU()
C
C Check to see if we only require the diagonal solution. For this special
C case, the solutions at each depth are independent.
C
	IF(SOL_TYPE(1:4) .EQ. 'DIAG')THEN
	  DO K=1,ND
	    DI=GET_DIAG(K)
C
C Perform the LU decomposition using DGETRF. We first equilibrize the matrix
C using DGEEQU sot the the maximum row and column values are approximately
C unity.
C
	    CALL DGEEQU(N,N,BA(1,1,DI,K),N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
              DO I=1,N
	        BA(I,J,DI,K)=BA(I,J,DI,K)*ROW_SF(I)*COL_SF(J)
	      END DO
	    END DO
	    CALL DGETRF(N,N,BA(1,1,DI,K),N,VEC,IFAIL)
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRF_DIAG'
	      GOTO 9999
	    END IF
C
C Now perform the solution.
C
	    CALL DGETRS(NO_TRANS,N,NSNG,BA(1,1,DI,K),N,VEC,STEQ(1,K),N,IFAIL)
	    DO J=1,N
	      STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	    END DO
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_DIAG'
	      GOTO 9999
	    END IF
	  END DO
	  FLAG=.TRUE.
	  RETURN
	END IF
C
C 
C
C General Case - Tridiagonal or Pentadiagonal. Set tridiagianl variable
C if tridiagonal.
C	
	IF(SOL_TYPE(1:3) .EQ. 'TRI')THEN
	  TRI=.TRUE.
	  IF(NUM_BNDS .LT. 3)THEN
	    WRITE(LUER,*)'Error in BLKBAND : NUM_BNDS too small ',
	1             'for solution type'
	    FLAG=.FALSE.
	    STOP
	  END IF
	ELSE IF(SOL_TYPE(1:3) .EQ. 'PEN')THEN
	  TRI=.FALSE.
	  IF(NUM_BNDS .LT. 5)THEN
	    WRITE(LUER,*)'Error in BLKBAND : NUM_BNDS too small ',
	1             'for solution type'
	    FLAG=.FALSE.
	    STOP
	  END IF
	ELSE
	  WRITE(LUER,*)'Error in BABLKBAND - invalid SOL_TYPE - ',
	1          'SOLTYPE= ',SOL_TYPE
	  STOP
	END IF
C
C Solve for C[1]^{-1} (Performs LU decomposition).
C
	K=1
	DI=GET_DIAG(K)
	CALL DGEEQU(N,N,BA(1,1,DI,K),N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	DO J=1,N
          DO I=1,N
	    BA(I,J,DI,K)=BA(I,J,DI,K)*ROW_SF(I)*COL_SF(J)
	  END DO
	END DO
	CALL DGETRF(N,N,BA(1,1,DI,K),N,VEC,IFAIL)
	IF(IFAIL .NE. 0)THEN
	  DESC='DGETRF_1'
	  GOTO 9999
	END IF
C
C Solve for phi(1) = C[1]^{-1}.L[1]
C
	DO J=1,N
	  STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	END DO
	CALL DGETRS(NO_TRANS,N,NSNG,BA(1,1,DI,K),N,VEC,STEQ(1,K),N,IFAIL)
	DO J=1,N
	  STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	END DO
	IF(IFAIL .NE. 0)THEN
	  DESC='DGETRS_1'
	  GOTO 9999
	END IF
C
C Solve for alpha(1) = C[1]^{-1}.D[1]
C
	DO J=1,N
          DO I=1,N
	    BA(I,J,DI+1,K)=BA(I,J,DI+1,K)*ROW_SF(I)
	  END DO
	END DO
	CALL DGETRS(NO_TRANS,N,N,BA(1,1,DI,K),N,VEC,BA(1,1,DI+1,K),N,IFAIL)
	DO J=1,N
          DO I=1,N
	    BA(I,J,DI+1,K)=BA(I,J,DI+1,K)*COL_SF(I)
	  END DO
	END DO
	IF(IFAIL .NE. 0)THEN
	  DESC='DGETRS_2'
	  GOTO 9999
	END IF
C
C Solve for gamma[1] = C[1]^{-1}.E[1]
C
	IF(.NOT. TRI)THEN
	  DO J=1,N
            DO I=1,N
	      BA(I,J,DI+2,K)=BA(I,J,DI+2,K)*ROW_SF(I)
	    END DO
	  END DO
	  CALL DGETRS(NO_TRANS,N,N,BA(1,1,DI,K),N,VEC,BA(1,1,DI+2,K),N,IFAIL)
	  DO J=1,N
            DO I=1,N
	      BA(I,J,DI+2,K)=BA(I,J,DI+2,K)*COL_SF(I)
	    END DO
	  END DO
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGETRS_3'
	  END IF
	END IF
C
C 
C
C Now do the general case. NB - the diagonal index is a function of the
C second index, hence the calls to GET_DIAG for DIKM1 and DIKM2.
C
	DO K=2,ND-1
C
	  DI=GET_DIAG(K)
	  DIKM1=GET_DIAG(K-1)
	  ANONZERO=.TRUE.
	  KM2=K-2
	  DIM2=DI-2
	  IF(K.EQ. 2 .OR. TRI)THEN
	    ANONZERO=.FALSE.
	    KM2=1		!Arrays not used anyway
	    DIM2=1
	  END IF
	  DIKM2=GET_DIAG(KM2)
C
C Computes phi[k]' (i.e. ' denotes not yet multiplied by m[k]^{-1} )
C Stored in L[k].
C
	  ENONZERO=.TRUE.
	  CALL MAT5PEN(STEQ(1,K),BA(1,1,DI-1,K),
	1              STEQ(1,K-1),BA(1,1,DIM2,K),
	1              STEQ(1,KM2),BA(1,1,DIKM2+1,KM2),
	1               WRKMAT,VEC,N,NSNG,ENONZERO,ANONZERO,IFAIL)
	  IF(IFAIL .NE. 0)THEN
	    DESC='MAT5PEN_1'
	    GOTO 9999
	  END IF
C
C Computes alpha[k]' (i.e. ' denotes not yet multiplied by m[k]^{-1} )
C Stored in D[k].
C
	  IF(.NOT. TRI)THEN
	    ENONZERO=.FALSE.
	    CALL MAT5PEN(BA(1,1,DI+1,K),BA(1,1,DI-1,K),
	1                  BA(1,1,DIKM1+2,K-1),BA(1,1,DIM2,K),
	1                  BA(1,1,DIKM2+1,KM2),BA(1,1,DIKM2+1,KM2),
	1                  WRKMAT,VEC,N,N,ENONZERO,ANONZERO,IFAIL)
	    IF(IFAIL .NE. 0)THEN
	      DESC='MAT5PEN_2'
	      GOTO 9999
	    END IF
	  END IF
C
C Computes m[k] - Stored in C[k]
C
	  ENONZERO=.TRUE.
   	  CALL MAT5PEN(BA(1,1,DI,K),BA(1,1,DI-1,K),
	1             BA(1,1,DIKM1+1,K-1),BA(1,1,DIM2,K),
	1             BA(1,1,DIKM2+2,KM2),BA(1,1,DIKM2+1,KM2),
	1               WRKMAT,VEC,N,N,ENONZERO,ANONZERO,IFAIL)
	  IF(IFAIL .NE. 0)THEN
	    DESC='MAT5PEN_3'
	    GOTO 9999
	  END IF
C
C Do LU decompostion of m[k]. We first equilibrilze BA.
C
	  CALL DGEEQU(N,N,BA(1,1,DI,K),N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	  DO J=1,N
            DO I=1,N
	      BA(I,J,DI,K)=BA(I,J,DI,K)*ROW_SF(I)*COL_SF(J)
	    END DO
	  END DO
	  CALL DGETRF(N,N,BA(1,1,DI,K),N,VEC,IFAIL)
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGETRF_2'
	    GOTO 9999
	  END IF
C
C Computes phi[k] (Stored in L[k])
C
	  DO J=1,N
	    STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	  END DO
	  CALL DGETRS(NO_TRANS,N,NSNG,BA(1,1,DI,K),N,VEC,STEQ(1,K),N,IFAIL)
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGETRS_4'
	    GOTO 9999
	  END IF
	  DO J=1,N
	    STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	  END DO
C
C Computes alpha[k] (Stored in D[k)).
C
	  DO J=1,N
            DO I=1,N
	      BA(I,J,DI+1,K)=BA(I,J,DI+1,K)*ROW_SF(I)
	    END DO
	  END DO
	  CALL DGETRS(NO_TRANS,N,N,BA(1,1,DI,K),N,VEC,BA(1,1,DI+1,K),N,IFAIL)
	  DO J=1,N
            DO I=1,N
	      BA(I,J,DI+1,K)=BA(I,J,DI+1,K)*COL_SF(I)
	    END DO
	   END DO
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGETRS_5'
	    GOTO 9999
	  END IF
C
C Computes gamma[k] (Stored in E[k]).
C
	  IF(K .NE. ND-1 .AND. .NOT. TRI)THEN
	    DO J=1,N
              DO I=1,N
	        BA(I,J,DI+2,K)=BA(I,J,DI+2,K)*ROW_SF(I)
	      END DO
	    END DO
	    CALL DGETRS(NO_TRANS,N,N,BA(1,1,DI,K),N,VEC,BA(1,1,DI+2,K),N,IFAIL)
	    DO J=1,N
              DO I=1,N
	        BA(I,J,DI+2,K)=BA(I,J,DI+2,K)*COL_SF(I)
	      END DO
	    END DO
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGETRS_6'
	      GOTO 9999
	    END IF
	  END IF
C
	END DO
C
C Perform the elimination on the last row. Alpha and Beta for
C the last row are zero.
C
	K=ND
	DI=GET_DIAG(K)
	DIKM1=GET_DIAG(K-1)
	DIKM2=GET_DIAG(K-2)
	IF(TRI)THEN
	  ANONZERO=.FALSE.
	ELSE
	  ANONZERO=.TRUE.
	END IF
C
C Computes phi[nd]' (i.e. ' denotes not yet multiplied by m[k]^{-1} )
C Stored in D[ND].
C
	ENONZERO=.TRUE.
	CALL MAT5PEN(STEQ(1,K),BA(1,1,DI-1,K),
	1           STEQ(1,K-1),BA(1,1,DI-2,K),
	1           STEQ(1,K-2),BA(1,1,DIKM2+1,K-2),
	1             WRKMAT,VEC,N,NSNG,ENONZERO,ANONZERO,IFAIL)
	IF(IFAIL .NE. 0)THEN
	  DESC='MAT5PEN_1'
	  GOTO 9999
	END IF
C
C Computes m[ND] - Stored in C[ND].
C
	CALL MAT5PEN(BA(1,1,DI,K),BA(1,1,DI-1,K),
	1            BA(1,1,DIKM1+1,K-1),BA(1,1,DI-2,K),
	1            BA(1,1,DIKM2+2,K-2),BA(1,1,DIKM2+1,K-2),
	1              WRKMAT,VEC,N,N,ENONZERO,ANONZERO,IFAIL)
	IF(IFAIL .NE. 0)THEN
	  DESC='MAT5PEN_2'
	  GOTO 9999
	END IF
C
	CALL DGEEQU(N,N,BA(1,1,DI,K),N,ROW_SF,COL_SF,
	1               ROW_CND,COL_CND,MAX_VAL,IFAIL)
	DO J=1,N
           DO I=1,N
	    BA(I,J,DI,K)=BA(I,J,DI,K)*ROW_SF(I)*COL_SF(J)
	  END DO
	END DO
	CALL DGETRF(N,N,BA(1,1,DI,K),N,VEC,IFAIL)
	IF(IFAIL .NE. 0)THEN
	  DESC='DGETRF_3'
	  GOTO 9999
	END IF
C
	DO J=1,N
	  STEQ(J,K)=STEQ(J,K)*ROW_SF(J)
	END DO
	CALL DGETRS(NO_TRANS,N,NSNG,BA(1,1,DI,K),N,VEC,STEQ(1,K),N,IFAIL)
	DO J=1,N
	  STEQ(J,K)=STEQ(J,K)*COL_SF(J)
	END DO
	IF(IFAIL .NE. 0)THEN
	  DESC='DGETRS_7'
	  GOTO 9999
	END IF
C
C
C**************************************************************************
C
C Now we start the back substitution.
C
	DO K=ND-1,1,-1
C
	  DI=GET_DIAG(K)
C
	  IF(K .NE. ND-1 .AND. .NOT. TRI)THEN
	    CALL DGEMV(NO_TRANS,N,N,DP_NEG_ONE,BA(1,1,DI+2,K),N,STEQ(1,K+2),
	1                 INT_ONE,DP_ONE,STEQ(1,K),INT_ONE)
	    IF(IFAIL .NE. 0)THEN
	      DESC='DGEMV_1'
	      GOTO 9999
	    END IF
	  END IF
C
	  CALL DGEMV(NO_TRANS,N,N,DP_NEG_ONE,BA(1,1,DI+1,K),N,STEQ(1,K+1),
	1                 INT_ONE,DP_ONE,STEQ(1,K),INT_ONE)
	  IF(IFAIL .NE. 0)THEN
	    DESC='DGEMV_2'
	    GOTO 9999
	  END IF
C
	END DO
C
C Successfule solution obtained.
C
	FLAG=.TRUE.
C
	RETURN
C
C Error handling routine.
C
9999	CONTINUE
	WRITE(LUER,*)'Error in LINPAC (or BLAS) routine',
	1               DESC,' in BABLKBAND'
	WRITE(LUER,100)K,IFAIL
100	FORMAT(1x,'depth=',I3,10x,'IFAIL=',I3)
	FLAG=.FALSE.
	RETURN
C
	END
