C
C Routine to output collisional data in 132 column format. Only 3 significant
C digits are output. OPTION is presently not utilized, while DESC is
C used to detrmine the file name.
C
	SUBROUTINE WR_COL(OMEGA,LEV_NAME,N,DESC,LU,OPTION)
	IMPLICIT NONE
C
C Created 11-Aug-1997
C
	INTEGER N,LU
	REAL*8 OMEGA(N,N)
	CHARACTER*(*) LEV_NAME(N),OPTION,DESC
C
C Internal variables.
C
	INTEGER I,J,K,L,M,IOS
	INTEGER N_PER_LINE,LMAX,LIM,ST_POS
	CHARACTER*80 FORM
	CHARACTER*132 STRING
	CHARACTER*30 TMP_NAME
C
C Number of digits to output OMEGA. Must include spaces used to separate 
C number from proceeding number. Must be the same as the number X in
C EX.3                                 
c
	INTEGER, PARAMETER :: DIG_PER_NUM=11
	INTEGER, PARAMETER :: IZERO=0
C
	FORM=TRIM(DESC)//'_COL'
	CALL GEN_ASCI_OPEN(LU,FORM,'UNKNOWN',' ',' ',IZERO,IOS)
C
C Determine maximum name length. Used for formatting. Then determine the
C maximum numbers to output in one line, allowing for the level name.
C
	LMAX=0
	DO I=1,N
	  LMAX=MAX(LEN_TRIM(LEV_NAME(I)),LMAX)
	END DO
	N_PER_LINE=MIN(N,(132-LMAX-2)/DIG_PER_NUM)
C
	DO M=1,(N+N_PER_LINE-1)/N_PER_LINE
	  LIM=MIN(M*N_PER_LINE,N)
C
C Determine header string. Only the last DIG_PER_NUM-2 characters in the 
C name are output.
C
	  ST_POS=LMAX+1
	  STRING=' '
	  DO J=(N_PER_LINE*(M-1)+1),LIM
            K=LEN_TRIM(LEV_NAME(J))
	    IF(K .LT. DIG_PER_NUM-2)THEN
	      TMP_NAME=' '
	      TMP_NAME(DIG_PER_NUM+1-K:DIG_PER_NUM)=LEV_NAME(J)
	    ELSE
	      TMP_NAME='  '//LEV_NAME(J)(K-DIG_PER_NUM+3:K)
	    END IF
	    STRING(ST_POS:)=TMP_NAME
	    ST_POS=ST_POS+DIG_PER_NUM
	  END DO
	  WRITE(LU,*)' '
	  WRITE(LU,'(A)')TRIM(STRING)
C
C We only ouput OMEGA(I,J) for J>=I. K and L help set up the appropriate
C format to do this.  For the INTEL compiler, K must be at least 1
C (altered 24-Jun-2003)
C
	  K=1
	  L=N_PER_LINE
	  DO I=1,MIN(N_PER_LINE*M,N)
	    IF(N_PER_LINE*M-I .LT. N_PER_LINE-1)THEN
              K=K+DIG_PER_NUM
	      L=L-1
	    END IF
	    WRITE(FORM,'(A,I4.4,A,I4.4,A,I4.4,A)')
	1                   '(A',LMAX,',',K,'X,1P,',L,'E11.3)'
	    WRITE(LU,FORM)
	1       LEV_NAME(I),(OMEGA(I,J),J=MAX((N_PER_LINE*(M-1)+1),I),LIM)
	  END DO
	END DO
C
	CLOSE(LU)
C
	RETURN
	END
