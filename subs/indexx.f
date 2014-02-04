C
C Routine creates an index array to sort a sequence of DOUBLE precision
C number into increasing (NUMER true) or decreasing (NUMER false)
C numerical order. Thus to sort use
C
C         ARRAY_OUT(I)=ARRAY_IN(INDX(I))
C
C Modified 22-may-1997 - SORTLOG installed.
C Modified 29-Jun-1989 - NUMER option installed. RANK routine created.
C Created -
C
C Related routines:
C                  SORTINT
C                  SORTDP
C                  SORTLOG
C                  SORTCHAR
C                  RANK  - Creates a rank vector from the index vector.
C                          Thus RANK(I) gives the location of the Ith
C                          variable in the newly sorted array.
C
	SUBROUTINE INDEXX(N,ARRIN,INDX,NUMER)
	IMPLICIT NONE
C
	INTEGER N
	INTEGER INDX(N)
	REAL*8 ARRIN(N)
	LOGICAL NUMER
C
	REAL*8 Q
	INTEGER L,IR,I,J,INDXT,ISAV
C
	DO J=1,N
	  INDX(J)=J
	END DO
	L=N/2+1
	IR=N
C
10	CONTINUE
	  IF(L .GT. 1)THEN
	    L=L-1
	    INDXT=INDX(L)
	    Q=ARRIN(INDXT)
	  ELSE
	    INDXT=INDX(IR)
	    Q=ARRIN(INDXT)
	    INDX(IR)=INDX(1)
	    IR=IR-1
	     IF(IR .EQ. 1)THEN
	       INDX(1)=INDXT
C
C For historical reasons, have sorted numbers into decreasing numerical
C order. Need to convert to numerical order if NUMER is true.
C Note that if N is odd, the middle INDX value is correct.
C
	       IF(NUMER)THEN
	         DO I=1,N/2
	           ISAV=INDX(I)
	           J=N-I+1
	           INDX(I)=INDX(J)
	           INDX(J)=ISAV
	         END DO
	       END IF
	       RETURN
	     END IF
	  END IF
	  I=L
	  J=L+L
20	  IF(J .LE. IR)THEN
	    IF(J .LT. IR)THEN
	      IF( ARRIN(INDX(J)) .GT. ARRIN(INDX(J+1)) )J=J+1
	    END IF
	    IF( Q .GT. ARRIN(INDX(J)) )THEN
	      INDX(I)=INDX(J)
	      I=J
	      J=J+J
	    ELSE
	      J=IR+1
	    END IF
	  GOTO 20
	END IF
	INDX(I)=INDXT
	GOTO 10
C
	END
C
	SUBROUTINE SORTINT(N,ARRIN,INDX,WORK)
	IMPLICIT NONE
	INTEGER N,I,INDX(N)
	INTEGER ARRIN(N),WORK(N)
C
	DO I=1,N
	  WORK(I)=ARRIN(I)
	END DO
	DO I=1,N
	  ARRIN(I)=WORK(INDX(I))
	END DO
C
	RETURN
	END
C
	SUBROUTINE SORTDP(N,ARRIN,INDX,WORK)
	IMPLICIT NONE
	INTEGER N,I,INDX(N)
	REAL*8 ARRIN(N),WORK(N)
C
	DO I=1,N
	  WORK(I)=ARRIN(I)
	END DO
	DO I=1,N
	  ARRIN(I)=WORK(INDX(I))
	END DO
C
	RETURN
	END
C
	SUBROUTINE SORTLOG(N,ARRIN,INDX,WORK)
	IMPLICIT NONE
	INTEGER N,I,INDX(N)
	LOGICAL ARRIN(N),WORK(N)
C
	DO I=1,N
	  WORK(I)=ARRIN(I)
	END DO
	DO I=1,N
	  ARRIN(I)=WORK(INDX(I))
	END DO
C
	RETURN
	END
C
	SUBROUTINE SORTCHAR(N,ARRIN,INDX,WORK)
	IMPLICIT NONE
	INTEGER N,I,INDX(N)
	CHARACTER*(*) ARRIN(N),WORK(N)
C
	DO I=1,N
	  WORK(I)=ARRIN(I)
	END DO
	DO I=1,N
	  ARRIN(I)=WORK(INDX(I))
	END DO
C
	RETURN
	END
C
	SUBROUTINE RANK(N,INDX,IRANK)
	IMPLICIT NONE
	INTEGER N,INDX(N),IRANK(N)
	INTEGER J
C
	DO J=1,N
	  IRANK(INDX(J))=J
	END DO
C
	RETURN
	END
