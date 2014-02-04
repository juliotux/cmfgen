	PROGRAM TST_1
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: IONE=1
	INTEGER  N,NR,I,J,K,L
	REAL*8, ALLOCATABLE :: A(:,:)
	REAL*8, ALLOCATABLE :: B(:,:)
	REAL*8, ALLOCATABLE :: C(:,:)
	REAL*8, ALLOCATABLE :: D(:,:)
!
	CHARACTER*1 NO_TRANS
	REAL*8 DP_ONE,DP_ZERO,DP_NEG_ONE
	PARAMETER (DP_ZERO=0.0)
	PARAMETER (DP_ONE=1.0)
	PARAMETER (DP_NEG_ONE=-1.0)
	PARAMETER (NO_TRANS='N')
!
	WRITE(6,*)'Input Matrix(N,N) dimension  N'
	READ(5,*)N
	NR=N
	ALLOCATE (A(N,N))
	ALLOCATE (B(N,N))
	ALLOCATE (C(N,N))
	ALLOCATE (D(N,N))
!
	DO I=1,N
	  DO J=1,N
	    A(J,I)=J+(I-1)*N
	  END DO
	END DO
	B=2*A
	C=2*A+5
C
C Computes A= A-B. C
C
	CALL TUNE(1,'DGEMM')
	CALL DGEMM(NO_TRANS,NO_TRANS,N,NR,N,DP_NEG_ONE,B,N,C,N,DP_ONE,
	1            A,N)
	CALL TUNE(2,'DGEMM')
	WRITE(6,*)A(5,5),A(N,N)
!
	DO I=1,N
	  DO J=1,N
	    A(J,I)=J+(I-1)*N
	  END DO
	END DO
	B=2*A
	C=2*A+5
!
	CALL TUNE(1,'M2')
	  DO J=1,N
	    DO I=1,N
	      DO K=1,N
	        A(I,J)=A(I,J)-B(I,K)*C(K,J)
	      END DO
	    END DO
	  END DO
	CALL TUNE(2,'M2')
	WRITE(6,*)A(5,5),A(N,N)
!
	DO I=1,N
	  DO J=1,N
	    A(J,I)=J+(I-1)*N
	  END DO
	END DO
	B=2*A
	C=2*A+5
!
	CALL TUNE(1,'TRANS')
	  DO J=1,N
	    DO I=1,N
	      D(J,I)=B(I,J)
	    END DO
	  END DO
	CALL TUNE(2,'TRANS')
!
	CALL TUNE(1,'M3')
	  DO J=1,N
	    DO I=1,N
	      DO K=1,N
	        A(I,J)=A(I,J)-D(K,I)*C(K,J)
	      END DO
	    END DO
	  END DO
	CALL TUNE(2,'M3')
	WRITE(6,*)A(5,5),A(N,N)
!
C
	CALL TUNE(3,' ')
C
	STOP
	END
C
C Routine to replace VMS TUNE routine for collecting tuning statistics for 
C a program section.
C
C Usage:
C          CALL TUNE(1,'Unique ID')
C              .......
C              Section of code.
C              .......
C          CALL TUNE(2,'Same Unique ID')
C
C To print results:
C
C           CALL TUNE(3,' ')
C
C Output from print:
C
C           CPUTOT() : CPU time accumulated in each IDENT
C           WALLTOT() : Total wall time bewteen each IDENT.
C
C NB: Total time for each code section is accumulated on successive calls.
C i.e. TUNE(1,'Unique ID') does not initialize the counters.
C
        SUBROUTINE TUNE(LRUN,IDENT)
        IMPLICIT NONE
!
! Altered 11-Nov-2000 : Call to F90 SYSTEM_CLOCK routine implemented.
!                       Wall time now returened as in original VMS routine.
!                       Counters now initialized if LRUN=0 is passed.
!
	INTEGER LRUN
	CHARACTER*(*) IDENT
!
	INTEGER, PARAMETER :: MAX_IDS=50
!	INTEGER, PARAMETER :: LUOUT=55
	INTEGER, PARAMETER :: LUOUT=6
!
        REAL*8 T0,OVERHEAD
        REAL*8 ST_CPU(MAX_IDS)
	REAL*8 CPUTOT(MAX_IDS)
        REAL*8 WALLTOT(MAX_IDS)
        CHARACTER*30 IDLIST(MAX_IDS)
        INTEGER I
!
	INTEGER*8 IEND_WALL
	INTEGER*8 IC0,IR0,IM0,IT1
	INTEGER*8 IST_WALL(MAX_IDS)
	REAL*8 CLK_PERIOD,RR0
!
	REAL*4 ETIME,TARRY(2)
!
	LOGICAL*4 FIRSTTIME
	DATA FIRSTTIME/.TRUE./
        SAVE FIRSTTIME,OVERHEAD
        SAVE ST_CPU,IST_WALL,CPUTOT,WALLTOT
	SAVE IC0,IR0,IM0,RR0
        SAVE IDLIST

        IF (FIRSTTIME)THEN
          FIRSTTIME=.FALSE.
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.D0
            IST_WALL(I)=0.D0
            CPUTOT(I)=0.D0
            WALLTOT(I)=0.D0
	    IDLIST(I)=' '
          END DO
	  CALL SYSTEM_CLOCK(IC0,IR0,IM0);    RR0=IR0
          T0=ETIME(TARRY)
          OVERHEAD=2.0D0*(ETIME(TARRY)-T0)
!	  OPEN(UNIT=LUOUT,STATUS='REPLACE',FILE='TIMING')
        ENDIF
!
! If LRUN =1, we are beginning the TIME bracket. Therefore we find the
! correct storage location first. 
!
	IF (LRUN .EQ. 1) THEN
	  DO I=1,MAX_IDS
            IF (IDENT .EQ. IDLIST(I))THEN
	      CALL SYSTEM_CLOCK(IST_WALL(I))
	      ST_CPU(I)=ETIME(TARRY)
	      RETURN	
	    END IF
	    IF (IDLIST(I) .EQ. ' ') THEN
	      IDLIST(I)=IDENT
	      CALL SYSTEM_CLOCK(IST_WALL(I))
	      ST_CPU(I)=ETIME(TARRY)
	      RETURN	
	    END IF
	  END DO
	  WRITE (LUOUT,'(A)')' ***** TOO MANY TUNING POINTS '
	  RETURN
!
	ELSE IF (LRUN .EQ. 2) THEN
!
! If LRUN=2, we are ending the TIME bracket. Therefore we call the timing 
! routine first.
!
          T0=ETIME(TARRY)
	  CALL SYSTEM_CLOCK(IEND_WALL)
          DO I=1,MAX_IDS
	    IF (IDENT.EQ.IDLIST(I))THEN
	      CPUTOT(I)=CPUTOT(I)+(T0-ST_CPU(I)-OVERHEAD)
	      IT1=IEND_WALL-IST_WALL(I)
	      IF(IT1 .LT. 0)IT1=IT1+IM0
	      WALLTOT(I)=WALLTOT(I)+IT1/RR0
	      RETURN
	    END IF 
	    IF (IDLIST(I).EQ.' ')EXIT
	  END DO
	  WRITE(LUOUT,*)' ***** UNMATCHED TUNING POINT '
	  RETURN
C
	ELSE IF (LRUN .EQ. 3) THEN
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,*)'Overhead is ',OVERHEAD
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,'(8X,''Identifier'',11x,''Elapsed'',11x,''  CPU'')')
	  WRITE(LUOUT,'(29x,''  Time '',11x,''  Time'')')
	  DO I=1,MAX_IDS
	    IF (IDLIST(I).EQ.' ') EXIT
	    WRITE(LUOUT,'(1X,a24,f15.6,2X,f15.6)')
	1   IDLIST(I),WALLTOT(I),CPUTOT(I)
          END DO
	ELSE IF(LRUN .EQ. 0) THEN
          DO  I=1,MAX_IDS
            ST_CPU(I)=0.D0
            IST_WALL(I)=0.D0
            CPUTOT(I)=0.D0
            WALLTOT(I)=0.D0
	    IDLIST(I)=' '
          END DO
	ELSE
	  WRITE (LUOUT,'(A)')' ***** ILLEGAL VALUE OF LRUN IN CALL TO TUNE '
	  WRITE(LUOUT,*)' LRUN=',LRUN
	  STOP
	ENDIF
        
	RETURN
	END
