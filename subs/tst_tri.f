	PROGRAM TST_TRI
	IMPLICIT NONE
C
	INTEGER N,NRHS,I,J
	PARAMETER (N=60)
	PARAMETER (NRHS=10000)
	REAL*8 A(N),B(N),C(N),D(N),F(N)
	REAL*8 X(N),W(N,NRHS),GAM(N),T1
!
	INTEGER INFO,IPIV(N)
!
C
	DO I=1,N
	  T1=SQRT(FLOAT(I))
	  A(I)=-0.48*T1
	  C(I)=-0.49*T1
	  B(I)=T1-A(I)-C(I)
	  X(I)=SQRT(2*I+0.5)
	  D(I)=0.0D0
	END DO
	A(1)=0.0
	C(N)=0.0
C
	D(1)=D(1)+B(1)*X(1)+C(1)*X(2)
	DO I=2,N-1
	  D(I)=D(I)+A(I)*X(I-1)+B(I)*X(I)+C(I)*X(I+1)
	END DO
	D(N)=D(N)+A(N)*X(N-1)+B(N)*X(N)
!
	DO J=1,NRHS
	  W(:,J)=D(:)+0.1*(J-1)
	END DO
!
	CALL TUNE(1,'LAP')
	  CALL DGTTRF(N,A(2),B,C,F,IPIV,INFO)
	  CALL DGTTRS('N',N,NRHS,A(2),B,C,F,IPIV,W,N,INFO)
	CALL TUNE(2,'LAP')
!
	DO I=1,N
	  WRITE(2,*)X(I),W(I,1),W(I,NRHS)
	END DO
C
C
	DO I=1,N
	  T1=SQRT(FLOAT(I))
	  A(I)=-0.48*T1
	  C(I)=-0.49*T1
	  B(I)=T1-A(I)-C(I)
	  X(I)=SQRT(2*I+0.5)
	  D(I)=0.0D0
	END DO
	A(1)=0.0
	C(N)=0.0
C
	D(1)=D(1)+B(1)*X(1)+C(1)*X(2)
	DO I=2,N-1
	  D(I)=D(I)+A(I)*X(I-1)+B(I)*X(I)+C(I)*X(I+1)
	END DO
	D(N)=D(N)+A(N)*X(N-1)+B(N)*X(N)
!
	DO J=1,NRHS
	  W(:,J)=D(:)+0.1*(J-1)
	END DO
!
	CALL TUNE(1,'TH')
	  CALL THOMAS(A,B,C,W,N,NRHS)
	CALL TUNE(2,'TH')
	CALL TUNE(1,'NO_ROUT')
	CALL TUNE(2,'NO_ROUT')
	CALL TUNE(3,' ')
!C
!C For Rybicki Hummer algorithm.
!C
!	DO I=1,N
!	  A(I)=-A(I)
!	  C(I)=-C(I)
!	  B(I)=B(I)-A(I)-C(I)
!	END DO
!	I=1
!	CALL TRI_RH_1V(A,B,C,D,N,I)
C
	DO I=1,N
	  WRITE(2,*)X(I),W(I,1),W(I,NRHS)
	END DO
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
	  WRITE(LUOUT,*)' '
	  WRITE(LUOUT,*)'Overhead is ',OVERHEAD
	  WRITE(LUOUT,*)' '
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
