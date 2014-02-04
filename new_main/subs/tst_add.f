	PROGRAM TST_ADD
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: IONE=1
!
	INTEGER NT,NUM_BNDS,ND,NVEC
	INTEGER I,J,K,L
!
	REAL*8, ALLOCATABLE :: VJ_R(:,:,:)
	REAL*8, ALLOCATABLE :: VJ_P(:,:,:)
	REAL*8, ALLOCATABLE :: VJ(:,:,:)
!
	LOGICAL, ALLOCATABLE :: NON_ZERO_VJ(:)
!
	REAL*8 T1,T2
!
	NT=1000; ND=60; NUM_BNDS=3
	CALL GEN_IN(NT,'NT')	
	CALL GEN_IN(NUM_BNDS,'NUM_BNDS')	
	CALL GEN_IN(ND,'ND')
!
	ALLOCATE (VJ_R(NT,NUM_BNDS,ND))	
	ALLOCATE (VJ_P(NT,NUM_BNDS,ND))	
	ALLOCATE (VJ(NT,NUM_BNDS,ND))
	ALLOCATE (NON_ZERO_VJ(NT))
!
	DO L=1,ND
	  DO J=1,NUM_BNDS
	    DO I=1,NT
	      VJ_R(I,J,L)=I+(J-1)*(L-1)
	      VJ_P(I,J,L)=I+3*(J-1)*(L-1)
	      VJ(I,J,L)=I+4*(J-1)*(L-1)
	    END DO
	  END DO
	END DO
	NON_ZERO_VJ(:)=.FALSE.
	NON_ZERO_VJ(50)=.TRUE.
	NON_ZERO_VJ(200)=.TRUE.
	NON_ZERO_VJ(300)=.TRUE.
!
	T2=2.345D00
	T1=4.345D00
	CALL TUNE(1,'STANDARD')
	DO K=1,100
	  VJ(10,2,30)=1.1*VJ(10,2,30)
	  DO L=1,ND
	    T1=4.345D00+0.001*(L-1)
	     DO J=1,NUM_BNDS
	       DO I=1,NT
	           VJ_P(I,J,L)=VJ_P(I,J,L)+T2*VJ(I,J,L)
	           VJ_R(I,J,L)=VJ_R(I,J,L)+T1*VJ(I,J,L)
	       END DO
	     END DO
	  END DO
 	END DO
	CALL TUNE(2,'STANDARD')
!
	T2=2.345D00
	T1=4.345D00
	CALL TUNE(1,'STAND_SPLIT')
	DO K=1,100
	  VJ(10,2,30)=1.1*VJ(10,2,30)
	  DO L=1,ND
	     DO J=1,NUM_BNDS
	       DO I=1,NT
	          VJ_P(I,J,L)=VJ_P(I,J,L)+T2*VJ(I,J,L)
	       END DO
	     END DO
	  END DO
 	END DO
!
	DO K=1,100
	  VJ(10,2,30)=1.1*VJ(10,2,30)
	  DO L=1,ND
	    T1=4.345D00+0.001*(L-1)
	     DO J=1,NUM_BNDS
	       DO I=1,NT
	           VJ_R(I,J,L)=VJ_R(I,J,L)+T1*VJ(I,J,L)
	       END DO
	     END DO
	  END DO
 	END DO
	CALL TUNE(2,'STAND_SPLIT')
	DO L=1,ND
	  DO J=1,NUM_BNDS
	    DO I=1,NT
	      VJ_R(I,J,L)=I+(J-1)*(L-1)
	      VJ_P(I,J,L)=I+3*(J-1)*(L-1)
	      VJ(I,J,L)=I+4*(J-1)*(L-1)
	    END DO
	  END DO
	END DO
!
	DO L=1,ND
	  DO J=1,NUM_BNDS
	    DO I=1,NT
	      VJ_R(I,J,L)=I+(J-1)*(L-1)
	      VJ_P(I,J,L)=I+3*(J-1)*(L-1)
	      VJ(I,J,L)=I+4*(J-1)*(L-1)
	    END DO
	  END DO
	END DO
!
	T2=2.345D00
	T1=4.345D00
	CALL TUNE(1,'STAND3')
	DO K=1,100
	  VJ(10,2,30)=1.1*VJ(10,2,30)
	  DO L=1,ND
	    T1=4.345D00+0.001*(L-1)
	    VJ_P(:,:,L)=VJ_P(:,:,L)+T2*VJ(:,:,L)
	    VJ_R(:,:,L)=VJ_R(:,:,L)+T1*VJ(:,:,L)
	  END DO
 	END DO
	CALL TUNE(2,'STAND3')
!
	DO L=1,ND
	  DO J=1,NUM_BNDS
	    DO I=1,NT
	      VJ_R(I,J,L)=I+(J-1)*(L-1)
	      VJ_P(I,J,L)=I+3*(J-1)*(L-1)
	      VJ(I,J,L)=I+4*(J-1)*(L-1)
	    END DO
	  END DO
	END DO
!
	T2=2.345D00
	T1=4.345D00
	CALL TUNE(1,'BLAS')
	NVEC=NT*ND*NUM_BNDS
	DO K=1,100
	   CALL DAXPY(NVEC,T2,VJ,IONE,VJ_P,IONE)
	   CALL DAXPY(NVEC,T1,VJ,IONE,VJ_R,IONE)
	END DO
	CALL TUNE(2,'BLAS')
!
	T2=2.345D00
	CALL TUNE(1,'STAND2')
	DO K=1,100
	  DO L=1,ND
	     T2=2.345D0+(L-1)*0.001D0
	     DO J=1,NUM_BNDS
	       DO I=1,NT
	         VJ_P(I,J,L)=VJ_P(I,J,L)+T2*VJ(I,J,L)
	       END DO
	     END DO
	  END DO
	END DO
	CALL TUNE(2,'STAND2')
!
	DO L=1,ND
	  DO J=1,NUM_BNDS
	    DO I=1,NT
	      VJ_R(I,J,L)=I+(J-1)*(L-1)
	      VJ_P(I,J,L)=I+3*(J-1)*(L-1)
	      VJ(I,J,L)=I+4*(J-1)*(L-1)
	    END DO
	  END DO
	END DO
!
	T2=2.345D00
	CALL TUNE(1,'BLAS2')
	NVEC=NT*NUM_BNDS
	DO K=1,100
	   DO L=1,ND
	     T2=2.345D0+(L-1)*0.001D0
	     CALL DAXPY(NVEC,T2,VJ(1,1,L),IONE,VJ_P(1,1,L),IONE)
	  END DO
	END DO
	CALL TUNE(2,'BLAS2')

	CALL TUNE(3,' ')
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
	  CLOSE(UNIT=LUOUT)
	  OPEN(UNIT=LUOUT,STATUS='OLD',FILE='TIMING',POSITION='APPEND')
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
