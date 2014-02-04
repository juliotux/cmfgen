	PROGRAM TST_TUNE
	IMPLICIT NONE
	INTEGER, PARAMETER :: NMAX=1000000
!
	REAL*8 A(NMAX)
	REAL*8 B(NMAX)
	REAL*8 SUM
!
	INTEGER I,J
!
	DO I=1,NMAX
	  A(I)=5.0D0*I/FLOAT(NMAX)
	  B(I)=2.0D0*I/FLOAT(NMAX)
	END DO
!
	CALL TUNE(1,'OUTER')
!
	CALL TUNE(1,'NULL')
	CALL TUNE(2,'NULL')
!
	CALL TUNE(1,'HOPE')
	  DO I=1,NMAX
	    A(I)=A(I)*B(I)
	  END DO
	CALL TUNE(2,'HOPE')
!
	WRITE(6,*)'Pause a little so offset WALL and CPU (1 to cont.)'
	READ(5,*)J
!
        CALL TUNE(1,'SEC')
        SUM=0.0D0
        DO J=1,10
          DO I=1,NMAX
            A(I)=SQRT(A(I))
          END DO
          DO I=1,NMAX
            SUM=SUM+A(I)
          END DO
        END DO
        CALL TUNE(2,'SEC')
        CALL TUNE(2,'OUTER')
        WRITE(6,*)SUM
        CALL TUNE(3,' ')
!
	STOP
	END
