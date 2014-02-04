	PROGRAM PLT_NRR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: ND_MAX=100
!
	REAL*8 R(ND_MAX)
	REAL*8 T(ND_MAX)
	REAL*8 ED(ND_MAX)
	REAL*8 DI(ND_MAX)
!
	REAL*8 NET(ND_MAX)
	REAL*8 RATE(ND_MAX)
!
	CHARACTER*80 FILENAME
	CHARACTER*80 STRING
!
	INTEGER I,IBEG
	INTEGER ND
	LOGICAL PLT_RATE
!
	PLT_RATE=.FALSE.
	CALL GEN_IN(ND,'Number of depth points in model')
	CALL GEN_IN(PLT_RATE,'Plot Recombination rates?')
!
	DO WHILE(1 .EQ. 1)
 	  FILENAME=' '
	  CALL GEN_IN(FILENAME,'File with data to be plotted')
	  IF(FILENAME .EQ. ' ')GOTO 1000
!
	  OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	    DO IBEG=1,ND,10
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(R(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(T(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(ED(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(DI(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING	!Photoiization
	      DO WHILE(1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(STRING .EQ. ' ')EXIT
	      END DO
	     
!
	      READ(11,'(A)')STRING      !Colisional ionization
	      DO WHILE(1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(STRING .EQ. ' ')EXIT
	      END DO
!
	      READ(11,'(A)')STRING      !Recombination
	      WRITE(6,*)TRIM(STRING)
	      DO WHILE(1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(STRING .EQ. ' ')EXIT
	      END DO
!
	      READ(11,'(A)')STRING      !Collisional recombination
	      DO WHILE(1 .EQ. 1)
	        READ(11,'(A)')STRING
	        IF(STRING .EQ. ' ')EXIT
	      END DO
!
	      READ(11,'(A)')STRING
	      READ(11,*)(NET(I),I=IBEG,MIN(IBEG+9,ND))
!
	      READ(11,'(A)')STRING
	      READ(11,'(A)')STRING
	      READ(11,*)(RATE(I),I=IBEG,MIN(IBEG+9,ND))
	      IF(IBEG+9 .LT. ND)READ(11,'(A)')STRING
	    END DO
	  CLOSE(UNIT=11)
!
	  IF(PLT_RATE)THEN
	    ED(1:ND)=DLOG10(ED(1:ND))
	    RATE(1:ND)=DLOG10(RATE(1:ND))+12.0D0
	    CALL DP_CURVE(ND,ED,RATE)
	  ELSE
	    ED(1:ND)=DLOG10(ED(1:ND))
	    NET(1:ND)=DLOG10( ABS(NET(1:ND))+1.0E-07)
	    CALL DP_CURVE(ND,ED,NET)
	  END IF
	END DO
!
1000	CONTINUE
	IF(PLT_RATE)THEN
	  CALL GRAMON_PGPLOT('Log(Ne)','Log(10\u12\d RR)',' ',' ')
	ELSE
	  CALL GRAMON_PGPLOT('Log(Ne)','Log(Net % RR)',' ',' ')
	END IF
!
	STOP
	END
