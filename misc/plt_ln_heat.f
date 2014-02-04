	PROGRAM PLT_LN_HEAT
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER, PARAMETER :: NMAX=400000
!
	REAL*8 NU(NMAX)
	REAL*8 LAM(NMAX)
	REAL*8 XV(NMAX)
	REAL*8 Y(NMAX)
	CHARACTER*40 NAME(NMAX)
!
	REAL*8, ALLOCATABLE :: LH(:,:)
	REAL*8, ALLOCATABLE :: SE_SCL(:,:)
	REAL*8, ALLOCATABLE :: SE_NOSCL(:,:)
!
	INTEGER ND
	INTEGER N_LINES
	INTEGER COUNT
	INTEGER I,K,ML,IBEG
	INTEGER IOS
	LOGICAL FILE_OPEN
!
	CHARACTER*80 FILENAME
	CHARACTER*80 STRING
        CHARACTER*30 UC
        EXTERNAL UC
	CHARACTER*20 PLT_OPT
!
	WRITE(6,'(A)')' '
	WRITE(6,'(A)')' Program to plot the radiative equilibrium equation.'
	WRITE(6,'(A)')' Can be used to show how RE changes as we integrate from blue to red.'
	WRITE(6,'(A)')' Designed to see effect of scaling the heating rates.'
	WRITE(6,'(A)')' Can be used to identify SL assignments that might be changed.'
	WRITE(6,'(A)')' '
!
100	CONTINUE
 	FILENAME='LINEHEAT'
	CALL GEN_IN(FILENAME,'File with data to be plotted')
	IF(FILENAME .EQ. ' ')GOTO 100
!
	OPEN(UNIT=20,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(20,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Number of depth points') .NE. 0)THEN
	         READ(STRING,*)ND
	         WRITE(6,'(A,I4)')' Number of depth points in the model is:',ND
	         EXIT
	      END IF
	    END DO
	  END IF
	INQUIRE(UNIT=20,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=20)
!
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(ND,'Number of data points (must be exact)')
	END IF
	N_LINES=NMAX
	CALL GEN_IN(N_LINES,'Maximum number of lines to be read')
!
	ALLOCATE (LH(ND,N_LINES))
	ALLOCATE (SE_SCL(ND,N_LINES))
	ALLOCATE (SE_NOSCL(ND,N_LINES))
!
	OPEN(UNIT=11,FILE=FILENAME,STATUS='OLD',ACTION='READ')
!
	COUNT=0
	DO ML=1,N_LINES
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ')
	    READ(11,'(A)',END=5000)STRING
	  END DO
	  IF(INDEX(STRING,'error in L due to') .NE. 0)GOTO 5000
	  STRING=ADJUSTL(STRING)
	  K=INDEX(STRING,' ')
	  STRING=ADJUSTL(STRING(K:))
	  K=INDEX(STRING,' ')
	  NAME(ML)=STRING(1:K-1) 
	  K=INDEX(STRING,')')
	  DO WHILE(K .NE. 0)
	    STRING(1:)=STRING(K+1:)
	    K=INDEX(STRING,')')
	  END DO
	  K=INDEX(STRING,'  ')
	  STRING(1:)=STRING(K:)
	  READ(STRING,*,ERR=200)NU(ML)
	  READ(11,*)LH(1:ND,ML)
	  READ(11,'(A)')STRING
	  IF(STRING .NE. ' ')THEN
	    WRITE(6,'(A)')' '
	    WRITE(6,*)'Error: invalid data format'
	    WRITE(6,*)'Check ND values'
	    WRITE(6,*)'Current line count is',COUNT
	    WRITE(6,'(A)')' '
	    STOP
	  END IF
	  READ(11,*)SE_SCL(1:ND,ML)
	  READ(11,*)SE_NOSCL(1:ND,ML)
	  COUNT=COUNT+1
	END DO
5000	CONTINUE
	N_LINES=COUNT
!
	LAM(1:N_LINES)=2.99702458D+03/NU(1:N_LINES)
!
	PLT_OPT='SS'
	DO WHILE(1 .EQ. 1)
	  WRITE(6,*)' '
	  WRITE(6,*)'Plot options are:'
	  WRITE(6,*)' LH:    Plot LH  at given depth'
	  WRITE(6,*)' SS:    Plot STEQ (scaling) at given depth'
	  WRITE(6,*)' SN:    Plot STEQ (no scaling) at a given depth'
	  WRITE(6,*)' FV:    Plot final values (SCL, NO SCL) as a function of depth'
	  WRITE(6,*)' WR:    Write SS & SN data at a single depth to file'
	  WRITE(6,*)' E(X):  Exit routine'
	  CALL GEN_IN(PLT_OPT,'Plot option')
	  IF(UC(PLT_OPT) .EQ. 'LH')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=LH(K,1:N_LINES)
	    CALL DP_CURVE(N_LINES,NU,Y)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','LH',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'SS')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SE_SCL(K,1:N_LINES)
	    WRITE(6,*)Y(N_LINES)
	    CALL DP_CURVE(N_LINES,NU,Y)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','STEQ(scaled)',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'SN')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    Y(1:N_LINES)=SE_NOSCL(K,1:N_LINES)
	    CALL DP_CURVE(N_LINES,NU,Y)
	    CALL GRAMON_PGPLOT('\gn(10\u15 \dHz)','STEQ(scaled)',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'FV')THEN
	    DO I=1,ND
	      XV(I)=I
	    END DO
	    Y(1:ND)=SE_SCL(1:ND,N_LINES)
	    CALL DP_CURVE(ND,XV,Y)
	    Y(1:ND)=SE_NOSCL(1:ND,N_LINES)
	    CALL DP_CURVE(ND,XV,Y)
	    CALL GRAMON_PGPLOT('Depth Index','STEQ(sc,scl)',' ',' ')
	  ELSE IF(UC(PLT_OPT) .EQ. 'WR')THEN
	    K=ND
	    CALL GEN_IN(K,'Depth for plotting')
	    DO I=1,N_LINES
	      WRITE(100,'(I7,4ES14.4,3X,A)')I,NU(I),SE_NOSCL(K,I)-SE_SCL(K,I),
	1                SE_NOSCL(K,I),SE_SCL(K,I),TRIM(NAME(I))
	    END DO
	  ELSE IF(UC(PLT_OPT) .EQ. 'EX' .OR. UC(PLT_OPT) .EQ. 'E')THEN
	    STOP
	  END IF
	END DO
!
200	WRITE(6,*)STRING
	STOP
!
	END
