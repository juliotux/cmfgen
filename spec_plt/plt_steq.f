!
! Simple program to plot data from the STEQ_VALS file. This can be used to check
!   convergence of RATE equations etc. It can also plot the corrections (which 
!   could also be plotted using PLT_SCR).
!
	PROGRAM PLT_STEQ
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 02-Jul-2014 - Bug fix -- last iteration was read was not being counted as an 
!                        iteration even though it was successfully read in.
! Created 14-Jan-2014
!
	TYPE STEQ_DATA
	  REAL*8, ALLOCATABLE :: STEQ(:,:)
	  REAL*8, ALLOCATABLE :: SOL(:,:)
	  REAL*8, ALLOCATABLE :: RE(:)
	END TYPE STEQ_DATA
!
	TYPE (STEQ_DATA) ST(200)
!
	REAL*8, ALLOCATABLE :: XV(:)
	REAL*8, ALLOCATABLE :: YV(:)
	REAL*8, ALLOCATABLE :: ZV(:)
!
	INTEGER, PARAMETER :: T_OUT=6
!
	INTEGER I,J,K,L
	INTEGER IT,ID,IV
	INTEGER ND,NT,NIT
	INTEGER IOS
!
	LOGICAL LOG_Y_AXIS
	LOGICAL DO_ABS
	LOGICAL NORMALIZE
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	CHARACTER*10 PLT_OPT
	CHARACTER*80 XLABEL
	CHARACTER*80 YLABEL
	CHARACTER(LEN=132) STRING
!
	REAL*8 T1,T2,T3
!
	XLABEL=' '
	YLABEL=' '
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the model directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       STEQ_VALS'
	WRITE(T_OUT,*)'                                           MODEL'
	WRITE(T_OUT,*)' '
!
	OPEN(UNIT=12,FILE='MODEL',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of depth') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)ND
	  DO WHILE(INDEX(STRING,'!Total number of variables') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)NT
C
100	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to read MODEL file'
	  CALL GEN_IN(NT,'Total number of levels')
	  CALL GEN_IN(ND,'Number of depth points')
	END IF 
	CLOSE(UNIT=12)
!
	OPEN(UNIT=12,FILE='STEQ_VALS',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  DO IT=1,200
	    IF(MOD(IT,10) .EQ. 1)THEN
	      IF(IT .NE. 1)WRITE(6,*)' '
	      WRITE(6,'(A,X,I4)',ADVANCE='NO')'Reading iteration',IT
	    ELSE
	      WRITE(6,'(2X,I4)',ADVANCE='NO')IT
	    END IF
	    STRING=' '
	    DO WHILE( INDEX(STRING,'Radiative') .EQ. 0)
	      READ(12,'(A)',ERR=5000,END=5000)STRING
	    END DO
	    ALLOCATE (ST(IT)%RE(ND))
	    READ(12,*,ERR=5000,END=5000)ST(IT)%RE
!
	    STRING=' '
	    DO WHILE( INDEX(STRING,'STEQ_ARRAY') .EQ. 0)
	      READ(12,'(A)',ERR=5000,END=5000)STRING
	    END DO
!
	    ALLOCATE (ST(IT)%STEQ(NT,ND))
	    READ(12,'(A)',ERR=5000,END=5000)STRING
	    DO L=1,ND,10
	      DO I=1,NT
	        READ(12,'(A)',ERR=5000,END=5000)STRING
	        K=INDEX(STRING,'*')+1
	        READ(STRING(K:),*,END=5000)(ST(IT)%STEQ(I,J),J=L,MIN(L+9,ND))
	      END DO
	      READ(12,'(A)',ERR=5000,END=5000)STRING
	    END DO
!
	    STRING=' '
	    DO WHILE( INDEX(STRING,'SOLUTION') .EQ. 0)
	      READ(12,'(A)',ERR=5000,END=5000)STRING
	    END DO
!
	    ALLOCATE (ST(IT)%SOL(NT,ND))
	    READ(12,'(A)',ERR=5000,END=5000)STRING
	    DO L=1,ND,10
	      DO I=1,NT
	        READ(12,'(A)',ERR=5000,END=5000)STRING
	        K=INDEX(STRING,'#')+1
	        READ(STRING(K:),*,END=5000)(ST(IT)%SOL(I,J),J=L,MIN(L+9,ND))
	      END DO
	      READ(12,'(A)',ERR=5000,END=6000)STRING
	    END DO
	    NIT=IT
	  END DO
6000	NIT=NIT+1		!We add one as we reach this point when trying to read end of file
5000	CONTINUE
	CLOSE(UNIT=12)
	WRITE(6,*)' '
	IF(NIT .EQ. 0)THEN
	  WRITE(6,*)'No full iteration successfully read from STEQ FILE'
	  STOP
	END IF
!
	J=MAX(NIT,ND); J=MAX(NT,J)
	ALLOCATE (XV(J))
	ALLOCATE (YV(J))
	ALLOCATE (ZV(J))
!
	WRITE(6,*)'   Number of depth points is:',ND
	WRITE(6,*)'Number of variables/depth is:',NT
	WRITE(6,*)'     Number of iterations is:',NIT
!
200	CONTINUE
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'ND  :: ',ND
	WRITE(T_OUT,*)'NT  :: ',NT
	WRITE(T_OUT,*)'NIT :: ',NIT
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)' '
        WRITE(T_OUT,*)'CIV_D  :: Plot correction for a given iteraton and variable as a function of depth index.'
!
        WRITE(T_OUT,*)'RIV_D  :: Plot rate for a given iteraton and variable as a function of depth index.'
        WRITE(T_OUT,*)'RDI_V  :: Plot rate for a given iteraton and depth as a function of variable.'
        WRITE(T_OUT,*)'RDV_I  :: Plot rate for a given variable and depth as a function of iteration.'
!
        WRITE(T_OUT,*)'RATD_V  :: Plot ratio of rate eqution at a given depth as a function of variable'
        WRITE(T_OUT,*)'RATV_D  :: Plot ratio of rate eqution for a given variable as a function of depth'
!
	WRITE(T_OUT,*)' '
!	WRITE(T_OUT,*)'LY  :: Switch to/from Log(Y) for options where appropriate (not full implemented)'
!	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'E   :: EXIT'
	WRITE(T_OUT,*)' '
	PLT_OPT='P'
	CALL GEN_IN(PLT_OPT,'Plot option: (E=exit)')
	CALL SET_CASE_UP(PLT_OPT,0,0)
!
	IF(PLT_OPT(1:2) .EQ. 'E ' .OR. PLT_OPT(1:2) .EQ. 'EX')STOP
!
	IF(PLT_OPT(1:2) .EQ. 'LY')THEN
	   WRITE(6,*)'LY option has not been implemented'
!	   IF(LOG_Y_AXIS)WRITE(6,*)'Switching to linear Y axis'
!	   IF(.NOT. LOG_Y_AXIS)WRITE(6,*)'Switching to logarithmic Y axis'
!	   LOG_Y_AXIS=.NOT. LOG_Y_AXIS
!	   GOTO 200
!
	ELSE IF(PLT_OPT(1:5) .EQ. 'CIV_D')THEN
	  CALL GEN_IN(IT,'Iteration number')
	  CALL GEN_IN(IV,'Variable')
	  YV(1:ND)=ST(IT)%SOL(IV,:)
	  DO I=1,ND
	    XV(I)=I
	  END DO
	  CALL DP_CURVE(ND,XV,YV) 
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:5) .EQ. 'RIV_D')THEN
	  CALL GEN_IN(IT,'Iteration number')
	  CALL GEN_IN(IV,'Variable')
	  YV(1:ND)=SIGN(ABS(ST(IT)%STEQ(IV,:))**0.25,ST(IT)%STEQ(IV,:))
	  DO I=1,ND
	    XV(I)=I
	  END DO
	  CALL DP_CURVE(ND,XV,YV) 
	  XLABEL='Depth index'
	  YLABEL='dN/dt\u1/4\d'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:5) .EQ. 'RDV_I')THEN
	  ID=5; IV=NT; NORMALIZE=.TRUE.
	  CALL GEN_IN(ID,'Depth')
	  CALL GEN_IN(IV,'Variable')
	  CALL GEN_IN(NORMALIZE,'Normalize plot by maximum value')
	  IF(NORMALIZE)THEN
	   T1=0
	   IF(IV .EQ. NT)THEN
	      DO I=1,NIT; T1=MAX(T1,ABS(ST(I)%RE(ID))); END DO
	      DO I=1,NIT; YV(I)=ST(I)%RE(ID)/T1;        END DO
	    ELSE
	      DO I=1,NIT; T1=MAX(T1,ABS(ST(I)%STEQ(IV,ID))); END DO
	      DO I=1,NIT; YV(I)=ST(I)%STEQ(IV,ID)/T1;        END DO
	    END IF
	    YLABEL='(dN/dt)/Max(|dN/dt|)'
	  ELSE
	   IF(IV .EQ. NT)THEN
	      DO I=1,NIT
	        YV(I)=SIGN(ABS(ST(I)%RE(ID))**0.25,ST(I)%RE(ID))
	      END DO
	    ELSE
	      DO I=1,NIT
	        YV(I)=SIGN(ABS(ST(I)%STEQ(IV,ID))**0.25,ST(I)%STEQ(IV,ID))
	      END DO
	    END IF
	    YLABEL='dN/dt\u1/4\d'
	  END IF
	  DO I=1,NIT
	    XV(I)=I
	  END DO
	  XLABEL='Iteration'
	  CALL DP_CURVE(NIT,XV,YV) 
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:5) .EQ. 'RDI_V')THEN
	  ID=1; IT=NIT
	  CALL GEN_IN(ID,'Depth')
	  CALL GEN_IN(IT,'Iteration number')
	  YV(:)=SIGN(ABS(ST(IT)%STEQ(:,ID))**0.25,ST(IT)%STEQ(:,ID))
	  DO I=1,NT
	    XV(I)=I
	  END DO
	  CALL DP_CURVE(NT,XV,YV) 
	  XLABEL='Variable'
	  YLABEL='dN/dt\u1/4\d'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:7) .EQ. 'RATD_V')THEN
	  ID=1; IT=1; K=NIT
	  CALL GEN_IN(ID,'Depth')
	  CALL GEN_IN(IT,'Iteration number of numerator')
	  CALL GEN_IN(K, 'Iteration number of divisor')
	  DO I=1,NT-1
	    IF(ST(K)%STEQ(I,ID) .NE. 0)THEN
	      YV(I)=ST(IT)%STEQ(I,ID)/ST(K)%STEQ(I,ID)
	    ELSE
	      YV(I)=1.0D+10
	    END IF
	  END DO
	  YV(NT)=ST(IT)%RE(ID)/ST(K)%RE(ID)
	  WRITE(6,*)'Any divisions by zero have been set to 10^10'
	  DO I=1,NT
	    XV(I)=I
	  END DO
	  CALL DP_CURVE(NT,XV,YV) 
	  XLABEL='Variable'
	  YLABEL='Rate equation ratio'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:6) .EQ. 'RATV_D')THEN
	  IV=1; IT=1; K=NIT
	  CALL GEN_IN(IV,'Variable')
	  CALL GEN_IN(IT,'Iteration number of numerator')
	  CALL GEN_IN(K, 'Iteration number of divisor')
	  IF(IV .EQ. NT)THEN
	    DO I=1,ND
	      IF(ST(K)%STEQ(IV,I) .NE. 0)THEN
	        YV(I)=ST(IT)%STEQ(IV,I)/ST(K)%STEQ(IV,I)
	      ELSE
	        YV(I)=1.0D+10
	      END IF
	    END DO
	  ELSE
	    DO I=1,ND
	      YV(I)=ST(IT)%RE(I)/ST(K)%RE(I)
	    END DO
	  END IF
	  WRITE(6,*)'Any divisions by zero have been set to 10^10'
	  DO I=1,ND
	    XV(I)=I
	  END DO
	  CALL DP_CURVE(ND,XV,YV) 
	  XLABEL='Depth index'
	  YLABEL='Rate equation ratio'
	  GOTO 200
	ELSE
	  CALL GRAMON_PGPLOT(XLABEL,YLABEL,' ',' ')
	  GOTO 200
	END IF
!
	END
