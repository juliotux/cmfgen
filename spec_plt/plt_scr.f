!
! Simple program to plot data from the SCRTEMP file. This can be used to check
!   Convergence of a program.
!   Convergence as a function of depth etc.
!
	PROGRAM PLT_SCR
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
	INTEGER ND,NT,NIT
!
! Altered 26-Feb-2015: Fixed error with FDG ouput if called multiple times.
! Altered 30-Mar-2014: Altered FDG option so as to print out adjacent values.
! Altered 11-Mar-2014: Installed INT option.
! Altered 25-Feb-2014: Modified WRST option.
! Altered 12-Dec-2013: Installed LY option so that we can plot very small populations on
!                        a logarithmic scale.
!
	REAL*8, ALLOCATABLE :: POPS(:,:,:)		!NT,ND,NIT
	REAL*8, ALLOCATABLE :: R_MAT(:,:)		!ND,NIT
	REAL*8, ALLOCATABLE :: V_MAT(:,:)		!ND,NIT
	REAL*8, ALLOCATABLE :: RAT(:,:)			!NT,ND
!
	REAL*8, ALLOCATABLE :: R(:)			!ND
	REAL*8, ALLOCATABLE :: V(:)			!ND
	REAL*8, ALLOCATABLE :: SIGMA(:)			!ND
C
	REAL*8, ALLOCATABLE :: X(:)			!NIT
	REAL*8, ALLOCATABLE :: Y(:)			!NIT
	REAL*8, ALLOCATABLE :: Z(:)			!NIT
!
	INTEGER, ALLOCATABLE :: I_BIG(:)		!NT
	REAL*8, ALLOCATABLE :: Z_BIG(:)			!NT
C
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: T_OUT=6
C
	INTEGER, SAVE :: FDG_COUNTER=0
	INTEGER, PARAMETER :: NLIM_MAX=10
	INTEGER LIMITS(NLIM_MAX)
	INTEGER NPLTS
	INTEGER IREC
	INTEGER IVAR
	INTEGER I
	INTEGER J
	INTEGER ID
	INTEGER IT,IT2
	INTEGER NY
	INTEGER NITSF
	INTEGER LST_NG
	INTEGER RITE_N_TIMES
	INTEGER LUSCR
	INTEGER K
	INTEGER IOS
C
	LOGICAL LOG_Y_AXIS
	LOGICAL NEWMOD
	LOGICAL WRITE_RVSIG
	LOGICAL DO_ABS
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	CHARACTER*10 TMP_STR
	CHARACTER*10 PLT_OPT
	CHARACTER*80 YLABEL
	CHARACTER*132 STRING
C
	REAL*8 T1,T2,T3
C
	LUSCR=26
	RITE_N_TIMES=1
	NEWMOD=.TRUE.
C
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine should be run from the data directory'
	WRITE(T_OUT,*)'It expects to find the following files:'
	WRITE(T_OUT,*)'                                       POINT1.DAT'
	WRITE(T_OUT,*)'                                       SCRTEMP.DAT'
	WRITE(T_OUT,*)'                                       MODEL.DAT'
	WRITE(T_OUT,*)' '
C
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
C
	OPEN(UNIT=12,FILE='POINT1',STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .EQ. 0)READ(12,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .EQ. 0)THEN
            IF(INDEX(STRING,'!Format date') .EQ. 0)THEN
              READ(STRING,*,IOSTAT=IOS)K,NIT
	    ELSE
              READ(12,*,IOSTAT=IOS)K,NIT
	    END IF
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Possible error reading POINT1'
	    CALL GEN_IN(NIT,'Number of iterations')
	  END IF
	CLOSE(UNIT=12)
C
	ALLOCATE (POPS(NT,ND,NIT))
	ALLOCATE (R_MAT(ND,NIT))
	ALLOCATE (V_MAT(ND,NIT))
	ALLOCATE (RAT(NT,ND))
!
	ALLOCATE (R(ND))
	ALLOCATE (V(ND))
	ALLOCATE (SIGMA(ND))
C
	J=MAX(NIT,ND); J=MAX(NT,J)
	ALLOCATE (X(J))
	ALLOCATE (Y(J))
	ALLOCATE (Z(J))
!
	WRITE(6,*)'   Number of depth points is:',ND
	WRITE(6,*)'Number of variables/depth is:',NT
	WRITE(6,*)'     Number of iterations is:',NIT
!
	ALLOCATE (I_BIG(NT))
	ALLOCATE (Z_BIG(NT))
C
	DO IREC=1,NIT
	  CALL SCR_READ_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  R_MAT(:,IREC)=R(:) 
	  V_MAT(:,IREC)=V(:) 
	END DO
C
200	CONTINUE
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'ND  :: ',ND
	WRITE(T_OUT,*)'NT  :: ',NT
	WRITE(T_OUT,*)'NIT :: ',NIT
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'F   :: Z(K)=100.0D0*(Y(K+1)-Y(K))/Y(K+1)'
	WRITE(T_OUT,*)'R   :: [Y(K+2)-Y(K+1)]/[Y(K+1)-Y(K)]'
	WRITE(T_OUT,*)'D   :: Z(K)=100.0D0*(Y(K)-Y(NIT))/Y(NIT)'
	WRITE(T_OUT,*)'Y   :: Z(K)=Y(K)'
	WRITE(T_OUT,*)' '
        WRITE(T_OUT,*)'PD  :: Plot 100.0D0*(Y(K)-Y(K-1))/Y(K) as a function of depth index.'
        WRITE(T_OUT,*)'PN  :: Plot a variable as a function of depth index.'
        WRITE(T_OUT,*)'PV  :: Plot a variable as a function of velocity.'
	WRITE(T_OUT,*)'PF  :: Plot 100.0D0*(Y(K+1)-Y(K))/Y(K+1) for all variables at a given depth.'
        WRITE(T_OUT,*)'VR  :: Plot velocity as a function of radius.'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'MED_R  :: Median corection as a function of depth'
	WRITE(T_OUT,*)'MR     :: Z(K)=100.0D0*(MEAN[Y(K-1)-Y(K-2)]/[Y(K)-Y(K-1)] - 1.0)'
	WRITE(T_OUT,*)'IR     :: Z(ID)=100.0D0*(MEAN[Y(K-1)-Y(K-2)]/[Y(K)-Y(K-1)] - 1.0)'
	WRITE(T_OUT,*)'WRST   :: Writes fractional corections to file (FRAC_COR -- same format as STEQ_VALS'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'FDG    :: Fudge individual values at a single depth and output to SCRTEMP'
	WRITE(T_OUT,*)'FDGV   :: Fudge values over a ranges of depths (% change) and output to SCRTEMP'
	WRITE(T_OUT,*)'INT    :: Interpolate values whose corrections are above a certain % limit'
	WRITE(T_OUT,*)'UNDO   :: Undo corrections over a range of depths'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'LY  :: Switch to/from Log(Y) for options where appropriate (not full implemented)'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'E   :: EXIT'
	WRITE(T_OUT,*)' '
	PLT_OPT='R'
	CALL GEN_IN(PLT_OPT,'Plot option: R(atio), F(rac) or D(elta), Y, IR, MR or E(xit)')
	CALL SET_CASE_UP(PLT_OPT,0,0)
!
	IF(PLT_OPT(1:2) .EQ. 'E ' .OR. PLT_OPT(1:2) .EQ. 'EX')STOP
	IF(PLT_OPT(1:2) .EQ. 'LY')THEN
	   WRITE(6,*)RED_PEN
	  IF(LOG_Y_AXIS)WRITE(6,*)'Switching to linear Y axis'
	  IF(.NOT. LOG_Y_AXIS)WRITE(6,*)'Switching to logarithmic Y axis'
	  WRITE(6,'(A)')DEF_PEN
	  TMP_STR=' '; CALL GEN_IN(TMP_STR,'Hit any character to continue')
	  LOG_Y_AXIS=.NOT. LOG_Y_AXIS
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:5) .EQ. 'CHK_R')THEN
	  DO IT=1,NIT
	    WRITE(6,'(2ES16.6)')R_MAT(1,IT),R_MAT(ND,IT)
	  END DO
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'WRST')THEN
	  I=1; K=NIT
	  CALL GEN_IN(K,'Iteration to check')
	  CALL GEN_IN(I,'Comparison iteration')
	  RAT(:,:)=(POPS(:,:,I)-POPS(:,:,K))/POPS(:,:,I)
	  OPEN(UNIT=11,FILE='FRAC_COR',STATUS='UNKNOWN')
	    WRITE(TMP_STR,'(I4)')K; TMP_STR=ADJUSTL(TMP_STR)
	    STRING='Fractional changes: 1- POP(IT='//TRIM(TMP_STR)
	    WRITE(TMP_STR,'(I4)')I; TMP_STR=ADJUSTL(TMP_STR)
	    STRING=TRIM(STRING)//')/POP(IT='//TRIM(TMP_STR)//')'
	    CALL WR2D_V2(RAT,NT,ND,STRING,'#',L_TRUE,11)
	  CLOSE(UNIT=11)
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'CHNG')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT
	    CALL GEN_IN(IT,'Iteraton # of lower iteration')
	  END DO
	  ID=ND
	  CALL GEN_IN(ID,'Depth index')
	  DO IVAR=1,NT
	    X(IVAR)=IVAR
	    T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	    T1=POPS(IVAR,ID,NIT)-POPS(IVAR,ID,IT-1)
	    IF(ABS(T2). GT. 1.0D-50*ABS(T1))THEN
	     Y(IVAR)=T1/T2
	    ELSE
	     Y(IVAR)=20000
	    END IF
	  END DO
	  CALL DP_CURVE(NT,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID','Ratio',' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:5) .EQ. 'MED_R')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT; CALL GEN_IN(IT,'Iteraton #')
	  END DO
	  DO ID=1,ND
	    Y(ID)=0.0D0
	    K=0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(ABS(T2). GT. 1.0D-50*ABS(T1))THEN
	        Z_BIG(IVAR)=T1/T2
	      ELSE
	        Z_BIG(IVAR)=20000
	      END IF
	    END DO
!
! Now find median value. First order.
!
	    CALL INDEXX(NT,Z_BIG,I_BIG,.TRUE.)
	    J=(NT+1)/2
	    Y(ID)=100.0D0*(Z_BIG(I_BIG(J))-1.0D0)
	    X(ID)=FLOAT(ID)
	    Ylabel='100(r\dmed\u-1)'
!	    IF(ID .EQ. 30)THEN
	      WRITE(74,*)'ID=',ID
	      DO J=1,NT
	        K=I_BIG(J)
!	        WRITE(74,*)J,I_BIG(J),Z_BIG(I_BIG(J))
	        WRITE(74,'(2I4,4ES16.8)')J,K,Z_BIG(K),POPS(K,ID,IT),POPS(K,ID,IT-1),POPS(K,ID,IT-2)
	      END DO
!	    END IF
	  END DO
	  CALL DP_CURVE(ND,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'MR')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT; CALL GEN_IN(IT,'Iteraton #')
	  END DO
	  DO_ABS=.TRUE.; CALL GEN_IN(DO_ABS,'Use absolute value')
	  DO ID=1,ND
	    Y(ID)=0.0D0
	    K=0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(T2 .NE. 0)THEN
	        T1=T1/T2
	        IF(DO_ABS)T1=ABS(T1)
	        Y(ID)=Y(ID)+ T1
	        K=K+1
	      ELSE
	        WRITE(6,*)'Problem with variable:',IVAR
	      END IF
	    END DO
	    Y(ID)=100.0D0*(Y(ID)/K-1.0D0)
	    X(ID)=FLOAT(ID)
	    Ylabel='100(AVE[r]-1)'
	  END DO
	  CALL DP_CURVE(ND,X,Y)
	  CALL GRAMON_PGPLOT('Depth ID',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'HR')THEN
	  IT=0
	  DO WHILE(IT .LT. 1 .OR. IT .GT. NIT)
	    IT=NIT
	    CALL GEN_IN(IT,'Iteration #')
	  END DO
	  ID=0
	  DO WHILE(ID .LT. 1 .OR. ID .GT. ND)
	    ID=ND/2; CALL GEN_IN(ID,'Depth index')
	  END DO
	  Y(ID)=0.0D0
!
	  K=MIN(201,NT)
	  T3=200.0D0/(K-1)
	  DO J=1,K
	    X(J)=-50.0D0+(J-1)*T3
	  END DO
	  Y(1:K)=0.0D0
!
	  DO IVAR=1,NT
	    T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	    T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	    IF(T2 .NE. 0)THEN
	      T1=(100*(T1/T2-1.0D0)+51)/T3
	      IF(T1 .LT. 1)T1=1
	      IF(T1 .GT. K)T1=K
	      Y(NINT(T1))=Y(NINT(T1))+1
	    ELSE
	      Y(1)=Y(1)+1
	    END IF
	    Ylabel='N(r)'
	  END DO
	  CALL DP_CURVE(K,X,Y)
	  CALL GRAMON_PGPLOT('100(r-1)',Ylabel,' ',' ')
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'IR')THEN
	  ID=0
	  DO WHILE(ID .LT. 1 .OR. ID .GT. ND)
	    ID=ND; CALL GEN_IN(ID,'Depth index')
	  END DO
	  DO IT=3,NIT
	    Y(IT)=0.0D0
	    DO IVAR=1,NT
	      T2=POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1)
	      T1=POPS(IVAR,ID,IT-1)-POPS(IVAR,ID,IT-2)
	      IF(T2 .NE. 0)THEN
	        Y(IT)=Y(IT)+ T1/T2
	        K=K+1
	      ELSE
	        WRITE(6,*)'Problem with variable:',IVAR,'for iteration',IT
	      END IF
	    END DO
	    Y(IT)=100.0D0*(Y(IT)/NT-1.0D0)
	    X(IT)=FLOAT(IT)
	    Ylabel='100(AVE[r]-1)'
	  END DO
	  IT=NIT-2
	  CALL DP_CURVE(IT,X(3),Y(3))
	  CALL GRAMON_PGPLOT('Iteration',Ylabel,' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PF')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    CALL GEN_IN(ID,'Depth (zero to exit)')
	    IF(ID .EQ. 0)EXIT
	    DO IVAR=1,NT
	      Y(IVAR)=100.0D0*(POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1))/POPS(IVAR,ID,IT)
	      X(IVAR)=IVAR
	    END DO
	    CALL DP_CURVE(NT,X,Y)
	  END DO
	  Ylabel=''
	  CALL GRAMON_PGPLOT('Depth',Ylabel,' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PD')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
	    IF(IVAR .EQ. 0)EXIT
	    DO ID=1,ND
	      Y(ID)=100.0D0*(POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1))/POPS(IVAR,ID,IT)
	      X(ID)=ID
	    END DO
	    CALL DP_CURVE(ND,X,Y)
	  END DO
	  YLABEL='[Y(I-1)-Y(I)]/Y(I) [%]'
	  CALL GRAMON_PGPLOT('Depth',Ylabel,' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PN')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      Z(ID)=LOG10(POPS(IVAR,ID,IT))
	      X(ID)=ID
	    END DO
	    IF(LOG_Y_AXIS)THEN
	      CALL DP_CURVE(ND,X,Z)
	      Ylabel='Log'
	    ELSE
	      CALL DP_CURVE(ND,X,Y)
	      Ylabel=''
	    END IF
	  END DO
	  CALL GRAMON_PGPLOT('Depth',Ylabel,' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PV')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    T1=1.0D0
	    IF(V(1) .GT. 10000.0D0)T1=1.0D-03
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      Z(ID)=LOG10(POPS(IVAR,ID,IT))
	      X(ID)=T1*V_MAT(ID,IT)
	    END DO
	    IF(LOG_Y_AXIS)THEN
	      CALL DP_CURVE(ND,X,Z)
	      Ylabel='Log'
	    ELSE
	      CALL DP_CURVE(ND,X,Y)
	      Ylabel=''
	    END IF
	  END DO
	  IF(V(1) .GT. 10000.0D0)THEN
	    CALL GRAMON_PGPLOT('V(Mm/s)',Ylabel,' ',' ')
	  ELSE
	    CALL GRAMON_PGPLOT('V(km/s)',Ylabel,' ',' ')
	  END IF
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PR')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    T1=1.0D0
	    IF(R(1) .GT. 1.0D+04)T1=1.0D-04
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      Z(ID)=LOG10(POPS(IVAR,ID,IT))
	      X(ID)=1.0D-04*R_MAT(ID,IT)
	    END DO
	    IF(LOG_Y_AXIS)THEN
	      CALL DP_CURVE(ND,X,Z)
	      Ylabel='Log'
	    ELSE
	      CALL DP_CURVE(ND,X,Y)
	      Ylabel=''
	    END IF
	  END DO
	  IF(R(1) .GT. 1.0D+04)THEN
	    CALL GRAMON_PGPLOT('R(10\u14 \dcm)',Ylabel,' ',' ')
	  ELSE
	    CALL GRAMON_PGPLOT('R(10\u10 \dcm)',Ylabel,' ',' ')
	  END IF
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'VR')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NIT)
	    IF(IT .EQ. 0)EXIT
	    T1=1.0D0; T2=1.0D0
	    IF(R(1) .GT. 1.0D+04)T1=1.0D-04
	    IF(V(1) .GT. 1.0D+04)T2=1.0D-04
	    DO ID=1,ND
	      X(ID)=T1*R_MAT(ID,IT)
	      Y(ID)=T2*V_MAT(ID,IT)
	    END DO
	    CALL DP_CURVE(ND,X,Y)
	  END DO
	  Ylabel='V(km/s)'
	  IF(V(1) .GT. 1.0D+04)Ylabel='V(Mm/s)'
	  IF(R(1) .GT. 1.0D+04)THEN
	    CALL GRAMON_PGPLOT('R(10\u14 \dcm)',Ylabel,' ',' ')
	  ELSE
	    CALL GRAMON_PGPLOT('R(10\u10 \dcm)',Ylabel,' ',' ')
	  END IF
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'UNDO')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT; LIMITS(:)=0; LIMITS(1)=1; LIMITS(2)=ND
	  CALL GEN_IN(IT,'Primary Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  IT2=IT-1
	  CALL GEN_IN(IT2,'Iteration to be merged (replaces values)',LOW_LIM=IZERO,UP_LIM=NIT)
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(LIMITS,J,NLIM_MAX,'Depths (L1 to L2, L3 to L4 etc)')
	    DO I=1,J,2
	      IF(LIMITS(I) .EQ. 0)EXIT
	      DO ID=LIMITS(I),LIMITS(I+1)
	        POPS(:,ID,IT)=POPS(:,ID,IT2)
	      END DO
	    END DO
	    IVAR=0
	  END DO
          NITSF=NITSF+1; IREC=NIT+FDG_COUNTER
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'A new record is writted every time FDG or FDGV is called'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:4) .EQ. 'FDGV')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT; LIMITS(:)=0; LIMITS(1)=1; LIMITS(2)=ND
	  CALL GEN_IN(IT,'Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  IF(IT .EQ. 0)GOTO 200
	  T1=0.0D0; T2=0.0D0
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    CALL GEN_IN(T1,'% change in variable')
	    IF(T1 .EQ. 0.0D0)CALL GEN_IN(T2,'Change in variable')
	    CALL GEN_IN(LIMITS,J,NLIM_MAX,'Depths (L1 to L2, L3 to L4 etc)')
	    DO I=1,J,2
	      IF(LIMITS(I) .EQ. 0)EXIT
	      DO ID=LIMITS(I),LIMITS(I+1)
	        T3=POPS(IVAR,ID,IT)
	        POPS(IVAR,ID,IT)=T3*(1.0D0+T1/100.0D0)+T2
	      END DO
	    END DO
	    IVAR=0
	  END DO
          NITSF=NITSF+1; IREC=NIT+FDG_COUNTER
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'A new record is writted every time FDG or FDGV is called'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'FDG')THEN
	  FDG_COUNTER=FDG_COUNTER+1
	  IT=NIT; ID=ND; IVAR=NT
	  CALL GEN_IN(IT,'Iteration # (zero to exit) - default is last iteration',
	1                   LOW_LIM=IZERO,UP_LIM=NIT)
	  IF(IT .EQ. 0)GOTO 200
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)',LOW_LIM=IZERO,UP_LIM=NT)
	    IF(IVAR .EQ. 0)EXIT
	    CALL GEN_IN(ID,'Depth of variable')
	    WRITE(6,'(7(9X,I5))')(I,I=MAX(ID-3,1),MIN(ID+3,ND))
	    WRITE(6,'(7ES14.4)')(POPS(IVAR,I,IT),I=MAX(ID-3,1),MIN(ID+3,ND))
	    T1=POPS(IVAR,ID,IT)
	    CALL GEN_IN(T1,'New value of variable')
	    POPS(IVAR,ID,IT)=T1
	  END DO
          NITSF=NITSF+1; IREC=NIT+FDG_COUNTER
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IT),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new (and last) iteration.'
	  WRITE(6,*)'A new record is writted every time FDG or FDGV is called'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'INT')THEN
	  IT=NIT; ID=ND; T2=100.0D0
	  CALL GEN_IN(IT,'Iteration # (zero to exit)')
	  CALL GEN_IN(ID,'Depth of variable')
	  CALL GEN_IN(T2,'Interpolate values with correction > >%')
	  DO IVAR=1,NT-1
	    T1=100.0D0*ABS(POPS(IVAR,ID,IT)-POPS(IVAR,ID,IT-1))/POPS(IVAR,ID,IT)
	    IF(T1 .GT. T2)THEN
	      WRITE(6,*)'Replacing population for variable',IVAR
	      IF(ID .EQ. 2 .OR. ID .EQ. ND)THEN
	        POPS(IVAR,ID,IT)=POPS(IVAR,ID-1,IT)
	      ELSE IF(ID .EQ. 1)THEN
	        POPS(IVAR,ID,IT)=POPS(IVAR,ID+1,IT)
	      ELSE
	        T1=LOG(R(ID)/R(ID-1))/LOG(R(ID+1)/R(ID-1))
	        POPS(IVAR,ID,IT)=EXP( T1*LOG(POPS(IVAR,ID+1,IT)) +
	1                     (1.0D0-T1)*LOG(POPS(IVAR,ID-1,1)) )
	      END IF
	    END IF
	  END DO
          NITSF=NITSF+1; IREC=IT
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new iteration.'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:3) .EQ. 'FIX_OSC')THEN
	  CALL FIX_POP_OSCILLATIONS(POPS(1,1,NIT),R,V,SIGMA,LUSCR,ND,NT)
          NITSF=NITSF+1; IREC=NIT
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	  WRITE(6,*)'Corrections written to SCRTEMP as new iteration.'
	  WRITE(6,*)'Restart program if you wish to compare to with pops from last iteration.'
	  WRITE(6,*)'Populations can be compared with older iterations.'
	  GOTO 200
!
	ELSE IF(PLT_OPT .EQ. 'R' .OR.
	1       PLT_OPT .EQ. 'F' .OR.
	1       PLT_OPT .EQ. 'D' .OR.
	1       PLT_OPT .EQ. 'Y')THEN
	  ID=ND
	  DO WHILE(0 .EQ. 0)	
	    NPLTS=0
	    IVAR=1
	    DO WHILE(IVAR .NE. 0)
500	      IVAR=0
	      WRITE(STRING,'(I5,A)')NT,'](0 to plot)'
	      DO WHILE(STRING(1:1) .EQ. ' ') ; STRING(1:)=STRING(2:) ; END DO
	      STRING='Variable to be plotted ['//STRING
	      CALL GEN_IN(IVAR,STRING)
	      IF(IVAR .LT. 0 .OR. IVAR .GT. NT)GO TO 500
	      IF(IVAR .EQ. 0)GOTO 1000
C
600	      CONTINUE
	      WRITE(STRING,'(I5,A)')ND,'](0 to plot)'
	      DO WHILE(STRING .EQ. ' ') ; STRING(1:)=STRING(2:); END DO
	      STRING='Depth of variable to be plotted ['//STRING
	      CALL GEN_IN(ID,STRING)
	      IF(ID .LE. 0 .OR. ID .GT. ND)GO TO 600
C
	      Y(1:NIT)=POPS(IVAR,ID,1:NIT)
C
  	      IF(PLT_OPT .EQ. 'F')THEN
  	        DO K=1,NIT-1
  	          Z(K)=100.0D0*(Y(K+1)-Y(K))/Y(K+1)
  	          X(K)=FLOAT(K)
  	        END DO
  	        NY=NIT-1
  	        T1=MAXVAL(ABS(Z(1:NY)))
  	        IF(T1 .LT. 1.0D-02)THEN
  	          Z(1:NY)=Z(1:NY)*1.0D+03
  	          YLABEL='\gDY/Y(%)\d \ux10\u3\d'
  	          WRITE(T_OUT,*)'Correction scaled by factor of 10^3'
  	        ELSE
	          YLABEL='\gDY/Y(%)'
	        END IF
!	      
	      ELSE IF(PLT_OPT .EQ. 'R')THEN
	        DO K=1,NIT-2
	          T1=Y(K+2)-Y(K+1)
	          T2=Y(K+1)-Y(K)
	          IF(T2 .NE. 0)THEN
	             Z(K)=T1/T2
	          ELSE
	             Z(K)=10.0
	          END IF
	          X(K)=FLOAT(K)+2
	        END DO
	        NY=NIT-2
	        YLABEL='\gDY(K+1)/\gDY(K)'
	      ELSE IF(PLT_OPT .EQ. 'D')THEN
	        DO K=1,NIT
	          Z(K)=100.0D0*(Y(K)-Y(NIT))/Y(NIT)
	          X(K)=FLOAT(K)
	        END DO
	        NY=NIT-1
	        T1=MAXVAL(ABS(Z(1:NY)))
	        IF(T1 .LT. 1.0D-02)THEN
	          Z(1:NY)=Z(1:NY)*1.0D+03
	          YLABEL='[Y(K)-Y(NIT)]/Y(NIT) [%]\d \ux10\u3\d'
	          WRITE(T_OUT,*)'Correction scaled by factor of 10^3'
	        ELSE
	        YLABEL='[Y(K)-Y(NIT)]/Y(NIT) [%]'
	        END IF
	      ELSE IF(PLT_OPT .EQ. 'Y')THEN
	        DO K=1,NIT
	          Z(K)=Y(K)
	          IF(LOG_Y_AXIS)Z(K)=LOG10(Z(K))
	          X(K)=FLOAT(K)
	        END DO
	        NY=NIT
	        YLABEL='Y(K)'
	        IF(LOG_Y_AXIS)YLABEL='Log Y(K)'
	      END IF
	      CALL DP_CURVE(NY,X,Z)
	      NPLTS=NPLTS+1
	    END DO
1000	    CONTINUE
	    IF(NPLTS .NE. 0)THEN
	      CALL GRAMON_PGPLOT('Iteration number K',Ylabel,' ',' ')
	    END IF
	    GOTO 200
	  END DO
	ELSE
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'Unrecognized command'
	  WRITE(6,'(A)')DEF_PEN
	  TMP_STR=' '; CALL GEN_IN(TMP_STR,'Hit any character to continue')
	  GOTO 200
	END IF
!
	END
