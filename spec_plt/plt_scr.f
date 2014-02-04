!
! Simple program to plot data from the SCRTEMP file. This can be used to check
!   Convergence of a program.
!   Convergence as a function of depth etc.
!
	PROGRAM PLT_SCR
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
	INTEGER ND,NT,NIT
C
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
	INTEGER, PARAMETER :: T_OUT=6
C
	INTEGER, PARAMETER :: NLIM_MAX=10
	INTEGER LIMITS(NLIM_MAX)
	INTEGER NPLTS
	INTEGER IREC
	INTEGER IVAR
	INTEGER I
	INTEGER J
	INTEGER ID
	INTEGER IT
	INTEGER NY
	INTEGER NITSF
	INTEGER LST_NG
	INTEGER RITE_N_TIMES
	INTEGER LUSCR
	INTEGER K
	INTEGER IOS
C
	LOGICAL NEWMOD
	LOGICAL WRITE_RVSIG
	LOGICAL DO_ABS
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
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
	WRITE(T_OUT,*)'WRST   :: Writes fractional corections to file (same format as STEQ_VALS'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'FDG    :: Fudge individual values at a single depth and output to SCRTEMP'
	WRITE(T_OUT,*)'FDGV   :: Fudge values over a ranges of depths (% change) and output to SCRTEMP'
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'E   :: EXIT'
	WRITE(T_OUT,*)' '
	PLT_OPT='R'
	CALL GEN_IN(PLT_OPT,'Plot option: R(atio), F(rac) or D(elta), Y, IR, MR or E(xit)')
	CALL SET_CASE_UP(PLT_OPT,0,0)
!
	IF(PLT_OPT(1:2) .EQ. 'E ' .OR. PLT_OPT(1:2) .EQ. 'EX')STOP
	IF(PLT_OPT(1:5) .EQ. 'CHK_R')THEN
	  DO IT=1,NIT
	    WRITE(6,'(2ES16.6)')R_MAT(1,IT),R_MAT(ND,IT)
	  END DO
	ELSE IF(PLT_OPT(1:4) .EQ. 'WRST')THEN
	  RAT(:,:)=(POPS(:,:,1)-POPS(:,:,NIT))/POPS(:,:,1)
	  OPEN(UNIT=11,FILE='POP_RATIOS',STATUS='UNKNOWN')
          CALL WR2D_V2(RAT,NT,ND,'Fractional changes from starting solution','#',L_TRUE,11)
	  CLOSE(UNIT=11)
	  GOTO 200
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
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
	    IF(IVAR .EQ. 0)EXIT
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      X(ID)=ID
	    END DO
	    CALL DP_CURVE(ND,X,Y)
	  END DO
	  Ylabel=''
	  CALL GRAMON_PGPLOT('Depth',Ylabel,' ',' ')
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PV')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
	    IF(IVAR .EQ. 0)EXIT
	    T1=1.0D0
	    IF(V(1) .GT. 10000.0D0)T1=1.0D-03
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      X(ID)=T1*V_MAT(ID,IT)
	    END DO
	    CALL DP_CURVE(ND,X,Y)
	  END DO
	  Ylabel=''
	  IF(V(1) .GT. 10000.0D0)THEN
	    CALL GRAMON_PGPLOT('V(Mm/s)',Ylabel,' ',' ')
	  ELSE
	    CALL GRAMON_PGPLOT('V(km/s)',Ylabel,' ',' ')
	  END IF
	  GOTO 200
	ELSE IF(PLT_OPT(1:2) .EQ. 'PR')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
	    IF(IVAR .EQ. 0)EXIT
	    T1=1.0D0
	    IF(R(1) .GT. 1.0D+04)T1=1.0D-04
	    DO ID=1,ND
	      Y(ID)=POPS(IVAR,ID,IT)
	      X(ID)=1.0D-04*R_MAT(ID,IT)
	    END DO
	    CALL DP_CURVE(ND,X,Y)
	  END DO
	  IF(R(1) .GT. 1.0D+04)THEN
	    CALL GRAMON_PGPLOT('R(10\u14 \dcm)',Ylabel,' ',' ')
	  ELSE
	    CALL GRAMON_PGPLOT('R(10\u10 \dcm)',Ylabel,' ',' ')
	  END IF
	  Ylabel=''
	  GOTO 200
!
	ELSE IF(PLT_OPT(1:2) .EQ. 'VR')THEN
	  IT=NIT; ID=ND
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IT,'Iteration # (zero to exit)')
	    IF(IT .EQ. 0)EXIT
	    IVAR=NT
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
	    IF(IVAR .EQ. 0)EXIT
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
	ELSE IF(PLT_OPT(1:4) .EQ. 'FDGV')THEN
	  IT=NIT; ID=ND; IVAR=NT; LIMITS(:)=0; LIMITS(1)=1; LIMITS(2)=ND
	  CALL GEN_IN(IT,'Iteration # (zero to exit)')
	  T1=0.0D0; T2=0.0D0
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
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
          NITSF=NITSF+1; IREC=IT
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	ELSE IF(PLT_OPT(1:3) .EQ. 'FDG')THEN
	  IT=NIT; ID=ND; IVAR=NT
	  CALL GEN_IN(IT,'Iteration # (zero to exit)')
	  DO WHILE(1 .EQ. 1)
	    CALL GEN_IN(IVAR,'Variable # (zero to exit)')
	    IF(IVAR .EQ. 0)EXIT
	    CALL GEN_IN(ID,'Depth of variable')
	    T1=POPS(IVAR,ID,IT)
	    CALL GEN_IN(T1,'New value of variable')
	    POPS(IVAR,ID,IT)=T1
	    IVAR=0
	  END DO
          NITSF=NITSF+1; IREC=IT
	  CALL SCR_RITE_V2(R,V,SIGMA,POPS(1,1,IREC),IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              NT,ND,LUSCR,NEWMOD)
	END IF
C
	ID=ND
	DO WHILE(0 .EQ. 0)	
	  NPLTS=0
	  IVAR=1
	  DO WHILE(IVAR .NE. 0)
500	    IVAR=0
	    WRITE(STRING,'(I5,A)')NT,'](0 to plot)'
	    DO WHILE(STRING(1:1) .EQ. ' ') ; STRING(1:)=STRING(2:) ; END DO
	    STRING='Variable to be plotted ['//STRING
	    CALL GEN_IN(IVAR,STRING)
	    IF(IVAR .LT. 0 .OR. IVAR .GT. NT)GO TO 500
	    IF(IVAR .EQ. 0)GOTO 1000
C
600	    CONTINUE
	    WRITE(STRING,'(I5,A)')ND,'](0 to plot)'
	    DO WHILE(STRING(1:1) .EQ. ' ') ; STRING(1:)=STRING(2:); END DO
	    STRING='Depth of variable to be plotted ['//STRING
	    CALL GEN_IN(ID,STRING)
	    IF(ID .LE. 0 .OR. ID .GT. ND)GO TO 600
C
	    Y(1:NIT)=POPS(IVAR,ID,1:NIT)
C
	    IF(PLT_OPT(1:1) .EQ. 'F')THEN
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
	    ELSE IF(PLT_OPT(1:1) .EQ. 'R')THEN
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
	    ELSE IF(PLT_OPT(1:1) .EQ. 'D')THEN
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
	    ELSE IF(PLT_OPT(1:1) .EQ. 'Y')THEN
	      DO K=1,NIT
	        Z(K)=Y(K)
	        X(K)=FLOAT(K)
	      END DO
	      NY=NIT
	      YLABEL='Y(K)'
	    ELSE
	      WRITE(T_OUT,*)' Option not recognized: try again'
	      GOTO 1000
	    END IF
	    CALL DP_CURVE(NY,X,Z)
	    NPLTS=NPLTS+1
	  END DO
1000	  CONTINUE
	  IF(NPLTS .NE. 0)THEN
	    CALL GRAMON_PGPLOT('Iteration number K',Ylabel,' ',' ')
	  ELSE
	    GOTO 200
	  END IF
	END DO
C
	END
