!
! Simple program to extrapolate 2 SN models in time to generate start conditions for
! a new SN model. At present, the program cannot handle models with
! different numbers of SL's.
!
	PROGRAM EXTRAP_SCR
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	USE READ_KEYWORD_INTERFACE
	IMPLICIT NONE
!
! Altered 25-Apr-2016: LST_NG rest to -100o as for new model.
! Altered 10-Feb-2016: Bug fix (for non identucal velocity grids).
! Altered 26-Feb-2016: Bug fix (for case of identical velocity grids).
! Created 19-Feb-2016: Based intially on PLT_SCR.
!
! We deine 3 memebers of the class. Two for the input models, and one
! for the extrapolated model.
!
	TYPE SCRATCH_DATA
	  REAL*8, ALLOCATABLE :: WRK_POPS(:,:)		!NT,ND
	  REAL*8, ALLOCATABLE :: POPS(:,:)		!NT,ND
	  REAL*8, ALLOCATABLE :: R(:)			!ND
	  REAL*8, ALLOCATABLE :: V(:)			!ND
	  REAL*8, ALLOCATABLE :: SIGMA(:)		!ND
	  REAL*8 SN_AGE
	  INTEGER ND
	  INTEGER NT
	END TYPE SCRATCH_DATA
	TYPE (SCRATCH_DATA) SCR(3)
!
	REAL*8, ALLOCATABLE :: X1(:),X2(:)
	REAL*8, ALLOCATABLE :: Y1(:),Y2(:)
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_OUT=6
!
	INTEGER I,J,K
	INTEGER IS
	INTEGER IREC
	INTEGER NIT
	INTEGER NITSF
	INTEGER LST_NG
	INTEGER RITE_N_TIMES
	INTEGER LUSCR
	INTEGER IOS
!
	LOGICAL NEWMOD
	LOGICAL WRITE_RVSIG
	LOGICAL SAME_GRID
	LOGICAL FILE_1_EXISTS
	LOGICAL FILE_2_EXISTS
	LOGICAL FILE_3_EXISTS
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=132) DIRECTORY
	CHARACTER(LEN=132) FILENAME
!
	REAL*8 T1,T2,T3
!
	LUSCR=26
	RITE_N_TIMES=1
	NEWMOD=.TRUE.
!
	WRITE(T_OUT,*)' '
	WRITE(T_OUT,*)'This routine requires the following files'
	WRITE(T_OUT,*)'     POINT1.DAT'
	WRITE(T_OUT,*)'     SCRTEMP.DAT'
	WRITE(T_OUT,*)'     MODEL'
	WRITE(T_OUT,*)'     VADAT'
!
	INQUIRE(FILE='POINT1', EXIST=FILE_1_EXISTS)
	INQUIRE(FILE='POINT2', EXIST=FILE_2_EXISTS)
	INQUIRE(FILE='SCRTEMP',EXIST=FILE_3_EXISTS)
	IF(FILE_1_EXISTS .OR. FILE_2_EXISTS .OR. FILE_3_EXISTS)THEN
	  WRITE(6,*)RED_PEN
	  WRITE(6,*)'Error: POINT1, POINT2 and SCRTEMP cannot exist in this woking directory'
	  WRITE(6,*)DEF_PEN
	  STOP
	END IF
!
	DO IS=1,2
!
! We add a slash to the directory name if it was not included.
!
	  WRITE(T_OUT,*)' '
	  DIRECTORY=' '
	  IF(IS .EQ. 1)THEN
	    CALL GEN_IN(DIRECTORY,'Directory path to SCRTEMP file at earliest time')
	  ELSE
	    CALL GEN_IN(DIRECTORY,'Directory path to SCRTEMP file at latest time')
	  END IF
	  J=LEN_TRIM(DIRECTORY)
	  IF(DIRECTORY(J:J) .NE. '/')DIRECTORY(J+1:J+1)='/'
!
	  FILENAME=TRIM(DIRECTORY)//'MODEL'
	  OPEN(UNIT=12,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)'Error unable to open '//TRIM(FILENAME)
	    WRITE(6,*)'Stopping execution'
	    WRITE(6,*)DEF_PEN
	    WRITE(6,*)' '
	    STOP
	  END IF
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of depth') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)SCR(IS)%ND
	  DO WHILE(INDEX(STRING,'!Total number of variables') .EQ. 0)
	    READ(12,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)GOTO 100
	  END DO
	  READ(STRING,*)SCR(IS)%NT
!
100	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to read MODEL file'
	    CALL GEN_IN(SCR(IS)%NT,'Total number of levels')
	    CALL GEN_IN(SCR(IS)%ND,'Number of depth points')
	  END IF 
	  CLOSE(UNIT=12)
!
	  FILENAME=TRIM(DIRECTORY)//'VADAT'
	  CALL READ_KEYWORD(SCR(IS)%SN_AGE,'SN_AGE',L_TRUE,FILENAME,L_TRUE,L_TRUE,LUSCR)
!
	  FILENAME=TRIM(DIRECTORY)//'POINT1'
	  OPEN(UNIT=12,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
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
	      STOP
	    END IF
	  CLOSE(UNIT=12)
!
	  ALLOCATE (SCR(IS)%POPS(SCR(IS)%NT,SCR(IS)%ND))
	  ALLOCATE (SCR(IS)%R(SCR(IS)%ND))
	  ALLOCATE (SCR(IS)%V(SCR(IS)%ND))
	  ALLOCATE (SCR(IS)%SIGMA(SCR(IS)%ND))
!
	  WRITE(6,*)'   Number of depth points is:',SCR(IS)%ND
	  WRITE(6,*)'Number of variables/depth is:',SCR(IS)%NT
	  WRITE(6,*)'     Number of iterations is:',NIT
!
! As we don't pass the filename to SCRTEMP, we need to executa system
! command to link POINT1 etc to thre write directory.
!
	  STRING='ln -sf '//TRIM(DIRECTORY)//'POINT1     POINT1'
	  CALL SYSTEM(STRING)
	  STRING='ln -sf '//TRIM(DIRECTORY)//'POINT2     POINT2'
	  CALL SYSTEM(STRING)
	  STRING='ln -sf '//TRIM(DIRECTORY)//'SCRTEMP    SCRTEMP'
	  CALL SYSTEM(STRING)
!
	  CALL SCR_READ_V2(SCR(IS)%R,SCR(IS)%V,SCR(IS)%SIGMA,SCR(IS)%POPS,
	1              NIT,NITSF,RITE_N_TIMES,LST_NG,WRITE_RVSIG,
	1              SCR(IS)%NT,SCR(IS)%ND,LUSCR,NEWMOD)
	  CALL SYSTEM('unlink POINT1')
	  CALL SYSTEM('unlink POINT2')
	  CALL SYSTEM('unlink SCRTEMP')
!
	END DO
!
	IF(SCR(1)%NT .NE. SCR(2)%NT)THEN
	  WRITE(6,*)' '
	  WRITE(6,*)'Error - at present the extrapolation cannot handle'//
	1                 ' a different number of super levels.'
	  WRITE(6,*)' '
	  STOP
	END IF
!
	IS=3
	SCR(IS)%ND=SCR(2)%ND
	SCR(IS)%NT=SCR(2)%NT
	ALLOCATE (SCR(IS)%POPS(SCR(IS)%NT,SCR(IS)%ND))
	ALLOCATE (SCR(IS)%R(SCR(IS)%ND))
	ALLOCATE (SCR(IS)%V(SCR(IS)%ND))
	ALLOCATE (SCR(IS)%SIGMA(SCR(IS)%ND))
!
	WRITE(T_OUT,*)' '
	SCR(3)%SN_AGE=SCR(2)%SN_AGE*(SCR(2)%SN_AGE/SCR(1)%SN_AGE)
	STRING='Age of new SN model in days '//RED_PEN//'[Must be exact]'//DEF_PEN
	CALL GEN_IN(SCR(3)%SN_AGE,STRING)
	WRITE(T_OUT,*)' '
!
! 1.0D-05 since V is in km/s, and R in units of 10^10 cm.
!
	T1=1.0D-05*24.0D0*3600.0D0*(SCR(3)%SN_AGE-SCR(2)%SN_AGE)
	DO I=1,SCR(3)%ND
	  SCR(3)%V(I)=SCR(2)%V(I)
	  SCR(3)%SIGMA(I)=SCR(2)%SIGMA(I)
	  SCR(3)%R(I)=SCR(2)%R(I)+T1*SCR(2)%V(I)
	END DO
	WRITE(6,*)'Set new R, V and SIGMA values'
!
! Since the density scales as 1/t^3, we scale all populations (but not
! the temperature) by the cube of the SN age.
!
	T1=SCR(1)%SN_AGE**3
	T2=SCR(2)%SN_AGE**3
	DO I=1,SCR(3)%NT-1
	  SCR(1)%POPS(I,:)=SCR(1)%POPS(I,:)*T1
	END DO
	DO I=1,SCR(3)%NT-1
	  SCR(2)%POPS(I,:)=SCR(2)%POPS(I,:)*T2
	END DO
!
! Check to see if both models have the same velocity grid.
!
	IF(SCR(1)%ND .EQ. SCR(2)%ND)THEN
	  SAME_GRID=.TRUE.
	  T2=0.0D0
	  DO I=1,SCR(2)%ND
	    T1=ABS(1.0D0-SCR(2)%V(I)/SCR(1)%V(I))
	    T2=MAX(T2,T1)
	  END DO
	  IF(T2 .GT. 1.0D-08)SAME_GRID=.FALSE.
	ELSE
	 SAME_GRID=.FALSE.
	END IF
	IF(SAME_GRID)WRITE(6,*)'Both models have an identical velocity grid'
!
	IF(SAME_GRID)THEN
	  T1=LOG(SCR(3)%SN_AGE/SCR(1)%SN_AGE)/LOG(SCR(2)%SN_AGE/SCR(1)%SN_AGE)
	  SCR(3)%POPS=T1*LOG(SCR(2)%POPS)+(1.0D0-T1)*LOG(SCR(1)%POPS)
	  SCR(3)%POPS=EXP(SCR(3)%POPS)
	  T3=SCR(3)%SN_AGE**3
	  DO I=1,SCR(3)%NT-1
	    SCR(3)%POPS(I,:)=SCR(3)%POPS(I,:)/T3
	  END DO
	ELSE
	  ALLOCATE (X1(SCR(1)%ND))
	  ALLOCATE (Y1(SCR(1)%ND))
	  ALLOCATE (X2(SCR(2)%ND))
	  ALLOCATE (Y2(SCR(2)%ND))
	  ALLOCATE (SCR(1)%WRK_POPS(SCR(3)%NT,SCR(3)%ND))
	  X1=DLOG(SCR(1)%V)
	  X2=DLOG(SCR(2)%V); X1(1)=X2(1); X1(SCR(1)%ND)=X2(SCR(2)%ND)
	  WRITE(6,*)'Starting to interpolate Model 1 onto Model 2 grid'
	  DO I=1,SCR(3)%NT
	    Y1=DLOG(SCR(1)%POPS(I,:))
	    CALL MON_INTERP(Y2,SCR(2)%ND,IONE,X2,SCR(2)%ND,Y1,SCR(1)%ND,X1,SCR(1)%ND)
	    SCR(1)%WRK_POPS(I,:)=EXP(Y2)
	  END DO 
	  WRITE(6,*)'Interpolated Model 2 onto Model 1 grid'
!
! We extrapolate in log space as this advoids -ve entries.
!
	  T1=LOG(SCR(3)%SN_AGE/SCR(1)%SN_AGE)/LOG(SCR(2)%SN_AGE/SCR(1)%SN_AGE)
	  SCR(3)%POPS=T1*LOG(SCR(2)%POPS)+(1.0D0-T1)*LOG(SCR(1)%WRK_POPS)
	  SCR(3)%POPS=EXP(SCR(3)%POPS)
	  T3=SCR(3)%SN_AGE**3
	  DO I=1,SCR(3)%NT-1
	    SCR(3)%POPS(I,:)=SCR(3)%POPS(I,:)/T3
	  END DO
	END IF
!
	STRING='ln -sf NEW_POINT1     POINT1'
	CALL SYSTEM(STRING)
	STRING='ln -sf NEW_POINT2     POINT2'
	CALL SYSTEM(STRING)
	STRING='ln -sf NEW_SCRTEMP    SCRTEMP'
	CALL SYSTEM(STRING)
!
	NITSF=1; IREC=0; IS=3; LST_NG=-1000
	CALL SCR_RITE_V2(SCR(IS)%R,SCR(IS)%V,SCR(IS)%SIGMA,SCR(IS)%POPS,IREC,NITSF,
	1              RITE_N_TIMES,LST_NG,L_TRUE,
	1              SCR(IS)%NT,SCR(IS)%ND,LUSCR,NEWMOD)
!
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')'You should now do: '//RED_PEN
	WRITE(6,'(A)')'        $cmfdist/com/mvscr.sh'//BLUE_PEN
	WRITE(6,'(A)')'to rename NEW_SCRTEMP to SCRTEMP etc'//DEF_PEN
	WRITE(6,'(A)')
!
	STOP
	END
