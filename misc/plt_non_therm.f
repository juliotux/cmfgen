!
! Program requires the following CMFGEN files:
!                                              NON_THERM_DEGRADATION_SPEC 
!                                              RVTJ
!                                              MODEL
!
	PROGRAM PLT_NON_THERM
	USE MOD_COLOR_PEN_DEF
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Cleaned: 06-Nov-2011
!
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: RAD_ENERGY(:)
!
	REAL*8, ALLOCATABLE :: XAXIS(:)
!
	REAL*8, ALLOCATABLE :: XKT(:)
	REAL*8, ALLOCATABLE :: SOURCE(:)
	REAL*8, ALLOCATABLE :: LELEC(:)
	REAL*8, ALLOCATABLE :: YE(:,:)
!
	INTEGER, PARAMETER :: NX=100000
	REAL*8 E_EV(NX)
	REAL*8 ED_EV(NX)
	REAL*8 dX
!
	REAL*8 SPEED_OF_LIGHT
	REAL*8 T1
	REAL*8 PI
	REAL*8 THERMAL_ED
	REAL*8 NON_THERMAL_ED
!
	EXTERNAL SPEED_OF_LIGHT
	CHARACTER(LEN=132) STRING
	CHARACTER(LEN=20)  XLABEL
!
	INTEGER ND
	INTEGER NKT
	INTEGER IOS
	INTEGER I,J
	INTEGER DPTH_INDX
	LOGICAL FILE_OPEN
!
	PI=4.0D0*ATAN(1.0D0)
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
          CALL GEN_IN(ND,'Number of depth points')
        END IF
!
	IOS=0
	ALLOCATE (R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ED(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (T(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RAD_ENERGY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (XAXIS(ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error -- unable to allocate vectors in PLT_NON_THERM'
	  WRITE(6,*)'Error is ',IOS
	  STOP
	END IF
!
	OPEN(UNIT=20,FILE='RVTJ',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Radius') .EQ. 0)
	    READ(20,'(A)')STRING
	  END DO
	  READ(20,*)(R(I),I=1,ND)
	  READ(20,'(A)')STRING
	  READ(20,*)(V(I),I=1,ND)
	  DO WHILE(INDEX(STRING,'Electron density') .EQ. 0)
	    READ(20,'(A)')STRING
	  END DO
	  READ(20,*)(ED(I),I=1,ND)
	  READ(20,'(A)')STRING
	  READ(20,*)(T(I),I=1,ND)
	CLOSE(UNIT=20)

	DO I=1,ND
	  R(I)=R(I)/R(ND)
	END DO
	XAXIS=V(1:ND)
	XLABEL='V(km/s)'
!
	OPEN(UNIT=20,FILE='NON_THERM_DEGRADATION_SPEC',STATUS='OLD',ACTION='READ')
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Number of energy grid:') .EQ. 0)
	    READ(20,'(A)')STRING
	  END DO
	  I=INDEX(STRING,':')
	  READ(STRING(I+1:),*)NKT
          WRITE(6,'(A,I4)')' Number of energy points in the model is: ',NKT
!
	  IOS=0
	  IF(IOS .EQ. 0)ALLOCATE (XKT(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (SOURCE(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LELEC(NKT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (YE(NKT,ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error -- unable to allocate XKT etc in PLT_NON_THERM'
	    WRITE(6,*)'Error is ',IOS
	    STOP
	  END IF
	  WRITE(6,*)'Allocated memory'
!
	  DPTH_INDX=0
	  DO WHILE(1 .EQ. 1)
	    DO WHILE(STRING .EQ. ' ')
	      READ(20,'(A)')STRING
	    END DO
	    IF(INDEX(STRING,'Energy grids') .NE. 0)THEN
	      READ(20,*)(XKT(I),I=1,NKT)
	      WRITE(6,*)'Read XKT'
	    ELSE IF(INDEX(STRING,'Source') .NE. 0)THEN
	      READ(20,*)(SOURCE(I),I=1,NKT)
	      WRITE(6,*)'SOURCE'
	    ELSE IF(INDEX(STRING,'Lelec') .NE. 0)THEN
	      READ(20,*)(LELEC(I),I=1,NKT)
	      WRITE(6,*)'Lelec'
	    ELSE IF(INDEX(STRING,'Depth:') .NE. 0)THEN
	      DPTH_INDX=DPTH_INDX+1
	      READ(20,*)(YE(I,DPTH_INDX),I=1,NKT)
	      WRITE(6,*)'Read depth:',DPTH_INDX
	      IF(DPTH_INDX .EQ. ND)EXIT
	    END IF
	    READ(20,'(A)')STRING
	  END DO
	CLOSE(UNIT=20)
!
	OPEN(UNIT=20,FILE='GENCOOL',STATUS='OLD',ACTION='READ')
	DO I=1,ND,10
	  STRING=' '
	  DO WHILE(INDEX(STRING,'Radiative decay heating') .EQ. 0)
	    READ(20,'(A)')STRING
	  END DO
	  READ(20,*)(RAD_ENERGY(J),J=I,MIN(I+9,ND))
	END DO
	CLOSE(UNIT=20)
!
	WRITE(6,*)' '
	WRITE(6,'(3A)')PG_PEN(2),' Plotting the non-thermal degradation spectrum in red',DEF_PEN
!
	DPTH_INDX=ND/2
10	CONTINUE
	CALL GEN_IN(DPTH_INDX,'Depth index for plotting electron degradation spectrum')
	IF(DPTH_INDX .LE. 0)GOTO 20
	CALL DP_CURVE(NKT,XKT,YE(1,DPTH_INDX))
	DO I=1,NKT
	  E_EV(I)=YE(I,DPTH_INDX)*LOG(XKT(I))/XKT(I)
	END DO
	CALL DP_CURVE(NKT,XKT,E_EV)
!	CALL DP_CURVE(NKT,XKT,LELEC)
	CALL GRAMON_PGPLOT('kT(ev)','YE',' ',' ')
	GOTO 10
20	CONTINUE
!
	DO I=1,ND
	  T1=SQRT(9.109389D-28/2.0D0/1.602177D-12)
	  DO J=1,NKT
	    YE(J,I)=YE(J,I)*T1/SQRT(XKT(J))
	  END DO
	END DO
!
	dX=100.0D0/NX
	DO I=1,NX
	  E_EV(I)=I*dX
	END DO
!
	DPTH_INDX=ND/2
30	CONTINUE
	CALL GEN_IN(DPTH_INDX,'Depth index for plotting electron and non-thermal electron distributions')
	IF(DPTH_INDX .LE. 0)STOP
!
	T1=0.86173324D0*T(DPTH_INDX)
	DO I=1,NX
	  ED_EV(I)=2.0D0*ED(DPTH_INDX)*SQRT(E_EV(I)/PI/T1)*EXP(-E_EV(I)/T1)
	END DO
	THERMAL_ED=0.0D0
	DO I=1,NX-1
	  THERMAL_ED=THERMAL_ED+(ED_EV(I)+ED_EV(I+1))*(E_EV(I+1)-E_EV(I))
	END DO
	THERMAL_ED=THERMAL_ED*0.5D0/T1
!
	NON_THERMAL_ED=0.0D0
	DO I=1,NKT-1
	  NON_THERMAL_ED=NON_THERMAL_ED+(YE(I,DPTH_INDX)+YE(I+1,DPTH_INDX))*(XKT(I+1)-XKT(I))
	END DO
	NON_THERMAL_ED=NON_THERMAL_ED*0.5D0*(RAD_ENERGY(DPTH_INDX)*6.24150974D+11)
	WRITE(6,*)' '
	WRITE(6,'(3A)')PG_PEN(2),' The non-thermal electron distribution is plotted in red',DEF_PEN
	WRITE(6,'(3A)')PG_PEN(3),' The thermal electron distribution is plotted in blue',DEF_PEN
	WRITE(6,*)'The temperature (10^4K) is ',T(DPTH_INDX)
	WRITE(6,*)ED(DPTH_INDX),THERMAL_ED,NON_THERMAL_ED
	WRITE(6,*)' '
	CALL DP_CURVE(NKT,XKT,YE(1,DPTH_INDX))
	CALL DP_CURVE(NX,E_EV,ED_EV)
	CALL GRAMON_PGPLOT('kT(ev)','YE',' ',' ')
	GOTO 30
!
	END
