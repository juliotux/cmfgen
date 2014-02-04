	SUBROUTINE PLT_ARN(XSPEC,ND,XV,YV,NYV)
	USE MOD_DISP
	USE MOD_NON_THERM
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NYV
	REAL*8 XV(NYV)
	REAL*8 YV(NYV)
	CHARACTER(LEN=*) XSPEC
!
	INTEGER ID
	INTEGER IT
	LOGICAL, SAVE :: NOT_READ_ARNAUD=.TRUE.
	CHARACTER(LEN=3), PARAMETER :: XKT_METHOD='lin'
!
	CHARACTER(LEN=30) UC
	EXTERNAL UC
!
! We now set these values by values passed from CMFGEN.
!
	INTEGER IOS
	INTEGER, SAVE :: NT_NKT=1000
	REAL*8, SAVE :: XKT_MIN=1.0D0
	REAL*8, SAVE :: XKT_MAX=1000.0D0
!
	IF(NOT_READ_ARNAUD)THEN
	   CALL READ_ARNAUD_ION_DATA_DISP(ND)
	   WRITE(6,*)'Succesfull read in non-thermal ioization cross-sections'
	END IF
!
	IF(.NOT. ALLOCATED(XKT))THEN
	  NKT=NT_NKT
	  WRITE(6,*)'Allocating memory to describe non-thermal electron distribution'
	  ALLOCATE (XKT(NKT),STAT=IOS)
	  ALLOCATE (dXKT(NKT),STAT=IOS)
          XKT(1:NKT)=0.0D0; dXKT(1:NKT)=0.0D0
	  CALL SET_XKT_ARRAY(XKT_MIN,XKT_MAX,NKT,XKT,dXKT,XKT_METHOD)
	  WRITE(6,*)'Succesfull allocated memory and set XKT'
!
!
! If needed, allocate memory for collisional ioinzation cross-sections.
!
          IOS=0
          DO IT=1,NUM_THD
            IF(THD(IT)%PRES)ALLOCATE (THD(IT)%CROSS_SEC(NKT),STAT=IOS)
              IF(IOS.NE.0) THEN
                 WRITE(6,*) 'Error allocating THD%CROSS_SEC',IOS
                 STOP
              END IF
           END DO
	   WRITE(6,*)'Successfully allocated memory for cross-sections'
!
	   CALL ARNAUD_CROSS_V3_DISP()
	END IF
!
	WRITE(6,'(3X,A,5X,A,5X,A,3X,A,3X,A)')'ID','Ion','IT','# of routes'
	DO ID=1,NUM_IONS
	  IF(XSPEC .EQ. UC(ION_ID(ID)))THEN
	    DO IT=1,MAX_NUM_THD
	      IF(THD(IT)%LNK_TO_ION .EQ. ID)THEN
	        WRITE(6,'(2X,I3,1X,A7,3X,I3,3X,I2,3X)')ID,TRIM(ION_ID(ID)),IT,THD(IT)%N_ION_ROUTES
	        YV(1:NKT)=THD(IT)%CROSS_SEC(1:NKT)
	        CALL DP_CURVE(NKT,XKT,YV)
	      END IF
	    END DO
	  END IF
	END DO
!
	RETURN
	END
