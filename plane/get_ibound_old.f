	SUBROUTINE GET_IBOUND(IBOUND,MU,NP,FREQ,INIT)
!
	INTEGER NP
	REAL*8 MU(NP)
	REAL*8 IBOUND(NP)
	REAL*8 FREQ
	LOGICAL INIT
!
	REAL*8, SAVE, ALLOCATABLE :: INTENSITY(:,:)
	REAL*8, SAVE, ALLOCATABLE :: NU_GRID(:)
	REAL*8, SAVE, ALLOCATABLE :: ANG_GRID(:)
	INTEGER, SAVE, ALLOCATABLE :: MU_INDX(:)
!
	INTEGER, SAVE :: NCF		!Number of incident frequencies
	INTEGER, SAVE :: NANG		!Number of incident angles
	INTEGER, SAVE:: CUR_FREQ_INDX
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
! Local variables
!
	REAL*8 IMU1,IMU2 
	REAL*8 T1,T2
	INTEGER I,J,ML,LU,LS
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER*80 STRING

	IF(FIRST_TIME)THEN
	  LU=7
	  OPEN(FILE='INCIDENT_INTENSITY',UNIT=LU,STATUS='OLD',ACTION='READ')
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	      READ(LU,'(A)')STRING
	    END DO
	    READ(STRING,*)NCF,NANG
	    ALLOCATE (NU_GRID(NCF))
	    ALLOCATE (ANG_GRID(NANG))
	    ALLOCATE (INTENSITY(NANG,NCF))
	    READ(LU,*)(ANG_GRID(J),J=1,NANG)
	    DO I=1,NCF
	      READ(LU,*)(NU_GRID(I),INTENSITY(J,I),J=1,NANG)
	    END DO
	  CLOSE(UNIT=LU)
	END IF
!
! Find locations for MU interpolation. Note that this must be done
! independent of how MU is ordered.
!
	IF(INIT)CUR_NU_INDX=1
	IF(FIRST_TIME .OR. NP .NE. NP_SAV)THEN
	  IF(ALLOCATED(MU_INDX))DEALLOCATE(MU_INDX)
	  ALLOCATE (MU_INDX(NP))
	  DO LS=1,NP
	    J=1
	    DO WHILE ( (MU(LS)-ANG_GRID(J))*(NU_GRID(J+1)-MU(LS)) .LT. 0.0D0 )
	      J=J+1
	      IF(J .GE. NANG)THEN
	        LUER=ERROR_LU()
	        WRITE(LUER,*)'Invalid incident MU range'
	        WRITE(LUER,*)'MU(LS)=',MU(LS)
	        WRITE(LUER,*)'ANG_GRID(1)=',ANG_GRID(1)
	        WRITE(LUER,*)'ANG_GRID(NUM_GRID)=',ANG_GRID(NANG)
	        STOP
	      END IF
            END DO
	    MU_INDX(LS)=J
	  END DO
	END IF
!
! Find frequencey location, and do simple linear interpolation.
! As we come into this routine, frequencies should be ordered. Thus
! we always search from last frequency index, unless starting again.
!
	IF(FREQ .LT. INTEN_NU(1) .OR. FREQ .GT. INTEN_NU(NCF))THEN
	  IBOUND=0.0D0
	ELSE
	  DO WHILE(FREQ .LT. INTEN_NU(CUR_NU_INDX+1))
	    CUR_NU_INDX=CUR_NU_INDX+1
	  END DO
	  ML=CUR_NU_INDX
	  T1=(INTEN_NU(ML)-FREQ)/(INTEN_NU(ML)-INTEN_NU(ML+1))
	  DO LS=1,NP	
	    J=MU_INDX(LS)
	    IMU1=(1.0D0-T1)*INTENSITY(J,ML)+T1*INTENSITY(J,ML+1)
	    IMU2=(1.0D0-T1)*INTENSITY(J+1,ML)+T1*INTENSITY(J+1,ML+1)
	    T2=(INTEN_MU(J)-MU(LS))/(INTEN_MU(J)-INTEN_MU(J+1))
	    IBOUND(LS)=(1.0D0-T2)*IMU1+T2*IMU2
	  END DO
	END IF
!
	RETURN
	END
