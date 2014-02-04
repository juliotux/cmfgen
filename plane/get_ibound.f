	SUBROUTINE GET_IBOUND(IBOUND,MU,NP,FREQ,INIT)
	IMPLICIT NONE
!
	INTEGER NP
	REAL*8 MU(NP)
	REAL*8 IBOUND(NP)
	REAL*8 FREQ
	LOGICAL INIT
!
	REAL*8, SAVE, ALLOCATABLE :: INTENSITY(:)
	REAL*8, SAVE, ALLOCATABLE :: NU_GRID(:)
	REAL*8, SAVE, ALLOCATABLE :: ANG_GRID(:)
	REAL*8, SAVE, ALLOCATABLE :: ANG_VARIATION(:)
	REAL*8, SAVE, ALLOCATABLE :: DIST(:)
!
	REAL*8, SAVE:: FREQ_SAV=0.0D0           !Previous frequency
	INTEGER, SAVE :: NCF			!Number of incident frequencies
	INTEGER, SAVE :: NANG			!Number of incident angles
	INTEGER, SAVE :: NP_SAV=0		!Number of angles in model
	INTEGER, SAVE :: CUR_FREQ_INDX
	LOGICAL, SAVE :: FIRST_TIME=.TRUE.
!
! Local variables
!
	REAL*8 T1,T2
	REAL*8 DILUTION_FACTOR
	INTEGER I,J,ML,LU,LS,IOS
	INTEGER CUR_NU_INDX
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER*80 STRING
!
! Read in incident intensity data. We only do this once. File has to have the
! correct format, although an arbitraty number of comments/blanklines can
! be provided at the top of the file. "!Format date" is used to signify the
! beginning of the data.
!
	IF(FIRST_TIME)THEN
	  LU=7
	  LUER=ERROR_LU()
	  OPEN(FILE='INCIDENT_INTENSITY',UNIT=LU,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening INCIDENT_INTENSITY in GET_IBOUND'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      STOP
	    END IF
!
	    STRING=' '
	    DO WHILE(INDEX(STRING,'!Format date') .EQ. 0)
	      READ(LU,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error getting format date in INCIDENT_INTENSITY in GET_IBOUND'
	        WRITE(LUER,*)'IOSTAT=',IOS
	        STOP
	      END IF
	    END DO
!
	    READ(LU,'(A)')STRING
	    IF(INDEX(STRING,'!Number of continuum frequencies') .EQ. 0)THEN
	      WRITE(LUER,*)'Error reading INCIDENT_INTENSITY in GET_IBOUND'
	      WRITE(LUER,*)'Data has incorrect format: reading NCF'
	      STOP
	    END IF
	    READ(STRING,*)NCF
!
	    READ(LU,'(A)')STRING
	    IF(INDEX(STRING,'!Number of angles') .EQ. 0)THEN
	      WRITE(LUER,*)'Error reading INCIDENT_INTENSITY in GET_IBOUND'
	      WRITE(LUER,*)'Data has incorrect format: reading NANG'
	      STOP
	    END IF
	    READ(STRING,*)NANG
!
	    READ(LU,'(A)')STRING
	    IF(INDEX(STRING,'!Dilution factor used for scaling') .EQ. 0)THEN
	      WRITE(LUER,*)'Error reading INCIDENT_INTENSITY in GET_IBOUND'
	      WRITE(LUER,*)'Data has incorrect format: reading dilution factor'
	      STOP
	    END IF
	    READ(STRING,*)DILUTION_FACTOR
!
	    ALLOCATE (NU_GRID(NCF))
	    ALLOCATE (INTENSITY(NCF))
	    ALLOCATE (ANG_GRID(NANG))
	    ALLOCATE (ANG_VARIATION(NANG))
!
	    STRING=' '
	    DO WHILE(INDEX(STRING,'Angle grid') .EQ. 0)
	      READ(LU,'(A)')STRING
	    END DO
	    READ(LU,*)(ANG_GRID(J),J=1,NANG)
!
	    STRING=' '
	    DO WHILE(INDEX(STRING,'Intensity variation with angle') .EQ. 0)
	      READ(LU,'(A)')STRING
	    END DO
	    READ(LU,*)(ANG_VARIATION(J),J=1,NANG)
!
	    DO I=1,NCF
	      READ(LU,*)NU_GRID(I),INTENSITY(I)
	      INTENSITY(I)=DILUTION_FACTOR*INTENSITY(I)
	    END DO
	  CLOSE(UNIT=LU)
	END IF
!
! Find locations for MU interpolation. Note that this must be done
! independent of how MU is ordered.
!
	IF(INIT)THEN
	  CUR_NU_INDX=1
	  FREQ_SAV=1.0D+30			!Big value
	END IF
	IF(FREQ .GT. FREQ_SAV)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Frequencies out of order in GET_IBOUND'
	  WRITE(LUER,*)'Old frequency is',FREQ_SAV
	  WRITE(LUER,*)'New frequency is',FREQ
	  STOP
	END IF
	FREQ_SAV=FREQ
!
	IF(FIRST_TIME .OR. NP .NE. NP_SAV)THEN
	  IF(ALLOCATED(DIST))DEALLOCATE(DIST)
	  ALLOCATE (DIST(NP))
	  DO LS=1,NP
	    J=1
	    DO WHILE ( (MU(LS)-ANG_GRID(J))*(ANG_GRID(J+1)-MU(LS)) .LT. 0.0D0 )
	      J=J+1
	      IF(J .GE. NANG)THEN
	        LUER=ERROR_LU()
	        WRITE(LUER,*)'Invalid incident MU range in GET_IBOUND'
	        WRITE(LUER,*)'MU(LS)=',MU(LS)
	        WRITE(LUER,*)'ANG_GRID(1)=',ANG_GRID(1)
	        WRITE(LUER,*)'ANG_GRID(NUM_GRID)=',ANG_GRID(NANG)
	        STOP
	      END IF
            END DO
	    T2=(ANG_GRID(J)-MU(LS))/(ANG_GRID(J)-ANG_GRID(J+1))
	    DIST(LS)=(1.0D0-T2)*ANG_VARIATION(J)+T2*ANG_VARIATION(J+1)
	  END DO
	END IF
!
! Find frequencey location, and do simple linear interpolation.
! As we come into this routine, frequencies should be ordered. Thus
! we always search from last frequency index, unless starting again.
!
	IF(FREQ .GT. NU_GRID(1) .OR. FREQ .LT. NU_GRID(NCF))THEN
	  IBOUND=0.0D0
	ELSE
	  DO WHILE(FREQ .LT. NU_GRID(CUR_NU_INDX+1))
	    CUR_NU_INDX=CUR_NU_INDX+1
	  END DO
	  ML=CUR_NU_INDX
	  T1=(NU_GRID(ML)-FREQ)/(NU_GRID(ML)-NU_GRID(ML+1))
	  T2=(1.0D0-T1)*INTENSITY(ML)+T1*INTENSITY(ML+1)
	  DO LS=1,NP	
	    IBOUND(LS)=T2*DIST(LS)
	  END DO
	END IF
!
	FIRST_TIME=.FALSE.
	RETURN
	END
