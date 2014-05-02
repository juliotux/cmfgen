	MODULE MOD_XRAY_FLUXES
	IMPLICIT NONE
!
	REAL*8, ALLOCATABLE :: X_NU(:)
	REAL*8, ALLOCATABLE :: LOG_X_TEMP(:)
	REAL*8, ALLOCATABLE :: LOG_X_ED(:)
	REAL*8, ALLOCATABLE :: XRAY_FLUXES(:,:,:)
!
	REAL*8, ALLOCATABLE :: X_EMISS1(:)
	REAL*8, ALLOCATABLE :: X_EMISS2(:)
!               
	REAL*8 BIN_MIN
	REAL*8 BIN_SIZE
	REAL*8 LOG_T_MIN
	REAL*8 DEL_LOG_T
	REAL*8 LOG_ED_MIN
	REAL*8 DEL_LOG_ED
!
	REAL*8 T_SHOCK1_SAV
	REAL*8 T_SHOCK2_SAV
!
	INTEGER N_BINS
	INTEGER N_TEMP
	INTEGER N_ED
!
	END MODULE MOD_XRAY_FLUXES
!
!
! Routine to return X-ray EMISSIVITIES for a set of NFREQ frequencies,
! and for 2 different shock temperatures (in units of 10^4 K).
! The returned emissivities have units ergs/cm^3/s/Hz/steradian.
! At present, we assume that the X-ray emissivity is independent of density.
!
! On the first call the tabulated RS data is read in from a file 
! RS_XRAY_FLUXES. FREQ may be monotonically increasing, or decreasing.
!
	SUBROUTINE RD_XRAY_SPEC(T_SHOCK1,T_SHOCK2,LU_IN)
	USE MOD_XRAY_FLUXES
	IMPLICIT NONE
!
	INTEGER LU_IN
	REAL*8 T_SHOCK1		!In units of 10^4 K
	REAL*8 T_SHOCK2
!
! Local variables:
!
	REAL*8, PARAMETER :: EV_TO_HZ=0.241838D0
	REAL*8 T1
!
	REAL*8 LOG_T_SHOCK1
	REAL*8 LOG_T_SHOCK2
!
	INTEGER I,J
	INTEGER IOS
	INTEGER T_INDX
	INTEGER ED_INDX
!
	CHARACTER*132 STRING
!
	REAL*8 FUN_PI,PI
	INTEGER ERROR_LU,LU_ER
	EXTERNAL ERROR_LU,FUN_PI
!
	LOGICAL FIRST_TIME
	DATA FIRST_TIME/.TRUE./
	SAVE FIRST_TIME
!
	LU_ER=ERROR_LU()
	IF(FIRST_TIME)THEN
	  FIRST_TIME=.FALSE.
!
! Read in RS data table. The fluxes are assumed to be units of
! 10^{-23} ergs/cm3/sec. They represent the cooling function per
! unit density of electrons and H ions.
!
	  WRITE(LU_ER,'(A)')' '
	  OPEN(UNIT=LU_IN,FILE='RS_XRAY_FLUXES',ACTION='READ',
	1       STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Error opening RS_XRAY_FLUXES'
	    WRITE(LU_ER,*)'Error occurred in RD_XRAY_SPEC'
	    STOP
	  END IF
!
! RS data assumed to be tabulated in bins equally spaced in eV, and ordered
! in increasing eV.
!
	  CALL RD_INT(N_BINS,'N_BINS',LU_IN,LU_ER,'# freq bins')
	  CALL RD_DBLE(BIN_MIN,'BIN_MIN',LU_IN,LU_ER,' ')
	  CALL RD_DBLE(BIN_SIZE,'BIN_SIZE',LU_IN,LU_ER,' ')
	  BIN_MIN=BIN_MIN*EV_TO_HZ	!Convert to units of 10^15 HZ
	  BIN_SIZE=BIN_SIZE*EV_TO_HZ
! 
! Temperature tabulated in equal increments of Log T.
!
	  CALL RD_INT(N_TEMP,'N_TEMP',LU_IN,LU_ER,'# freq bins')
	  CALL RD_DBLE(LOG_T_MIN,'LOG_T_MIN',LU_IN,LU_ER,' ')
	  CALL RD_DBLE(DEL_LOG_T,'DEL_LOG_T',LU_IN,LU_ER,' ')
	  LOG_T_MIN=LOG_T_MIN-4.0D0		!Convert from K to units of 10^4 K
! 
! Temperature tabulated in equal increments of Log Ne.
!
	  CALL RD_INT(N_ED,'N_ED',LU_IN,LU_ER,'# freq bins')
	  CALL RD_DBLE(LOG_ED_MIN,'LOG_ED_MIN',LU_IN,LU_ER,' ')
	  CALL RD_DBLE(DEL_LOG_ED,'DEL_LOG_ED',LU_IN,LU_ER,' ')
	  WRITE(LU_ER,'(A)')' '
!
! Now that we have the vector sizes, we cab allocate memory for the X-ray 
! table and vectors.
!
	  ALLOCATE (X_NU(N_BINS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (X_EMISS1(N_BINS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (X_EMISS2(N_BINS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LOG_X_TEMP(N_TEMP),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LOG_X_ED(N_ED),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (XRAY_FLUXES(N_BINS,N_TEMP,N_ED),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Error allocating memory (1)'
	    WRITE(LU_ER,*)'Error occurred in RD_XRAY_SPEC'
	    WRITE(LU_ER,*)'STAT=',IOS
	    STOP
	  END IF
!
! Axes computation
!
	  DO I=1,N_BINS
	     X_NU(I)=BIN_MIN+(I-1)*BIN_SIZE
	  END DO
	  DO I=1,N_TEMP
	    LOG_X_TEMP(I)=LOG_T_MIN+(I-1)*DEL_LOG_T
	  END DO	
	  DO I=1,N_ED
       	    LOG_X_ED(I)=LOG_T_MIN+(I-1)*DEL_LOG_ED
	  END DO
!
! Read in fluxes for each parameter set. Blank lines, comments,
! and a single header record are ignored.
!
	  DO J=1,N_ED
	    DO I=1,N_TEMP
	      STRING='!'
	      DO WHILE (STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!'
	1           .OR. INDEX(STRING,'Temperature(K)') .NE. 0)
	        READ(LU_IN,'(A)')STRING
	       END DO
	       BACKSPACE(LU_IN)
	       READ(LU_IN,*)XRAY_FLUXES(:,I,J)
	    END DO
!
! Convert fluxes from units of 10^{-23} ergs/cm^3/s/bin to units of
! ergs/cm^3/s/steradian/Hz. NB: Because of the corrections earlier, BIN_SIZE 
! is in units of 10^15 Hz.
!
! To put into pogram units, we multiply by an additional factor of 10^10.
!
! NB: 1.0D-28 = 1.0E+10*1.0E-23/1.0E+15
!
	    PI=FUN_PI()
	    XRAY_FLUXES(:,:,:)=1.0D-28*XRAY_FLUXES(:,:,:)/(PI*4.0D0)/BIN_SIZE
	  END DO
!
	  CLOSE(LU_IN)
	END IF
!            
	ED_INDX=1
	IF(N_ED .NE. 1)THEN
	  WRITE(LU_ER,'(70A)')('*',I=1,70)
	  WRITE(LU_ER,'(70A)')('*',I=1,70)
	  WRITE(LU_ER,*)'Warning in RD_XRAY_SPEC'
	  WRITE(LU_ER,*)'Data presently assumes no electron density dependence'
	  WRITE(LU_ER,*)'Using RS data for the lowest electron density'
	  WRITE(LU_ER,'(70A)')('*',I=1,70)
	END IF
!
! Compute the X-ray emissivities at the two shock temperatures of interest.
! These emissivities are computed on the X-ray frequency grid. The data is 
! stored for subsequent calls.
!
	IF(T_SHOCK1 .EQ. 0.0D0)THEN
	  X_EMISS1(1:N_BINS)=0.0D0
	ELSE IF(T_SHOCK1 .NE. T_SHOCK1_SAV)THEN
	  LOG_T_SHOCK1=LOG10(T_SHOCK1)
	  T_INDX=INT( (LOG_T_SHOCK1-LOG_T_MIN)/DEL_LOG_T )+1
	  IF(T_INDX .LT. 1 .OR. T_INDX .GT. N_TEMP-1)THEN
	    WRITE(LU_ER,*)'Error: T_SHOCK outside range'
	    WRITE(LU_ER,*)'T_SHOCK=',T_SHOCK1
	    WRITE(LU_ER,*)'Error occurred in RD_XRAY_SPEC'
	    STOP             
	  END IF
	  T1=(LOG_T_SHOCK1-LOG_X_TEMP(T_INDX))/
	1            (LOG_X_TEMP(T_INDX+1)-LOG_X_TEMP(T_INDX))
	  DO I=1,N_BINS
	    X_EMISS1(I)=(1.0D0-T1)*XRAY_FLUXES(I,T_INDX,ED_INDX)
	1               +T1*XRAY_FLUXES(I,T_INDX+1,ED_INDX)
	  END DO
	  T_SHOCK1_SAV=T_SHOCK1
	END IF
!
	IF(T_SHOCK2 .EQ. 0.0D0)THEN
	  X_EMISS2(1:N_BINS)=0.0D0
	ELSE IF(T_SHOCK2 .NE. T_SHOCK2_SAV)THEN
	  LOG_T_SHOCK2=LOG10(T_SHOCK2)
	  T_INDX=INT( (LOG_T_SHOCK2-LOG_T_MIN)/DEL_LOG_T )+1
	  IF(T_INDX .LT. 1 .OR. T_INDX .GT. N_TEMP-1)THEN
	    WRITE(LU_ER,*)'Error: T_SHOCK outside range'
	    WRITE(LU_ER,*)'T_SHOCK=',T_SHOCK2
	    WRITE(LU_ER,*)'Error occurred in RD_XRAY_SPEC'
	    STOP
	  END IF
	  T1=(LOG_T_SHOCK2-LOG_X_TEMP(T_INDX))/
	1            (LOG_X_TEMP(T_INDX+1)-LOG_X_TEMP(T_INDX))
	  DO I=1,N_BINS
	    X_EMISS2(I)=(1.0D0-T1)*XRAY_FLUXES(I,T_INDX,ED_INDX)
	1               +T1*XRAY_FLUXES(I,T_INDX+1,ED_INDX)
	  END DO
	  T_SHOCK2_SAV=T_SHOCK2
	END IF
!
	RETURN
	END
