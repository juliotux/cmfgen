!
! Read in collsion strengths from a file. The values are tabulated as a 
! function of temperature.
!
! OMEGA_SCALE is used by SUBCOL_MULTI to scale all collsion strengths by a
! constant value.
!
	PROGRAM SCALE_OMEGA_V1
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created 25-Mar-2011 : Based on gen_omega_rd_v2.f
!
	REAL*8 OMEGA_SCALE
	REAL*8 OMEGA_SET
	REAL*8 SCALE_FACTOR
!
! T_TABLE: Tempearure values (10^4 K) at which collsion strengths tabulated.
!
	INTEGER, PARAMETER :: MAX_TVALS=100
	REAL*8 T_TABLE(MAX_TVALS)
	REAL*8 OMEGA(MAX_TVALS)
!
	INTEGER NUM_TRANS
	INTEGER NUM_TVALS
	INTEGER TRANS_CNT
!
	INTEGER, PARAMETER :: LUER=6
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUOUT=10
!
! Local variables.
!
	INTEGER I		!Used for lower level.
	INTEGER K
	INTEGER L
	INTEGER IOS
!
	CHARACTER(LEN=80) IN_FILE
	CHARACTER(LEN=80) OUT_FILE
	CHARACTER(LEN=500) STRING
!
! Open file with collisional data.
!
	IN_FILE=' '
	CALL GEN_IN(IN_FILE,'File with collisional data')
	OPEN(LUIN,FILE=IN_FILE,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',IN_FILE
	  WRITE(LUER,*)'IOS=',IOS	
	  STOP
	END IF
!
	OUT_FILE=TRIM(IN_FILE)//'_FUDGED'
	CALL GEN_IN(OUT_FILE,'Output file for revided collisional data')
	OPEN(LUOUT,FILE=OUT_FILE,STATUS='UNKNOWN',ACTION='WRITE',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening ',TRIM(OUT_FILE)
	  WRITE(LUER,*)'IOS=',IOS	
	  STOP
	END IF
!
	SCALE_FACTOR=1.0D0
	CALL GEN_IN(SCALE_FACTOR,'Factor to fudge collisional data')
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Number of transitions')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',IN_FILE
	    WRITE(LUER,*)'!Numer of transitions string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	  WRITE(LUOUT,'(A)')TRIM(STRING)
	END DO
!
	READ(STRING,*)NUM_TRANS
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Number of T values OMEGA tabulated at')  
	1               .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',IN_FILE
	    WRITE(LUER,*)'!Numer of T values string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	  WRITE(LUOUT,'(A)')TRIM(STRING)
	END DO
	READ(STRING,*)NUM_TVALS
	IF(NUM_TVALS .GT. MAX_TVALS)THEN
	  WRITE(LUER,*)'Error reading collisonal data from '//IN_FILE
	  WRITE(LUER,*)'NU_TVALS too small'
	  STOP
	END IF
!         
! OMEGA_SCALE provides a means of scaling the OMEGA that have not been read in.
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Scaling factor for OMEGA') .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',IN_FILE
	    WRITE(LUER,*)'!Scaling factor for Omega string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	  WRITE(LUOUT,'(A)')TRIM(STRING)
	END DO
	READ(STRING,*)OMEGA_SCALE
!         
! OMEGA_SET provides a means to set OMEGA for transition for which
! no atomic data is availaable, and for which the oscilator strength
! is zero.
!
	STRING=' '
	DO WHILE( INDEX(STRING,'!Value for OMEGA if f=0') .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',IN_FILE
	    WRITE(LUER,*)'!Value for OMEGA if f=0 string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	  WRITE(LUOUT,'(A)')TRIM(STRING)
	END DO
	READ(STRING,*)OMEGA_SET
!
! We do this check here so that OMEGA_SCALE and OMEGA_SET are defined.
!
	IF(NUM_TRANS .EQ. 0)THEN
	  CLOSE(LUIN)
	  STOP			!i.e use approximate formulae only.
	END IF
!
	STRING=' '
	DO WHILE( INDEX(STRING,'Transition\T')  .EQ. 0 )
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading ',IN_FILE
	    WRITE(LUER,*)'Transition\T string not found'
	    WRITE(LUER,*)'IOS=',IOS
	    STOP
	  END IF
	  WRITE(LUOUT,'(A)')TRIM(STRING)
	END DO
	L=INDEX(STRING,'ion\T')
	READ(STRING(L+5:),*)(T_TABLE(I),I=1,NUM_TVALS )
!
	STRING=' '
	DO WHILE(STRING .EQ. ' ')
	  READ(LUIN,'(A)')STRING
	  IF(STRING .EQ. ' ')WRITE(LUOUT,'(A)')TRIM(STRING)
	END DO             
	BACKSPACE(LUIN)
!
! TRANS_CNT is the index used to read in all collisional data. 
!
	DO TRANS_CNT=1,NUM_TRANS
!
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ')
	    READ(LUIN,'(A)')STRING
	    IF(STRING .EQ. ' ')WRITE(LUOUT,'(A)')TRIM(STRING)
	  END DO
	  K=INDEX(STRING,'-')+1
	  DO I=K,K+50
	    L=I
	    IF(STRING(L:L) .NE. ' ')EXIT
	  END DO
	  K=INDEX(STRING(L:),' ')+L
	  DO I=K,K+50
	    L=I
	    IF(STRING(L:L) .NE. ' ')EXIT
	  END DO
	  READ(STRING(L:),*,IOSTAT=IOS)(OMEGA(I),I=1,NUM_TVALS)
	  DO I=1,NUM_TVALS
	    OMEGA(I)=SCALE_FACTOR*OMEGA(I)
	  END DO
!
	  IF(MINVAL(OMEGA(1:NUM_TVALS)) .GT. 0.01 .AND. MAXVAL(OMEGA(1:NUM_TVALS)) .LT. 9999.0D0)THEN
	    WRITE(LUOUT,'(A,30F10.4)')STRING(1:L-1),(OMEGA(I),I=1,NUM_TVALS)
	  ELSE
	    WRITE(LUOUT,'(A,30ES10.2)')STRING(1:L-1),(OMEGA(I),I=1,NUM_TVALS)
	  END IF
	END DO
!
	STOP
	END
