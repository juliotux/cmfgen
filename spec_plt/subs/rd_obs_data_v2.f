C
C Routine reads in a SIT (or other observational) data from a file 
C which has the format:
C
C **********************************
C **********************************
C  Several lines of description.
C **********************************
C **********************************
C
C FLUX_UNIT=					!Must be first
C WAVE_UNIT=
C AIR_LAM=
C SCALE_FACTOR=
C
C Where:
C      The Flux Unit can Jy, or ergs/cm^2/s/Ang, or Norm.
C      The WAVE_UNIT can be Anstroms
C      AIR_LAM can be TRUE (air) or FALSE (vac)
C
C If WAVE_UNIT is not specifed, Angstoms is assumed
C If AIR_LAM is not specifed, AIR is assumed if L(MAX) > 3500Ang.
C
C WAVE_UNIT, AIR_LAM, and SCAKE_FACTOR can occur in any order, but must
C occur (if present) after FLUX_UNIT. No blank lines are allowed between
c specifiers.
C
C And the a list of wavelenghts (in Angstroms) and Fluxes (in Jy or 
C ergs/cm^2/s/Ang). One pair of data is assumed to be on each
C line.
C 
C The routine returns the wavelngth (in Angstroms) and the flux in 
C Jansky's. NORM is treated like Janskies, and the flux is not altered.
C
C More than one data set can be include in the file. Such a data set is
C separated from the previous data set by a row of * (At least 20). The format
C is then identical to the initial data set, with INDEPENDENT specifiers.
C The column format must be identical.
C
	SUBROUTINE RD_OBS_DATA_V2(LAM,FLUX,NMAX,NPTS,
	1                         FILENAME,COLS,IOS)
	IMPLICIT NONE
C
C Altered 07-Jul-2011 : Improved error message when reading bad data.
C Altered 11-May-2008 : Altered to handle blank lines/comments at end of file.
C Altered 18-Nov-1999 : Altered to allow for multiple data sets in the same
C                        file.
C Altered 18-Jun-1995 : COLS inserted to allow reading of data in different
C                        columns. Changed to _V2
C Altered 29-May-1997 : CONFUSE_CNT variable installed.
C Altered 11-Apr-1997 : Optional WAVE_UNIT and AIR_LAM units installed.
C                         Now allow conversion from AIR to VAC wavelengths.
C Altered 02-May-1996 : IOS installed in CALL.
C
	INTEGER NMAX
	INTEGER NPTS
	INTEGER IOS
C
C COLS(1) indicates which column the wavelength information is in.
C COLS(2) indicates which column the flux information is in.
C 
	INTEGER COLS(2)
	REAL*4 LAM(NMAX),FLUX(NMAX)
	CHARACTER*(*) FILENAME
!
! Local variables and arrays.
!
	REAL*8 TEMP_STORE( MAX(COLS(1),COLS(2)) )
	REAL*8 TO_JANSKY
	REAL*8 LAM_ST
	REAL*8 DEL_LAM
	REAL*8 T1
	INTEGER I,J,K
	INTEGER CONFUSE_CNT
	INTEGER N_STR
	INTEGER NLST
	CHARACTER*80 STRING(10)
	CHARACTER*80 WAVE_UNIT
	CHARACTER*80 AIR_LAM
	CHARACTER*80 FLUX_UNIT
	CHARACTER*20 FLUX_KEY_WORD
	CHARACTER*20 DATA_FORM
	LOGICAL FINISHED
C
	INTEGER, PARAMETER :: T_OUT=6
	INTEGER, PARAMETER :: IZERO=0
C
	REAL*8 VAL_SCALE_FAC
	REAL*8 LAM_VAC
	EXTERNAL LAM_VAC
C
	IOS=0
	OPEN(UNIT=10,FILE=FILENAME,ACTION='READ',STATUS='OLD',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Unable to open file in RD_OBS_DATA'
	  RETURN
	END IF
C
	NPTS=0
	TO_JANSKY=1.0D+23*1.0E-08/2.998E+10
C
	NLST=0
5000	CONTINUE
	FINISHED=.FALSE.
C
C Get default flux unit.
C
	STRING(1)=' '
	DO WHILE(INDEX(STRING(1),'FLUX_UNIT=') .EQ. 0 .AND.
	1                INDEX(STRING(1),'FLUX_UNIT_2=') .EQ. 0)
	  WRITE(T_OUT,'(A)')STRING(1)
	  READ(10,'(A)',IOSTAT=IOS)STRING(1)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'FLUX_UNIT not found in file for RD_OBS_DATA'
	    WRITE(T_OUT,'(A)')STRING(1)
	    IOS=1
	    CLOSE(UNIT=10)
	    RETURN
	  END IF
	  I=I+1
	END DO
	FLUX_UNIT=STRING(1)
	WRITE(T_OUT,'(A)')' '
C
C Read in all keywords.
C
	N_STR=0
	DO WHILE(INDEX(STRING(N_STR+1),'=') .NE. 0)
	  N_STR=N_STR+1
	  READ(10,'(A)')STRING(N_STR+1)
	END DO
	BACKSPACE(UNIT=10)
C
C If the flux data is not in column 2 (the old default), we do a search
C for the unit associated with the data column.
C
	IF(COLS(2) .NE. 2)THEN
	  FLUX_KEY_WORD='FLUX_UNIT_'
	  IF(COLS(2) .LT. 10)THEN
	      WRITE(FLUX_KEY_WORD(11:12),'(I1,A1)')COLS(2),'='
	  ELSE IF(COLS(2) .LT. 100)THEN
	      WRITE(FLUX_KEY_WORD(11:13),'(I2,A1)')COLS(2),'='
	  END IF
	  K=2
	  DO WHILE (K .LE. N_STR)
	    IF(INDEX(STRING(K),TRIM(FLUX_KEY_WORD)) .NE. 0)THEN
	      FLUX_UNIT=STRING(K)
	      K=N_STR+1			!Finish loop
	    END IF
	    K=K+1
	  END DO
	END IF
C
	WAVE_UNIT='ANGSTROMS'				!Default
	DO K=2,N_STR
	  IF(INDEX(STRING(K),'WAVE_UNIT=') .NE. 0)THEN
	    I=INDEX(STRING(K),'=')
	    WAVE_UNIT=STRING(K)(I+1:)
	    DO WHILE(WAVE_UNIT(1:1) .EQ. ' ')
	      WAVE_UNIT(2:)=WAVE_UNIT(1:)
	    END DO
	    I=INDEX(WAVE_UNIT,' ')
	    WAVE_UNIT=WAVE_UNIT(1:I-1)
	    CALL SET_CASE_UP(WAVE_UNIT,IZERO,IZERO)
	    IF(WAVE_UNIT .EQ. 'ANGSTROMS' )THEN
	    ELSE IF(WAVE_UNIT .EQ. 'MICROMETERS')THEN
	    ELSE IF(WAVE_UNIT .EQ. 'UM')THEN
	        WAVE_UNIT='MICROMETERS'
	    ELSE IF(WAVE_UNIT .EQ. 'HZ')THEN
	    ELSE
	      WRITE(T_OUT,*)'Error -- only Angstroms, micrometers or Hz ',
	1                   ' handled in RD_OBS_DATA'
	      WRITE(T_OUT,*)'Edit RD_OBS_DATA to handle other units'
	      CLOSE(UNIT=10)
	      RETURN
	    END IF
	  END IF
	END DO
C
	AIR_LAM='UNKNOWN'
	DO K=2,N_STR
	  IF(INDEX(STRING(K),'AIR_LAM=') .NE. 0)THEN
	    I=INDEX(STRING(K),'=')
	    AIR_LAM=STRING(K)(I+1:)
	    DO WHILE(AIR_LAM(1:1) .EQ. ' ')
	      AIR_LAM(2:)=AIR_LAM(1:)
	    END DO
	    I=INDEX(AIR_LAM,' ')
	    AIR_LAM=AIR_LAM(1:I-1)
	    CALL SET_CASE_UP(AIR_LAM,IZERO,IZERO)
	    IF(INDEX(AIR_LAM,'FALSE') .NE. 0)THEN
	      AIR_LAM='FALSE'
	    ELSE IF(INDEX(AIR_LAM,'TRUE') .NE. 0)THEN
	      AIR_LAM='TRUE'
	    ELSE
	      WRITE(T_OUT,*)'Error -- dont recognize AIR_LAM option in',
	1                    ' RD_OBS_DATA'
	      CLOSE(UNIT=10)
	      RETURN
	    END IF
	  END IF
	END DO
!
	DATA_FORM=' '
	DO K=2,N_STR
	  IF(INDEX(STRING(K),'DATA_FORM=') .NE. 0)THEN
	    I=INDEX(STRING(K),'=')
	    DATA_FORM=STRING(K)(I+1:)
	    DO WHILE(DATA_FORM(1:1) .EQ. ' ')
	      DATA_FORM(2:)=DATA_FORM(1:)
	    END DO
	    I=INDEX(DATA_FORM,' ')
	    DATA_FORM=DATA_FORM(1:I-1)
	    CALL SET_CASE_UP(DATA_FORM,IZERO,IZERO)
	    IF(INDEX(DATA_FORM,'HR_IUE') .NE. 0)THEN
	      DATA_FORM='HR_IUE'
	    ELSE
	      WRITE(T_OUT,*)'Error -- dont recognize DATA_FORM option in',
	1                    ' RD_OBS_DATA'
	      CLOSE(UNIT=10)
	      RETURN
	    END IF
	  END IF
	END DO
!
	IF(DATA_FORM .EQ. 'HR_IUE')THEN
	  LAM_ST=1
          DO K=2,N_STR
	    IF(INDEX(STRING(K),'LAM_ST=') .NE. 0)THEN
	      I=INDEX(STRING(K),'=')
	      READ(STRING(K)(I+1:),*)LAM_ST
	      EXIT
	    END IF
	  END DO
!
	  DEL_LAM=1
          DO K=2,N_STR
	    IF(INDEX(STRING(K),'DEL_LAM=') .NE. 0)THEN
	      I=INDEX(STRING(K),'=')
	      READ(STRING(K)(I+1:),*)DEL_LAM
	      EXIT
	    END IF
	  END DO
!
	  NPTS=1
          DO K=2,N_STR
	    IF(INDEX(STRING(K),'NPIX=') .NE. 0)THEN
	      I=INDEX(STRING(K),'=')
	      READ(STRING(K)(I+1:),*)NPTS
	      EXIT
	    END IF
	  END DO
	END IF
C
C Get scale factor if present.
C
	VAL_SCALE_FAC=1.0D0
	DO K=2,N_STR
	  IF(INDEX(STRING(K),'SCALE_FACTOR=') .NE. 0)THEN
	    I=INDEX(STRING(K),'=')
	    READ(STRING(K)(I+1:),*)VAL_SCALE_FAC
	  END IF
	END DO
C
	IF(DATA_FORM .EQ. 'HR_IUE')THEN
	  DO I=NLST+1,NLST+NPTS
	    READ(10,*)FLUX(I)
	    IF(FLUX(I) .LT. 0)FLUX(I)=0.0D0
	    LAM(I)=LAM_ST+(I-NLST-1)*DEL_LAM
	  END DO
	ELSE IF(COLS(1) .EQ. 1 .AND. COLS(2) .EQ. 2)THEN
	  DO I=NLST+1,NMAX
10	    READ(10,'(A)',END=1000)STRING(1)
	    IF(INDEX(STRING(1),'********************') .NE. 0)GOTO 2000
	    IF(INDEX(STRING(1)(1:1),'!').NE. 0)GOTO 10
	    IF(STRING(1) .EQ. ' ')GOTO 10
	    READ(STRING(1),*,IOSTAT=IOS)LAM(I),FLUX(I)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Error reading data. I index=',I
	      EXIT
	    END IF
	    NPTS=I
	  END DO
	ELSE
	  K=MAX(COLS(1),COLS(2))
	  DO I=NLST+1,NMAX
20	    READ(10,'(A)',END=1000)STRING(1)
	    IF(INDEX(STRING(1),'********************') .NE. 0)GOTO 2000
	    IF(INDEX(STRING(1)(1:1),'!').NE. 0)GOTO 20
	    IF(STRING(1) .EQ. ' ')GOTO 20
	    READ(STRING(1),*,END=1000,IOSTAT=IOS)(TEMP_STORE(J),J=1,K)
	    IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error reading observational data. String is:'
	     WRITE(T_OUT,*)TRIM(STRING(1))
	     RETURN
	    END IF
	    LAM(I)=TEMP_STORE(COLS(1))
	    FLUX(I)=TEMP_STORE(COLS(2))
	    NPTS=I
	  END DO
	END IF
!
! If we reach here we have read in NMAX points. If we can still do another
! successfull read we have more data to read in.
!
	READ(10,*,END=1000)T1,T1
	WRITE(T_OUT,*)'Warning: ALL the data has not been read in'
1000	CONTINUE
	FINISHED=.TRUE.
2000    CONTINUE
C
C Perform scaling.
C
	FLUX(NLST+1:NPTS)=FLUX(NLST+1:NPTS)*VAL_SCALE_FAC
C
	IF(WAVE_UNIT .EQ. 'ANGSTROMS')THEN
C
C Do nothing as unit we want.
C
	ELSE IF(WAVE_UNIT .EQ. 'MICROMETERS')THEN
	  LAM(NLST+1:NPTS)=1.0D+04*LAM(NLST+1:NPTS)
	ELSE IF(WAVE_UNIT .EQ. 'HZ')THEN
	  LAM(NLST+1:NPTS)=2.99794E+18/LAM(NLST+1:NPTS)
	  AIR_LAM='FALSE'
	END IF
C
	IF(INDEX(FLUX_UNIT,'ergs/cm^2/s/Ang') .NE. 0)THEN
	  DO I=NLST+1,NPTS
	    FLUX(I)=TO_JANSKY*FLUX(I)*(LAM(I)**2)
	  END DO
C
	ELSE IF(INDEX(FLUX_UNIT,'ergs/cm^2/s/Hz') .NE. 0)THEN
	  DO I=NLST+1,NPTS
	    FLUX(I)=1.0D+23*FLUX(I)
	  END DO
	ELSE IF(INDEX(FLUX_UNIT,'mJy') .NE. 0 .OR.
	1          INDEX(FLUX_UNIT,'milli-Jansky') .NE. 0)THEN
	  FLUX(NLST+1:NPTS)=1.0D-03*FLUX(NLST+1:NPTS)
C
	ELSE IF(INDEX(FLUX_UNIT,'Jy') .NE. 0 .OR.
	1          INDEX(FLUX_UNIT,'Jansky') .NE. 0)THEN
C
C Do nothing as unit we want.
C
	ELSE IF(INDEX(FLUX_UNIT,'Norm') .NE. 0 .OR.
	1          INDEX(FLUX_UNIT,'Jansky') .NE. 0)THEN
	  WRITE(T_OUT,*)'Warning --- normalized data'
C
	ELSE
	  WRITE(T_OUT,*)'Invalid flux unit in RD_OBS_DAT'
	  WRITE(T_OUT,*)'DATA MAY BE GARBAGE'
	END IF
C
	IF(AIR_LAM .EQ. 'UNKNOWN')THEN
	  IF(MAX(LAM(NLST+1),LAM(NPTS)) .GT. 3500.0)THEN
	    AIR_LAM='TRUE'
	    WRITE(T_OUT,*)' Wavelengths assumed to be in AIR'
	  ELSE
	    AIR_LAM='FALSE'
	    WRITE(T_OUT,*)' Wavelengths assumed to be in VACUUM'
	  END IF
	END IF
	IF(AIR_LAM .EQ. 'TRUE')THEN
	  CONFUSE_CNT=0
	  DO I=NLST+1,NPTS
	    IF(LAM(I) .GT. 2000)THEN
	      T1=LAM(I)
	      LAM(I)=LAM_VAC(T1)		!Double precis arg.
	    ELSE
	      CONFUSE_CNT=CONFUSE_CNT+1
            END IF
	  END DO
	  IF(CONFUSE_CNT .NE. 0)THEN
	      WRITE(T_OUT,*)' '
	      WRITE(T_OUT,*)('*',I=1,50)
	      WRITE(T_OUT,*)('*',I=1,50)
	      WRITE(T_OUT,*)'Warning: possible confused AIR and VAC',
	1                   ' wavelengths in RD_OBS_DATA'
	      WRITE(T_OUT,*)CONFUSE_CNT,' wavelengths less than 2000Ang',
	1                     ' assumed to be vacuum'
	      WRITE(T_OUT,*)('*',I=1,50)
	      WRITE(T_OUT,*)('*',I=1,50)
	      WRITE(T_OUT,*)' '
	  END IF
	END IF
!
! Check to see if we have finished reading all the data. If not we continue.
!
	IF( .NOT. FINISHED)THEN
	  NLST=NPTS
	  GOTO 5000
	END IF
!
	CLOSE(UNIT=10)
!
	RETURN
	END
