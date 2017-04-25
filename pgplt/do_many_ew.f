!
! Altered:  26-May-2015 : Include the EWG option. Measure the EW for a group of lines from a file.
!
	SUBROUTINE DO_MANY_EW(TYPE_CURVE,NPLTS,MAX_PLTS)
	USE NEW_GEN_IN_INTERFACE
	IMPLICIT NONE
!
	INTEGER NPLTS
	INTEGER MAX_PLTS
	CHARACTER(LEN=2) TYPE_CURVE(MAX_PLTS)
	CHARACTER(LEN=132), SAVE :: ID_FILENAME=' '
!
! Variables used by EWG option	
!
	REAL*4, ALLOCATABLE :: WLINEO(:)		!Array for initial wavelength for each line
	REAL*4, ALLOCATABLE :: WLINEF(:)		!Array for final wavelength for each line
	REAL*4, ALLOCATABLE :: EWGR(:)
	REAL*4, ALLOCATABLE :: CENT(:)
	REAL*4 LRANGE(2)
	REAL*4 DOP_VEL
	REAL*8 TAU_CUT
!
        INTEGER N_LINE_IDS
        LOGICAL OBSERVED_WAVE(5000)
        CHARACTER*10 LINE_ID(5000)
        REAL*4 ID_WAVE(5000)
        REAL*4 ID_WAVE_OFF(5000)
        REAL*4 ID_Y_OFF(5000)
        REAL*4 TAU(5000)
	REAL*4 XPAR(2),YPAR(2),XT(2),YT(2)
	REAL*4 T1,T2,T3
!
	INTEGER I,J
	INTEGER CNT
	INTEGER IOS
	INTEGER N_IONNAME
	INTEGER NLINES
!
	LOGICAL GEN_FILE
	LOGICAL FILE_PRES
	CHARACTER(LEN=10) IONNAME(10)
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=80) FILNAME
	CHARACTER(LEN=80) ID_FILNAME
	CHARACTER(LEN=80) TMP_STR
	CHARACTER(LEN=80) FILE_OPTION
	CHARACTER(LEN=20), ALLOCATABLE :: INAMES(:)	!Array of ion names of each line considered for ew
!
	INTEGER, PARAMETER :: T_OUT=6
!	
! EWG Option compute the EW for a group of lines from a file. Assuming cont=1 !!!!
! They are two options:
! 	1) It's possible to generate the file using the LINE_IN (or other file with 
!	   same format) from DISPGEN, the ion name (e.g. SkIII), the wavelength 
!	   range, minimum optical depth at the center (intensity for plnid) and
!	   a doppler width in km/s that define the spectral width taken account
!	   to calculate the EW.
!	2) Only to input the file name. The file must have the following format:
!	   !
!	   ! Ion name    Wo   Wf   Lambda
!	   !
!	   FeIII       2212    2232  2220
!           ....        ...    ...   ...
!	   where Wo is the initial wavelength to measure the EW, Wf is the final one
!	   and Lambda is the central VACUUM wavelength of the line. 
!
!	Raul E. Puebla (27 - May - 2015)
!
	IF(ID_FILENAME .EQ. ' ')ID_FILENAME='LINE_ID'
	GEN_FILE=.TRUE.
	FILNAME='LINE_RANGES'
	CALL NEW_GEN_IN(GEN_FILE,'Do you want to create a file from LINE_ID?')
!	
! This section write the file that will be read below to calculate EWs for each line.
! Doppler width (km/s) is an ad hoc value. Is better to explore the regions 
! in order ensure that the whole line and wigns are taken account. 
! THe file containing the lines id MUST have the same format of LINE_ID from dispgen.
! The number of ions must be in the same format used in CMFGEN (e.g. Nk2 for Ni II).
! (Puebla 27-May-2015) 
!
	IF (GEN_FILE) THEN
	  CALL NEW_GEN_IN(ID_FILNAME,'File with line IDs')
	  CALL SET_CASE_UP(ID_FILNAME,1,0)
	  OPEN(UNIT=33,FILE=TRIM(ID_FILNAME),STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file',ID_FILNAME
	    GOTO 1000
	  END IF
	  DOP_VEL=20.0
	  CALL NEW_GEN_IN(DOP_VEL,'Doppler width (km/s)')
	  LRANGE(1)=4000.0
	  LRANGE(2)=5000.0
	  CALL NEW_GEN_IN(LRANGE(1),'Initial wavelength (A):')
	  CALL NEW_GEN_IN(LRANGE(2),'Final wavelength (A):')
	  TAU_CUT=0.1
	  CALL NEW_GEN_IN(TAU_CUT,'Omit lines with central optical depth <')
	  DO I=1,10
	    IONNAME(I)='NO MORE'
	    CALL NEW_GEN_IN(IONNAME(I),'Lines from ion (ALL or up to 10 lines)?:')
	    IF(IONNAME(I) .EQ. 'NO MORE')EXIT
	    N_IONNAME=I
	  END DO
!
	  TMP_STR='!'
	  DO WHILE(TMP_STR(1:1) .EQ. '!')
	    READ(33,'(A)')TMP_STR
	  END DO
	  BACKSPACE(33)
	  J=0
	  DO WHILE(J+1 .LE. 5000)
	    READ(33,*,END=2500)LINE_ID(J+1),ID_WAVE(J+1),TAU(J+1)
	    DO I=1,N_IONNAME
	      IF ((ID_WAVE(J+1)-LRANGE(1))*(LRANGE(2)-ID_WAVE(J+1)) .GT. 0.0
	1         .AND. TAU(J+1) .GE. TAU_CUT) THEN
	        IF ( TRIM(LINE_ID(J+1)) .EQ. TRIM(IONNAME(I)) .OR. TRIM(IONNAME(I)) .EQ. 'ALL') THEN
	          J=J+1
	        END IF
	     END IF 
	    END DO
	  END DO
2500	  CONTINUE
	  NLINES=J
!
500	  CALL NEW_GEN_IN(FILNAME,'Output file which can be later to read:')
	  INQUIRE(FILE=FILNAME,EXIST=FILE_PRES)
	  IF(FILE_PRES)THEN
	    WRITE(6,*)TRIM(FILNAME)//' file exists'
	    CALL NEW_GEN_IN(FILE_OPTION,'NEW name, OVERwrite, or APPEND')
	    IF(FILE_OPTION(1:3) .EQ. 'NEW')THEN
	      GOTO 500
	    ELSE IF(FILE_OPTION(1:4) .EQ. 'OVER')THEN
	      OPEN(UNIT=31,FILE=FILNAME,STATUS='OLD',ACTION='WRITE',IOSTAT=IOS)
	    ELSE IF(FILE_OPTION(1:4) .EQ. 'APPEND')THEN
	      OPEN(UNIT=31,FILE=FILNAME,STATUS='OLD',POSITION='APPEND',ACTION='WRITE',IOSTAT=IOS)
	    ELSE
	    END IF
	  ELSE
	    OPEN(UNIT=31,FILE=FILNAME,STATUS='NEW',ACTION='WRITE',IOSTAT=IOS)
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file',FILNAME
	    GOTO 1000
	  END IF
!
	  WRITE(31,'(A16,I5)')'Number of lines:',NLINES
	  WRITE(31,'(A1)')'!'
	  WRITE(31,'(A32)')'! Ion   Wo   Wf   Central Lambda'
	  WRITE(31,'(A1)')'!'
!
! Here the file is written in the proper format that is read below
!
	  REWIND(33)
	  TMP_STR='!'
	  DO WHILE(TMP_STR(1:1) .EQ. '!')
	    READ(33,'(A)')TMP_STR
	  END DO
	  BACKSPACE(33)
	  J=0
	  DO WHILE(J+1 .LE. 5000)
 	    READ(33,*,END=3500)LINE_ID(J+1),ID_WAVE(J+1),TAU(J+1),ID_WAVE_OFF(J+1),ID_Y_OFF(J+1)
 	    WRITE(17,*)LINE_ID(J+1),ID_WAVE(J+1),J+1
 	    DO I=1,N_IONNAME
	      IF ((ID_WAVE(J+1)-LRANGE(1))*(LRANGE(2)-ID_WAVE(J+1)) .GT. 0.0
	1          .AND. TAU(J+1) .GE. TAU_CUT) THEN
	        IF ( TRIM(LINE_ID(J+1)) .EQ. TRIM(IONNAME(I)) .OR. IONNAME(I) .EQ. 'ALL') THEN
 	           WRITE(16,*)ID_WAVE(J+1),TRIM(LINE_ID(J+1)),' ',TRIM(IONNAME(I))
	           T1=ID_WAVE(J+1)+ID_WAVE(J+1)*(DOP_VEL/3.0E5)
	           T2=ID_WAVE(J+1)-ID_WAVE(J+1)*(DOP_VEL/3.0E5)
	           WRITE(31,*)TRIM(LINE_ID(J+1)),T2,T1,ID_WAVE(J+1)
	           J=J+1
	          END IF
	        ENDIF
	      END DO
	    END DO
3500	    CONTINUE
	    CLOSE(UNIT=33)
	    REWIND(UNIT=31)
	ELSE
!
! This section open a file that already exists in the directory 
! with the information of the lines ranges and names.
!
	  FILNAME='ew_file1'
	  CALL NEW_GEN_IN(FILNAME,'FILE=')
	  OPEN(UNIT=31,FILE=FILNAME,STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error opening file'//TRIM(FILNAME)
	    GOTO 1000
	  END IF
	END IF
!
! Read list of lines with intervals or center and Vdop
!
	CNT=0
	STRING='!'
	i=0
	DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	  READ(31,'(A)',IOSTAT=IOS)STRING
	  WRITE(T_OUT,*)TRIM(STRING)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error reading file'
	    CLOSE(UNIT=31)
	    GOTO 1000
	  END IF
	  IF(INDEX(STRING,'Number of lines:') .NE. 0)THEN
	    I=INDEX(STRING,'Number of lines:')
	    READ(STRING(I+16:),*)CNT
	    STRING='!'
	  END IF
	  i=i+1
	END DO
	WRITE(T_OUT,*)'Number of lines',CNT
!	REWIND(UNIT=31)
!
! Allocate the necessary vectors to calculate the EWs
!
	ALLOCATE(INAMES(CNT),STAT=IOS)
	ALLOCATE(WLINEO(CNT),STAT=IOS)
	ALLOCATE(WLINEF(CNT),STAT=IOS)
	ALLOCATE(EWGR(CNT),STAT=IOS)
	ALLOCATE(CENT(CNT),STAT=IOS)
	T3=1.0D0
	EWGR(:)=0.0D0
	CENT(:)=0.0D0
	XT(1)=XPAR(1)
	XT(2)=XPAR(2)
	YT(1)=T3
	YT(2)=T3
!
! Plot a continuum line cont = 1.0
!
	CALL PGLINE(2,XT,YT)
	BACKSPACE(UNIT=31)
!
! Read the lines data: ion name, initial wavelength and final wavelength.
!
	DO I=1,CNT
	  READ(31,*,IOSTAT=IOS)INAMES(I),WLINEO(I),WLINEF(I)
	ENDDO
!
! ew_group compute the EWs values and the centroids.
!
	CALL EW_GROUP(CNT,INAMES,WLINEO,WLINEF,EWGR,CENT)
!
! Makes the vertical lines for the ranges around each lines
!
	TYPE_CURVE(NPLTS)='V'
!	 
! Deallocate the used vectors
!
1000	CONTINUE
	IF(ALLOCATED(INAMES))DEALLOCATE(INAMES)
	IF(ALLOCATED(WLINEO))DEALLOCATE(WLINEO)
	IF(ALLOCATED(WLINEF))DEALLOCATE(WLINEF)
	IF(ALLOCATED(EWGR))DEALLOCATE(EWGR)
	IF(ALLOCATED(CENT))DEALLOCATE(CENT)
	CLOSE(UNIT=31)
!
	RETURN
	END
