!
! Subroutine to read the energy levels and oscillator strengths from a
! data file. The Einstein A coefficients are also returned, and are
! computed from the oscillator strengths.
!
! Oscillator strengths (and A values) are returned only for transitions which
! meet certain criteria [GF_ACTION='SET_ZERO']. Alternatively, the transitions 
! which don't meet the selection criteria, can have their f values returned 
! as negative values [GF_ACTION='SET_NEG']. 
!
! Selection criteria.
!
! (1) A minimum number (MIN_NUM_TRANS) of downward transitions is selected
!     from each level, regardless of the gf values.
!
! (2) All transitions with a lower level less than LEV_CUT are retained,
!     regardless of their f-value.
!
! (3) Transitions with gf < GF_CUT, and which don't satisfy the other
!     selection criteria, are not selected (or set negative).
!
! In addition to the above criteria, we can omit lines for which at
! least one level has not been observed.
!
	SUBROUTINE GENOSC_V8(EINA,FEDGE,STAT_WT,LEVNAME,
	1              ARAD,GAM2,GAM4,OBSERVED_LEVEL,
	1              IONIZATION_EN,ZION,
	1              OSCDATE,N,NTRET,
	1              GF_ACTION,GF_CUT,LEV_CUT,MIN_NUM_TRANS,ONLY_OBS_LINES,
	1              LUIN,LUOUT,FILNAME)
	IMPLICIT NONE
!
! Altered 19-May-2015 - Added '17-Jun-2014' as a valid format date (9-Jun-2015). 
! Altered 23-Dec-2004 - Check on ENERGY level ordering, and name matching when
!                          reading in transitions included.
! Altered 21-Nov-2000 - Changed to V8
!                       Option ONLY_OBS_LINES inserted.
! Altered 19-Oct-2000 - Changed to V7
!                       ARAD,GAM2,GAM4, and OBSERVED level installed in call.
!                       FORMAT date is now search for in data file.
! Created 02-Feb-1998 - GF_ACTION and MIN_NUM_TRANS installed. Renamed from V5.
!
	INTEGER N			!Number of levels to be returned.
	INTEGER LUIN,LUOUT
!
! The following data values are returned.
!
	INTEGER NTRET			!Number of transitions returned
	REAL*8 EINA(N,N)		!f(i,j)/A(j,i) values (i < j)
	REAL*8 FEDGE(N)			!Level ionization energy (10^15 Hz)
	REAL*8 STAT_WT(N)		!Statistical weight of level
	REAL*8 ARAD(N)			!Inverse radiative lifetime of upper level.
	REAL*8 GAM2(N)			!Collisional broadening parameter
	REAL*8 GAM4(N)			!Collisional broadening parameter
	LOGICAL OBSERVED_LEVEL(N)	!If true, level energy is known.
!
	REAL*8 IONIZATION_EN		!Ionization energy of ion (cm^-1}
	REAL*8 ZION			!Charge on core (i.e. ion charge +1)
	CHARACTER*(*) OSCDATE		!Date oscillator file was written.
!
	CHARACTER*(*) LEVNAME(N)	!
	CHARACTER*(*) FILNAME		!
!
! The following allows us to ignore transitions with gf < GF_CUT and whose
! lower level is > LEV_CUT. A minimum of MIN_NUM_TRANS downward transitions is 
! kept for each level, independent of there gf value. GF_ACTION indicates 
! whether to set the f values of the weak transitions to zero, or negative.
!
	REAL*8    GF_CUT
	INTEGER LEV_CUT
	INTEGER MIN_NUM_TRANS
	CHARACTER*(*) GF_ACTION
	LOGICAL   ONLY_OBS_LINES
!
! External functions.
!
	EXTERNAL SPEED_OF_LIGHT,ICHRLEN,RD_FREE_VAL,ERROR_LU,WARNING_LU
	REAL*8 SPEED_OF_LIGHT,RD_FREE_VAL
	INTEGER ICHRLEN,ERROR_LU,WARNING_LU
!
! Local variables
!
	CHARACTER*40 LOCNAME(N)
	CHARACTER*40 LOW_NAME,UP_NAME
	CHARACTER*11 FORMAT_DATE
	REAL*8 T1,T2,T3,SPEED_LIGHT
	INTEGER I,J,K,NW,L1,L2,IOS,LEV_ID
	INTEGER MAXLEN,CUT_CNT
	INTEGER LUER,LUWARN
	CHARACTER(LEN=200) STRING
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
!
! Variables for deleting weak transitions.
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	REAL*8 DOWN(N)
	INTEGER INDX(N)
!
! Variables for free-format internal reads.
!
	INTEGER NEXT,STR_LEN
	CHARACTER*80 DESC
	DATA STR_LEN/80/
!
! Constants for opacity etc. These must be set in the calling program. 
! They are used to keep absolute consistency between the f and A values.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
	SPEED_LIGHT=SPEED_OF_LIGHT()			!cm/s^-1
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
!
! Initialize arrays.
!
	EINA(1:N,1:N)=0.0D0
	GAM2(1:N)=0.0D0
	GAM4(1:N)=0.0D0
	ARAD(1:N)=0.0D0
	OBSERVED_LEVEL(1:N)=.TRUE.
!
! CUT_CNT is used to determine the number of transitions cut from the
! calculations due to low oscillator values.
!
	CUT_CNT=0
!
	CALL GEN_ASCI_OPEN(LUIN,FILNAME,'OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening '//FILNAME//' in GENOSCIL'
	  STOP
	END IF
!
! We keep reading the file until we come upon the format date or
! date (old format). From then on the file has to have a fixed 
! format. All file header information is written to unit LUOUT. 
! This can be the MODEL file, and hence will contain descriptions
! of the atomic model data used.
!
	  FORMAT_DATE='OLD'
	  L1=0; L2=0
	  DO WHILE(L1 .EQ. 0 .AND. L2 .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error reading in Oscillator Date from '//
	1                  FILNAME
	      WRITE(LUER,*)'IOSTAT=',IOS
	      STOP
	    END IF
	    WRITE(LUOUT,'(A)')TRIM(STRING)
	    L1=INDEX(STRING,'!Format date')
	    IF(L1 .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      FORMAT_DATE(1:11)=STRING(1:11)
	      IF(FORMAT_DATE .NE. '17-Oct-2000' .AND. FORMAT_DATE .NE. '17-Jun-2014')THEN
	        WRITE(LUER,*)'Error in GENOSC_V8: reading ',TRIM(FILNAME)
	        WRITE(LUER,*)'Invalid format date: ',TRIM(FORMAT_DATE)
	        STOP
	      END IF
	      READ(LUIN,'(A)')STRING		!Now get data date
	      L2=INDEX(STRING,'!Date')
	      IF(L2 .EQ. 0)THEN
	        WRITE(LUER,*)'Error in GENOSC_V8: reading ',TRIM(FILNAME)
	        WRITE(LUER,*)'!Date record not found, or incorrectly located'
	        WRITE(LUER,*)TRIM(STRING)
	        STOP
	      END IF
	    ELSE
	      L2=INDEX(STRING,'!Date')
	    END IF
	  END DO
	  STRING=ADJUSTL(STRING)
	  OSCDATE(1:11)=STRING(1:11)
	  IF(TRIM(FORMAT_DATE) .EQ. 'OLD')THEN
	    WRITE(LUWARN,*)'Possible error -- FORMAT Date not found in oscilator file'
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Number of energy levels')
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='NW Read in GENOSCIL-'//FILNAME
	    NW=RD_FREE_VAL(STRING,IONE,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1     'Error reading in # of energy levels from '//FILNAME
	    STOP
	  END IF
	  IF(N .GT. NW)THEN
	    WRITE(LUER,*)'Error reading '//FILNAME//' in GENOSCIL '//
	1         ' - incorrect number energy levels'
	    STOP
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Ionization energy') 
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='Ionization energy read in GENOSCIL-'//FILNAME
	    IONIZATION_EN=RD_FREE_VAL(STRING,IONE,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1    'Error reading in Ionization Energy from '//FILNAME
	    STOP
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Screened nuclear charge') 
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='NW Read in GENOSCIL-'//FILNAME
	    ZION=RD_FREE_VAL(STRING,IONE,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1     'Error reading in Screened Charge from '//FILNAME
	    STOP
	  END IF
!
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Number of transitions')
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='NTRET Read in GENOSCIL-'//FILNAME
	    NTRET=RD_FREE_VAL(STRING,IONE,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1     'Error reading in # of transitions from '//FILNAME
	    STOP
	  END IF
!
! Skip blank record.
!
	  READ(LUIN,'(A)')STRING
	  IF(STRING .NE. ' ')THEN
	    WRITE(LUER,*)'Error reading blank(1) from '//FILNAME
	    STOP
	  END IF	    
!
! We first read the record into a character string so that we can do
! an unformatted read on the real variables. LEVNAME name (or transition)
! must be separated by at least 2 spaces from the real data. 
! Level names need not be the same length. Note that FEDGE is initially 
! the excitation energy in cm^-1.
!
	  MAXLEN=MIN( LEN(LEVNAME(1)),LEN(LOCNAME(1)) )
	  DO I=1,N
	    READ(LUIN,'(A)')STRING
	    STRING=ADJUSTL(STRING)
	    L1=INDEX(STRING,'  ')-1
	    IF( L1 .LE. 0 .OR. L1 .GT. MAXLEN )THEN
	      WRITE(LUER,*)
	1            'Error reading in Level Names from '//FILNAME
	      WRITE(LUER,*)'Invalid length of name'
	      STOP
	    END IF
	    LEVNAME(I)=STRING(1:L1)
	    LOCNAME(I)=LEVNAME(I)
!
	    IOS=0
	    IF(FORMAT_DATE .EQ. '17-Oct-2000' .OR. FORMAT_DATE .EQ.  '17-Jun-2014')THEN
	      READ(STRING(L1+1:),*,IOSTAT=IOS)STAT_WT(I),FEDGE(I),T1,T2,T3,LEV_ID,
	1                         ARAD(I),GAM2(I),GAM4(I)
	      IF(ABS(LEV_ID) .NE. I)THEN
	         WRITE(LUER,*)'Error reading levels in GENOSC_V8'
	         WRITE(LUER,*)'Invalid level ID'
	         WRITE(LUER,*)TRIM(STRING)
	         STOP
	      ELSE IF(IOS .NE. 0)THEN
	         WRITE(LUER,*)'Error reading levels in GENOSC_V8'
	         WRITE(LUER,*)'IOS=',IOS
	         WRITE(LUER,*)TRIM(STRING)
	         STOP
	      END IF
	      IF(LEV_ID .LT. 0)OBSERVED_LEVEL(I)=.FALSE.
	    ELSE
	      READ(STRING(L1+1:),*)STAT_WT(I),FEDGE(I)
	    END IF
	    FEDGE(I)=(IONIZATION_EN-FEDGE(I))*SPEED_LIGHT*1.0D-15
	  END DO
!
! Read in(i.e. skip) additional level names.
!
	  IF(N .LT. NW)THEN
	    DO I=N+1,NW
	      READ(LUIN,'(A)')STRING
	      L1=INDEX(STRING,'  ')-1
	      IF( L1 .LE. 0 .OR. L1 .GT. MAXLEN )THEN
	        WRITE(LUER,*)
	1             'Error reading in Level Names from '//FILNAME
	        STOP
	      END IF
	    END DO
	  END IF
!
! Check the energy levels are in order.
!
	  DO I=2,N
	   IF(FEDGE(I) .GT. FEDGE(I-1))THEN
	      WRITE(LUER,*)' '
	      WRITE(LUER,*)'Warning/error reading in Level Names from '//TRIM(FILNAME)
	      WRITE(LUER,*)'Energy levels are out of order: levels are:',I-1,I
	    END IF
	  END DO
!
! If GF_CUT is large, we assume we just want the ENERGY levels.
!
	IF(GF_CUT .GT. 1000.0D0 .OR. NTRET .EQ. 0)THEN
	  CLOSE(LUIN)
	  NTRET=0
	  RETURN
	END IF
!
! Reads in all transitions, but notes only those important for the
! current model atom. No error checking is performed unless N=NW.
! (This could be done but would need to introduce a local character
!  array because of the additional names).
!
!
! Skip all records until we come across the header.
!
	L1=0
	DO WHILE(L1 .EQ. 0 .OR. STRING(1:1) .EQ. '!')
	  READ(LUIN,'(A)')STRING
	  L1=INDEX(STRING,'Transition')
	END DO
	L1=INDEX(STRING,'Lam(A)')
	IF(L1 .EQ. 0)THEN
	  WRITE(LUER,*)
	1    'Error reading oscillator header in '//FILNAME
	  STOP
	END IF
!
! Skip one blank line
!
	READ(LUIN,'(A)')STRING
	IF(STRING .NE. ' ')THEN
	  WRITE(LUER,*)'Error reading blank(2) from '//FILNAME
	  STOP
	END IF
!
! Now read in the transition data. We use the level indices, rather than the
! names to assign the oscillator strengths.
!
	DO K=1,NTRET
	  READ(LUIN,'(A)')STRING
	  L1=INDEX(STRING,'-')			!In case 2 spaces at end of name
	  LOW_NAME=ADJUSTL(STRING(1:L1-1))
	  STRING(1:)=STRING(L1+1:)
	  L1=INDEX(STRING,'  ')			!skip over upper level
          UP_NAME=ADJUSTL(STRING(1:L1-1))
	  STRING(1:)=STRING(L1:)
!
	  READ(STRING,*)T1			!Oscillator value
	  DO L2=1,3				!Skip over f, A, and Lambda
	    STRING=ADJUSTL(STRING)
	    L1=INDEX(STRING,'  ')
	    STRING(1:)=STRING(L1:)
	  END DO
	  L1=INDEX(STRING,'-')
	  READ(STRING(1:L1-1),*)I
	  READ(STRING(L1+1:),*)J
	  IF(I .LE. N .AND. J .LE. N)THEN
	    EINA(I,J)=ABS(T1)
	    IF(LOW_NAME .NE. LEVNAME(I) .OR.
	1       UP_NAME .NE. LEVNAME(J) .OR. I .GE. J)THEN
	      WRITE(LUER,*)'Invalid transition format in '//FILNAME
	      WRITE(LUER,*)'Indices and level names don''t match'
	      WRITE(LUER,*)'I,J=',I,J
	      STOP
	    END IF
	  END IF
	END DO			!K (Reading loop)
!
	CLOSE(UNIT=LUIN)
!
! Compute the Einstein A coefficients.
!
	T1=OPLIN/EMLIN*TWOHCSQ
	DO I=1,N-1
	  DO J=I+1,N                
	    EINA(J,I)=T1*EINA(I,J)*STAT_WT(I)/STAT_WT(J)
	1     *( (FEDGE(I)-FEDGE(J))**2 )
	  END DO
	END DO
!
	IF(FORMAT_DATE .EQ. 'OLD')THEN
	  GAM2(1:N)=0.0D0
	  GAM4(1:N)=0.0D0
	  ARAD(1:N)=0.0D0
	  DO I=1,N-1
	    DO J=I+1,N
	      ARAD(J)=ARAD(J)+EINA(J,I)
	    END DO
	  END DO
	END IF
!	  
!
! Delete very weak transitions from the line list. We first ensure that all
! parameters are reasonable.
!
	IF(GF_CUT .GT. 0 .AND. MAX(MIN_NUM_TRANS,LEV_CUT) .LT. N)THEN
!
	  IF(GF_ACTION .NE. 'SET_NEG' .AND. GF_ACTION .NE. 'SET_ZERO')THEN
	    WRITE(LUER,'(1X,A)')'Error in GENOSC_V8'
	    WRITE(LUER,'(1X,A,A)')'Invalid GF_ACTION ',GF_ACTION
	    WRITE(LUER,'(1X,A,A)')'Currently reading: ',FILNAME
	    STOP
	  END IF
!
	  K=MAX(LEV_CUT,1)
	  K=MAX(K,MIN_NUM_TRANS)
	  DO J=K+1,N				!Upper index
	    DOWN(1:J-1)=EINA(J,1:J-1)		!Einstein A values for level J
	    I=J-1
	    CALL INDEXX(I,DOWN,INDX,L_FALSE)	!Sort into reverse numerical order
	    T1=DOWN(INDX(MIN_NUM_TRANS))	! 
	    DO I=1,J-1				!Lower index
	      IF( STAT_WT(I)*EINA(I,J) .LT. GF_CUT .AND.
	1              I .GT. LEV_CUT .AND. EINA(J,I) .LT. T1)THEN
	        IF(EINA(I,J) .NE. 0)CUT_CNT=CUT_CNT+1
	        IF(GF_ACTION .EQ. 'SET_NEG')THEN
	          EINA(I,J)=-EINA(I,J)
	        ELSE
	          EINA(I,J)=0.0D0
	          EINA(J,I)=0.0D0
	        END IF
	      END IF
	    END DO
	  END DO
!
	  IF(CUT_CNT .NE. 0)THEN
	    WRITE(LUWARN,'(1X,A,I5,A,A)')'***Warning**** --- ',CUT_CNT,
	1     ' weak transitions cut in GENOSC_V8 --- ',TRIM(FILNAME)
 	  END IF
	END IF
!
	CUT_CNT=0
	IF(ONLY_OBS_LINES)THEN
	  DO I=1,N
	    DO J=I+1,N
	      IF(.NOT. OBSERVED_LEVEL(I) .OR. .NOT. OBSERVED_LEVEL(J))THEN
	        IF(EINA(I,J) .NE. 0)THEN
	          EINA(I,J)=0
	          EINA(J,I)=0
	          CUT_CNT=CUT_CNT+1
	        END IF
	      END IF
	    END DO
	  END DO
	END IF
	IF(CUT_CNT .NE. 0)THEN
	    WRITE(LUWARN,'(1X,A,I5,A,A)')'***Warning**** --- ',CUT_CNT,
	1     ' unobserved transitions cut in GENOSC_V8 --- ',TRIM(FILNAME)
 	END IF
!
!
	NTRET=0.0D0			!Number of principal transitions
	DO J=2,N
	  DO I=1,J-1
	    IF(EINA(I,J) .GT. 0)NTRET=NTRET+1
	  END DO
	END DO
	WRITE(LUOUT,'('' Final number of transitions found is: '',I6)')NTRET
!
	RETURN
	END
