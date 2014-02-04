C
C Subroutine to compute the oscillator strengths and Einstein A coefficients
C for an arbitrary species. Two auxilary arrays are returned. One of 
C these indicates if the transition is allowed and contains a transition 
C number. The second array is used to indicate the transitions of secondary 
C importance.
C
	SUBROUTINE GENOSC_V5(EINA,FEDGE,STAT_WT,LEVNAME,
	1                   IONIZATION_EN,ZION,
	1                   OSCDATE,N,NTRET,
	1                   GF_CUT,LEV_CUT,
	1                   LUIN,LUOUT,FILNAME)
	IMPLICIT NONE
C
C Altered 05-Feb-1998 - ROutine allways cut out "negative" oscillator strengths.
C Altered 05-Dec-1996 - GEN_ASCI_OPEN installed.
C Altered 25-May-1996 - NLOC removed - LOCANME dimensioned using N (F90)
C                       CUT_CNT instituted.
C                       GEN_SEQ_OPEN removed.
C                       Addittional levelnames are NOT stored in LOCNAME.
C                       IONE inserted in RD_FREE_VAL calls.
C Altered 24-OCt-1995 - GION deleted fromc all an now readin with 
C                       photoionization cross section.
C                       Blanks suppresed from beginning of level names.
C
C Altered 05-Oct-1995 - GION inserted in call. GION represents the statistical
C                       weight of the next ionization stage and is used in 
C                       the absence of that ionization stage.
C Altered  1-Aug-1995 - GF_CUT and LEV_CUT installed.
C Altered 27-Jul-1995 - Bug fix. Oscilator strengths were not read in when
C                         the first name in the transition ID had 2 or more
C                         spaces and the end .
C Created 30-May-1995 - Based on GENOSCIL --- called V2 for this reason.
C                       Arrays SECND and TRANS were deleted from call.
C
	INTEGER N,NTRET,LUIN,LUOUT
	REAL*8 EINA(N,N),FEDGE(N),STAT_WT(N)
	REAL*8 IONIZATION_EN,ZION
	CHARACTER*(*) LEVNAME(N)
	CHARACTER*(*) FILNAME,OSCDATE
C
C The following allows us to ignore transitions with gf < GF_CUT and whose
C lower level is > LEV_CUT.
C
	REAL*8 GF_CUT
	INTEGER LEV_CUT
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
C External functions.
C
	EXTERNAL SPEED_OF_LIGHT,ICHRLEN,RD_FREE_VAL,ERROR_LU
	REAL*8 SPEED_OF_LIGHT,RD_FREE_VAL
	INTEGER ICHRLEN,ERROR_LU
C
C Local variables
C
	CHARACTER*30 LOCNAME(N)
	REAL*8 T1,SPEED_LIGHT
	INTEGER I,J,K,NW,L1,L2,IOS,BIGLEN,MAXLEN,LUER,CUT_CNT
	CHARACTER*132 STRING
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
C
C Variables for free-format internal reads.
C
	INTEGER NEXT,STR_LEN
	CHARACTER*80 DESC
	DATA STR_LEN/80/
C
C Constants for opacity etc.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
C
	SPEED_LIGHT=SPEED_OF_LIGHT()			!cm/s^-1
	LUER=ERROR_LU()
C
C Initialize arrays.
C
	DO I=1,N
	  DO J=1,N
	     EINA(J,I)=0.0D0
	  END DO
	END DO
C
C CUT_CNT is used to determine the number of transitions cut from the
C calculations due to low oscilator values.
C
	CUT_CNT=0
C
	CALL GEN_ASCI_OPEN(LUIN,FILNAME,'OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening '//FILNAME//' in GENOSCIL'
	  STOP
	END IF
C
C We keep reading the file until we come upon the date. From
C then on the file has to have a fixed format. All file header
C information is written to unit LUOUT. This can be the MODEL file,
C and hence will contain relevant model data.
C
50	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Date')
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading in Oscilator Date from '//
	1                  FILNAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
	  IF(L1 .NE. 0)THEN
	    OSCDATE(1:11)=STRING(1:11)
	  ELSE
	    GOTO 50
	  END IF
C
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
C
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Ionization energy') 
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='IOnization energy read in GENOSCIL-'//FILNAME
	    IONIZATION_EN=RD_FREE_VAL(STRING,IONE,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1    'Error reading in Ionization Energy from '//FILNAME
	    STOP
	  END IF
C
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
C
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
C
C Skip blank record.
C
	  READ(LUIN,'(A)')STRING
	  IF(STRING .NE. ' ')THEN
	    WRITE(LUER,*)'Error reading blank(1) from '//FILNAME
	    STOP
	  END IF	    
C
C We first read the record into a character string so that we can do
C an unformatted read on the real variables. LEVNAME name (or transition)
C must be separated by at least 2 spaces from the real data. 
C Level names need not be the same length. Note that FEDGE is initially 
C the excitation energy in cm^-1.
C
	  BIGLEN=0
	  MAXLEN=MIN( LEN(LEVNAME(1)),LEN(LOCNAME(1)) )
	  DO I=1,N
	    READ(LUIN,'(A)')STRING
	    L1=INDEX(STRING,'  ')-1
	    IF( L1 .LE. 0 .OR. L1 .GT. MAXLEN )THEN
	      WRITE(LUER,*)
	1            'Error reading in Level Names from '//FILNAME
	      WRITE(LUER,*)'Invalid length of name'
	      STOP
	    END IF
	    LEVNAME(I)=STRING(1:L1)
	    DO WHILE(LEVNAME(I)(1:1) .EQ. ' ')		!Strip leading blanks.
	      LEVNAME(I)(1:)=LEVNAME(I)(2:)
	    END DO
	    LOCNAME(I)=LEVNAME(I)
	    DESC='STAT_WT read in GENOSCIL-'//FILNAME
	    STAT_WT(I)=RD_FREE_VAL(STRING,L1+1,STR_LEN,NEXT,DESC)
	    DESC='FEDGE read in GENOSCIL-'//FILNAME
	    FEDGE(I)=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,DESC)
	    FEDGE(I)=(IONIZATION_EN-FEDGE(I))*SPEED_LIGHT*1.0D-15
	    BIGLEN=MAX(BIGLEN,L1)
	  END DO
C
C Read in additional level names.
C
	  IF(N .LT. NW)THEN
	    DO I=N+1,NW
	      READ(LUIN,'(A)')STRING
	      L1=INDEX(STRING,'  ')-1
	      IF( L1 .LE. 0 .OR. L1 .GT. MAXLEN )THEN
	        WRITE(LUER,*)
	1             'Error eading in Level Names from '//FILNAME
	        STOP
	      END IF
	    END DO
	  END IF
C
C Reads in all transitions, but notes only those important for the
C current model atom. No error checking is performed unless N=NW.
C (This could be done but would need to introduce a local character
C  array because of the additional names).
C
C
C Skip all records until we come across the header.
C
999	READ(LUIN,'(A)')STRING
	L1=INDEX(STRING,'Transition')
	IF(L1 .EQ. 0)GOTO 999
	L1=INDEX(STRING,'Lam(A)')
	IF(L1 .EQ. 0)THEN
	  WRITE(LUER,*)
	1    'Error reading oscilator header in '//FILNAME
	  STOP
	END IF
C
C Skip one blank line
C
	READ(LUIN,'(A)')STRING
	IF(STRING .NE. ' ')THEN
	  WRITE(LUER,*)'Error reading blank(2) from '//FILNAME
	  STOP
	END IF
C
	DO K=1,NTRET
	  READ(LUIN,'(A)')STRING
C
	  L1=INDEX(STRING,'-')			!In case 2 spaces at end of name
	  STRING(1:)=STRING(L1:)
	  L1=INDEX(STRING,'  ')
	  STRING(1:)=STRING(L1:)
	  READ(STRING,*)T1		!Oscilator value
	  DO L2=1,3
	    DO WHILE(STRING(1:1) .EQ. ' ')
	      STRING=STRING(2:)
	    END DO
	    L1=INDEX(STRING,'  ')
	    STRING(1:)=STRING(L1:)
	  END DO
	  L1=INDEX(STRING,'-')
	  READ(STRING(1:L1-1),*)I
	  READ(STRING(L1+1:),*)J
	  IF(I .LE. N .AND. J .LE. N)THEN
	    IF(STAT_WT(I)*ABS(T1) .LT. GF_CUT .AND. I .GT. LEV_CUT)THEN
	      T1=0.0D0
	      CUT_CNT=CUT_CNT+1
	    END IF
	    EINA(I,J)=T1
	  END IF
	END DO			!K (Reading loop)
	IF(CUT_CNT .NE. 0)THEN
	  WRITE(LUER,'(1X,A,I5,A)')'***Warning**** --- ',CUT_CNT,
	1       'transitions cut in GENOSC_V5 --- ',FILNAME
	END IF
C
	CLOSE(UNIT=LUIN)
C
C Check for transitions of secondary importance.
C
	NTRET=0.0D0			!Actual number of transitions
	DO I=2,N
	  DO J=1,I
	    IF(EINA(J,I) .LT. 0)EINA(J,I)=-EINA(J,I)
	    IF(EINA(J,I) .NE. 0)NTRET=NTRET+1
	  END DO
	END DO
C
C Compute the Einstein A coefficients.
C
	T1=OPLIN/EMLIN*TWOHCSQ
	DO I=1,N-1
	  DO J=I+1,N                
	    EINA(J,I)=T1*EINA(I,J)*STAT_WT(I)/STAT_WT(J)
	1     *( (FEDGE(I)-FEDGE(J))**2 )
	  END DO
	END DO
C
	WRITE(LUOUT,'('' Final number of transitions found is: '',I6)')
	1         NTRET
C
	RETURN
	END
