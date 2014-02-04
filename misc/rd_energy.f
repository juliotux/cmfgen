C
C Subroutine to read in the levale names and energies for use by 
C WR_F_TO_S.
C
	SUBROUTINE RD_ENERGY(LEVNAME,STAT_WT,ENERGY,FEDGE,N,NMAX,
	1                   IONIZATION_EN,ZION,
	1                   OSCDATE,FILNAME,LUIN,LUOUT,IOS)
	IMPLICIT NONE
C
	INTEGER N,NMAX
	INTEGER LUIN,LUOUT,IOS
	REAL*8 STAT_WT(NMAX)
	REAL*8 ENERGY(NMAX)
	REAL*8 FEDGE(NMAX)
	REAL*8 IONIZATION_EN,ZION
	CHARACTER*(*) LEVNAME(NMAX)
	CHARACTER*(*) FILNAME,OSCDATE
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
	INTEGER, PARAMETER :: IZERO=0
	REAL*8 SPEED_LIGHT
	INTEGER I,J,L1,BIGLEN,MAXLEN,LUER
	CHARACTER*132 STRING
	CHARACTER*40 LOCNAME(NMAX)
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
	LEVNAME(:)=' '
	FEDGE(:)=0
	ENERGY(:)=0
	STAT_WT(:)=0
C
	CALL GEN_ASCI_OPEN(LUIN,FILNAME,'OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening '//FILNAME//' in GENOSCIL'
	  IOS=2
	  RETURN
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
	    RETURN
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
	    DESC='N Read in GENOSCIL-'//FILNAME
	    N=RD_FREE_VAL(STRING,1,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1     'Error reading in # of energy levels from '//FILNAME
	    IOS=3
	    RETURN
	  END IF
	  IF(N .GT. NMAX)THEN
	    WRITE(LUER,*)'Error reading '//FILNAME//' in GENOSCIL '//
	1         ' - insufficient storage'
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
	    IONIZATION_EN=RD_FREE_VAL(STRING,1,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1    'Error reading in Ionization Energy from '//FILNAME
	    IOS=4
	    RETURN
	  END IF
C
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Screened nuclear charge') 
	  IF(L1 .NE. 0)THEN
	    L1=INDEX(STRING,'  ')
	    DESC='NW Read in GENOSCIL-'//FILNAME
	    ZION=RD_FREE_VAL(STRING,1,L1-1,NEXT,DESC)
	  ELSE
	    WRITE(LUER,*)
	1     'Error reading in Screened Charge from '//FILNAME
	    IOS=5
	    RETURN
	  END IF
C
	  READ(LUIN,'(A)')STRING
	  J=ICHRLEN(STRING)
	  WRITE(LUOUT,'(A)')STRING(1:J)
	  L1=INDEX(STRING,'!Number of transitions')
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)
	1     'Error reading in # of energy levels from '//FILNAME
	    IOS=6
	    RETURN
	  END IF
C
C Skip blank record.
C
	  READ(LUIN,'(A)')STRING
	  IF(STRING .NE. ' ')THEN
	    WRITE(LUER,*)'Error reading blank(1) from '//FILNAME
	    IOS=7
	    RETURN
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
	    ENERGY(I)=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,DESC)
	    FEDGE(I)=(IONIZATION_EN-ENERGY(I))*SPEED_LIGHT*1.0D-15
	    BIGLEN=MAX(BIGLEN,L1)
	  END DO
C
	RETURN
	END
