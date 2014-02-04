C
C Program to rewrite a new Departure Coefficient file to accomodate an
C increase in the number of levels due to level splitting.
C
C Required: OLD Departure coefficient file.
C Required: OLD oscilator file
C
C Required: NEW oscilator file
C
C The number of levels for the OLD case is  determined by the number of 
C levels in the DC file.
C
C The number of levels for the NEW case is  determined by the last level
C corresponding to the HIGHEST level in the OLD data set.
C
C The name of any new levels can be individually reset if there is a 
C known mismatch. For example, 4z2Z renamed to 4f2Fo 
C (provided not last level).
C
	PROGRAM REWRITE_DC
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
C
C Altered: 01-Jun-1998 Big fix --- STRING length changed from 80 to 132
C                         so that clump factor output.
C Altered; 15-Mar-1997 !Cleaned; GEN_IN installed etc
C                        T_OUT installed ect.
C Altered: 22-OCT-1996 !Cleaned.
C Created: 16-Sep-1996
C
C Storage for oscilators etc from old oscilator file.
C A_O is needed for the call, but is not used in the main routine.
C
	INTEGER N_O
	REAL*8, ALLOCATABLE :: A_O(:,:)
	REAL*8, ALLOCATABLE :: EDGE_O(:)
	REAL*8, ALLOCATABLE :: G_O(:)
	CHARACTER*30, ALLOCATABLE ::LEVNAME_O(:)
	CHARACTER*132 OLD_OSC_FILE
C
C Storage for oscilators etc from new oscilator file.
C
	INTEGER N_N
	INTEGER, ALLOCATABLE :: INDX(:)
	REAL*8, ALLOCATABLE :: A_N(:,:)
	REAL*8, ALLOCATABLE :: EDGE_N(:)
	REAL*8, ALLOCATABLE :: G_N(:)
	CHARACTER*30, ALLOCATABLE ::LEVNAME_N(:)
	CHARACTER*132 NEW_OSC_FILE
C
	INTEGER, PARAMETER :: T_IN=5		!Terminal input
	INTEGER, PARAMETER :: T_OUT=6		!Terminal output
	INTEGER, PARAMETER :: LUIN=10		!File input
	INTEGER, PARAMETER :: LUSCR=11	!
	INTEGER, PARAMETER :: DCIN=15
	INTEGER, PARAMETER :: DCOUT=16
C
	REAL*8, ALLOCATABLE :: DC(:)
C
	INTEGER ND
	REAL*8 LSTAR
	REAL*8 RSTAR
	CHARACTER*132 DC_FILE
	CHARACTER*30 TMP_NAME
C
	REAL*8 GF_LEV_CUT
	REAL*8 EN_LEV_CUT
	REAL*8 Z
	REAL*8 T1
	CHARACTER*30 OSCDATE
	CHARACTER*132 STRING,STRING_N
	CHARACTER*132 LNK_FILE
	CHARACTER*132 DCOUT_FILE
	LOGICAL OLD_OSC_EXIST,CREATE_LNK_FILE
C
	INTEGER I,J,K,KJ
	INTEGER IZERO,IOS
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
C
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08		!pi*e*e/m/c*1.0E+10
	EMLIN=5.27296E-03		!pc*1.0E+025/4.0/pi
	DC_FILE=' '
C
	IOS=1
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(DC_FILE,'Name of file with old depart. coef.')
	  CALL GEN_ASCI_OPEN(DCIN,DC_FILE,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)WRITE(T_OUT,*)' Error opening DC file: Try again'
	END DO                                                    
C
C Check whether the file has a record containing 'Format date'. Its presence
C effects the way we read the file. If it has, we save it to output. Note that
C the header to each depth (i.e. that contianing R, Ne, etc) output to the 
C final file has EXACTLY the same format as the main input file.
C
	I=0
	STRING_N=' '
	DO WHILE(INDEX(STRING_N,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(DCIN,'(A)')STRING_N
	END DO
	IF( INDEX(STRING_N,'!Format date') .EQ. 0)THEN
	   REWIND(DCIN)
	   STRING_N=' '
	END IF
	READ(DCIN,*)RSTAR,LSTAR,N_O,ND
C
	ALLOCATE (A_O(N_O,N_O))
	ALLOCATE (EDGE_O(N_O))
	ALLOCATE (G_O(N_O))
	ALLOCATE (LEVNAME_O(N_O))
C
	OLD_OSC_EXIST=.FALSE.
	OLD_OSC_FILE=' '
	DO WHILE(.NOT. OLD_OSC_EXIST)
	  CALL GEN_IN(OLD_OSC_FILE,
	1       'Oscillator file assoc. with old D.C. file')
	  INQUIRE(FILE=OLD_OSC_FILE,EXIST=OLD_OSC_EXIST)
	  IF(.NOT. OLD_OSC_EXIST)
	1       WRITE(T_OUT,*)'Error opening OLD_OSC_FILE'
	END DO
	CALL GENOSC_V5(A_O,EDGE_O,G_O,
	1                   LEVNAME_O,T1,Z,
	1                   OSCDATE,N_O,I,EN_LEV_CUT,GF_LEV_CUT,
	1                   LUIN,LUSCR,OLD_OSC_FILE)
C
	DO I=1,N_O
	  DO WHILE(LEVNAME_O(I)(1:1) .EQ. ' ')
	    LEVNAME_O(I)(1:)=LEVNAME_O(I)(2:)
	  END DO
	END DO
C                                        
200	WRITE(T_OUT,'(A)',ADVANCE='NO')
	IOS=1
	NEW_OSC_FILE=' '
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(NEW_OSC_FILE,
	1    'Oscillator file to be assoc. with NEW D.C. file')
C
C Do an initial open to get the total number of available energy levels.
C
	  CALL GEN_ASCI_OPEN(LUIN,NEW_OSC_FILE,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)
	1    WRITE(T_OUT,*)'Error opening NEW OSC. file: Try again'
	END DO
C
	J=0
	DO WHILE(J .EQ. 0)
	  READ(LUIN,'(A)')STRING
	  J=INDEX(STRING,'!Number of energy levels')
	END DO
	READ(STRING,*)N_N
	CLOSE(UNIT=LUIN)
C
	ALLOCATE (DC(N_N))
	ALLOCATE (INDX(N_N))
	ALLOCATE (A_N(N_N,N_N))
	ALLOCATE (EDGE_N(N_N))
	ALLOCATE (G_N(N_N))
	ALLOCATE (LEVNAME_N(N_N))
C
C Now re-open to use standard routine to get levelnames etc.
C
	CALL GENOSC_V5(A_N,EDGE_N,G_N,
	1                   LEVNAME_N,T1,Z,
	1                   OSCDATE,N_N,I,EN_LEV_CUT,GF_LEV_CUT,
	1                   LUIN,LUSCR,NEW_OSC_FILE)
C
	DO I=1,N_N
	  DO WHILE(LEVNAME_N(I)(1:1) .EQ. ' ')
	    LEVNAME_N(I)(1:)=LEVNAME_N(I)(2:)
	  END DO
	END DO
C
C Find level in NEW data set corresponding to last level in OLD DATA set.
C It is assumed that the level ordering is the same.
C
	K=INDEX(LEVNAME_O(N_O),'[')
	IF(K .EQ. 0)THEN
	  J=N_N
	  KJ=INDEX(LEVNAME_N(N_N),'[')-1
	  IF(KJ .LE. 0)KJ=LEN(TRIM(LEVNAME_N(N_N)))
	  DO WHILE(LEVNAME_N(J)(1:KJ) .NE. LEVNAME_O(N_O))
	    J=J-1
	    IF(J .LT. 1)THEN
	      WRITE(T_OUT,*)' No final states match'
	      WRITE(T_OUT,*)' This may be due to a level-name change'
	      WRITE(T_OUT,*)' Use link routine'
	      STOP
	    END IF
	    KJ=INDEX(LEVNAME_N(J),'[')-1
	    IF(KJ .LT. 0)KJ=LEN(TRIM(LEVNAME_N(J)))
	  END DO
	ELSE
	  J=N_N
	  DO WHILE(LEVNAME_N(J) .NE. LEVNAME_O(N_O))
	    J=J-1           
	    IF(J .LT. 1)THEN
	      WRITE(T_OUT,*)' No final state match'
	      WRITE(T_OUT,*)' This may be due to a level-name change'
	      WRITE(T_OUT,*)' Use link routine'
	      STOP
	    END IF
	  END DO
	END IF
C
C NB: Setting N_N=J will make A_N inaccessible, but we don't need it anyway.
C
	N_N=J
C	
	INDX(1:N_N)=0.0D0
C
	DO I=1,N_N
C
C Identify those levels whose level names are IDENTICAL in both data sets.
C
500	  CONTINUE
	  DO J=1,N_O
	    IF(LEVNAME_N(I) .EQ. LEVNAME_O(J))THEN
	      IF(INDX(I) .NE. 0)THEN
	        WRITE(T_OUT,*)' Error: More than 1 level match'
	        WRITE(T_OUT,*)' New level name is ',LEVNAME_N(I)
	        WRITE(T_OUT,*)' Current old level name is ',LEVNAME_O(J)
	        WRITE(T_OUT,*)' Previous match was ',LEVNAME_O(INDX(I))
	        STOP
	      ELSE 
	        INDX(I)=J
	        IF(G_N(I) .NE. G_O(J))THEN
	          WRITE(T_OUT,*)' Error: inconsistent statistical weights'
	          WRITE(T_OUT,*)' G new=',G_N(I)
	          WRITE(T_OUT,*)' G old=',G_O(J)
	          WRITE(T_OUT,*)' New level name is ',LEVNAME_N(I)
	          WRITE(T_OUT,*)' Current old level name is ',LEVNAME_O(J)
	          STOP
	        END IF
	      END IF
	    END IF
	  END DO
C
C Now identify those levels which must have been split into individual
C j states.
C
	  IF(INDX(I) .EQ. 0)THEN
	    K=INDEX(LEVNAME_N(I),'[')
	    IF(K .EQ. 0)THEN
	      WRITE(T_OUT,*)' Error no match for state'
	      WRITE(T_OUT,*)' New level name is ',LEVNAME_N(I)
	      STOP
	    END IF
	    K=K-1		!don't want [ included.
	    DO J=1,N_O
	      IF(LEVNAME_N(I)(1:K) .EQ. LEVNAME_O(J))THEN
	        IF(INDX(I) .NE. 0)THEN   
	          WRITE(T_OUT,*)' Error: More than 1 level match'
	          WRITE(T_OUT,*)' New level name is ',LEVNAME_N(I)
	          WRITE(T_OUT,*)' Current old level name is ',LEVNAME_N(J)
	          WRITE(T_OUT,*)' Previous match was ',LEVNAME_N(INDX(I))
	          STOP
	        ELSE 
	          INDX(I)=J
	        END IF
	      END IF
	    END DO
	    IF(INDX(I) .EQ. 0)THEN
	      WRITE(T_OUT,*)' Error: No match for level ',LEVNAME_N(I)
C
C Allows for the possible renaming of a level.
C
	      TMP_NAME=' '
	      CALL GEN_IN(TMP_NAME,'Input matching old level name')
	      J=INDEX(LEVNAME_N(I),'[')
	      LEVNAME_N(I)=TRIM(TMP_NAME)//LEVNAME_N(I)(J:)
	      GOTO 500				!Begin search again.
	    END IF
	  END IF
	END DO
C
	DCOUT_FILE='DCOUT'
	CALL GEN_IN(DCOUT_FILE,'Output file for departure coefficient')
	CALL GEN_ASCI_OPEN(DCOUT,DCOUT_FILE,'UNKNOWN',' ',
	1                                    'WRITE',IZERO,IOS)

C
C Check whether the file has a record containing 'Format date'. Its presence
C effects the way we read the file. If it has, we save it to output. Note that
C the header to each depth (i.e. that contianing R, Ne, etc) output to the 
C final file has EXACTLY the same format as the main input file.
C
	IF(STRING_N .NE. ' ')THEN
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')TRIM(STRING_N)
	END IF
	WRITE(DCOUT,'(A)')' '
	WRITE(DCOUT,'(F9.4,3X,1P,E12.4,5X,0P,I4,5X,I4)')RSTAR,LSTAR,N_N,ND
C
C Can now read in DC file an rewrite out.
C
	DO K=1,ND
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ')
	    READ(DCIN,'(A)')STRING
	  END DO
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')STRING(1:LEN(STRING))
	  READ(DCIN,*)(DC(J),J=1,N_O)
	  WRITE(DCOUT,'(1X,1P,5E15.5)')(DC(INDX(I)),I=1,N_N)
	END DO
C
C Create a one-to-one link file which may be simpler to use.
C
	CREATE_LNK_FILE=.TRUE.
	CALL GEN_IN(CREATE_LNK_FILE,'Create a link file')
	IF(CREATE_LNK_FILE)THEN
	  LNK_FILE='LINK'
	  CALL GEN_IN(LNK_FILE,'Name for link file')
	  CALL GEN_ASCI_OPEN(DCOUT,LNK_FILE,
	1              'UNKNOWN',' ','WRITE',IZERO,IOS)
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A,A)')'New oscillator file was ',
	1         TRIM(NEW_OSC_FILE)
	  WRITE(DCOUT,'(A,A)')'Old oscillator file was ',
	1         TRIM(OLD_OSC_FILE)
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(75A)')('*',I=1,75)
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(2X,A,6X,A,T25,A,T45,A)')'New','Old',
	1              'New name','Old name'
	  WRITE(DCOUT,'(A)')'!'		!Signifies beginning of data.
	  DO I=1,N_N
	    WRITE(DCOUT,'(1X,I4,5X,I4,T25,A,T45,A)')I,INDX(I),
	1       TRIM(LEVNAME_N(I)),TRIM(LEVNAME_O(INDX(I)))
	  END DO
	END IF
C
	STOP
	END
