C
C Program to rewrite a new Departure Coefficient file to accommodate an
C increase in the number of levels which is due to an increase in the size
C of the model atom.
C
C Required: Departure coefficient file with departure coefficients of
C              lower levels. This file will normally be used for input to
C              the new model.
C Required: Departure coefficient file with departure coeffiecients of
C              upper levels. This file determines the number of levels
C              output. This file is from some othe model in which an
C              extended atom was used.
C
	PROGRAM APPEND_DC
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
C
C Altered: 22-Jul-1997 - Format date now handled.
C Created: 19-May-1997
C
C Storage for oscillators etc from new oscilator file.
C
	INTEGER, PARAMETER :: T_IN=5		!Terminal input
	INTEGER, PARAMETER :: T_OUT=6		!Terminal output
	INTEGER, PARAMETER :: LUIN=10		!File input
	INTEGER, PARAMETER :: LUSCR=11	!
	INTEGER, PARAMETER :: DCIN_O=15
	INTEGER, PARAMETER :: DCIN_N=16
	INTEGER, PARAMETER :: DCOUT=17
	INTEGER, PARAMETER :: IONE=1
C
	REAL*8, ALLOCATABLE :: R_N(:)
	REAL*8, ALLOCATABLE :: DI_N(:)
	REAL*8, ALLOCATABLE :: ED_N(:)
	REAL*8, ALLOCATABLE :: DC_N(:,:)
	REAL*8, ALLOCATABLE :: TA(:)
C
	REAL*8, ALLOCATABLE :: R_O(:)
	REAL*8, ALLOCATABLE :: DI_O(:)
	REAL*8, ALLOCATABLE :: ED_O(:)
	REAL*8, ALLOCATABLE :: DC_O(:,:)
	CHARACTER*132, ALLOCATABLE :: HEAD(:)
C
	INTEGER N_N
	INTEGER ND_N
	REAL*8 LSTAR_N
	REAL*8 RSTAR_N
C
	INTEGER N_O
	INTEGER ND_O
	REAL*8 LSTAR_O
	REAL*8 RSTAR_O
C
	CHARACTER*80 DC_FILE
C
	REAL*8 T1
	CHARACTER*132 STRING_O,STRING_N
	CHARACTER*80 DCOUT_FILE
C
	INTEGER I,J,K,KJ
	INTEGER IZERO,IOS
C
C Constants for opacity etc.
C
	IOS=1
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(DC_FILE,'Name of file with depart. coef.')
	  CALL GEN_ASCI_OPEN(DCIN_N,DC_FILE,'OLD',' ','READ',IZERO,IOS)
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
	  READ(DCIN_N,'(A)')STRING_N
	END DO
	IF( INDEX(STRING_N,'!Format date') .EQ. 0)THEN
	   REWIND(DCIN_N)
	   STRING_N=' '
	END IF
	READ(DCIN_N,*)RSTAR_N,LSTAR_N,N_N,ND_N
C
	ALLOCATE (DC_N(N_N,ND_N))
	ALLOCATE (R_N(ND_N))
	ALLOCATE (DI_N(ND_N))
	ALLOCATE (ED_N(ND_N))
	ALLOCATE (HEAD(ND_N))
C
	DO K=1,ND_N
	  HEAD(K)=' '
	  DO WHILE(HEAD(K) .EQ. ' ')
	    READ(DCIN_N,'(A)')HEAD(K)
	  END DO
	  READ(HEAD(K),*)R_N(K),DI_N(K),ED_N(K)
	  READ(DCIN_N,*)(DC_N(J,K),J=1,N_N)
	END DO
C
	IOS=1
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(DC_FILE,'Name of file with old depart. coef.')
	  CALL GEN_ASCI_OPEN(DCIN_O,DC_FILE,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)WRITE(T_OUT,*)
	1      ' Error opening old DC file: Try again'
	END DO
C                                                  
C Check whether the file has a record containing 'Format date'. Its presence
C effects the way we read the file. 
C
	I=0
	STRING_O=' '
	DO WHILE(INDEX(STRING_O,'!Format date') .EQ. 0 .AND. I .LE. 10)
	  I=I+1
	  READ(DCIN_O,'(A)')STRING_O
	END DO
	IF( INDEX(STRING_O,'!Format date') .EQ. 0)REWIND(DCIN_O)
	READ(DCIN_O,*)RSTAR_O,LSTAR_O,N_O,ND_O
	CLOSE(DCIN_O)
C
	IF(N_O .LE. N_N)THEN
	  WRITE(6,*)'Error --  N_O must be larger than N_N'
	  STOP
	END IF
	IF( (STRING_O .EQ. ' ' .AND. STRING_N .NE. ' ')
	1      .OR. (STRING_O .NE. ' ' .AND. STRING_N .EQ. ' ') )THEN
	  WRITE(6,*)' '
	  WRITE(6,'(80A)')('*',I=1,80)
	  WRITE(6,*)'Warning -- program presently assumes both',
	1           ' data sets have the same clumping factor'
	  WRITE(6,'(80A)')('*',I=1,80)
	  WRITE(6,*)' '
	END IF
        
C
	ALLOCATE (DC_O(N_O,ND_N))
	ALLOCATE (TA(ND_N))
C
	CALL REGRID_B_ON_NE(DC_O,ED_N,DI_N,TA,
	1                N_O,IONE,N_O,ND_N,DCIN_O,DC_FILE)
C
	DCOUT_FILE='DCOUT'
	CALL GEN_IN(DCOUT_FILE,'Output file for departure coefficient')
	CALL GEN_ASCI_OPEN(DCOUT,DCOUT_FILE,'UNKNOWN',' ',
	1                                    'WRITE',IZERO,IOS)
C
C Output format string date if present.
C
	IF(STRING_N .NE. ' ')THEN
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')STRING_N
	END IF
	WRITE(DCOUT,'(A)')' '
	WRITE(DCOUT,'(F9.4,3X,1P,E12.4,5X,0P,I4,5X,I4)')
	1             RSTAR_N,LSTAR_N,N_O,ND_N
C
	DO K=1,ND_N
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')TRIM(HEAD(K))
	  WRITE(DCOUT,'(1X,1P,5E15.5)')(DC_N(I,K),I=1,N_N),
	1                           (DC_O(I,K),I=N_N+1,N_O)
	END DO
C
	STOP
	END
