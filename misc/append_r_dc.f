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
	PROGRAM APPEND_R_DC
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
	CHARACTER*132, ALLOCATABLE :: HEAD(:)
C
	REAL*8, ALLOCATABLE :: R_O(:)
	REAL*8, ALLOCATABLE :: DI_O(:)
	REAL*8, ALLOCATABLE :: ED_O(:)
	REAL*8, ALLOCATABLE :: DC_O(:,:)
	CHARACTER*132, ALLOCATABLE :: HEAD_O(:)
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
!
	INTEGER N_ADD,N_NEW
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
	  CALL GEN_IN(DC_FILE,'Name of main file with departure coeficents.')
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
	  CALL GEN_IN(DC_FILE,'Name of file with extend R grid.')
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
!
        ALLOCATE (DC_O(N_O,ND_O))
        ALLOCATE (R_O(ND_O))
        ALLOCATE (DI_O(ND_O))
        ALLOCATE (ED_O(ND_O))
        ALLOCATE (HEAD_O(ND_O))
!
        DO K=1,ND_O
          HEAD_O(K)=' '
          DO WHILE(HEAD_O(K) .EQ. ' ')
            READ(DCIN_O,'(A)')HEAD_O(K)
          END DO
          READ(HEAD_O(K),*)R_O(K),DI_O(K),ED_O(K)
          READ(DCIN_O,*)(DC_O(J,K),J=1,N_O)
        END DO
!
	N_ADD=12
	CALL GEN_IN(N_ADD,'Number of depths from old file')
!	DCOUT_FILE='DCOUT'
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
	1             RSTAR_N,LSTAR_N,N_O,ND_N+N_ADD
C
	N_NEW=MIN(N_O,N_N)
	DO K=1,N_ADD
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')TRIM(HEAD_O(K))
	  WRITE(DCOUT,'(1X,1P,5E15.5)')(DC_O(I,K),I=1,N_NEW)
	END DO
!
	DO K=1,ND_N
	  WRITE(DCOUT,'(A)')' '
	  WRITE(DCOUT,'(A)')TRIM(HEAD(K))
	  WRITE(DCOUT,'(1X,1P,5E15.5)')(DC_N(I,K),I=1,N_NEW)
	END DO
C
	STOP
	END
