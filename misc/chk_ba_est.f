!
! Small program to read in two BA_ASI_N_D? files from successive iterations. 
! Ideally the populations for the iterations should be identical, except for one
! value. This code allows a comparision between the actual changes, and those predicted
! from solution of the linearized equations. Useful for testing.
!
	PROGRAM CHK_BA_EST
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Created: 18-Nov-2009
!
	REAL*8, ALLOCATABLE :: CMAT_A(:,:)
	REAL*8, ALLOCATABLE :: POPS_A(:)
	REAL*8, ALLOCATABLE :: STEQ_A(:)
!
	REAL*8, ALLOCATABLE :: CMAT_B(:,:)
	REAL*8, ALLOCATABLE :: POPS_B(:)
	REAL*8, ALLOCATABLE :: STEQ_B(:)
!
	REAL*8, ALLOCATABLE :: dSTEQ(:)
	REAL*8 T1
!
	INTEGER NT
	INTEGER I
	INTEGER IVAR
	INTEGER IOS
	LOGICAL FILE_OPEN
	CHARACTER(LEN=80) FILENAME
	CHARACTER(LEN=132) STRING
!
	INTEGER, PARAMETER :: LU=20
!
	OPEN(UNIT=LU,FILE='MODEL',STATUS='OLD',IOSTAT=IOS)
	  IF(IOS .EQ. 0)THEN
	    DO WHILE(1 .EQ. 1)
	      READ(LU,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)EXIT
	      IF(INDEX(STRING,'!Total number of variables') .NE. 0)THEN
	        READ(STRING,*)NT
	        WRITE(6,'(A,I4)')' Number of variables in the model is:',NT
	        EXIT
	      END IF
	    END DO
	  END IF
	INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	IF(FILE_OPEN)CLOSE(UNIT=LU)
!
	 IF(IOS .NE. 0)THEN
	  WRITE(6,*)' Unable to open MODEL file to get # of depth points'
	  CALL GEN_IN(NT,'Number of variables (NT)')
	END IF
!
	ALLOCATE (CMAT_A(NT,NT),STEQ_A(NT),POPS_A(NT))
	ALLOCATE (CMAT_B(NT,NT),STEQ_B(NT),POPS_B(NT))
	ALLOCATE (dSTEQ(NT))
!
	FILENAME='BA_ASCI_N_D84'
	CALL GEN_IN(FILENAME,'Current file with BA and STEQ data')
	CALL RD_BA_MAT(CMAT_A,STEQ_A,POPS_A,NT,FILENAME,LU)
!
	CALL GEN_IN(FILENAME,'Old file with BA and STEQ data')
	CALL RD_BA_MAT(CMAT_B,STEQ_B,POPS_B,NT,FILENAME,LU)
!
	WRITE(6,*)'The following variables were changed -- ideally only one should have changed'
	IVAR=0
	DO I=1,NT
	  IF(POPS_A(I) .NE. POPS_B(I))THEN
	     WRITE(6,*)I,POPS_A(I),POPS_B(I)
	     IVAR=I
	  END IF
	END DO
	IF(IVAR .EQ. 0)THEN
	  WRITE(6,*)'Error -- all level populations are identical'
	  STOP
	END IF
!
	DO I=1,NT
	  dSTEQ(I)=CMAT_B(I,IVAR)*(POPS_A(IVAR)-POPS_B(IVAR))/POPS_B(IVAR)
	END DO
!
	WRITE(45,'(/,A,5(5X,A))')'     I','   STEQ_A','   STEQ_B','STEQ(A-B)','PRED(A-B)','  % Diff.'
	DO I=1,NT
	  T1=STEQ_A(I)-STEQ_B(I)
	  IF(T1 .NE. 0)T1=100.0D0*(dSTEQ(I)/T1-1.0D0)
	  IF(ABS(T1) .LT. 100000)THEN
	     WRITE(45,'(2X,I4,4ES14.4,F14.4)')I,STEQ_A(I),STEQ_B(I),STEQ_A(I)-STEQ_B(I),
	1            dSTEQ(I),T1
	  ELSE
	     WRITE(45,'(2X,I4,5ES14.4)')I,STEQ_A(I),STEQ_B(I),STEQ_A(I)-STEQ_B(I),
	1            dSTEQ(I),T1
	  END IF
	END DO
!
	STOP
	END
