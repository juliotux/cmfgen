!
! Simple program to tst GET_IBOUND which is to be used to read in
! incident intesnities on a plane parallel atmosphere.
!
	PROGRAM TST_GET_IB
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	REAL*8, ALLOCATABLE :: MU(:)
	REAL*8, ALLOCATABLE :: IBOUND(:)
	REAL*8 FREQ
	INTEGER NP,LS
	LOGICAL INIT
!
	NP=5
	CALL GEN_IN(NP,'Number of angles')
	ALLOCATE (MU(NP))
	ALLOCATE (IBOUND(NP))
!
	MU(1)=0.0D0; MU(NP)=1.0D0
	DO LS=2,NP-1
	  MU(LS)=(LS-1.0D0)/(NP-1.0D0)
	END DO
!
	FREQ=20.0D0
	INIT=.TRUE.
	DO WHILE(1 .EQ. 1)
	  FREQ=FREQ/2.0D0
	  CALL GEN_IN(FREQ,'Next frequency (<0 to exit)')
	  IF(FREQ .LE. 0)EXIT
	  CALL GET_IBOUND(IBOUND,MU,NP,FREQ,INIT)
	  INIT=.FALSE.
	  WRITE(20,*)' '
	  WRITE(20,'(A,ES14.6)')'Freq=',FREQ
	  WRITE(20,'(5ES14.6)')(IBOUND(LS),LS=1,NP)
	END DO
!
	STOP
	END
