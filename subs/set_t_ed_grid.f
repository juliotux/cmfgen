!
! Subroutine to provide simple T and ED estimates to allow
! the computation of profile sacing and profile limits.
!

	SUBROUTINE SET_T_ED_GRID(T,ED,ND,LUIN)
	IMPLICIT NONE
!
! Created 06-Jan-2014
!
	INTEGER ND
	INTEGER LUIN
	REAL*8 T(ND)
	REAL*8 ED(ND)
!
	REAL*8 DELTA_T
	REAL*8 DELTA_ED
	REAL*8 ATOM_MIN_VAL,ATOM_MAX_VAL
	REAL*8 T_MIN_VAL,T_MAX_VAL
!
	INTEGER N_ED
	INTEGER N_T
	INTEGER I,J,L
	INTEGER, PARAMETER :: LUER=6			!Unit for error messages
!
	OPEN(UNIT=LUIN,STATUS='OLD',ACTION='READ',FILE='GRID_PARAMS')
	  READ(LUIN,*)N_T,N_ED
	  IF(N_ED*N_T .NE. ND)THEN
	    WRITE(LUER,*)'Error --- invalid table size specified in GRID_PARAMS'
	    WRITE(LUER,*)'Inconsistent with ND'
	    WRITE(LUER,*)'If this is meant to be a CMFGEN model delete the file GRID_PARAMS'
	    STOP
	  END IF
	  READ(LUIN,*)T_MIN_VAL,T_MAX_VAL
	  READ(LUIN,*)ATOM_MIN_VAL,ATOM_MAX_VAL
	CLOSE(UNIT=LUIN)
!
	DELTA_T=0.0D0
	DELTA_ED=0.0D0
	IF(N_ED .NE. 1)DELTA_ED=LOG(ATOM_MAX_VAL/ATOM_MIN_VAL)/(N_ED-1)
	IF(N_T .NE. 1)DELTA_T=LOG(T_MAX_VAL/T_MIN_VAL)/(N_T-1)
!
! Since only need for profiles, assume once ionized.
!
	L=0
	DO J=1,N_ED
	  DO I=1,N_T
	    L=L+1
	    ED(L)=EXP(LOG(ATOM_MIN_VAL)+(J-1)*DELTA_ED)
	    T(L)=EXP(LOG(T_MIN_VAL)+(I-1)*DELTA_T)
	  END DO
	END DO
!
	RETURN
	END
