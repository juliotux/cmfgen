!
! Subroutine designed to allow the type of iteration to be specified.
! Options should be placed in IT_SPECIFIER.
! Options are:
!
! LAM              -- Lambda iteration.
! FULL             -- Regular full iteration.
! FULL(FIXT)       -- Hold T fixed.
! FULL(FIXBA)      -- Hold BA matrix fixed.
! FULL(FIXT,FIXBA) -- Hole BA matrix and T fixed.
! NORM             -- Use CMFGEN decision
!
! The number after the itration specified indicates the number of times that
! option is done. The iteration cycle is looped, unless the NORM command
! is found at the top of the file. In this case, the remaining specifcations
! are ignored. Only 20 specifications are considered.
!
! When using the FIXBA option, users should make sure that the FIXT options is
! consistent.
!
	SUBROUTINE SPECIFY_IT_CYCLE_V2(MAIN_COUNTER,COMPUTE_BA,LAMBDA_ITERATION,FIXED_T)
	USE UPDATE_KEYWORD_INTERFACE
	IMPLICIT NONE
!
! Altered 10-Jul-2015: Changed to V2. Added MAIN_COUNTER to call.
!                        Added REVISE_R option (comments added 16-Aug-2015; cur_hmi)
! Altered 06-Apr-2014: Modified for allow NONE for NORM, and FIX_T and FIX_BA
!                        (common typo errors).
! Created 11-Mar-2014
!
	INTEGER MAIN_COUNTER
	LOGICAL COMPUTE_BA
	LOGICAL FIXED_T
	LOGICAL LAMBDA_ITERATION
!
	INTEGER, PARAMETER :: NSTR_MAX=20
	CHARACTER(LEN=80) STORE(NSTR_MAX)
	CHARACTER(LEN=80) TEMP_STR
	CHARACTER(LEN=20) IT_OPTION
!
	INTEGER CYCLE_LENGTH
	INTEGER ITS_TO_DO
	INTEGER ITS_DONE
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER IOS
	INTEGER I,K,L
	INTEGER LU_IN
	INTEGER CNT
	LOGICAL, SAVE :: FIXED_T_SAVED=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	LU_IN=7
	OPEN(UNIT=LU_IN,FILE='IT_SPECIFIER',STATUS='OLD',IOSTAT=IOS)
	IF(IOS .NE. 0)RETURN
!
	K=0
	DO WHILE(K .LT. NSTR_MAX)
	  READ(LU_IN,'(A)',END=100)STORE(K+1)
	  K=K+1
	END DO
100	CONTINUE
	CYCLE_LENGTH=K
	REWIND(LU_IN)
!
! Skip over comments -- place them at end of option list.
! Check if we want to force a Revsion of the R grid. The ADJUST_R_DEFAULTS
! file must exist.
!
	CNT=0
	DO WHILE(CNT .LT. CYCLE_LENGTH)
	  CNT=CNT+1
	  STORE(1)=ADJUSTL(STORE(1))
	  IF(STORE(1)(1:1) .EQ. ' ' .OR. STORE(1)(1:1) .EQ. '!')THEN
	    TEMP_STR=STORE(1)
	    STORE(1:CYCLE_LENGTH-1)=STORE(2:CYCLE_LENGTH)
	    STORE(CYCLE_LENGTH)=TEMP_STR
	  ELSE IF(STORE(1)(1:8) .EQ. 'REVISE_R')THEN
	    TEMP_STR=STORE(1)
	    STORE(1:CYCLE_LENGTH-1)=STORE(2:CYCLE_LENGTH)
	    STORE(CYCLE_LENGTH)=TEMP_STR
	    CALL GET_LU(K,'SPECIFY_IT_CYCLE_V2')
	    I=MAIN_COUNTER+1
	    CALL UPDATE_KEYWORD(I,'[STRT_ITS]','ADJUST_R_DEFAULTS',L_TRUE,L_FALSE,K)
	    I=1
	    CALL UPDATE_KEYWORD(I,'[N_ITS]','ADJUST_R_DEFAULTS',L_FALSE,L_TRUE,K)
	  ELSE
	    EXIT
	  END IF
	END DO
!
	STORE(1)=ADJUSTL(STORE(1))
	L=INDEX(STORE(1),' ')
	IT_OPTION=STORE(1)(1:L-1)
	CALL SET_CASE_UP(IT_OPTION,IZERO,IZERO)
!
! We allow for the use of FIX_T and FIX_BA instead of FIXT and FIXBA.
! We will also treat NORM and NONE as equivalent.
!
	DO WHILE(INDEX(IT_OPTION,'FIX_') .NE. 0)
	  L=INDEX(IT_OPTION,'FIX_')
	  IT_OPTION(L+3:)=IT_OPTION(L+4:)
	END DO
!
	IF(IT_OPTION(1:4) .EQ. 'NORM' .OR.  IT_OPTION(1:4) .EQ. 'NONE')THEN
	  CLOSE(LU_IN)
	  RETURN
	ELSE IF(IT_OPTION(1:3) .EQ. 'LAM')THEN
	  COMPUTE_BA=.TRUE.
	  LAMBDA_ITERATION=.TRUE.
	  FIXED_T=.TRUE.
	ELSE IF(IT_OPTION .EQ. 'FULL')THEN
	  LAMBDA_ITERATION=.FALSE.
	  COMPUTE_BA=.TRUE.
	  FIXED_T=.FALSE.
	ELSE IF(IT_OPTION .EQ. 'FULL(FIXT)')THEN
	  LAMBDA_ITERATION=.FALSE.
	  COMPUTE_BA=.TRUE.
	  FIXED_T=.TRUE.
	ELSE IF(IT_OPTION .EQ. 'FULL(FIXBA)')THEN
	  LAMBDA_ITERATION=.FALSE.
	  COMPUTE_BA=.FALSE.
	  FIXED_T=FIXED_T_SAVED
	ELSE IF(IT_OPTION .EQ. 'FULL(FIXT,FIXBA)' .OR. IT_OPTION .EQ. 'FULL(FIXBA,FIXT)')THEN
	  LAMBDA_ITERATION=.FALSE.
	  COMPUTE_BA=.FALSE.
	  FIXED_T=.TRUE.
	ELSE
	  WRITE(6,*)'Unrecognized option in IT_SPECIFICATION'
	  WRITE(6,*)'Option is ',TRIM(IT_OPTION),' and will be omitted from file'
	  DO I=2,CYCLE_LENGTH
	    WRITE(LU_IN,'(A)')TRIM(STORE(I))
	  END DO
	  CLOSE(LU_IN)
	  RETURN
	END IF
	FIXED_T_SAVED=FIXED_T
!
	TEMP_STR=ADJUSTL(STORE(1)(L:))
	WRITE(6,*)TEMP_STR
	WRITE(6,*)STORE(1)
	IF(TEMP_STR .EQ. ' ')THEN
	  ITS_TO_DO=1
	  ITS_DONE=0
	ELSE
	  READ(TEMP_STR,*,IOSTAT=IOS)ITS_TO_DO
	  IF(IOS .NE. 0 .OR. ITS_TO_DO .LE. 0)THEN
	    ITS_TO_DO=1
	    ITS_DONE=0
	  ELSE
	    L=INDEX(TEMP_STR,' ')
	    TEMP_STR=ADJUSTL(TEMP_STR(L:))
	    IF(L .EQ. 0)THEN
	      ITS_DONE=0
	    ELSE
	      READ(TEMP_STR,*,IOSTAT=IOS)ITS_DONE
	      IF(IOS .NE. 0)ITS_DONE=0
	    END IF
	  END IF
	END IF
!
	ITS_DONE=ITS_DONE+1
	WRITE(6,*)TRIM(IT_OPTION),ITS_TO_DO,ITS_DONE
	REWIND(LU_IN)
	IF(ITS_TO_DO .EQ. ITS_DONE)THEN
	  DO I=2,CYCLE_LENGTH
	    WRITE(LU_IN,'(A)')TRIM(STORE(I))
	  END DO
	  WRITE(LU_IN,'(A,T21,I4)')TRIM(IT_OPTION),ITS_TO_DO
	ELSE
	  WRITE(LU_IN,'(A,T21,I4,3X,I4)')TRIM(IT_OPTION),ITS_TO_DO,ITS_DONE
	  DO I=2,CYCLE_LENGTH
	    WRITE(LU_IN,'(A)')TRIM(STORE(I))
	  END DO
	END IF
!
	CLOSE(LU_IN)
	RETURN
	END
