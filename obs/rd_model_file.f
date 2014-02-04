!
! Set or routines to read the MODEL.DAT file created by CMFGEN. Can be
! used to get variables not available in the POP files.
!
	MODULE MOD_MODEL_FILE
	IMPLICIT NONE
!
! Finalized 5-Jan-1999
!
	INTEGER N_REC
	INTEGER LU_ER
	INTEGER REC_SPEC_HD
	CHARACTER*132, ALLOCATABLE :: STR_STORE(:)
!
	END MODULE MOD_MODEL_FILE
! 
!
! Subroutines read in entire model file, and stores in STR_STORE for
! access by specic reading routines.
!
	SUBROUTINE RD_MODEL_FILE(FILE_IN,LU_IN,IOS)
	USE MOD_MODEL_FILE
	IMPLICIT NONE
!
! Passed 
!
	INTEGER IOS
	INTEGER LU_IN
	CHARACTER*(*) FILE_IN
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
!Local
!
	INTEGER I
!
	IOS=0
	LU_ER=ERROR_LU()
!
	OPEN(UNIT=LU_IN,FILE=FILE_IN,STATUS='OLD',ACTION='READ',
	1    IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LU_ER,*)'Error opening file in RD_MODEL_FILE'
	    WRITE(LU_ER,*)'IOS=',IOS
	    WRITE(LU_ER,*)'FILE=',FILE_IN
	    RETURN
	  END IF
!
	  N_REC=0
	  DO WHILE (IOS .EQ. 0)
	    READ(LU_IN,'(A)',IOSTAT=IOS)
	    N_REC=N_REC+1
	  END DO
	  N_REC=N_REC-1		!As last read unsuccessful
	  IOS=0
!
	  ALLOCATE (STR_STORE(N_REC))
	  REWIND(LU_IN)
	  DO I=1,N_REC
	    READ(LU_IN,'(A)')STR_STORE(I)
	    IF(INDEX(STR_STORE(I),'Species') .NE. 0 .AND.
	1      INDEX(STR_STORE(I),'N_F') .NE. 0 .AND.
	1      INDEX(STR_STORE(I),'N_S') .NE. 0)REC_SPEC_HD=I
	  END DO
	CLOSE(LU_IN)
!
	RETURN
	END
!
!
!
	SUBROUTINE CLEAN_MODEL_STORE
	USE MOD_MODEL_FILE
	IF(ALLOCATED(STR_STORE))DEALLOCATE (STR_STORE)
	RETURN
	END
!
! 
!
! Read in nuber of super levels for each species, and the dielectronic
! options.
!
	SUBROUTINE RD_MODEL_SPEC_INFO(SPECIES,NC2_F,C2_PRES,
	1                    NC2_S,DIE_AUTO_C2,DIE_WI_C2)
	USE MOD_MODEL_FILE
	IMPLICIT NONE
!
! Passed
!
	INTEGER NC2_F
	LOGICAL C2_PRES
	CHARACTER*(*) SPECIES
!
! Returned
!
	INTEGER NC2_S
	LOGICAL DIE_AUTO_C2
	LOGICAL DIE_WI_C2
!
! Local vaiables
!
	INTEGER I,K
	INTEGER LNG
	CHARACTER*30 UC
	EXTERNAL UC
!
	IF(.NOT. C2_PRES)THEN
	  NC2_S=1
	  DIE_AUTO_C2=.FALSE.
	  DIE_WI_C2=.FALSE.
	  RETURN
	END IF
!
! We only set NC2_S this way if it wasn't read in from the POP file.
!
	IF(NC2_S .EQ. 0)THEN
	  LNG=LEN_TRIM(SPECIES)
	  DO I=REC_SPEC_HD+2,N_REC
	    IF( STR_STORE(I)(2:LNG+1) .EQ. TRIM(SPECIES) .OR.
	1         STR_STORE(I)(2:LNG+1) .EQ. TRIM(UC(SPECIES)) )THEN
	      READ(STR_STORE(I)(LNG+2:),*)K,NC2_S
	      IF(K .NE. NC2_F)THEN
	        WRITE(LU_ER,*)'Error in RD_MODEL_SPEC_INFO: Inconsistenmt N_F'
	        WRITE(LU_ER,*)'Species=',Species
	        WRITE(LU_ER,*)'Passed N_F=',NC2_F
	        WRITE(LU_ER,*)'N_F in model file=',K
	        STOP
	      END IF
	      GOTO 1000
	    END IF
	  END DO
	  WRITE(LU_ER,*)'Unable to read N_S from MODEl file'
	  WRITE(LU_ER,*)'Error in RD_MODEL_SPEC_INFO: Species: ',TRIM(SPECIES)
	  STOP
	END IF
!
1000	CONTINUE
	LNG=LEN_TRIM(SPECIES)
	DO I=1,N_REC
	  K=INDEX(STR_STORE(I),'[')
	  IF(K .NE. 0)THEN
	    IF( STR_STORE(I)(K:K+LNG+5) .EQ. '[DIE_'//TRIM(SPECIES)//']' .OR.
	1       STR_STORE(I)(K:K+LNG+5) 
	1                    .EQ. '[DIE_'//TRIM(UC(SPECIES))//']' )THEN
	      READ(STR_STORE(I),*)DIE_AUTO_C2,DIE_WI_C2
	      GOTO 2000
	    END IF
	  END IF
	END DO
	WRITE(LU_ER,*)'Unable to read DIE options from MODEL file'
	WRITE(LU_ER,*)'Error in RD_MODEL_SPEC_INFO: Species: ',TRIM(SPECIES)
	STOP
!
2000	CONTINUE
	RETURN
	END
!
!
!
! The remaining routines alow individual options to be input.
!
	SUBROUTINE RD_MODEL_LOG(VAR,DESC)
	USE MOD_MODEL_FILE
	IMPLICIT NONE
!
	LOGICAL VAR
	CHARACTER*(*) DESC
	INTEGER I,K
!
	DO I=1,N_REC
	  K=INDEX(STR_STORE(I),'['//TRIM(DESC)//']')
	  IF(K .NE. 0)THEN
	    READ(STR_STORE(I),*)VAR
	    GOTO 2000
	  END IF
	END DO
	WRITE(LU_ER,*)'Unable to read option from MODEL file'
	WRITE(LU_ER,*)'Error in RD_MODEL_LOG: Key word: ',TRIM(DESC)
	STOP
!
2000	CONTINUE
	RETURN
	END
!
!
!
	SUBROUTINE RD_MODEL_DBLE(VAR,DESC)
	USE MOD_MODEL_FILE
	IMPLICIT NONE
!
	REAL*8 VAR
	CHARACTER*(*) DESC
	INTEGER I,K
!
	DO I=1,N_REC
	  K=INDEX(STR_STORE(I),'['//TRIM(DESC)//']')
	  IF(K .NE. 0)THEN
	    READ(STR_STORE(I),*)VAR
	    GOTO 2000
	  END IF
	END DO
	WRITE(LU_ER,*)'Unable to read option from MODEL file'
	WRITE(LU_ER,*)'Error in RD_MODEL_DBLE: Key word: ',TRIM(DESC)
	STOP
!
2000	CONTINUE
	RETURN
	END
!
	SUBROUTINE RD_MODEL_INT(VAR,DESC)
	USE MOD_MODEL_FILE
	IMPLICIT NONE
!
	INTEGER VAR
	CHARACTER*(*) DESC
	INTEGER I,K
!
	DO I=1,N_REC
	  K=INDEX(STR_STORE(I),'['//TRIM(DESC)//']')
	  IF(K .NE. 0)THEN
	    READ(STR_STORE(I),*)VAR
	    GOTO 2000
	  END IF
	END DO
	WRITE(LU_ER,*)'Unable to read option from MODEL file'
	WRITE(LU_ER,*)'Error in RD_MODEL_INT: Key Word: ',TRIM(DESC)
	STOP
!
2000	CONTINUE
	RETURN
	END
