!
! Routine to read in charge exchange reactions and cross-sections.
! The first data line in the file should contain the string:
!
!	"N		!Number of charge exchange reactions"
!
! where N is the nubmer of chage exchage reactions. For a reaction of the form
!
!        X++ + Y+  -->  X+ + Y++
!
! the SPECIES and LEVEL of reactant should be specified IN ORDER on the
! same line:
!           X++ must be specfied before Y+
!           X+  must be specfied before Y++
!
	SUBROUTINE RD_CHG_EXCH_V3(LUIN,INCL_CHG_EXCH)
	USE CHG_EXCH_MOD_V3
	IMPLICIT NONE
!
! Altered 22-Sep-2016: Error reporting improved
! Altered 04-Dec-200:  Bug fix: Could enter infinite loop when left adjusting
!                      reaction string.
!                      Now use / to allow the specification of one alternate
!                      name (not {} since some names now contain these 
!                      brackets).
	INTEGER LUIN
	LOGICAL INCL_CHG_EXCH
!
! Local variables
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER IOS
	INTEGER N_COEF
	INTEGER I,K,L,LB,RB
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
	CHARACTER*132 STRING,OLD_STRING
	CHARACTER*11  FORMAT_DATE 
!
	IF(INCL_CHG_EXCH)THEN
	  DO_CHG_EXCH=.TRUE.
	ELSE
	  DO_CHG_EXCH=.FALSE.
	  N_CHG=0
	  RETURN
	END IF
!
	LUER=ERROR_LU()
	WRITE(LUER,*)'Beginning to read Charge data'
	CALL GEN_ASCI_OPEN(LUIN,'CHG_EXCH_DATA','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Unable to open CHG_EXCH_DATA in RD_CHG_EXCH'
	  WRITE(LUER,*)'IOS=',IOS
	  STOP
	END IF
!
! Read in number of reactions. We first skip over all comment lines. The
! first data line in the file should contain the string
!    "!Number of charge exchange reactions"
! Subsequent line which are between charge exchange reactions, and
! which begin with a ! or a blank, are ignored.
!
	IOS=0
        L=0
	DO WHILE(L .EQ. 0 .AND. IOS .EQ. 0)
	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	  L=INDEX(STRING,'!Number of charge exchange reactions')
	END DO
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RD_CHG_EXCH'
	  WRITE(LUER,*)'Number of charge exchange reactions sting not found'
	  STOP
	ELSE
	 READ(STRING,*)N_CHG_RD
	END IF
!
	FORMAT_DATE=' '
	READ(LUIN,'(A)')STRING
	IF( INDEX(STRING,'!Format date') .NE. 0)THEN
	  FORMAT_DATE=STRING(1:11)
	  READ(LUIN,'(A)')STRING
	  IF(INDEX(STRING,'!Modification date') .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RD_CHG_EXCH'
            WRITE(LUER,*)'Modification file date must follow',
	1                       ' format date'
	    STOP
	  END IF
        ELSE
!
! For consistency with old format
!
          BACKSPACE(LUIN)
	END IF
!
! We only allocate arrays necessary to store the atomic data which
! is read in.
!
	IF(.NOT. ALLOCATED(TYPE_CHG_RD))THEN
          ALLOCATE (TYPE_CHG_RD(N_CHG_RD),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (TLO_CHG_RD(N_CHG_RD),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (THI_CHG_RD(N_CHG_RD),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (CHG_INCLUDED_RD(N_CHG_RD),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (COEF_CHG_RD(N_CHG_RD,N_COEF_MAX),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (SPEC_ID_CHG_RD(N_CHG_RD,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (LEV_NAME_CHG_RD(N_CHG_RD,4),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ALT_LEV_NAME_CHG_RD(N_CHG_RD,4),STAT=IOS)
	END IF
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RD_CHG_EXCH'
	  WRITE(LUER,*)'Unable to allocate all vectors in RD_CHG_EXCH'
	  STOP
	END IF
        ALT_LEV_NAME_CHG_RD(1:N_CHG_RD,1:4)='No_alternative_level_name'
!
	DO I=1,N_CHG_RD
	  STRING=' '
	  DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,'(A,I3)')'Error reading charge exchange reaction # ',I
	      WRITE(LUER,'(A,I3)')'Last reaction read follows'
	      WRITE(LUER,'(A)')TRIM(OLD_STRING)
	      STOP
	    END IF
	  END DO
	  OLD_STRING=STRING
!
	  DO K=1,4
	    STRING=ADJUSTL(STRING)
	    IF(STRING .EQ. ' ')THEN
	      WRITE(LUER,*)'Insufficient information in reaction string (1).'
	      WRITE(LUER,*)' ',TRIM(OLD_STRING)
	      STOP
	    END IF
	    L=INDEX(STRING,' ')
	    IF( L .NE. INDEX(STRING,'  '))THEN
	      WRITE(LUER,*)'Error in RD_CHG_EXCH'
	      WRITE(LUER,*)'Use at least 2 spaces to separate reaction data'
	      WRITE(LUER,*)TRIM(OLD_STRING)
	      STOP
	    END IF
	    SPEC_ID_CHG_RD(I,K)=STRING(1:L-1)
	    STRING(1:)=STRING(L+1:)
	    IF(STRING .EQ. ' ')THEN
	      WRITE(LUER,*)'Insufficient information in reaction string (2).'
	      WRITE(LUER,*)' ',TRIM(OLD_STRING)
	      STOP
	    END IF
	    STRING=ADJUSTL(STRING)
!
	    L=INDEX(STRING,' ')
	    IF( L .NE. INDEX(STRING,'  '))THEN
	      WRITE(LUER,*)'Error in RD_CHG_EXCH'
	      WRITE(LUER,*)'Use at least 2 spaces to separate reaction data'
	      STOP
	    END IF
	    IF(FORMAT_DATE .EQ. '01-Oct-1999')THEN
	      LB=INDEX(STRING(1:L-1),'{')
	      IF(LB .EQ. 0)THEN
	        LEV_NAME_CHG_RD(I,K)=STRING(1:L-1)
	      ELSE
	        RB=INDEX(STRING(1:L-1),'}')
	        IF(RB .LE. LB+2)THEN
	           WRITE(LUER,*)'Error in RD_CHG_EXCH'
	           WRITE(LUER,*)
	1            'Inavlid number of } for alternative level name'
	          STOP
	        END IF
	        LEV_NAME_CHG_RD(I,K)=STRING(1:LB-1)
	        ALT_LEV_NAME_CHG_RD(I,K)=STRING(LB+1:RB-1)
	      END IF
	    ELSE
	      LB=INDEX(STRING(1:L-1),'/')
	      IF(LB .EQ. 0)THEN
	        LEV_NAME_CHG_RD(I,K)=STRING(1:L-1)
	      ELSE
	        LEV_NAME_CHG_RD(I,K)=STRING(1:LB-1)
	        ALT_LEV_NAME_CHG_RD(I,K)=STRING(LB+1:L-1)
	      END IF
	    END IF
	    STRING(1:)=STRING(L+1:)
	  END DO
!
	  READ(LUIN,*,IOSTAT=IOS)TYPE_CHG_RD(I)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading TYPE of charge reaction for the following reaction'
	    WRITE(LUER,*)' ',TRIM(OLD_STRING)
	    STOP
	  END IF
!
	  READ(LUIN,*,IOSTAT=IOS)N_COEF
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading # of charge coefficeints for the following reaction'
	    WRITE(LUER,*)' ',TRIM(OLD_STRING)
	    STOP
	  END IF
	  IF(N_COEF .GT. N_COEF_MAX)THEN
	    WRITE(LUER,*)'Error in RD_CHG_EXH --- N_COEF_MAX too small'
	    STOP
	  END IF
!
	  IF(FORMAT_DATE .EQ. ' ')THEN
	    READ(LUIN,*,IOSTAT=IOS)(COEF_CHG_RD(I,K),K=1,N_COEF)
	    TLO_CHG_RD(I)=0.0D0
	    THI_CHG_RD(I)=100.0D0
	  ELSE
	    READ(LUIN,*,IOSTAT=IOS)(COEF_CHG_RD(I,K),K=1,N_COEF),TLO_CHG_RD(I),THI_CHG_RD(I)
	  END IF
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error reading charge coefficients for the following reaction'
	    WRITE(LUER,*)' ',TRIM(OLD_STRING)
	    STOP
	  END IF
!
! Do some consistency checks on the charge exhcnage reactions. We check that
! all species are unique, and that the matced species (i.e. 1 & 3, 2 & 4) have
! the same Element identifier.
!
	  DO K=1,3
	    DO L=K+1,4
	      IF(SPEC_ID_CHG_RD(I,L) .EQ. SPEC_ID_CHG_RD(I,K))THEN
	        WRITE(LUER,*)'Error in RD_CHG_EXCH'
	        WRITE(LUER,*)'Duplication of species ID, reaction:',I
	        WRITE(LUER,'(A,3X,A)')SPEC_ID_CHG_RD(I,L),SPEC_ID_CHG_RD(I,K)
	        STOP
	      END IF
	    END DO
	  END DO
!
	  L=1
	  IF(SPEC_ID_CHG_RD(I,1)(2:2) .GE. 'a' .AND. 
	1             SPEC_ID_CHG_RD(I,1)(2:2) .LE. 'z')L=2
	  IF(SPEC_ID_CHG_RD(I,1)(1:L) .NE. SPEC_ID_CHG_RD(I,3)(1:L))THEN
	    WRITE(LUER,*)'Error in RD_CHG_EXCH',
	1        ' SPECIES 1 and 3 don''t match for reaction',I
	    STOP
	  END IF
!
	  L=1
	  IF(SPEC_ID_CHG_RD(I,2)(2:2) .GE. 'a' .AND. 
	1              SPEC_ID_CHG_RD(I,2)(2:2) .LE. 'z')L=2
	  IF(SPEC_ID_CHG_RD(I,2)(1:L) .NE. SPEC_ID_CHG_RD(I,4)(1:L))THEN
	    WRITE(LUER,*)'Error in RD_CHG_EXCH',
	1	' SPECIES 2 and 4 don''t match for reaction',I
	    WRITE(6,'(A,5X,A)')SPEC_ID_CHG_RD(I,2)(1:L),SPEC_ID_CHG_RD(I,4)(1:L)
	    STOP
	  END IF
!
	END DO
!
! This will be reset later.
!
	CHG_INCLUDED_RD(:)=.TRUE.
!
	WRITE(LUER,*)'Charge data successfully read'
	CLOSE(LUIN)
!
	RETURN
	END
