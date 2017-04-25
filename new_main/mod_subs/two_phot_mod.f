!
! Module containing 2-photon data. It is also used to
! contain rate data so that STEQ (statistical equilibrium 
! equations) and BA matrices can be updated.
!
	MODULE TWO_PHOT_MOD
	INTEGER N_TWO
!
! Altered 19-Aug-2015: Extra variables added to improve two photon opacity calculation (cur_hmi, 8-Jul-2015).
! Created 26-Jun-1998
!
! Reaction data
!
! TYPE_TWO       : Indentifies formula to evaluate photon emission as a 
!                     function of frequency.
! SPEC_ID        : Species identifer (e.g. HI)
! LOW_NAME_TWO   : Level name of lower level
! UP_NAME_TWO    : Level name of upper level
! A_LOW_NAME_TWO : Alterantive name for lower level.
! A_LOW_NAME_TWO : Alterantive name for upper level.
!
          INTEGER, ALLOCATABLE :: TYPE_TWO(:)
	  CHARACTER*8, ALLOCATABLE :: SPEC_ID_TWO(:)
	  CHARACTER*30, ALLOCATABLE :: LOW_NAME_TWO(:)
	  CHARACTER*30, ALLOCATABLE :: UP_NAME_TWO(:)
	  CHARACTER*30, ALLOCATABLE :: A_LOW_NAME_TWO(:)
	  CHARACTER*30, ALLOCATABLE :: A_UP_NAME_TWO(:)
!
! COEF_TWO       : Array containing fit information.
!
	  INTEGER, PARAMETER :: MAX_COEF=10
          INTEGER, ALLOCATABLE :: N_COEF_TWO(:)
	  REAL*8, ALLOCATABLE:: COEF_TWO(:,:)
	
!
! Arrays required to update STEQ (statistical equilibrium) and BA (variation
! of STEQ) arrays.
!
! Z_TWO          : Charge on core (or ion charge_
! G_LOW_TWO      : Statistical weight of lower level
! G_UP_TWO       : Statistical weight of upper level
! LOW_LEV_TWO    : Pointer to lower level in POPS and STEQ arrays
! UP_LEV_TWO     : Pointer to upper level in POPS and STEQ arrays
! FS_RAT_LOW     : LTE ratio of population in FULL ATOM to that in the super 
!                       level for the lower level
! FS_RAT_UP      : LTE ratio of population in FULL ATOM to that in the super 
!                       upper for the lower level
! FREQ_TWO       : Maximum frequency for 2-photon transition
!
! DOWN_RATE_TWO  : Rate of transitions into the lower level
! UP_RATE_TWO    : Rate of transitions into the upper level
!
!TWO_PHOT_AVAILABLE  : Inidcates transition has been set, and species available.
!
          REAL*8, ALLOCATABLE :: Z_TWO(:)
          REAL*8, ALLOCATABLE :: G_LOW_TWO(:)
          REAL*8, ALLOCATABLE :: G_UP_TWO(:)
          REAL*8, ALLOCATABLE :: FS_RAT_LOW(:,:)
          REAL*8, ALLOCATABLE :: FS_RAT_UP(:,:)
          REAL*8, ALLOCATABLE :: FREQ_TWO(:)
!
	  REAL*8, ALLOCATABLE :: DOWN_RATE_TWO(:,:)
	  REAL*8, ALLOCATABLE :: UP_RATE_TWO(:,:)
	  REAL*8, ALLOCATABLE :: PHOT_OC_TWO(:,:)
!
          INTEGER, ALLOCATABLE :: ION_ID_TWO(:)
          INTEGER, ALLOCATABLE :: ION_LOW_LEV_TWO(:)
          INTEGER, ALLOCATABLE :: ION_UP_LEV_TWO(:)
          INTEGER, ALLOCATABLE :: LOW_LEV_TWO(:)
          INTEGER, ALLOCATABLE :: UP_LEV_TWO(:)
          INTEGER, ALLOCATABLE :: LST_FREQ_INDX_TWO(:)
!	  
	  LOGICAL, ALLOCATABLE :: TWO_PHOT_AVAILABLE(:)
	  LOGICAL, ALLOCATABLE :: TWO_PHOT_COEF_FIXED(:)
!
! Use to indicate thate data arrays need to be initialized. This is
! reset to FALSE in STEQ_BA_TWO_PHOT.
!
	  LOGICAL INITIALIZE_TWO
	  LOGICAL DO_TWO_PHOT
!
	  CHARACTER(LEN=12) :: TWO_METHOD=' '
	  CHARACTER(LEN=11) TWO_PHOT_FORMAT_DATE
!
	END MODULE TWO_PHOT_MOD
!
! Routine to read in data for 2-photon emisison transitions. For example,
! transition for the 2s2Se state of H.
!
! The first data line in the file should contain the string:
!
!	"N		!Number of 2-photon transitions"
!
! where N is the nubmer of 2 photon transitions.
!
	SUBROUTINE RD_TWO_PHOT(LUIN,INCL_TWO_PHOT)
	USE TWO_PHOT_MOD	
	IMPLICIT NONE
!
! Created 26-Jun-1998
!
	INTEGER LUIN
	LOGICAL INCL_TWO_PHOT
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER IOS
	INTEGER I,J,L
	CHARACTER*132 STRING
!
	TWO_PHOT_FORMAT_DATE=' '
	IF(INCL_TWO_PHOT)THEN
	  DO_TWO_PHOT=.TRUE.
	ELSE
	  DO_TWO_PHOT=.FALSE.
	  N_TWO=0
	  RETURN
	END IF
!
	LUER=ERROR_LU()
	CALL GEN_ASCI_OPEN(LUIN,'TWO_PHOT_DATA','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RD_TWO_PHOT'
	  WRITE(LUER,*)'Unable to open TWO_PHOT_DATA'
	  STOP
	END IF
!
! Read in number of reactions. We first skip over all comment lines. The
! first data line in the file should contain the string
!    "!Number of 2-photon transitions"
! Subsequent line which are BETWEEN transitions, and
! which begin with a ! or a blank, are ignored.
!
	  IOS=0
          L=0
	  DO WHILE(L .EQ. 0 .AND. IOS .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    L=INDEX(STRING,'!Number of 2-photon transitions')
	    I=INDEX(STRING,'!Format date')
	    IF(I .NE. 0)THEN
	      STRING=ADJUSTL(STRING)
	      I=INDEX(STRING,'  ')
	      TWO_PHOT_FORMAT_DATE=STRING(1:I)
	    END IF
	  END DO
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in RD_TWO_PHOT'
	    WRITE(LUER,*)'Number of 2-photon transitions string not found'
	    STOP
	  ELSE
	   READ(STRING,*)N_TWO
	  END IF
!
! We only allocate arrays necessary to store the atomic data which
! is read in.
!
	  IF(.NOT. ALLOCATED(TYPE_TWO))THEN
            ALLOCATE (TYPE_TWO(N_TWO))
	    ALLOCATE (SPEC_ID_TWO(N_TWO))
	    ALLOCATE (LOW_NAME_TWO(N_TWO))
	    ALLOCATE (UP_NAME_TWO(N_TWO))
	    ALLOCATE (A_LOW_NAME_TWO(N_TWO))
	    ALLOCATE (A_UP_NAME_TWO(N_TWO))
	    ALLOCATE (N_COEF_TWO(N_TWO))
	    ALLOCATE (COEF_TWO(N_TWO,MAX_COEF))
	    ALLOCATE (TWO_PHOT_COEF_FIXED(N_TWO))
	    INITIALIZE_TWO=.TRUE.
	  END IF
!
	  IF(INITIALIZE_TWO)THEN
	    A_LOW_NAME_TWO(:)=' '
	    A_UP_NAME_TWO(:)=' '
	    COEF_TWO(:,:)=0.0D0
	  END IF
	  TWO_PHOT_COEF_FIXED(:)=.FALSE.
!
! Read in 2-photon transition data.
!
	  DO J=1,N_TWO
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    END DO
!
	    L=INDEX(STRING,'!2-photon identification')
	    IF(L .EQ. 0)THEN
	      WRITE(LUER,*)'Error in RD_TWO_PHOT'
	      WRITE(LUER,*)'Invalid format for 2-photon data file'
	      WRITE(LUER,*)STRING
	      STOP
	    END IF
	    STRING(L:)=' '
C
C Replace all TABS with 2 spaces. NB: CHAR(9) is TAB.
C
	    DO WHILE(INDEX(STRING,CHAR(9)) .NE. 0)
	      L=INDEX(STRING,CHAR(9)) 
	      STRING(L:)='  '//STRING(L+1:)
	    END DO
C
	    IF(STRING .EQ. ' ')THEN
	      WRITE(LUER,*)'Error in RD_TWO_PHOT'
	      WRITE(LUER,*)'Invalid 2-photon string - no species ID'
	      STOP
	    END IF
	    DO WHILE(STRING(1:1) .EQ. ' ')
	      STRING(1:)=STRING(2:)
	    END DO
	    L=INDEX(STRING,' ')
	    SPEC_ID_TWO(J)=STRING(1:L-1)
	    STRING(1:)=STRING(L+1:)
!
	    IF(STRING .EQ. ' ')THEN
	      WRITE(LUER,*)'Error in RD_TWO_PHOT'
	      WRITE(LUER,*)'Invalid 2-photon string -- no lower level ID'
	      STOP
	    END IF
	    DO WHILE(STRING(1:1) .EQ. ' ')
	      STRING(1:)=STRING(2:)
	    END DO   
	    L=INDEX(STRING,' ')
	    LOW_NAME_TWO(J)=STRING(1:L-1)

	    STRING(1:)=STRING(L+1:)
	    IF(STRING .EQ. ' ')THEN
	      WRITE(LUER,*)'Error in RD_TWO_PHOT'
	      WRITE(LUER,*)'Invalid 2-photn string -- no upper level ID'
	      STOP
	    END IF
	    DO WHILE(STRING(1:1) .EQ. '  ')
	      STRING(1:)=STRING(2:)
	    END DO   
	    L=INDEX(STRING,'  ')
	    UP_NAME_TWO(J)=STRING(1:L-1)
	    STRING(1:)=STRING(L+1:)
!
! Get alternate level names. These are useful for H and He2 in case
! we split the n levels into individual 'l' states.
!
	    IF(STRING .NE. ' ')THEN		!Get alternate level names
!
	      DO WHILE(STRING(1:1) .EQ. ' ')
	        STRING(1:)=STRING(2:)
	      END DO   
	      L=INDEX(STRING,'  ')
	      A_LOW_NAME_TWO(J)=STRING(1:L-1)
	      STRING(1:)=STRING(L+1:)
!
	      IF(STRING .EQ. ' ')THEN
	        WRITE(LUER,*)'Error in RD_TWO_PHOT'
	        WRITE(LUER,*)'Invalid 2-photon string'
	        WRITE(LUER,*)'No alternate upper level ID'
	        STOP
	      END IF
	      DO WHILE(STRING(1:1) .EQ. ' ')
	        STRING(1:)=STRING(2:)
	      END DO   
	      L=INDEX(STRING,'  ')
	      A_UP_NAME_TWO(J)=STRING(1:L-1)
	      STRING(1:)=STRING(L+1:)
	    END IF
!
	    READ(LUIN,*)TYPE_TWO(J)
	    READ(LUIN,*)N_COEF_TWO(J)
	    DO I=1,N_COEF_TWO(J)
	      READ(LUIN,*)COEF_TWO(J,I)
	    END DO
	  END DO
!
	CLOSE(LUIN)
!
	RETURN
	END
