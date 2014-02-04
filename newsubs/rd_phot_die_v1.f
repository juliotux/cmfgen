!
! Subroutine to read in Low temperature dielectronic data for a speciec XzV.
! These data will be combined directly with the photoionization cross-sections.
!
! The user indicates whether to include NORMAL (ie AUTO) and WI LTDR
! transitions. LTDR transition are returned only for those levels included
! in the model atom.
!
! NB : AUTO --- Levels permitted to autoioonize in pure LS coupling.
!      WI   --- Levels not permitted to autoionize in LS copuling.
!
	SUBROUTINE RD_PHOT_DIE_V1(ID,EDGE,LEVELNAME,NXzV,GION_GS,
	1            VSMOOTH_KMS,DO_AUTO,DO_WI,
	1            DESC,LUIN,LUOUT,FILENAME)
!
! Photoionization data module.
!
	USE PHOT_DATA_MOD
!
	IMPLICIT NONE
!
! Altered 18-Feb-2006 : Minor bug fix -- accessed outside array when checking CROSS_TYPE.
! Altered 02-Oct-2005 : Changed eallocation section and have to compile with -g for pgf95
! Altered 02-Sep-2004 : Reinserted writing out header infromation located before date.
! Altered 18-MAr-2002 : Bug fix for profile with negligble Doppler core.
! Altered 12-Mar-2002 : Warning written if DO_AUTO is true and if using
!                         Opacity Project cross-sections
! Altered 28-Jul-2001 : New format for dielectronic data installed.
!                         Free-free resonances now treated.
! Altered 04-Jul-2001 : Bux fixed for double IF clause. Second argument
!                         was being accessed outside range. Inserted
!                         MIN(J,CNT) instead of just J.
! Altered 21-Dec-1997 : Bug fix caused by handling split lower levels.
! Altered 05-Dec-1996 : GEN_ASCI_OPEN used to OPEN ASCI files.
! Altered 28-May-1996 : GEN_SEQ_OPEN replaced by F90 OPEN
! Created 18-Dec-1995
!
	INTEGER ID			!Species identifier (integer key)
	INTEGER NXzV			!Number of levels in FULL atom
	REAL*8 EDGE(NXzV)		!Ionization frequency (10^15 Hz)
	REAL*8 GION_GS			!St. Weight of g.s. of ion (eg CIII).
	INTEGER LUIN,LUOUT
	LOGICAL DO_AUTO
	LOGICAL DO_WI
	CHARACTER*(*) LEVELNAME(NXzV)
!
! VSMOOTH_KMS is the velocity (in km/s) used to smooth the Lorentzian/Voigt
! profile to ensure that it is not undersampled. In principal it should
! be a function of T and VTURB.
!
	REAL*8 VSMOOTH_KMS
	CHARACTER*(*) DESC
	CHARACTER*(*) FILENAME
!
! External functions.
!
	REAL*8 FUN_PI,SPEED_OF_LIGHT
	INTEGER ERROR_LU,ICHRLEN
	EXTERNAL ERROR_LU,FUN_PI,ICHRLEN,SPEED_OF_LIGHT
!
! Common block with opacity/emissivity constants.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
! Local variables
!
	REAL*8, PARAMETER :: IZERO=0
!
	REAL*8 EDGEDIE
	REAL*8 EINA
	REAL*8 GUPDIE
	REAL*8 C_KMS
	REAL*8 IONIZATION_ENERGY
	CHARACTER(LEN=80) TRANSDIE
	CHARACTER(LEN=30) LS_NAME
!
! NB: DIELEV is double precision so that it can be passed directly to
! INDEXX.
!
	REAL*8, ALLOCATABLE :: DIELEV(:)
	REAL*8, ALLOCATABLE :: VEC_DP_WRK(:)
	INTEGER, ALLOCATABLE :: VEC_INDX(:)
!
	REAL*8, ALLOCATABLE :: DIE_LEV_G(:)
	REAL*8, ALLOCATABLE :: DIE_LEV_ENERGY(:)
	REAL*8, ALLOCATABLE :: DIE_LEV_AUTO(:)
	CHARACTER(LEN=30), ALLOCATABLE :: DIE_LEV_NAME(:)
!
	INTEGER NUM_D_RD
	INTEGER NUM_FREE_FREE
	INTEGER NUM_LEVELS
!
	INTEGER I,J,LOOP,IOS,L1,L2,LUER
	INTEGER MNL,MNL_CNT,CNT,CNT_BEG
	INTEGER INC,WIINC,MIS,WIMIS,INDX_HASH
!
	INTEGER NUM_OF_MATCHES
	INTEGER LOW_PNT
	INTEGER UP_PNT 
	CHARACTER(LEN=30) LOW_NAME
	CHARACTER(LEN=30) UP_NAME
	LOGICAL NEW_FILE_FORMAT
!
	REAL*8 GSUM
	REAL*8 T1
	REAL*8 DEL_NU,NU_DOP
	CHARACTER*132 STRING
	LOGICAL SPLIT_J
!
	REAL*8 A10,A20,A30
	REAL*8 A1,A2,A3,M1,M2,M3		!Effective recombination rate.
	REAL*8 WIA1,WIA2,WIA3,WIM1,WIM2,WIM3    !(included and missing).
!
	LUER=ERROR_LU()
!
! Most of the dielectronic storage locations are initialized by
! RD_PHOT_XzV/
!
	PD(ID)%VSM_KMS=VSMOOTH_KMS
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
        NEW_FILE_FORMAT=.FALSE.
!
! Check whether we might be including permitted dielectronic transitions twice.
! This is not a sophisticated check --- it simply issues a warning to the user
! to double check atomic data assignments etc.
!
	IF(DO_AUTO)THEN
	  DO I=1,PD(ID)%MAX_TERMS
	    L1=PD(ID)%CROSS_TYPE(I,1)
	    IF(L1 .EQ. 20 .OR. L1 .EQ. 21)THEN
	      WRITE(LUER,'(A)')' '
	      WRITE(LUER,'(1X,79A)')('*',J=1,79)
	      WRITE(LUER,'(1X,79A)')('*',J=1,79)
	      WRITE(LUER,'(1X,A,A)')'WARNING when reading ',TRIM(FILENAME)
	      WRITE(LUER,'(1X,A)')'You have indicated that you wish to include permitted dielectronic transitions.'
	      WRITE(LUER,'(1X,A)')'However, the photoionization file contains OPACITY-PROJECT cross-sections.'
	      WRITE(LUER,'(1X,A)')'These may already contain the deilectronic transitions as resonances.'
	      WRITE(LUER,'(1X,79A)')('*',J=1,79)
	      WRITE(LUER,'(1X,79A)')('*',J=1,79)
	      EXIT
	    END IF
	  END DO
	END IF
!
! 
!
! Open file containing dielectronic data. We first check the data
! format (only 2 are available at present). We then read in
! heaer information (such as the number of transitions etc.).
!
	CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	   WRITE(LUER,*)'Error in RD_XzV_PHOT_DIE - cant open '//FILENAME
	   STOP
	  END IF
!
! Read until come across first header. This should be '!Format date' if file
! is in the new format, or simply 'Date' if in old format.
!
	  L1=0; L2=0
	  DO WHILE (L1 .EQ. 0 .AND. L2 .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    WRITE(LUOUT,'(A)')STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error reading ''DATE'' from '//FILENAME
	      WRITE(LUER,*)'IOSTAT=',IOS
	      STOP
	    END IF
	    L1=INDEX(STRING,'!Format date')
	    L2=INDEX(STRING,'!Date')
	  END DO
	  IF(L1 .NE. 0)THEN
	    NEW_FILE_FORMAT=.TRUE.
	    STRING=ADJUSTL(STRING)
	    L1=INDEX(STRING,' ')
	    IF(STRING(1:L1) .NE. '27-Jul-2001')THEN
	      WRITE(LUER,*)'Format date is invalid'
	      STOP
	    END IF
	  END IF
!
! Initialize variables. Latter we check they have been successfully read in.
!
	  NUM_LEVELS=-1000
	  NUM_D_RD=-1000
	  NUM_FREE_FREE=-1000
	  IONIZATION_ENERGY=-1000
	  SPLIT_J=.FALSE.
!
! Read in all keywords. We continue reading until we come across 
! a blank line. The data after 'Format date' (or 'Date') must be
! contiguous.
!
	  STRING='Not yet finished '
          DO WHILE( STRING .NE. ' ')
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error finding blank line after header date in '//FILENAME
	      STOP
	    END IF
	    WRITE(LUOUT,'(A)')STRING
	    IF(INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
	      READ(STRING,*)NUM_LEVELS
	    ELSE IF(INDEX(STRING,'!Number of dielectronic transitions') .NE. 0)THEN
	      READ(STRING,*)NUM_D_RD
	    ELSE IF(INDEX(STRING,'!Number of free-free transitions') .NE. 0)THEN
	      READ(STRING,*)NUM_FREE_FREE
	    ELSE IF(INDEX(STRING,'!Ionization energy') .NE. 0)THEN
	      READ(STRING,*)IONIZATION_ENERGY
	    ELSE IF(INDEX(STRING,'!Split J levels') .NE. 0)THEN
	      READ(STRING,*)SPLIT_J
	    ELSE IF(INDEX(STRING,'!Format date') .NE. 0)THEN
	      WRITE(LUER,*)'Error reading dielectronic date from '//FILENAME
	      WRITE(LUER,*)'Format date must be first date in file'
	      STOP
	    END IF
	  END DO
!
! Now check that the data read in was valid.
!
	  IF(NEW_FILE_FORMAT .AND. NUM_LEVELS .LE. 0)THEN
	    WRITE(LUER,*)'Error reading # of levels from '//FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  ELSE IF(.NOT. NEW_FILE_FORMAT)THEN
	   NUM_LEVELS=0
	  END IF
!
	  IF(NUM_D_RD .LT. 0)THEN
	    WRITE(LUER,*)'Error reading # of delectronic transitions from '//FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  END IF
!
	  IF(NEW_FILE_FORMAT .AND. NUM_FREE_FREE .LT. 0)THEN
	    WRITE(LUER,*)'Error reading # of free-fee transitions from '//FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  ELSE IF(.NOT. NEW_FILE_FORMAT)THEN
	    NUM_FREE_FREE=0
	  END IF
!
	  IF(NEW_FILE_FORMAT .AND. IONIZATION_ENERGY .LT. 0)THEN
	    WRITE(LUER,*)'Error reading ionization energy from '//FILENAME
	    WRITE(LUER,*)'IOSTAT=',IOS
	    STOP
	  ELSE IF(.NOT. NEW_FILE_FORMAT)THEN
	    IONIZATION_ENERGY=0
	  END IF
!
! 
! Now allocate all dynamic arrays.
!
! NB: We set NDIE_MAX=3*NUM_D_RD to allow for the extra dielectronic
!     transitions introudced when the data is NOT split into individual J
!     states, but the lower levels in the model atom might be. 
!     This is not very satisfactory.
!
	IF(SPLIT_J)THEN
	  PD(ID)%NDIE_MAX=NUM_D_RD
	ELSE
	  PD(ID)%NDIE_MAX=3*NUM_D_RD
	END IF
	ALLOCATE (PD(ID)%OSC(PD(ID)%NDIE_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (PD(ID)%GAMMA(PD(ID)%NDIE_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (PD(ID)%NU_ZERO(PD(ID)%NDIE_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (PD(ID)%NU_MIN(PD(ID)%NDIE_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (PD(ID)%NU_MAX(PD(ID)%NDIE_MAX),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE (DIELEV(PD(ID)%NDIE_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (VEC_DP_WRK(PD(ID)%NDIE_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (VEC_INDX(PD(ID)%NDIE_MAX),STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE (PD(ID)%ST_INDEX(NXzV),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (PD(ID)%END_INDEX(NXzV),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error allocating memory in RD_PHOT)DIE_V2'
	  WRITE(LUER,*)'Currently reading ',TRIM(FILENAME)
	  STOP
	END IF
!
! Allocate memory for free-free transitions.
!
	PD(ID)%NUM_FF=NUM_FREE_FREE
	IF(NUM_FREE_FREE .NE. 0)THEN
	  IF(IOS .EQ. 0)ALLOCATE (PD(ID)%FF_GF(NUM_FREE_FREE))
	  IF(IOS .EQ. 0)ALLOCATE (PD(ID)%FF_GAMMA(NUM_FREE_FREE),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PD(ID)%FF_NU_ZERO(NUM_FREE_FREE),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PD(ID)%FF_NU_EXCITE(NUM_FREE_FREE),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PD(ID)%FF_NU_MIN(NUM_FREE_FREE),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (PD(ID)%FF_NU_MAX(NUM_FREE_FREE),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating FF memory in RD_PHOT)DIE_V2'
	    WRITE(LUER,*)'Currently reading ',TRIM(FILENAME)
	    STOP
	  END IF
	END IF
!
!
! The following loop reads in, and saves, data from a file in the new
! format.
!
	IF(NEW_FILE_FORMAT)THEN
!
! Allocate memory for energy levels.
!
	  ALLOCATE (DIE_LEV_NAME(NUM_LEVELS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DIE_LEV_ENERGY(NUM_LEVELS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DIE_LEV_G(NUM_LEVELS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (DIE_LEV_AUTO(NUM_LEVELS),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocating FF memory in RD_PHOT)DIE_V2'
	    WRITE(LUER,*)'Currently reading ',TRIM(FILENAME)
	    STOP
	  END IF
!
! Now read in the energy level data. Blank lines and comments (i.e., lines beginning
! with a !) can be inserted between the energy levels.
!
	  DO I=1,NUM_LEVELS
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
              READ(LUIN,'(A)')STRING
	    END DO
	    STRING=ADJUSTL(STRING)
	    L1=INDEX(STRING,'  ')
	    DIE_LEV_NAME(I)=STRING(1:L1-1)
	    READ(STRING(L1:),*)DIE_LEV_G(I),DIE_LEV_ENERGY(I),DIE_LEV_AUTO(I)
!
	    WRITE(LUER,*)DIE_LEV_NAME(I),DIE_LEV_G(I)
	  END DO
!
! Find header line for the transition data. This signifies the beginning of the
! atomic data (except for blank lines and comments')
!
	  STRING=' '
          DO WHILE(INDEX(STRING,'Transition') .EQ. 0 .AND.
	1          INDEX(STRING,'Lam(A)') .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error reading Transition header in '//FILENAME
	      STOP
	    END IF
	  END DO
!
! Finally we can read in the dielectronic transition data. The first
! data set is for transition to BOUND levels. Blank lines and comments 
! (i.e., lines beginning with a !) can be inserted between the energy levels.
!
	  NUM_OF_MATCHES=0
	  CNT=0
	  DO LOOP=1,NUM_D_RD
	     STRING=' '
	     DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
               READ(LUIN,'(A)')STRING
	     END DO
!
! Check whether this LTDR transition can be accomodated.
!
	    IF(.NOT. SPLIT_J .AND. CNT+5 .GT. PD(ID)%NDIE_MAX)THEN
	      WRITE(LUER,*)'PD(ID)%NDIE_MAX not large enough in ',
	1        'RD_XzV_PHOT_DIE - '//FILENAME
	      WRITE(LUER,*)'PD(ID)%NDIE_MAX=',PD(ID)%NDIE_MAX,'NUM_DIE=',NUM_D_RD
	      STOP
	    END IF
	    STRING=ADJUSTL(STRING)
	    WRITE(LUER,*)TRIM(STRING)
!
! Allows for gaps between first level name and '-'.
!
	    L1=INDEX(STRING,'-'); LOW_NAME=STRING(1:L1-1)
	    STRING=STRING(L1+1:)
	    STRING=ADJUSTL(STRING)
	    L1=INDEX(STRING,'  '); UP_NAME=STRING(1:L1-1)
!
! Find match with upper-level name:
!
	    UP_PNT=0
	    DO I=1,NUM_LEVELS
	      IF(UP_NAME .EQ. DIE_LEV_NAME(I))THEN
	         UP_PNT=I
	         EXIT
	      END IF
	    END DO
	    IF(UP_PNT .EQ. 0)THEN
	      WRITE(LUER,*)'Error matching upper levelname in ',TRIM(FILENAME)
	      WRITE(LUER,*)'UP_NAME=',UP_NAME
	      STOP
	    END IF
!
! Now match-lower level name: If the dielectronic levels are split,
! and the model atom is not, we can get more than one match.
!
	    CNT_BEG=CNT+1; GSUM=0.0D0
	    DO I=1,NXzV
	      J=INDEX(LEVELNAME(I),'[')-1
	      IF(J .EQ. 0)J=LEN_TRIM(LEVELNAME(I))
	      IF(LEVELNAME(I) .EQ. LOW_NAME)THEN
	         CNT=CNT+1
                 DIELEV(CNT)=I
	         READ(STRING(L1:),*)PD(ID)%OSC(CNT)
	         IF(CNT .NE. CNT_BEG)THEN
	           WRITE(LUER,*)'Error matching lower-level names'
	           WRITE(LUER,*)'Currently reading',TRIM(FILENAME)
	           WRITE(LUER,*)'Too many matches'
	           STOP
	         END IF
	      ELSE IF(LEVELNAME(I)(1:J) .EQ. LOW_NAME)THEN
	         CNT=CNT+1
                 DIELEV(CNT)=I
	         READ(STRING(L1:),*)PD(ID)%OSC(CNT)
	      END IF
	    END DO
	    IF(INDEX(STRING,'#') .EQ. 0)THEN
	      IF(.NOT. DO_AUTO)CNT=CNT_BEG-1
	    ELSE IF(.NOT. DO_WI)THEN
	      CNT=CNT_BEG-1
	    END IF
!
	    IF(CNT .LT. CNT_BEG)EXIT
	    NUM_OF_MATCHES=NUM_OF_MATCHES+1
!
! Store data in the photoionization module. If lower levels
! are split, and upper are not, more than one match may be present.
! In that case we assume A is proportional to the statistical weight
! of the lower level, which is true in LS couplinga. This is equivalent
! to assuming that the f-values for the individual multiplet transitions
! are identical.
!
	     DO J=CNT_BEG,CNT
	       PD(ID)%GAMMA(J)=DIE_LEV_AUTO(UP_PNT)/4.0D0/FUN_PI()
	       PD(ID)%NU_ZERO(J)=1.0D-10*C_KMS*
	1             (DIE_LEV_ENERGY(UP_PNT)-IONIZATION_ENERGY)+EDGE(NINT(DIELEV(J)))
	       WRITE(LUER,*) PD(ID)%OSC(J),PD(ID)%GAMMA(J),PD(ID)%NU_ZERO(J),DIELEV(J)
	     END DO      
	  END DO 
!
	  WRITE(LUER,*)'Number of transitions read in is',NUM_OF_MATCHES,CNT     
!
! Now read in the free-free tranistions. These are always included in the
! CMFGEN calculation, if present.
!
	  DO LOOP=1,NUM_FREE_FREE
	     STRING=' '
	     DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
               READ(LUIN,'(A)')STRING
	     END DO
	     STRING=ADJUSTL(STRING)
!
! Get level names involved in this transition. We allow for gaps between 
! first level name and '-', and extra spaces.
!
	    L1=INDEX(STRING,'-'); LOW_NAME=STRING(1:L1-1)
	    STRING=STRING(L1+1:);  STRING=ADJUSTL(STRING)
	    L1=INDEX(STRING,'  '); UP_NAME=STRING(1:L1-1)
	    STRING=STRING(L1+1:);  STRING=ADJUSTL(STRING)
!
! Find match with lower-level name:
!
	    LOW_PNT=0
	    DO I=1,NUM_LEVELS
	      IF(LOW_NAME .EQ. DIE_LEV_NAME(I))THEN
	         LOW_PNT=I
	         EXIT
	      END IF
	    END DO
	    IF(LOW_PNT .EQ. 0)THEN
	      WRITE(LUER,*)'Error in RD_PHOT_DIE when reading free-free data'
	      WRITE(LUER,*)'ID=',ID,'  Filename=',TRIM(FILENAME)
	      WRITE(LUER,*)'No match for',LOW_NAME
	      STOP
	    END IF
!
! Find match with upper-level name:
!
	    UP_PNT=0
	    DO I=1,NUM_LEVELS
	      IF(UP_NAME .EQ. DIE_LEV_NAME(I))THEN
	         UP_PNT=I
	         EXIT
	      END IF
	    END DO
	    IF(UP_PNT .EQ. 0)THEN
	      WRITE(LUER,*)'Error in RD_PHOT_DIE when reading free-free data'
	      WRITE(LUER,*)'ID=',ID,'  Filename=',TRIM(FILENAME)
	      WRITE(LUER,*)'No match for',UP_NAME
	      STOP
	    END IF
!
	    IF(DIE_LEV_ENERGY(UP_PNT) .LT. DIE_LEV_ENERGY(LOW_PNT))THEN
	      WRITE(LUER,*)'Error in RD_PHOT_DIE when reading free-free data'
	      WRITE(LUER,*)'ID=',ID,'  Filename=',TRIM(FILENAME)
	      WRITE(LUER,*)'Energy levels reversed'
	      WRITE(LUER,*)LOW_NAME,UP_NAME
	      STOP
	    END IF
!
	    WRITE(LUER,*)'LOW_PNT=',LOW_PNT
	    PD(ID)%FF_NU_ZERO(LOOP)=1.0D-10*C_KMS*(DIE_LEV_ENERGY(UP_PNT)-DIE_LEV_ENERGY(LOW_PNT))
	    PD(ID)%FF_NU_EXCITE(LOOP)=1.0D-10*C_KMS*(DIE_LEV_ENERGY(LOW_PNT)-IONIZATION_ENERGY)
	    READ(STRING,*)PD(ID)%FF_GF(LOOP)
	    PD(ID)%FF_GF(LOOP)=PD(ID)%FF_GF(LOOP)*DIE_LEV_G(LOW_PNT)
	    PD(ID)%FF_GAMMA(LOOP)=(DIE_LEV_AUTO(LOW_PNT)+DIE_LEV_AUTO(UP_PNT))/4.0D0/FUN_PI()
	    WRITE(LUER,*)'UP_PNT=',UP_PNT
!
	  END DO
!
!
!
! Read in data file which is in the old 'OBSOLETE' data format.
!
	ELSE
!
! Skip all records until we come across the header.
!
	  STRING=' '
          DO WHILE(INDEX(STRING,'Transition') .EQ. 0)
	    READ(LUIN,'(A)')STRING
	  END DO
	  L1=INDEX(STRING,'Lam(A)')
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error reading Transition header in '//FILENAME
	    STOP
	  END IF
!
! Zero summation scalers for computing recombination rates.
!
	  A1=0.0D0
	  A2=0.0D0
	  A3=0.0D0
	  M1=0.0D0
	  M2=0.0D0
	  M3=0.0D0
	  WIA1=0.0D0
	  WIA2=0.0D0
	  WIA3=0.0D0
	  WIM1=0.0D0
	  WIM2=0.0D0
	  WIM3=0.0D0
	  INC=0				!LTDR counters.
	  WIINC=0
	  MIS=0
	  WIMIS=0
!
! Begin the reading of the data
	  CNT=0
	  DO LOOP=1,NUM_D_RD
!
	    CNT=CNT+1
!
! Check whether this LTDR transition can be accomodated.
!
	    IF(CNT+5 .GT. PD(ID)%NDIE_MAX)THEN
	      WRITE(LUER,*)'PD(ID)%NDIE_MAX not large enough in ',
	1        'RD_XzV_PHOT_DIE - '//FILENAME
	      WRITE(LUER,*)'PD(ID)%NDIE_MAX=',PD(ID)%NDIE_MAX,'NUM_DIE=',NUM_D_RD
	      STOP
	    END IF
!
! Skip blank records or comments  before reading in dielectronic transitions.
!
	    STRING=' '
	    DO WHILE(STRING .EQ. ' ' .OR. STRING(1:1) .EQ. '!')
	      READ(LUIN,'(A)')STRING
	    END DO
!
! Allows for gaps between first level name and '-'.
!
	    L1=INDEX(STRING,'-')
	    L1=INDEX(STRING(L1+1:),'  ')+L1-1
	    IF( L1 .LE. 0)THEN
	      WRITE(LUER,*)
	1       'Error reading in Transition Names from '//FILENAME
	      STOP
	    END IF
	    IF( LEN(TRANSDIE) .LT. L1+LEN_TRIM(DESC)+2 )THEN
	      WRITE(LUER,*)'Error - transition name in ',
	1                        'RD_XzV_PHOT_DIE too small'
	      WRITE(LUER,*)'STRING(1:L1)=',STRING(1:L1)
	      WRITE(LUER,*)'DESC=',DESC
	      STOP
	    END IF
	    TRANSDIE=TRIM(DESC)//'('//STRING(1:L1)//')'
!
! Perform a free-format read: f, A, Lam, Gu, and EDGE (energy above ionization limit).
!
! LEVEL is the name of the lower bound state, FL is the frequency of the
! stabalizing transition, GuP is the statistical weight of the autoionizing
! level, and EINA is the Einstein A of the stabalizing transition.
!
	    READ(STRING(L1+1:),*)PD(ID)%OSC(CNT),EINA,T1,GUPDIE,EDGEDIE
	    EDGEDIE=-EDGEDIE   	  !Since above ionization limit.
!
! Now compute the effective recombination coefficient (units of
! 10^{-12} for 10^4, 2 x 10^4, and 3 x 10^4 K.
!
	    T1=HDKT*EDGEDIE
	    A10=2.07D-10*GUPDIE*EINA/GION_GS
	    A20=A10*EXP(0.5D0*T1)/(2.0D0**1.5D0)
	    A30=A10*EXP(T1/3.0D0)/(3.0D0**1.5D0)
	    A10=A10*EXP(T1)
!
! Is this lower level in model atom ?
!
	    L2=INDEX(STRING,'-')-1
	    IF( L2 .LE. 0 .OR. L2 .GE. L1)THEN
	      WRITE(LUER,*)'Error reading in Level Names from '//FILENAME
	      WRITE(2,'(A)')STRING
	      WRITE(LUER,*)'L2=',L2,'L1=',L1
	      STOP
	    END IF
!
! Identify lower level of the transition. If the lower level is split,
! but the DIELECTRONI! rates are not, the photoionization cross-section
! (i.e oscillator strength) is taken to be independent of J. This is NOT
! CORRECT, but will not matter when the J levels are being combined as a
! single SUPER levels (normal case for CNO).
!
! NB: MNL_CNT refers to the number of model levels which the lower level
!               of the dielectronic transition was matched with.
!
	    MNL=0
	    MNL_CNT=0
!
	    SPLIT_J=.FALSE.
	    IF(INDEX(STRING(1:L2),'[') .NE. 0 .AND.
	1           INDEX(STRING(1:L2),']') .NE. 0)SPLIT_J=.TRUE.
!
! Search for level match only if transition to be included.
!
	    INDX_HASH=INDEX(STRING,'#')
	    IF( (DO_AUTO .AND. INDX_HASH .EQ. 0) .OR.
	1          (DO_WI .AND. INDX_HASH .NE. 0) )THEN
	      DO I=1,NXzV
	        IF(SPLIT_J)THEN
	          J=0
	        ELSE
	          J=INDEX(LEVELNAME(I),'[')-1
	        END IF
	        IF(J .LE. 0)J=ICHRLEN(LEVELNAME(I))
	        IF(STRING(1:L2) .EQ. LEVELNAME(I)(1:J))THEN
	           MNL=I
	           MNL_CNT=MNL_CNT+1
	           IF(MNL_CNT .NE. 1)THEN
	             CNT=CNT+1
	             IF(CNT .GT. PD(ID)%NDIE_MAX)THEN
	               WRITE(LUER,*)'NDIE_MAX not large enough in ',
	1                     'RD_XzV_PHOT_DIE: ',FILENAME
	               WRITE(LUER,*)'NDIE_MAX=',PD(ID)%NDIE_MAX,'NUM_DIE=',NUM_D_RD
	               STOP
	             END IF
	             PD(ID)%OSC(CNT)=PD(ID)%OSC(CNT-1)
	             IF(LS_NAME .NE. LEVELNAME(I)(1:J))THEN
	               WRITE(LUER,*)'Error in RD_XzV_PHOT_DIE'
	               WRITE(LUER,*)'Inconsistent lowere level names'
	               WRITE(LUER,*)'LS_NAME=',LS_NAME
	               WRITE(LUER,*)'LEVELNAME=',LEVELNAME(I)(1:J)
	               STOP
	             END IF
	           ELSE
	             LS_NAME=LEVELNAME(I)(1:J)
	           END IF
	           DIELEV(CNT)=MNL
	           PD(ID)%NU_ZERO(CNT)=EDGE(MNL)-EDGEDIE
!
! We already no this transition must be included. We simply set a different
! width depending on whether it is forbidden or allowed in pure LS coupling.
!
	           IF(INDX_HASH .EQ. 0)THEN
	             PD(ID)%GAMMA(CNT)=1.0D+13/4.0D0/FUN_PI()
	           ELSE IF(INDX_HASH .NE. 0)THEN
	             PD(ID)%GAMMA(CNT)=1.0D+12/4.0D0/FUN_PI()
	           END IF
	        END IF
	      END DO
	    END IF
!
	    IF(SPLIT_J .AND. MNL_CNT .GT. 1)THEN
	      WRITE(LUER,*)'Error in RD_XzV_PHOT_DIE'
	      WRITE(LUER,*)'Multiple levels matched with transition'
	      WRITE(LUER,*)TRANSDIE
	      STOP
	    END IF
!
! Check the total rates being included. If the transition was not included,
! we must subtract 1 from CNT, as CNT was updated on the READ.
!
	    IF(INDX_HASH .EQ. 0)THEN			!Not WI transition.
	      IF(DO_AUTO .AND. MNL_CNT .NE. 0)THEN
	        A1=A1+A10
	        A2=A2+A20
	        A3=A3+A30
	        INC=INC+1
	      ELSE
	        M1=M1+A10
	        M2=M2+A20
	        M3=M3+A30
	        CNT=CNT-1
	        MIS=MIS+1
	      END IF
	    ELSE
	      IF(DO_WI .AND. MNL_CNT .NE. 0)THEN
	        WIA1=WIA1+A10
	        WIA2=WIA2+A20
	        WIA3=WIA3+A30
	        WIINC=WIINC+1
	      ELSE
	        WIM1=WIM1+A10
	        WIM2=WIM2+A20
	        WIM3=WIM3+A30
	        CNT=CNT-1
	        WIMIS=WIMIS+1
	      END IF
	    END IF
	  END DO
!
	  IF( (INC+WIINC+MIS+WIMIS) .NE. NUM_D_RD)THEN
	     WRITE(LUER,*)'Error in RDGENDIE -'//FILENAME
	     WRITE(LUER,*)'Invalid summation of included transitions'
	     STOP
	  END IF
!
	  WRITE(LUOUT,900)
900	  FORMAT(/,' Summary of dielectronic transitions ',
	1                  'included (LS : WI) ',/,
	1 '[Units 10^-12] ( ) denotes percentage of LTDR NOT included')
	  IF( (A1+M1) .NE. 0)THEN
	    M1=100.0D0*M1/(A1+M1)
	    M2=100.0D0*M2/(A2+M2)
	    M3=100.0D0*M3/(A3+M3)
	    WRITE(LUOUT,1000)INC,(INC+MIS),A1,M1,A2,M2,A3,M3
1000	    FORMAT( X,I4,'(',I4,')',3( 2X,1PE10.3,'(',0PF6.2,')' )  )
	  END IF
!
	  IF( (WIA1+WIM1) .NE. 0)THEN
	    WIM1=100.0*WIM1/(WIA1+WIM1)
	    WIM2=100.0*WIM2/(WIA2+WIM2)
	    WIM3=100.0*WIM3/(WIA3+WIM3)
	    WRITE(LUOUT,1000)WIINC,(WIINC+WIMIS),WIA1,WIM1,WIA2,WIM2,
	1                  WIA3,WIM3
	  END IF
	END IF
!
!
!
! Close file, as all required data has been read in.
!
	CLOSE(UNIT=LUIN)
!
! Sort the dielectronic transitions according to their lower level (stored
! in DIELEV).
!
	PD(ID)%NUM_DIE=CNT			!Required in PHOT_XzV (Mod variable)
	CALL INDEXX(CNT,DIELEV,VEC_INDX,.TRUE.)
	CALL SORTDP(CNT,DIELEV,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(CNT,PD(ID)%OSC,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(CNT,PD(ID)%GAMMA,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(CNT,PD(ID)%NU_ZERO,VEC_INDX,VEC_DP_WRK)
!
! Now count the nuber of dielectronic transitions associated with each
! lower level, and store the start and end indexs for transitions to a
! particular level.
!
! Setting ST_INDEX=1 and END_INDEX=0 for all levels ensures that no dieletronic
! transitions are handled for a particular level unless some are available.
!
	PD(ID)%ST_INDEX(:)=1
	PD(ID)%END_INDEX(:)=0
	J=1
	DO I=1,NXzV
	  DO WHILE(J .LE. CNT .AND. NINT(DIELEV(MIN(J,CNT))) .EQ. I)
	    IF(PD(ID)%END_INDEX(I) .EQ. 0)THEN
	      PD(ID)%ST_INDEX(I)=J
	      PD(ID)%END_INDEX(I)=J
	    ELSE
	      PD(ID)%END_INDEX(I)=PD(ID)%END_INDEX(I)+1
	    END IF
	    J=J+1
          END DO
	END DO
!
! Throughout the program, Frequency is in units of 10^15 Hz. We therefore
! put GAMMA in the same units.
!
	DO I=1,PD(ID)%NUM_DIE
	  PD(ID)%GAMMA(I)=1.0D-15*PD(ID)%GAMMA(I)
	END DO
!
! Compute the frequency range associated with each transition. This
! range is chosen to give 0.1% accuracy in the integral across the full
! VOIGT profile.
!
	T1=VSMOOTH_KMS/2.998D+05
	DO I=1,PD(ID)%NUM_DIE
	  NU_DOP=PD(ID)%NU_ZERO(I)*T1
	  DEL_NU=64.0D0*PD(ID)%GAMMA(I)
	  IF(DEL_NU .LT. 6.0D0*NU_DOP)DEL_NU=6.0D0*NU_DOP
	  PD(ID)%NU_MIN(I)=PD(ID)%NU_ZERO(I)-DEL_NU
	  PD(ID)%NU_MIN(I)=MAX(PD(ID)%NU_MIN(I),EDGE( NINT(DIELEV(I)) ))
	  PD(ID)%NU_MAX(I)=PD(ID)%NU_ZERO(I)+DEL_NU
	END DO
!
	IF(NUM_FREE_FREE .NE. 0)THEN
!
! Throughout the program, Frequency is in units of 10^15 Hz. We therefore
! put GAMMA in the same units.
!
	  DO I=1,PD(ID)%NUM_FF
	    PD(ID)%FF_GAMMA(I)=1.0D-15*PD(ID)%FF_GAMMA(I)
	  END DO
	WRITE(LUER,*)'Done FF units'
!
! Compute the frequency range associated with each transition. This
! range is chosen to give 1.0% accuracy in the integral across the full
! (pure) VOIGT profile.
!
	  T1=VSMOOTH_KMS/2.998D+05
	  DO I=1,PD(ID)%NUM_FF
	    NU_DOP=PD(ID)%FF_NU_ZERO(I)*T1
	    DEL_NU=64.0D0*PD(ID)%FF_GAMMA(I)
	    IF(DEL_NU .LT. 6*NU_DOP)DEL_NU=6.0D0*NU_DOP
	    PD(ID)%FF_NU_MIN(I)=PD(ID)%FF_NU_ZERO(I)-DEL_NU
	    PD(ID)%FF_NU_MAX(I)=PD(ID)%FF_NU_ZERO(I)+DEL_NU
	    WRITE(LUER,*)PD(ID)%FF_GF(I),PD(ID)%FF_NU_ZERO(I),PD(ID)%FF_NU_EXCITE(I)
	    WRITE(LUER,*)PD(ID)%FF_NU_MIN(I),PD(ID)%FF_NU_MAX(I)
	  END DO
	  WRITE(LUER,*)'Done FF limits units'
	END IF
!
! Clean up allocated memory.
!
	IF(ALLOCATED(DIE_LEV_G))DEALLOCATE(DIE_LEV_G)
	IF(ALLOCATED(DIE_LEV_NAME))DEALLOCATE(DIE_LEV_NAME)
	IF(ALLOCATED(DIE_LEV_ENERGY))DEALLOCATE(DIE_LEV_ENERGY)
	IF(ALLOCATED(DIE_LEV_AUTO))DEALLOCATE(DIE_LEV_AUTO)
!
	RETURN
	END
