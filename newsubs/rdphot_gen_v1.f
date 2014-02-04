!
! General program to read in the photoionization cross-sections for XzV.
! Dielectronic transitions can also be read in if necessary, and they
! are then treated with the photoionization cross-sections.
!
! NB: EDGE, XzV_LEVELNAME, GION_GS, NXzV refer to the FULL atom.
!
	SUBROUTINE RDPHOT_GEN_V1(EDGE,XzV_LEVELNAME,
	1             GION_GS,AT_NO_DUM,ZXzV,NXzV,
	1             XzV_ION_LEV_ID,N_PHOT,MAX_N_PHOT,
	1             XzSIX_PRES,EDGEXzSIX_F,GXzSIX_F,F_TO_S_XzSIX,
	1             XzSIX_LEVNAME_F,NXzSIX_F,
	1             X_RAYS,ID,DESC,
	1             LUIN,LUOUT)
!
!  Data modules required.
!
	USE PHOT_DATA_MOD		!Contains all photoionization data
!
	IMPLICIT NONE
!
! Altered 28-Mar-2003 : Check ensuring that ID < NPD_MAX installed.
! Altered 03-Jul-2001 : Only compute NEF for levels below ionization limit.
!                         May need work.
! Altered 16-Dec-1999 : FILENAME nolonger set to upper case [for UNIX].
! Altered 01-Dec-1999 : Altered to handle alternate state names in PHOT file.
!                          Done becuase there is still some inconsistency in
!                          naming conventions. Alternate names are seperated by
!                          a "/". Spaces may appear before the slash. This
!                          was done to allow earlier versions of CMFGEN to read
!                          the same data file.
! Altered 12-Dec-1997 : Bug fix. DO_PHOT was not being properly set for
!                          B level photoionizations when B level not present.
!                          In this case require DO_PHOT(1,PHOT_ID)=.TRUE.
!                          [not PHOT_ID,1] Thus when computing level 1 cross-
!                          section, we include the cross-section for PHOT_ID.
!
! Altered  5-Dec-1996 : GEN_ASCI_OPEN used for opening ASCI files.
! Altered Sep-18-1996 : Dimension of LST_U, LST_LOC and LST_CROSS changed.
!                       Name of LST_U changed to LST_FREQ.
!
!
! Passed variables.
!
	INTEGER NXzV
	INTEGER N_PHOT
	INTEGER MAX_N_PHOT
	REAL*8 EDGE(NXzV)
	CHARACTER*30 XzV_LEVELNAME(NXzV)
!
! Statistical weight of resulting ion for PHOT_ID=1. Used only if
! ion is not present.
!
	REAL*8 GION_GS
	REAL*8 ZXzV
	REAL*8 AT_NO_DUM
	INTEGER XzV_ION_LEV_ID(MAX_N_PHOT)
	INTEGER ID
	CHARACTER*(*) DESC
!
	LOGICAL XzSIX_PRES
	INTEGER NXzSIX_F
	REAL*8 EDGEXzSIX_F(NXzSIX_F)
	REAL*8 GXzSIX_F(NXzSIX_F)
	INTEGER F_TO_S_XzSIX(NXzSIX_F)
	CHARACTER*(*) XzSIX_LEVNAME_F(NXzSIX_F)
!
	LOGICAL X_RAYS
	INTEGER LUIN,LUOUT
!
! Common block with opacity/emissivity constants.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
! External Functions.
!
	REAL*8 RD_FREE_VAL,SPEED_OF_LIGHT
	INTEGER ICHRLEN,ERROR_LU
	EXTERNAL RD_FREE_VAL,ERROR_LU,ICHRLEN,SPEED_OF_LIGHT
!
! Local variables.
!
	INTEGER, PARAMETER :: IZERO=0
!
	INTEGER NTERMS(MAX_N_PHOT)
	INTEGER N_DP(MAX_N_PHOT)
	CHARACTER*30, ALLOCATABLE :: LOC_TERM_NAME(:)
!
	REAL*8 LOC_EXC_FREQ_ION,TOTAL_WT
	REAL*8 T1,NU_INF
	CHARACTER*20 CROSS_UNIT
	INTEGER PHOT_ID,IOS
	INTEGER I,J,K
	INTEGER N_HD,L1,L2,LN
	INTEGER NEXT,LUER
	INTEGER END_CROSS
	INTEGER NPNTS
!
	LOGICAL SPLIT_J
!
	INTEGER MAX_HEAD
	PARAMETER (MAX_HEAD=15)
	CHARACTER*132 HEAD_STR(MAX_HEAD)
	CHARACTER*1 APPEND_CHAR
!
	CHARACTER*132 STRING
	CHARACTER*80 FILENAME
!
	INTEGER, PARAMETER :: N_FIN_MAX=5			!Up to 5 alternate names
	INTEGER N_FIN_STATES
	CHARACTER*30 FINAL_STATE(MAX_N_PHOT,N_FIN_MAX)
	LOGICAL FOUND_FINAL_STATE
!
	LUER=ERROR_LU()
	FINAL_STATE(:,:)=' '
!
! 
! Open the data files to determine the number of energy levels and the
! number of data pairs. This is done for all files belonging to the current
! species so that we can allocate the necessary memory.
!
	PHOT_ID=0
	N_PHOT=1
	DO WHILE (PHOT_ID .LT. N_PHOT)
	  PHOT_ID=PHOT_ID+1
!
	  FILENAME='PHOT'//TRIM(DESC)//'_'
	  I=ICHRLEN(FILENAME)
	  APPEND_CHAR=CHAR( ICHAR('A')+(PHOT_ID-1) )
	  WRITE(FILENAME(I+1:I+1),'(A)')APPEND_CHAR
!	  CALL SET_CASE_UP(FILENAME,IZERO,IZERO)
	  CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ','READ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Unable to open file',FILENAME
	    STOP
	  END IF
!
! Read in number of energy levels.
!
	  IOS=0
          L1=0.0
	  DO WHILE(L1 .EQ. 0 .AND. IOS .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    L1=INDEX(STRING,'!Number of energy levels')
	  END DO
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Number of energy levels string not found'
	    STOP
	  ELSE
	    READ(STRING,*,IOSTAT=IOS)NTERMS(PHOT_ID)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)' Error reading NTERMS in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'PHOT_ID=',PHOT_ID
	      STOP
	    END IF
	  END IF
!
! Read in number of photionization routes.
!
	  IF(PHOT_ID .EQ. 1)THEN
	    IOS=0
	    L1=0
	    DO WHILE(L1 .EQ. 0 .AND. IOS .EQ. 0)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      L1=INDEX(STRING,'!Number of photoionization routes')
	    END DO
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'Number of photoionization routes string not found'
	      STOP
	    ELSE
	      READ(STRING,*,IOSTAT=IOS)N_PHOT
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)' Error reading N_PHOT in RDPHOT_GEN_V1: ',DESC
	        WRITE(LUER,*)'PHOT_ID=',PHOT_ID
	        STOP
	      END IF
	      IF(N_PHOT .GT. MAX_N_PHOT)THEN
	        WRITE(LUER,*)' Error in RDPHOT_GEN_V1: ',DESC,
	1                       ' --- MAX_N_PHOT too small'
	        WRITE(LUER,*)' MaX_N_PHOT=',MAX_N_PHOT
	        WRITE(LUER,*)' N_PHOT=',N_PHOT
	        STOP
	      END IF
	    END IF
	  END IF
!
! Read in number of data pairs.
!
	  L1=0
	  IOS=0
	  DO WHILE(L1 .EQ. 0 .AND. IOS .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    L1=INDEX(STRING,'!Total number of data pairs')
	  END DO
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Total number of data pairs not found'
	    STOP
	  ELSE
	    READ(STRING,*,IOSTAT=IOS)N_DP(PHOT_ID)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)' Error reading number of data ',
	1                      ' pairs in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'PHOT_ID=',PHOT_ID
	      STOP
	    END IF
	  END IF
	  CLOSE(UNIT=LUIN)
	END DO
!
! Allocate arrays and vectors whose size is determined by the passed
! parameters.
!
	IF(ID .GT. NPD_MAX)THEN
	  WRITE(LUER,*)'Error in RDPHOT_GEN_V1: NPD_MAX too small'
	  WRITE(LUER,*)'ID=',ID,'NPD_MAX=',NPD_MAX
	  STOP
	END IF
	ALLOCATE (PD(ID)%NEF(NXzV,N_PHOT))
	ALLOCATE (PD(ID)%EXC_FREQ(N_PHOT))
	ALLOCATE (PD(ID)%GION(N_PHOT))
	ALLOCATE (PD(ID)%A_ID(NXzV,N_PHOT))
	ALLOCATE (PD(ID)%DO_PHOT(N_PHOT,N_PHOT))
!
! Now perform the dynamic allocation of memory for variables determined by the
! size of the input data.
!
	PD(ID)%MAX_CROSS=0
	PD(ID)%MAX_TERMS=0
	DO PHOT_ID=1,N_PHOT
	  PD(ID)%MAX_CROSS=PD(ID)%MAX_CROSS+N_DP(PHOT_ID)
	  PD(ID)%MAX_TERMS=MAX(PD(ID)%MAX_TERMS,NTERMS(PHOT_ID))
	END DO
!
	ALLOCATE (PD(ID)%NU_NORM(PD(ID)%MAX_CROSS))
	ALLOCATE (PD(ID)%CROSS_A(PD(ID)%MAX_CROSS))
	ALLOCATE (PD(ID)%ST_LOC(PD(ID)%MAX_TERMS,N_PHOT))
	ALLOCATE (PD(ID)%END_LOC(PD(ID)%MAX_TERMS,N_PHOT))
	ALLOCATE (PD(ID)%CROSS_TYPE(PD(ID)%MAX_TERMS,N_PHOT))
!
	ALLOCATE (LOC_TERM_NAME(PD(ID)%MAX_TERMS))
	ALLOCATE (PD(ID)%LST_CROSS(NXzV,N_PHOT))
	ALLOCATE (PD(ID)%LST_FREQ(NXzV,N_PHOT))
	ALLOCATE (PD(ID)%LST_LOC(NXzV,N_PHOT))
!
! Perform default initializations.
!
	PD(ID)%DO_PHOT(:,:)=.FALSE.
	PD(ID)%EXC_FREQ(:)=0.0D0
	XzV_ION_LEV_ID(:)=0
	PD(ID)%A_ID(:,:)=0
!
! Set values for ionizations to the XzSIX ground state.
!
	PD(ID)%DO_PHOT(1,1)=.TRUE.
	PD(ID)%EXC_FREQ(1)=0.0D0
	XzV_ION_LEV_ID(1)=1
	PD(ID)%AT_NO=AT_NO_DUM
	PD(ID)%ZION=ZXzV
!
! END_CROSS provides an indication of the last storage location used. It is
! initialized here, since all photoionization routes use the same storage
! locations.
!
	END_CROSS=0
!
! If the XzSIX ground state is not present (hence XzSIX_PRES is .FALSE.) we
! include the X-RAY opacity with the regular photoionization cross-section.
!
	IF(X_RAYS .AND. .NOT. XzSIX_PRES)THEN
	  PD(ID)%DO_KSHELL_W_GS=.TRUE.
	END IF
!
! 
! This section reads and stores the cross sections, looping over each final
! state. The final states should have logical name (files)
!
!               PHOTXzV_n where n=1 is the ground state;
!                          where n=2, 3 etc denote excited states.
!
	DO PHOT_ID=1,N_PHOT
!
	  WRITE(LUOUT,'(A)')'  '
!
! Open file containing cross sections. We first read in all header info.
!
	  FILENAME='PHOT'//TRIM(DESC)//'_'
	  I=ICHRLEN(FILENAME)
	  APPEND_CHAR=CHAR( ICHAR('A')+(PHOT_ID-1) )
	  WRITE(FILENAME(I+1:I+1),'(A)')APPEND_CHAR
!	  CALL SET_CASE_UP(FILENAME,IZERO,IZERO)
	  CALL GEN_ASCI_OPEN(LUIN,FILENAME,'OLD',' ','READ',IZERO,IOS)
!
! Find the record with the date the file was written. This begins the
! important information. NB: After the date the header records can be
! in any order. There must be NO BLANK lines however, as this signifies
! the start of the data.
!
	  L1=0
	  DO WHILE(L1 .EQ. 0)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    L1=INDEX(STRING,'!Date')
	    IF(IOS .NE. 0)THEN
 
	      WRITE(LUER,*)'Error reading date from ',FILENAME
	      WRITE(LUER,*)'PHOT_ID=',PHOT_ID,'IOSTAT=',IOS
	      STOP
	    END IF
	    WRITE(LUOUT,'(A)')STRING(1:ICHRLEN(STRING))
	  END DO
!
	  N_HD=1
	  READ(LUIN,'(A)')HEAD_STR(N_HD)
	  WRITE(LUOUT,'(A)')HEAD_STR(N_HD)
	  DO WHILE(HEAD_STR(N_HD) .NE. '  ')
	    N_HD=N_HD+1
	    IF(N_HD .GT. MAX_HEAD)THEN
	      WRITE(LUER,*)'Insufficient header space in RDPHOT_GEN_V1: ',DESC
	      STOP
	    END IF
	    READ(LUIN,'(A)')HEAD_STR(N_HD)
	    WRITE(LUOUT,'(A)')HEAD_STR(N_HD)(1:ICHRLEN(HEAD_STR(N_HD)))
	  END DO
!
! Read in number of energy levels.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Number of energy levels')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Number of energy levels string not found'
	    STOP
	  ELSE
	    READ(HEAD_STR(J),*,IOSTAT=IOS)NTERMS(PHOT_ID)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)' Error reading NTERMS in RDPHOT_GEN_V1 - 2:',DESC
	      WRITE(LUER,*)' PHOT_ID=',PHOT_ID
	      STOP
	    END IF
	  END IF
!
! Read in units for cross-sections.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Cross-section unit')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Unit for cross-sections not found'
	    STOP
	  ELSE
	    L1=INDEX(HEAD_STR(J),'  ')
	    CROSS_UNIT=HEAD_STR(J)(1:L1-1)
	    IF(CROSS_UNIT .NE. 'Megabarns')THEN
	      WRITE(LUER,*)' Error in RDPHOT_GEN_V1 - incorrect',
	1                        ' cross-section unit:',DESC
	      WRITE(LUER,*)'CROSS _UNIT=',CROSS_UNIT
	      WRITE(LUER,*)'UNIT should be Megabarns'
	      STOP
	    END IF
	  END IF
!
! Are photoioization cross-sections split for individual J states.
! This is done to faciliated determining level correspondance.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Split J levels')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Split J levels indication not FOUND!'
	    STOP
	  ELSE
	    IOS=0
	    L1=INDEX(HEAD_STR(J),'  ')
	    READ(HEAD_STR(J)(1:L1-1),*,IOSTAT=IOS)SPLIT_J
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'Can''t read SPLIT_J value'
	    END IF
	  END IF
!
! Now read in the final state. We allow for the possibility
! of alternate names, separetd by a '/'. No double spaces
! must be present in the names.
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Final state in ion')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Final ion state string not found.'
	    STOP
	  ELSE
!
! We can extract the Final level. To do so we frist remove
! unwanted garbage, TABS, and excess blanks.
!
	    HEAD_STR(J)(L1:)=' '
	    DO WHILE(INDEX(HEAD_STR(J),CHAR(9)) .NE. 0)
	      K=INDEX(HEAD_STR(J),CHAR(9))	!Remove tabs 
	      HEAD_STR(J)(K:K)=' '
	    END DO
!
	    L1=LEN_TRIM(HEAD_STR(J))
	    K=INDEX(HEAD_STR(J)(1:L1),' ')	!Remove blanks
	    DO WHILE(K .NE. 0) 
	      HEAD_STR(J)(K:)=HEAD_STR(J)(K+1:)
	      L1=L1-1
	      K=INDEX(HEAD_STR(J)(1:L1),' ')
	    END DO
!
! Can now get final level name
!	        
	    L2=INDEX(HEAD_STR(J),'  ')
	    L1=0
	    DO WHILE(L1 .LT. N_FIN_MAX)
	      L1=L1+1
	      K=INDEX(HEAD_STR(J)(1:L2),'/')-1
	      IF(K .LE. 0)K=L2
	      FINAL_STATE(PHOT_ID,L1)=HEAD_STR(J)(1:K)
	      IF(K .EQ. L2)EXIT
	      HEAD_STR(J)(1:)=HEAD_STR(J)(K+2:)
	      L2=INDEX(HEAD_STR(J),'  ')
	    END DO
	    N_FIN_STATES=L1
	  END IF
!
! Now read in the excitation energy of the final state above the
! ground state. Excitation energy should be in cm^-1
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Excitation energy of final state')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Excitation energy not found.'
	    STOP
	  ELSE
	    L2=INDEX(HEAD_STR(J),'  ')
	    LOC_EXC_FREQ_ION=RD_FREE_VAL(HEAD_STR(J),1,L2-1,NEXT,
	1                     ' EXC_FREQ RDPHOT_GEN_V1')
	    LOC_EXC_FREQ_ION=1.0D-15*LOC_EXC_FREQ_ION*SPEED_OF_LIGHT()
	  END IF
!
	  J=0
	  L1=0
	  DO WHILE(L1 .EQ. 0 .AND. J .LT. N_HD)
	    J=J+1
	    L1=INDEX(HEAD_STR(J),'!Statistical weight of ion')
	  END DO
	  IF(L1 .EQ. 0)THEN
	    WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	    WRITE(LUER,*)'Statistical weight of ion g.s. not found'
	    STOP
	  ELSE
	    L2=INDEX(HEAD_STR(J),'  ')
	    T1=RD_FREE_VAL(HEAD_STR(J),1,L2-1,NEXT,
	1                   ' EXC_FREQ RDPHOT_GEN_V1')
	    PD(ID)%GION(PHOT_ID)=T1
	    IF(PHOT_ID .EQ. 1)GION_GS=T1
	  END IF
! 
!
! We now determine and initialize variables so that the final state
! is correctly linked. We only check the name of the final state if that
! species is present.
!
	  IF(PHOT_ID .EQ. 1)THEN
	    IF(XzSIX_PRES)THEN
	      LN=INDEX(XzSIX_LEVNAME_F(1),'[')-1
	      IF(LN .LT. 0)LN=ICHRLEN(XzSIX_LEVNAME_F(1))
	      FOUND_FINAL_STATE=.FALSE.
	      DO L1=1,N_FIN_STATES
	        IF(XzSIX_LEVNAME_F(1)(1:LN) .EQ. FINAL_STATE(PHOT_ID,L1))THEN
	          FOUND_FINAL_STATE=.TRUE.
	          STRING=FINAL_STATE(PHOT_ID,1)		!Switch names for later output
	          FINAL_STATE(PHOT_ID,1)=FINAL_STATE(PHOT_ID,L1)
	          FINAL_STATE(PHOT_ID,L1)=STRING
	          EXIT
	        END IF
	      END DO
	      IF(.NOT. FOUND_FINAL_STATE)THEN
	        WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	        WRITE(LUER,*)'Invalid final state PHOT_ID=',PHOT_ID
	        WRITE(LUER,*)'XzSIX_LEVNAME(1)=',XzSIX_LEVNAME_F(1)
	        WRITE(LUER,*)FINAL_STATE
	        STOP
	      END IF
	    END IF
	  ELSE
!
! Determine the correspondence of the B_level with the super levels.
! It is assumed that PHOT_ID corresponds to the fround state.
!
	    PD(ID)%EXC_FREQ(PHOT_ID)=0.0D0
	    TOTAL_WT=0.0D0
	    IF(XzSIX_PRES .AND. NXzSIX_F .GE. 2)THEN
	      DO I=1,NXzSIX_F
	        LN=INDEX(XzSIX_LEVNAME_F(I),'[')-1
	        IF(LN .LT. 0)LN=ICHRLEN(XzSIX_LEVNAME_F(I))
	        FOUND_FINAL_STATE=.FALSE.
	        DO L1=1,N_FIN_STATES
	          IF(XzSIX_LEVNAME_F(I)(1:LN) .EQ. FINAL_STATE(PHOT_ID,L1))THEN
	            FOUND_FINAL_STATE=.TRUE.
	            STRING=FINAL_STATE(PHOT_ID,1)		!Switch names for later output
	            FINAL_STATE(PHOT_ID,1)=FINAL_STATE(PHOT_ID,L1)
	            FINAL_STATE(PHOT_ID,L1)=STRING
	            EXIT
	          END IF
	        END DO
	        IF(FOUND_FINAL_STATE)THEN
	          IF(XzV_ION_LEV_ID(PHOT_ID) .EQ. 0)THEN
	            XzV_ION_LEV_ID(PHOT_ID)=F_TO_S_XzSIX(I)
	          ELSE IF(XzV_ION_LEV_ID(PHOT_ID) .NE. F_TO_S_XzSIX(I))THEN
	            WRITE(2,*)'Warning in RDPHOT_GEN_V1: ',DESC
	            WRITE(LUER,*)'Super levels of final states do not match'
!
! Need to decide how to treat this better
!	            STOP
!
	          END IF
	          PD(ID)%EXC_FREQ(PHOT_ID)=PD(ID)%EXC_FREQ(PHOT_ID) +
	1                               GXzSIX_F(I)*EDGEXzSIX_F(I)
	          TOTAL_WT=TOTAL_WT+GXzSIX_F(I)
	          PD(ID)%DO_PHOT(PHOT_ID,PHOT_ID)=.TRUE.
	        END IF
	      END DO
	      IF(XzV_ION_LEV_ID(PHOT_ID) .NE. 0)THEN
	        PD(ID)%EXC_FREQ(PHOT_ID)=EDGEXzSIX_F(1)-PD(ID)%EXC_FREQ(PHOT_ID)/TOTAL_WT
	      ELSE
	        PD(ID)%DO_PHOT(1,PHOT_ID)=.TRUE.
	        PD(ID)%EXC_FREQ(PHOT_ID)=LOC_EXC_FREQ_ION
	      END IF
	    ELSE
	      PD(ID)%DO_PHOT(1,PHOT_ID)=.TRUE.
	      PD(ID)%EXC_FREQ(PHOT_ID)=LOC_EXC_FREQ_ION
	    END IF
	  END IF
!
! 
!
! Read in the cross-sections
!
	  DO J=1,NTERMS(PHOT_ID)
	    STRING=' '
	    DO WHILE( INDEX(STRING,'Configuration name') .EQ. 0)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error reading in configuration name from ',
	1                       FILENAME
	        WRITE(LUER,*)'Currently reading term',J
	        STOP
	      END IF
	    END DO
	    L1=INDEX(STRING,'!Configuration')-1
	    L1=ICHRLEN(STRING(1:L1))
	    LOC_TERM_NAME(J)=STRING(1:L1)
	    DO WHILE(LOC_TERM_NAME(J)(1:1) .EQ. ' ')
	      LOC_TERM_NAME(J)=LOC_TERM_NAME(J)(2:)
	    END DO
!
	    STRING=' '
	    DO WHILE( INDEX(STRING,'Type of cross-section') .EQ. 0)
	      READ(LUIN,'(A)',IOSTAT=IOS)STRING
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error reading Type of cross-section from ',
	1         FILENAME
	        WRITE(LUER,*)'Currently reading term',J
	        STOP
	      END IF
	    END DO
	    READ(STRING,*)PD(ID)%CROSS_TYPE(J,PHOT_ID)
!
	    READ(LUIN,*)NPNTS
	    PD(ID)%ST_LOC(J,PHOT_ID)=END_CROSS+1
	    PD(ID)%END_LOC(J,PHOT_ID)=END_CROSS+NPNTS
	    END_CROSS=PD(ID)%END_LOC(J,PHOT_ID)
	    IF(PD(ID)%END_LOC(J,PHOT_ID) .GT. PD(ID)%MAX_CROSS)THEN
	      WRITE(LUER,*)'Error in RDPHOT_GEN_V1: ',DESC
	      WRITE(LUER,*)'MAX_CROSS is too small'
	      STOP
	    END IF
!
! We can either read in a tabulated data set, or parameters for
! fitting a function. For a tabulated data set there must be the frequency
! and cross-section on each line --- for the other cases one parameter
! is specified per line. For compatibility with the tabular form, a zero
! must be written after each parameter. NB: In the second case PD(ID)%NU_NORM is
! not a frequency, but rather a parameter.
!
	    IF(PD(ID)%CROSS_TYPE(J,PHOT_ID) .LT. 20)THEN
	      READ(LUIN,*)(PD(ID)%CROSS_A(K),
	1                    K=PD(ID)%ST_LOC(J,PHOT_ID),PD(ID)%END_LOC(J,PHOT_ID))
	    ELSE
	      READ(LUIN,*)(PD(ID)%NU_NORM(K),PD(ID)%CROSS_A(K),
	1                    K=PD(ID)%ST_LOC(J,PHOT_ID),PD(ID)%END_LOC(J,PHOT_ID))
	    END IF
	  END DO
!
! Match and identify input cross sections to levels in code.
!
! Match each of the TERMS with the corresponding levels in the full atom.
! We ensure that each level has a photoionization cross-section.
!
! To get the Rydberg constant we assume that the atomic mass in AMU is just
! twice the atomic number. Only exeption is H.
!
	  NU_INF=1.0D-15*109737.31*SPEED_OF_LIGHT()
	  IF(PD(ID)%AT_NO .EQ. 1)THEN
	    NU_INF=NU_INF/(1+5.48597D-04)
	  ELSE
	    NU_INF=NU_INF/(1+5.48597D-04/(2*PD(ID)%AT_NO))
	  END IF
	  DO I=1,NXzV
	    IF(SPLIT_J)THEN
	      L1=ICHRLEN(XzV_LEVELNAME(I))
	    ELSE
	      L1=INDEX(XzV_LEVELNAME(I),'[')-1
	      IF(L1 .LE. 0)L1=ICHRLEN(XzV_LEVELNAME(I))
	    END IF
	    DO J=1,NTERMS(PHOT_ID)
	      IF(XzV_LEVELNAME(I)(1:L1) .EQ. LOC_TERM_NAME(J))THEN
	        IF(PD(ID)%A_ID(I,PHOT_ID) .EQ. 0)THEN
	          PD(ID)%A_ID(I,PHOT_ID)=J
	        ELSE
	          WRITE(LUER,*)'Error in handling photoionization ',
	1                'cross-sections for ',FILENAME
	          WRITE(LUER,*)'Non-unique cross-section for ',
	1                        XzV_LEVELNAME(I)
	          STOP
	        END IF
	      END IF
	    END DO
	    IF(PD(ID)%A_ID(I,PHOT_ID) .EQ. 0)THEN
	      WRITE(LUER,*)'Error in handling photoionization ',
	1                   'cross-sections for ',FILENAME
	      WRITE(LUER,*)'Cross section unavailable for level ',
	1                      XzV_LEVELNAME(I)
	      STOP
	    END IF
	    T1=NU_INF/(EDGE(I)+PD(ID)%EXC_FREQ(PHOT_ID))
	    PD(ID)%NEF(I,PHOT_ID)=0.0D0
	    IF(T1 .GT. 0)PD(ID)%NEF(I,PHOT_ID)=ZXzV*DSQRT(T1)
	  END DO
!
	  CLOSE(UNIT=LUIN)
	END DO			!PHOT_ID
!
! ALPHA_BF is the constant needed to evaluate the photoionization cross-section
! for a hydrogenic ion. It is in the correct units for R in units of
! 10^10 cm, and CHI.R dimensionless.
!
	PD(ID)%ALPHA_BF=2.815D-06*(PD(ID)%ZION)**4
!
! Output summary of photoionization routes.
!
	WRITE(LUOUT,*)'Summary of Photoionization routes for ',TRIM(DESC)
	WRITE(LUOUT,'(1X,A,T40,(2X,I2))')
	1            'Final State/Phot ID',(I,I=1,N_PHOT)
	DO PHOT_ID=1,N_PHOT
	  WRITE(LUOUT,'(1X,A,T40,(3X,L1))')
	1    FINAL_STATE(PHOT_ID,1),(PD(ID)%DO_PHOT(PHOT_ID,I),I=1,N_PHOT)
	END DO
!
! PD(ID)%NUM_PHOT_ROUTES is used to refer to the actual number of photoionization
! routes available. If the final state is not available, they go to the
! ground state.
!
	PD(ID)%NUM_PHOT_ROUTES=N_PHOT
!
! Set N_PHOT to the actual number of final states that are being used.
! This is returned in the call.
!
	N_PHOT=0
	DO I=1,PD(ID)%NUM_PHOT_ROUTES
	  IF(PD(ID)%DO_PHOT(I,I))N_PHOT=N_PHOT+1
	END DO
!
	DEALLOCATE (LOC_TERM_NAME)
!
! Initialize storage locations for use by PHOT_XzV. These are used to store
! photoionization cross-sections to speed up computation.
!
	PD(ID)%LST_CROSS(:,:)=0.D0
	PD(ID)%LST_FREQ(:,:)=0.D0
	DO J=1,PD(ID)%NUM_PHOT_ROUTES		!
	  DO I=1,NXzV			!Loop over levels
	    K=PD(ID)%A_ID(I,J)			!Get pointer to cross-section
	    PD(ID)%LST_LOC(I,J)=(PD(ID)%ST_LOC(K,J) + PD(ID)%END_LOC(K,J))/2
	  END DO
	END DO
!
	RETURN
	END
