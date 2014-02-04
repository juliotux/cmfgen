!
! Routine to plot OPACITY photoionization cross-sections. Data from 2 distinct
! routines may be plotted. Presently assumed level ordering is the same in
! both files. 
!
	PROGRAM PLT_PHOT_RAW
	USE GEN_IN_INTERFACE
	USE HYD_BF_PHOT_DATA
	USE MOD_COLOR_PEN_DEF
	IMPLICIT NONE
!
! Altered 20-Sep-2012 : Bug fix with incorrect link test for 2nd photoionization data set.
!                         Cleaned.
! Altered 16-Apr-2008 : Read in energy levels form oscilator file (if available).
!                       Deleted log plot section (since can be done in pgplot)
!                       Can now eneter level index for name. If unrecognized level,
!                          all names are dumped to terminal.
! Altered 15-Jan-2008 : TYPE_1 etc changed to allocatable arrays.
!                       Errors returned if files cannot be opened.
! Created 09-Jun-1999
!
	INTEGER NPAIRS_1
	INTEGER NLEV_1
	INTEGER, ALLOCATABLE :: TYPE_1(:)
	INTEGER, ALLOCATABLE ::  NUM_VALS_1(:)
	INTEGER, ALLOCATABLE ::  LOC_1(:)
	REAL*8, ALLOCATABLE :: NU_1(:)
	REAL*8, ALLOCATABLE :: CROSS_1(:)
	REAL*8, ALLOCATABLE :: ENERGY_1(:)
	REAL*8, ALLOCATABLE :: STAT_WT_1(:)
!
	CHARACTER(LEN=40), ALLOCATABLE :: NAME_1(:)
	REAL*8 ZION_1
	REAL*8 GION_1 
	REAL*8 EXC_EN_1
	REAL*8 AMASS
	LOGICAL SPLITJ_1
!
! Storage for 2nd photoionization data set.
!
	INTEGER NPAIRS_2
	INTEGER NLEV_2
	INTEGER, ALLOCATABLE :: TYPE_2(:)
	INTEGER, ALLOCATABLE ::  NUM_VALS_2(:)
	INTEGER, ALLOCATABLE ::  LOC_2(:)
	CHARACTER(LEN=40), ALLOCATABLE :: NAME_2(:)
	REAL*8, ALLOCATABLE :: NU_2(:)
	REAL*8, ALLOCATABLE :: CROSS_2(:)
	REAL*8 ZION_2
	REAL*8 GION_2 
	REAL*8 EXC_EN_2
	LOGICAL SPLITJ_2
!
! These are used to read the Oscillator file.
!
	INTEGER NELEV
	REAL*8, ALLOCATABLE :: ENERGY(:)
	REAL*8, ALLOCATABLE :: G(:)
	REAL*8 IONIZATION_ENERGY
	CHARACTER(LEN=40), ALLOCATABLE :: E_NAME(:)
!
	CHARACTER*40 LEVEL_NAME1
	CHARACTER*40 LEVEL_NAME2
!
	REAL*8 T1,T2
	REAL*8 FREQ_SCL_FAC
!
	INTEGER I,J,K
	INTEGER ICOUNT
	INTEGER NT
	INTEGER INDX_1,INDX_2,NV
	INTEGER NXT_LOC
	INTEGER IOS
	LOGICAL OKAY
	LOGICAL CREATE_SUMMARY
	CHARACTER*132 STRING,FILENAME
!
	INTEGER, PARAMETER :: NCF_MAX=10000
	REAL*8 XV(NCF_MAX),YV(NCF_MAX)
!
	INTEGER, PARAMETER :: NREC_MAX=5
	REAL*8 TEMP_VEC(NREC_MAX)
	REAL*8 TOTAL_REC_VEC(NREC_MAX)
	REAL*8 LEVEL_REC_VEC(NREC_MAX)
!
	REAL*8 STAT_WEIGHT,GION,EDGE
	REAL*8 LEVEL_REC,TOTAL_REC,TEMP
	REAL*8 ANG_TO_HZ,SPEED_OF_LIGHT
	REAL*8, PARAMETER :: RONE=1.0D0
	LOGICAL DO_WAVE_PLT
        LOGICAL PLOT_REL_TO_GS_EDGE
	LOGICAL DO_RECOM
	LOGICAL DO_ALL_RECOM
	LOGICAL DO_SEQ_PLTS
	LOGICAL OSCILLATOR_FILE_AVAIL
	LOGICAL OUT_PHOT
	EXTERNAL SPEED_OF_LIGHT
!
	CHARACTER(LEN=30) UC; EXTERNAL UC
	CHARACTER(LEN=30) XLAB
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER, PARAMETER :: LUER=6
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUSCR=8
	INTEGER, PARAMETER :: LUOUT=12
!
! Set constants.
!
        CHIBF=2.815D-06
        CHIFF=3.69D-29
        HDKT=4.7994145D0
        TWOHCSQ=0.0147452575D0
        OPLIN=2.6540081D+08
        EMLIN=5.27296D-03
	ANG_TO_HZ=1.0D-07*SPEED_OF_LIGHT()
	AMASS=40.0D0
!
! Read in bound-free gaunt factors for individual n states of hydrogen,
! and hydrogenic cross-sections for individual l states (n =0 to 30,
! l=0 to n-1)
!
        CALL RD_HYD_BF_DATA(LUIN,LUSCR,LUER)
!
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')' Program designed to plot/compare opacity project data from photoionization files'
	WRITE(6,'(A)')DEF_PEN
!
10	FILENAME='PHOT1'
	CALL GEN_IN(FILENAME,RED_PEN//'FIRST'//DEF_PEN//' photoionization file:')
	NPAIRS_1=0
	OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,'(A)')RED_PEN
	    WRITE(6,*)'Unable to open file: ',TRIM(FILENAME)
	    WRITE(6,*)'IOSTAT=',IOS
	    WRITE(6,'(A)')DEF_PEN
	    GOTO 10
	  END IF
	  DO WHILE(NPAIRS_1 .EQ. 0)
	    READ(10,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(6,'(A)')RED_PEN
	      WRITE(6,*)'Error reading file: ',TRIM(FILENAME)
	      WRITE(6,*)'IOSTAT=',IOS
	      WRITE(6,'(A)')DEF_PEN
	      STOP
	    END IF
	    IF( INDEX(STRING,'!Screened nuclear charge') .NE. 0)THEN
	      READ(STRING,*)ZION_1
	    ELSE IF( INDEX(STRING,'Statistical weight of ion') .NE. 0)THEN
	      READ(STRING,*)GION_1
	    ELSE IF( INDEX(STRING,'!Split J levels') .NE. 0)THEN
	      READ(STRING,*)SPLITJ_1
	    ELSE IF( INDEX(STRING,'!Excitation energy of final state') .NE. 0)THEN
	      READ(STRING,*)EXC_EN_1
	    ELSE IF( INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
	      READ(STRING,*)NLEV_1
	      ALLOCATE (TYPE_1(NLEV_1),NUM_VALS_1(NLEV_1),LOC_1(NLEV_1),
	1               NAME_1(NLEV_1),ENERGY_1(NLEV_1),STAT_WT_1(NLEV_1),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(6,'(A)')RED_PEN
	        WRITE(6,*)'Error in PLT_PHOT_RAW -- unable to allocate storage (1)'
	        WRITE(6,*)'STAT=',IOS
	        WRITE(6,*)'NLEV_1=',NLEV_1
	        WRITE(6,'(A)')DEF_PEN
	        STOP
	      END IF
	    END IF
	    IF( INDEX(STRING,'!Total number of data pairs') .NE. 0)THEN
	      READ(STRING,*)NPAIRS_1
	    END IF
	  END DO
!
	  ALLOCATE (NU_1(NPAIRS_1))
	  ALLOCATE (CROSS_1(NPAIRS_1))
!
	  NXT_LOC=1
	  DO J=1,NLEV_1
	    DO WHILE(INDEX(STRING,'!Configuration name') .EQ. 0)
	      READ(10,'(A)')STRING
	    END DO
	    K=INDEX(STRING,'  ')
	    NAME_1(J)=STRING(1:K-1)
	    READ(10,*)TYPE_1(J)
	    READ(10,*)NUM_VALS_1(J)
	    WRITE(16,*)TRIM(NAME_1(J)),TYPE_1(J),NUM_VALS_1(J)
	    IF(TYPE_1(J) .EQ. 9 .AND. MOD(NUM_VALS_1(J),8) .NE. 0)THEN
	      WRITE(6,*)'Invalid number of data points for cross-section 9'
	      STOP
	    END IF
	    IF(TYPE_1(J) .GE. 20 .AND. TYPE_1(J) .LE. 24)THEN
	      READ(10,*)(NU_1(I),CROSS_1(I), I=NXT_LOC,NXT_LOC+NUM_VALS_1(J)-1)
	    ELSE
	      READ(10,*)(CROSS_1(I), I=NXT_LOC,NXT_LOC+NUM_VALS_1(J)-1)
	    END IF
	    LOC_1(J)=NXT_LOC
	    NXT_LOC=NXT_LOC+NUM_VALS_1(J)
	    STRING=' '
	  END DO
!
	I=INDEX(FILENAME,'_A')
	IF(I .EQ.  0)I=INDEX(FILENAME,'_B')
	IF(I .GT. 6)THEN
	  FILENAME=FILENAME(5:I)//'F_OSCDAT'	
	ELSE
	  FILENAME=' '
	END IF
!
! Read in energy levels from oscilator file, if it is available.
! May need to change if oscilator file format changes.
!
	WRITE(6,'(A)')BLUE_PEN
	WRITE(6,'(A)')' Not all options available if no oscillator file '
	WRITE(6,'(A)')' Can only treat opacity data when no oscillator file'
	WRITE(6,'(A)')DEF_PEN
!
14	CONTINUE
	CALL GEN_IN(FILENAME,'File with oscillator data: "" for no file')
	IF(FILENAME .NE. " ")THEN
	  OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)'Unable to open oscillator file: ',TRIM(FILENAME)
	    WRITE(6,*)'IOSTAT=',IOS
	    WRITE(6,*)DEF_PEN
	    GOTO 14
	  END IF
	  STRING=' '
	  DO WHILE(INDEX(STRING,'!Number of energy levels') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  READ(STRING,*)NELEV
	  DO WHILE(INDEX(STRING,'!Ionization energy') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  READ(STRING,*)IONIZATION_ENERGY
	  DO WHILE(INDEX(STRING,'!Number of transitions') .EQ. 0)
	    READ(10,'(A)')STRING
	  END DO
	  READ(10,'(A)')STRING		!Get blank line
!	  CALL GEN_IN(NELEV,'Number of energy levels')
	  ALLOCATE(G(NELEV),E_NAME(NELEV),ENERGY(NELEV))
!
	  G(:)=0.0D0;ENERGY(:)=0.0D0
	  DO I=1,NELEV
	    READ(10,'(A)')STRING
	    K=INDEX(STRING,'  ')
	    E_NAME(I)=STRING(1:K-1)
	    READ(STRING(K:),*)G(I),ENERGY(I)
	    IF(.NOT. SPLITJ_1)THEN
	      IF(E_NAME(I)(K-1:K-1) .EQ.  ']')THEN
	        K=INDEX(E_NAME(I),'[')
	        E_NAME(I)(K:)=' '
	      END IF
	    END IF
	  END DO
	  CLOSE(UNIT=10)
!
! Now need to match names. We assume photoionization files do not
! have [].
!
	  ENERGY_1(1:NLEV_1)=0.0D0; STAT_WT_1(1:NLEV_1)=0.0D0
	  WRITE(6,*)BLUE_PEN
	  WRITE(6,*)'Entering name comparison loop'
	  DO I=1,NLEV_1
	    DO J=1,NELEV
	      IF(E_NAME(J) .EQ. NAME_1(I))THEN
	        ENERGY_1(I)=ENERGY_1(I)+G(J)*ENERGY(J)
	        STAT_WT_1(I)=STAT_WT_1(I)+G(J)
	        WRITE(126,'(I5,ES15.6,F7.1,4X,A)')J,ENERGY(J),G(J),E_NAME(J)
	      END IF
	    END DO
	    IF(STAT_WT_1(I) .EQ. 0)THEN
	      WRITE(6,*)RED_PEN,'Error - no match for level',NAME_1(I),I
	    ELSE
	      ENERGY_1(I)=ENERGY_1(I)/STAT_WT_1(I)
	      WRITE(126,'(I5,ES15.6,4X,A)')I,ENERGY_1(I),TRIM(NAME_1(I))
	      IF(STAT_WT_1(I) .NE. 0.0D0)ENERGY_1(I)=1.0D-15*SPEED_OF_LIGHT()*
	1        (IONIZATION_ENERGY-ENERGY_1(I))
	      WRITE(126,'(I5,ES15.6,4X,A)')I,ENERGY_1(I),TRIM(NAME_1(I))
	    END IF
	  END DO
	  OSCILLATOR_FILE_AVAIL=.TRUE.
	  WRITE(6,*)'Number of levels in oscillator file is',NELEV
	  WRITE(6,*)DEF_PEN
	ELSE
	  ENERGY_1(1:NLEV_1)=0.0D0; STAT_WT_1(1:NLEV_1)=0.0D0
	  OSCILLATOR_FILE_AVAIL=.FALSE.
	END IF
!
!
!
20	NPAIRS_2=0
	FILENAME='PHOT2'
	WRITE(6,*)RED_PEN
	CALL GEN_IN(FILENAME,'SECOND'//DEF_PEN//'photoionization file ("" for null):')
	IF(FILENAME .EQ. ' ')GOTO 1000
	OPEN(UNIT=10,STATUS='OLD',FILE=FILENAME,ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)RED_PEN
	    WRITE(6,*)'Unable to open file: ',TRIM(FILENAME)
	    WRITE(6,*)'IOSTAT=',IOS
	    GOTO 20
	  END IF
	  DO WHILE(NPAIRS_2 .EQ. 0)
	    READ(10,'(A)',IOSTAT=IOS)STRING
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)RED_PEN
	      WRITE(6,*)'Error reading file: ',TRIM(FILENAME)
	      WRITE(6,*)'IOSTAT=',IOS
	      WRITE(6,*)DEF_PEN
	      STOP
	    END IF
	    IF( INDEX(STRING,'!Screened nuclear charge') .NE. 0)THEN
	      READ(STRING,*)ZION_2
	    ELSE IF( INDEX(STRING,'Statistical weight of ion') .NE. 0)THEN
	      READ(STRING,*)GION_2
	    ELSE IF( INDEX(STRING,'!Excitation energy of final state') .NE. 0)THEN
	      READ(STRING,*)EXC_EN_2
	    ELSE IF( INDEX(STRING,'!Split J levels') .NE. 0)THEN
	      READ(STRING,*)SPLITJ_2
	    ELSE IF( INDEX(STRING,'!Number of energy levels') .NE. 0)THEN
	      READ(STRING,*)NLEV_2
	      ALLOCATE (TYPE_2(NLEV_2),NUM_VALS_2(NLEV_2),LOC_2(NLEV_2), NAME_2(NLEV_2),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)RED_PEN
	        WRITE(6,*)'Error in PLT_PHOT_RAW -- unable to allocate storage (2)'
	        WRITE(6,*)'STAT=',IOS
	        WRITE(6,*)'NLEV_2=',NLEV_2
	        WRITE(6,*)DEF_PEN
	        STOP
	      END IF
	    END IF
	    IF( INDEX(STRING,'!Total number of data pairs') .NE. 0)THEN
	      READ(STRING,*)NPAIRS_2
	    END IF
	  END DO
!
	  ALLOCATE (NU_2(NPAIRS_2))
	  ALLOCATE (CROSS_2(NPAIRS_2))
!
	  NXT_LOC=1
	  DO J=1,NLEV_2
	    DO WHILE(INDEX(STRING,'!Configuration name') .EQ. 0)
	      READ(10,'(A)')STRING
	    END DO
	    K=INDEX(STRING,'  ')
	    NAME_2(J)=STRING(1:K-1)
	    READ(10,*)TYPE_2(J)
	    READ(10,*)NUM_VALS_2(J)
	    IF(TYPE_2(J) .EQ. 9 .AND. MOD(NUM_VALS_2(J),8) .NE. 0)THEN
	      WRITE(6,*)RED_PEN
	      WRITE(6,*)'Invalid number of data points for cross-section 9'
	      WRITE(6,*)DEF_PEN
	      STOP
	    END IF
	    IF(TYPE_2(J) .GE. 20 .AND. TYPE_2(J) .LE. 24)THEN
	      READ(10,*)(NU_2(I),CROSS_2(I), I=NXT_LOC,NXT_LOC+NUM_VALS_2(J)-1)
	    ELSE
	      READ(10,*)(CROSS_2(I), I=NXT_LOC,NXT_LOC+NUM_VALS_2(J)-1)
	    END IF
	    LOC_2(J)=NXT_LOC
	    NXT_LOC=NXT_LOC+NUM_VALS_2(J)
	    STRING=' '
	  END DO
	CLOSE(UNIT=10)
1000	CONTINUE
!
! 
!
	CREATE_SUMMARY=.FALSE.
	WRITE(6,'(A)')' '
	CALL GEN_IN(CREATE_SUMMARY,'Create a summary of photoionization cross sections (1st file only)?')
	IF(CREATE_SUMMARY)THEN
	  OPEN(UNIT=11,FILE='Phot_summary',STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(11,'(A,T42,A,6X,A,3(10X,A))')'Level','Type',' Np','X1','X2','X2'
	    DO K=1,40
	      DO J=1,NLEV_1
	        IF(TYPE_1(J) .EQ. K)THEN
	          WRITE(11,'(A,T42,I4,2X,I6,5ES13.4)')TRIM(NAME_1(J)),TYPE_1(J),NUM_VALS_1(J),
	1            CROSS_1(LOC_1(J):LOC_1(J)+MIN(4,NUM_VALS_1(J)-1))
	        END IF
	      END DO
	    END DO
	  CLOSE(UNIT=11)
	END IF
!
	DO_ALL_RECOM=.FALSE.
	OUT_PHOT=.FALSE.
	IF(OSCILLATOR_FILE_AVAIL)CALL GEN_IN(DO_ALL_RECOM,'Compute recombination rates for all levels?')
	IF(OSCILLATOR_FILE_AVAIL)CALL GEN_IN(OUT_PHOT,'Output tabulated set of photoioization cross-sections')
	IF(DO_ALL_RECOM)THEN
	  WRITE(6,*)BLUE_PEN
	  WRITE(6,*)'Data will be written to RECOM_SUM'
	  OPEN(UNIT=LUOUT,FILE='RECOM_SUM',STATUS='UNKNOWN',ACTION='WRITE')
	  TEMP_VEC(1)=0.5D0; TEMP_VEC(2)=1.0D0; TEMP_VEC(3)=2.0; TEMP_VEC(4)=5.0D0; TEMP_VEC(5)=10.0D0
	  TOTAL_REC_VEC(:)=0.0D0
	  CALL GEN_IN(TEMP_VEC,NT,NREC_MAX,'Temperature in 10^4K (5 values max)')
	  WRITE(6,*)GREEN_PEN
	  WRITE(6,'(A,T30,5(5X,F6.2))')' Temperature (10^4 K)=',(TEMP_VEC(I),I=1,NT)
	  WRITE(LUOUT,'(A,T30,5(5X,F6.2))')' Level / Temperature (10^4 K)',(TEMP_VEC(I),I=1,NT)
	  DO INDX_1=1,NLEV_1
	    EDGE=ENERGY_1(INDX_1)
	    STAT_WEIGHT=STAT_WT_1(INDX_1)
	    IF(TYPE_1(INDX_1) .GE. 20 .AND. TYPE_1(INDX_1) .LE. 23)THEN
	      NV=NUM_VALS_1(INDX_1)
	      XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	      YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	      FREQ_SCL_FAC=EDGE+EXC_EN_1
	    ELSE IF(TYPE_1(INDX_1) .EQ. 24)THEN
	      NV=NUM_VALS_1(INDX_1)-1
	      XV(1:NV)=NU_1(LOC_1(INDX_1+1):LOC_1(INDX_1)+NV)
	      YV(1:NV)=CROSS_1(LOC_1(INDX_1+1):LOC_1(INDX_1)+NV)
	      FREQ_SCL_FAC=1.0D0
	    ELSE
	      NV=1000
	      CALL RAW_SUBPHOT_V2(YV,XV,CROSS_1(LOC_1(INDX_1)),TYPE_1(INDX_1),NUM_VALS_1(INDX_1),
	1                    EDGE,EXC_EN_1,ZION_1,AMASS,NAME_1(INDX_1),NV)
	      FREQ_SCL_FAC=XV(1)
	      XV(1:NV)=XV(1:NV)/FREQ_SCL_FAC
	    END IF
	    IF(OUT_PHOT)THEN
	       WRITE(70,'(A,T50,A)')NAME_1(INDX_1),'!Configuration name'
	       WRITE(70,'(I2,T50,A)')21,'!Type of cross-section'
	       IF(FREQ_SCL_FAC .GT. 1.0001D0*EDGE)THEN
	         WRITE(STRING,*)NV+2
	       ELSE
	         WRITE(STRING,*)NV
	       END IF
	       WRITE(70,'(A,T50,A)')TRIM(STRING),'!Number of cross-section points'
	       IF(FREQ_SCL_FAC .GT. 1.0001D0*EDGE)THEN
	         WRITE(70,'(2ES12.4)')1.0D0,0.0D0
	         WRITE(70,'(2ES12.4)')FREQ_SCL_FAC/EDGE/1.0001D0,0.0D0
	       END IF
	       DO I=1,NV
	         WRITE(70,'(2ES12.4)')XV(I)*FREQ_SCL_FAC/EDGE,YV(I)
	       END DO
	    ELSE
	      T1=EDGE+EXC_EN_1
	      DO I=1,NT
	        CALL RECOM_OPAC_V2(YV,XV,T1,FREQ_SCL_FAC,STAT_WEIGHT,GION_1,NV,NV,LEVEL_REC_VEC(I),TEMP_VEC(I))
	        TOTAL_REC_VEC(I)=TOTAL_REC_VEC(I)+LEVEL_REC_VEC(I)
	      END DO
	      WRITE(6,'(A,T30,5ES11.3)')TRIM(NAME_1(INDX_1)),(LEVEL_REC_VEC(I),I=1,NT)
	      WRITE(LUOUT,'(X,A,T30,5ES11.3)')TRIM(NAME_1(INDX_1)),(LEVEL_REC_VEC(I),I=1,NT)
	      FLUSH(LUOUT)
	    END IF
	  END DO
	  WRITE(6,'(A,T30,5ES11.3)')' Total Recom. Rate/ion=',(TOTAL_REC_VEC(I),I=1,NT)
	  WRITE(LUOUT,'(A,T30,5ES11.3)')' Total Recom. Rate/ion=',(TOTAL_REC_VEC(I),I=1,NT)
!
	  IF(EXC_EN_1 .NE. 0.0D0)THEN
	    WRITE(6,*)' '
	    WRITE(6,*)'These give the total recombination rate / ground state ion'
	    WRITE(6,*)' '
	    WRITE(LUOUT,*)' '
	    WRITE(LUOUT,*)'These give the total recombination rate / ground state ion'
	    WRITE(LUOUT,*)' '
	    T1=1.0D0
	    CALL GEN_IN(T1,'Statistical weight for ion ground level')
	    T1=GION_1/T1
	    T2=-HDKT*EXC_EN_1
	    WRITE(6,'(A,T30,5ES11.3)')' Recom. Rate/g.s. ion=',(T1*TOTAL_REC_VEC(I)*EXP(T2/TEMP_VEC(I)),I=1,NT)
	    WRITE(LUOUT,'(A,T30,5ES11.3)')' Recom. Rate/g.s. ion=',(T1*TOTAL_REC_VEC(I)*EXP(T2/TEMP_VEC(I)),I=1,NT)
	  END IF
	  WRITE(6,*)DEF_PEN
	  CLOSE(LUOUT)
	END IF
!
	WRITE(6,*)BLUE_PEN
	WRITE(6,*)'These reamining options refer to individual cross-sections'
	WRITE(6,*)DEF_PEN
	DO_WAVE_PLT=.FALSE.;  CALL GEN_IN(DO_WAVE_PLT,'Plot versus wavelength?')
	DO_RECOM=.FALSE.;     CALL GEN_IN(DO_RECOM,'Compute recombination rate?')
!
	WRITE(6,*)BLUE_PEN
	WRITE(6,*)'Plotting relative to the ground state photoianiozatiom limit only makes sense'
	WRITE(6,*)'for levels below the ionization limit.'
	WRITE(6,*)DEF_PEN
        PLOT_REL_TO_GS_EDGE=.FALSE.
	CALL GEN_IN(PLOT_REL_TO_GS_EDGE,'Use frequency in units of iomization energy to the ground state?')
!
	LEVEL_NAME1=' '
	DO_SEQ_PLTS=.FALSE.
	CALL GEN_IN(DO_SEQ_PLTS,'Plot photoionization cross-section for a sequence: eg., 6d (T or F)')
	IF(DO_SEQ_PLTS)THEN
	  DO WHILE(1 .EQ. 1)
	    WRITE(6,'(A)')' '
	    CALL GEN_IN(LEVEL_NAME1,'Sequence desciptor')
	    IF(LEVEL_NAME1 .EQ. ' ')STOP
	    WRITE(6,'(A)')' '
	    ICOUNT=0
	    DO J=1,NLEV_1
	      IF(INDEX(NAME_1(J),TRIM(LEVEL_NAME1)) .NE. 0)THEN
                ICOUNT=ICOUNT+1
	        INDX_1=J
	        EDGE=ENERGY_1(INDX_1)
	        WRITE(6,'(4X,A,A,T40,ES10.4)')'Level name/energy is: ',TRIM(NAME_1(INDX_1)),EDGE
	        IF(TYPE_1(INDX_1) .GE. 20 .AND. TYPE_1(INDX_1) .LE. 23)THEN
	          NV=NUM_VALS_1(INDX_1)
	          XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	          YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	        ELSE IF(TYPE_1(INDX_1) .EQ. 24)THEN
	          NV=NUM_VALS_1(INDX_1)-1
	          XV(1:NV)=NU_1(LOC_1(INDX_1)+1:LOC_1(INDX_1)+NV)
	          YV(1:NV)=CROSS_1(LOC_1(INDX_1)+1:LOC_1(INDX_1)+NV)
	        ELSE
	          NV=2000
	          CALL RAW_SUBPHOT_V2(YV,XV,CROSS_1(LOC_1(INDX_1)),TYPE_1(INDX_1),NUM_VALS_1(INDX_1),
	1                    EDGE,EXC_EN_1,ZION_1,AMASS,NAME_1(INDX_1),NV)
	          XV(1:NV)=XV(1:NV)/(EDGE+EXC_EN_1)
	        END IF
	        IF(DO_WAVE_PLT)THEN
	          XV(1:NV)=ANG_TO_HZ/XV(1:NV)
	        END IF
	        CALL DP_CURVE(NV,XV,YV)
	      END IF
	      IF(ICOUNT .EQ. 50)THEN
	        WRITE(6,*)RED_PEN//'Maximum number of plots exceeded'//DEF_PEN
	        EXIT
	      END IF
	    END DO
	    WRITE(6,'(A)')' '
	   CALL GRAMON_PGPLOT('\gn/\gn\do\u ','\gs(Mb)',TRIM(LEVEL_NAME1),' ')
	  END DO
	END IF
!
	XLAB='\gn/\gn\do\u'
	IF(DO_WAVE_PLT)XLAB='\gl(\A)'
	DO WHILE(1 .EQ. 1)
!
! We can now use a level name, or an index (as stored in PHOT file).
!
	  INDX_1=0
	  LEVEL_NAME1='1'
	  DO WHILE(INDX_1 .EQ. 0)
	    WRITE(6,'(A)')' '
100	    CALL GEN_IN(LEVEL_NAME1,'Level name [or index] for File 1 (P to plot, E to exit)')
	    IF(UC(LEVEL_NAME1) .EQ. 'P')THEN
	      CALL GRAMON_PGPLOT(XLAB,'\gs(Mb)',LEVEL_NAME1,' ')
	      CALL GEN_IN(LEVEL_NAME1,'Level name [or index] for File (E to exit)')
	    END IF
	    IF( UC(LEVEL_NAME1) .EQ. 'E' .OR. UC(LEVEL_NAME1(1:2)) .EQ. 'EX' 
	1         .OR. UC(LEVEL_NAME1) .EQ. ' ')STOP
	    DO J=1,NLEV_1
	      IF(LEVEL_NAME1 .EQ. NAME_1(J))THEN
                INDX_1=J
                EXIT
	      END IF
	    END DO
	    IF(INDX_1 .EQ. 0)THEN
	      READ(LEVEL_NAME1,*,IOSTAT=IOS)INDX_1
	      IF(INDX_1 .EQ. 0)THEN
	         WRITE(6,*)'Invalid index'
	         GOTO 100
	      END IF
	      IF(IOS .NE. 0)THEN
	        WRITE(6,*)'Error - level name not found'
	        DO I=1,NLEV_1,5
	          WRITE(6,'(5A14)')(TRIM(NAME_1(J)),J=I,MAX(I+4,NLEV_1))
	        END DO
	        INDX_1=0
	      ELSE IF(OSCILLATOR_FILE_AVAIL)THEN
	        IF(INDX_1 .GT. NELEV)THEN
	          WRITE(6,*)'Invalid index; index should be < ',NELEV
	          GOTO 100
	        END IF
	        LEVEL_NAME1=E_NAME(INDX_1)
	        DO I=1,NLEV_1
	          IF(LEVEL_NAME1 .EQ. NAME_1(I))THEN
	            INDX_1=I
	            EXIT
	          END IF
	        END DO
	      ELSE
	        LEVEL_NAME1=NAME_1(INDX_1)
	      END IF
	    END IF
	  END DO
	  EDGE=ENERGY_1(INDX_1)
	  WRITE(6,'(/,8X,A,A)')     '          Level_1 name is: ',NAME_1(INDX_1)
	  WRITE(6,'(8X,A,ES15.8,A)')'       Energy of level is: ',EDGE,' (10^15Hz)'
	  WRITE(6,'(8X,A,I2,/)')    ' Type of cross-section is: ',TYPE_1(INDX_1)
!
	  IF(TYPE_1(INDX_1) .GE. 20 .AND. TYPE_1(INDX_1) .LE. 23)THEN
	    NV=NUM_VALS_1(INDX_1)
	    XV(1:NV)=NU_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	    YV(1:NV)=CROSS_1(LOC_1(INDX_1):LOC_1(INDX_1)+NV-1)
	    FREQ_SCL_FAC=EDGE+EXC_EN_1
	  ELSE IF(TYPE_1(INDX_1) .EQ. 24)THEN
	    NV=NUM_VALS_1(INDX_1)-1
	    XV(1:NV)=NU_1(LOC_1(INDX_1)+1:LOC_1(INDX_1)+NV)
	    YV(1:NV)=CROSS_1(LOC_1(INDX_1)+1:LOC_1(INDX_1)+NV)
	    FREQ_SCL_FAC=1.0D0
	  ELSE
	    NV=1000
	    CALL RAW_SUBPHOT_V2(YV,XV,CROSS_1(LOC_1(INDX_1)),TYPE_1(INDX_1),NUM_VALS_1(INDX_1),
	1                    EDGE,EXC_EN_1,ZION_1,AMASS,NAME_1(INDX_1),NV)
	    FREQ_SCL_FAC=XV(1)
	    XV(1:NV)=XV(1:NV)/FREQ_SCL_FAC
	  END IF
	  IF(DO_WAVE_PLT)THEN
	    XV(1:NV)=ANG_TO_HZ/FREQ_SCL_FAC/XV(1:NV)
	  ELSE IF(PLOT_REL_TO_GS_EDGE)THEN
	    XV(1:NV)=XV(1:NV)*FREQ_SCL_FAC/EDGE
	  END IF
	  CALL DP_CURVE(NV,XV,YV)
!
	  STAT_WEIGHT=MAX(RONE,STAT_WT_1(INDX_1))
	  GION=1 
	  IF(EDGE .EQ. 0.0D0)THEN
	    WRITE(6,'(A)')' '
	    WRITE(6,'(A)')'A non-zero ionization energy is only required to compute'
	    WRITE(6,'(A)')'the recombination rate.'
	    WRITE(6,'(A)')' '
	    CALL GEN_IN(EDGE,'Ionization energy (10^15 Hz)')
	    DO_RECOM=.TRUE.
	    IF(EDGE .EQ. 0)DO_RECOM=.FALSE.
	  END IF
	  IF(DO_RECOM .AND. EDGE .NE. 0)THEN
	    TEMP=1.0D0
	    CALL GEN_IN(TEMP,'Temperature in 10^4K')
	    IF(TEMP .NE. 0)THEN
	      WRITE(6,'(A,ES10.4,3X,A,F5.1,3X,A,F4.1,3X,A,I6,3X,A,F6.2)')
	1        ' EDGE=',EDGE,'g=',STAT_WEIGHT,'gion=',GION_1,'NV=',NV,'T(10^K)=',TEMP
	      T1=EDGE+EXC_EN_1
	      CALL RECOM_OPAC_V2(YV,XV,T1,FREQ_SCL_FAC,STAT_WEIGHT,GION_1,NV,NV,TOTAL_REC,TEMP)
	      WRITE(6,'(A,ES11.4)')' Total Rec=',TOTAL_REC
	    END IF
	  END IF
!
	  IF(NPAIRS_2 .NE. 0)THEN
	    LEVEL_NAME2=NAME_1(INDX_1)
	    INDX_2=0
	    DO WHILE(INDX_2 .EQ. 0)
	      CALL GEN_IN(LEVEL_NAME2,'Level name for File 2')
              DO J=1,NLEV_2
	        IF(LEVEL_NAME2 .EQ. NAME_2(J))THEN
                  INDX_2=J
                  EXIT
	        END IF
	      END DO
	      IF(INDX_2 .EQ. 0)THEN
	        READ(LEVEL_NAME2,*,IOSTAT=IOS)INDX_2
	        IF(IOS .NE. 0)THEN
	          WRITE(6,*)'Error - level name not found'
	          DO I=1,NLEV_2,5
	            WRITE(6,'(5A14)')(TRIM(NAME_2(J)),J=I,MAX(I+4,NLEV_2))
	          END DO
	          INDX_2=0
	        ELSE
	          LEVEL_NAME2=NAME_2(INDX_2)
	        END IF
	      END IF
	    END DO
!
	    EDGE=ENERGY_1(INDX_1)
	    IF(INDX_2 .GT. 0)THEN
	      IF(TYPE_2(INDX_2) .EQ. 20 .OR. TYPE_2(INDX_2) .LE. 23)THEN
	        NV=NUM_VALS_2(INDX_2)
	        XV(1:NV)=NU_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	        YV(1:NV)=CROSS_2(LOC_2(INDX_2):LOC_2(INDX_2)+NV-1)
	        FREQ_SCL_FAC=EDGE+EXC_EN_2
	      ELSE
	        NV=1000
	        CALL RAW_SUBPHOT_V2(YV,XV,CROSS_2(LOC_2(INDX_2)),TYPE_2(INDX_2),NUM_VALS_2(INDX_2),
	1                    EDGE,EXC_EN_2,ZION_2,AMASS,NAME_2(INDX_2),NV)
	        FREQ_SCL_FAC=XV(1)
	        XV(1:NV)=XV(1:NV)/FREQ_SCL_FAC
	      END IF
	      IF(DO_RECOM .AND. OSCILLATOR_FILE_AVAIL)THEN
	        IF(TEMP .NE. 0)THEN
	          T1=EDGE+EXC_EN_2
	          CALL RECOM_OPAC_V2(YV,XV,T1,FREQ_SCL_FAC,STAT_WEIGHT,GION,NV,NV,TOTAL_REC,TEMP)
	          WRITE(6,*)'Total Rec=',TOTAL_REC
	        END IF
	      END IF
	    END IF
	    IF(DO_WAVE_PLT)THEN
	      XV(1:NV)=ANG_TO_HZ/FREQ_SCL_FAC/XV(1:NV)
	    ELSE IF(PLOT_REL_TO_GS_EDGE)THEN
	      XV(1:NV)=XV(1:NV)*FREQ_SCL_FAC/EDGE
	    END IF
	    CALL DP_CURVE(NV,XV,YV)
	  END IF
!
	END DO		!Loop over levels
!
	END
