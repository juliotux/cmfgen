!
! Routinine to compute the frequency grid for CMFGEN. Routine also
! allocates the vectors needed for the line data, sets the line data,
! and puts the line data into numerical order.
!
! The option FREQ_GRID_OPTION was introduced to allow an old frequency grid to be used.
! This could be useful when testing new options, and tesing the code against
! earlier versions. Its still possible that the grids may not be absolutely
! identical.
!
! If IOPT is 0, and anything else except 1, or 2 (only 1 is curently implmented), 
! the default frequncy grid will used.
!
	SUBROUTINE SET_FREQUENCY_GRID_V2(NU,FQW,LINES_THIS_FREQ,NU_EVAL_CONT,
	1               NCF,NCF_MAX,N_LINE_FREQ,ND,
	1               OBS_FREQ,OBS,N_OBS,LUIN,IMPURITY_CODE)
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE LINE_VEC_MOD
	IMPLICIT NONE
!
! Created 8-Jun-2004
!
	INTEGER NCF_MAX
	INTEGER NCF				!Total number of continuum points.
	INTEGER N_LINE_FREQ                     !Number of lines
	INTEGER N_OBS
	INTEGER ND
	INTEGER LUIN
!
	REAL*8 NU(NCF_MAX)
	REAL*8 FQW(NCF_MAX)
	REAL*8 NU_EVAL_CONT(NCF_MAX)
	INTEGER LINES_THIS_FREQ(NCF_MAX)
!
	REAL*8 OBS_FREQ(NCF_MAX)
	REAL*8 OBS(NCF_MAX)
	REAL*8 CHIL(ND)
	REAL*8 LOC_ED(ND)
	REAL*8 LOC_T(ND)
!
	LOGICAL IMPURITY_CODE
	LOGICAL GRID_EXISTS
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
	REAL*8 SPEED_OF_LIGHT
	REAL*8 LAMVACAIR
!
	REAL*8 NU_MAX_OBS
	REAL*8 NU_MIN_OBS
	REAL*8 T1,T2
	REAL*8 C_KMS
!
	INTEGER I,J,K
	INTEGER ID
	INTEGER ML
	INTEGER NL,NUP
	INTEGER MNL,MNUP
	INTEGER NLINES_PROF_STORE
	INTEGER NFREQ_PROF_STORE
	LOGICAL FIRST
!
	CHARACTER(LEN=132) TEMP_CHAR
!
! Constants for opacity etc. These are set in CMFGEN.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 OPLIN,EMLIN
!
	LUER=ERROR_LU()
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
!
!
! 
! Set up lines that will be treated with the continuum calculation.
! This section of code is also used by the code treating purely lines
! (either single transition Sobolev or CMF, or overlapping Sobolev).
!
! To define the line transitions we need to operate on the FULL atom models.
! We thus perform separate loops for each species. 
!
	ML=0			!Initialize line counter.
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNUP=2,ATM(ID)%NXzV_F
	      NUP= ATM(ID)%F_TO_S_XzV(MNUP)+ ATM(ID)%EQXzV-1
	      DO MNL=1,MNUP-1
	        NL= ATM(ID)%F_TO_S_XzV(MNL)+ ATM(ID)%EQXzV-1
	        IF(ATM(ID)%AXzV_F(MNL,MNUP) .NE. 0)ML=ML+1
	      END DO
	    END DO
	  END IF
	END DO
	N_LINE_FREQ=ML
!
! Now that we have the number of lines, we can allocate the needed memory.
! VEC_TRANS_NAME is allocated temporaruly so that we can output the full 
! transitions name to TRANS_INFO.
!
	ALLOCATE (VEC_TRANS_NAME(N_LINE_FREQ),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_FREQ(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_STRT_FREQ(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_OSCIL(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_EINA(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_ARAD(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_VDOP_MIN(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_DP_WRK(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_INDX(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_ID(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_NL(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_NUP(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_MNL_F(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_MNUP_F(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_INT_WRK(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( PROF_LIST_LOCATION(N_LINE_FREQ) ,STAT=IOS)
!
	IF(IOS .EQ. 0)ALLOCATE( VEC_SPEC(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_CHAR_WRK(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( PROF_TYPE(N_LINE_FREQ) ,STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE( VEC_TRANS_TYPE(N_LINE_FREQ) ,STAT=IOS)
!
        IF(IOS .EQ. 0)ALLOCATE( LINE_ST_INDX_IN_NU(N_LINE_FREQ) ,STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE( LINE_END_INDX_IN_NU(N_LINE_FREQ) ,STAT=IOS)
        IF(IOS .EQ. 0)ALLOCATE( LINE_LOC(N_LINE_FREQ) ,STAT=IOS)
!
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Unable to allocate memory in SET_FREQUENCY_GRID'
	  STOP
	END IF
!
! PROF_T_ED should have the same format as a DC input file, and is used to
! provide an estimate of the electron density for STARK and VOIGT profile
! limits.
!
	INQUIRE(FILE='GRID_PARAMS',EXIST=GRID_EXISTS)
	IF(FIX_DOP .OR. GLOBAL_LINE_PROF(1:8) .EQ. 'DOP_SPEC')THEN
	ELSE IF(GRID_EXISTS)THEN
	  CALL SET_T_ED_GRID(LOC_T,LOC_ED,ND,7)
	ELSE
	  CALL RD_T_ED(LOC_T,LOC_ED,ND,7,'PROF_T_ED')
	END IF
!
! Now get the lines
!
	ML=0			!Initialize line counter.
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNUP=2,ATM(ID)%NXzV_F
	      NUP= ATM(ID)%F_TO_S_XzV(MNUP)+ ATM(ID)%EQXzV-1
	      DO MNL=1,MNUP-1
	        NL= ATM(ID)%F_TO_S_XzV(MNL)+ ATM(ID)%EQXzV-1
	        IF( ATM(ID)%AXzV_F(MNL,MNUP) .NE. 0)THEN
	          ML=ML+1
	          VEC_FREQ(ML)= ATM(ID)%EDGEXzV_F(MNL)- ATM(ID)%EDGEXzV_F(MNUP)
	          VEC_SPEC(ML)=ION_ID(ID)
	          VEC_ID(ML)=ID
	          VEC_NL(ML)=NL
	          VEC_NUP(ML)=NUP     
	          VEC_MNL_F(ML)=MNL
	          VEC_MNUP_F(ML)=MNUP     
	          VEC_OSCIL(ML)=ATM(ID)%AXzV_F(MNL,MNUP)
	          VEC_EINA(ML)=ATM(ID)%AXzV_F(MNUP,MNL)
	          VEC_ARAD(ML)= ATM(ID)%ARAD(MNL)+ATM(ID)%ARAD(MNUP)
	          VEC_TRANS_TYPE(ML)=ATM(ID)%XzV_TRANS_TYPE
	          VEC_TRANS_NAME(ML)=TRIM(VEC_SPEC(ML))//
	1           '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP))//'-'//
	1           TRIM( ATM(ID)%XzVLEVNAME_F(MNL))//')'
!
	          T1=12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
	          VDOP_FIX=T1
!
	          T1=VEC_OSCIL(ML)*OPLIN
	          T2=ATM(ID)%GXzV_F(MNL)/ATM(ID)%GXzV_F(MNUP)
	          DO I=1,ND
	            CHIL(I)=ABS(T1*(ATM(ID)%XzV_F(MNL,I)-T2*ATM(ID)%XzV_F(MNUP,I)))
	          END DO
	          PROF_TYPE(ML)=ATM(ID)%XzV_PROF_TYPE
	          IF(GLOBAL_LINE_PROF .NE. 'NONE')PROF_TYPE(ML)=GLOBAL_LINE_PROF
!
	          T2=0.0D0
	          CALL SET_PROF_LIMITS_V4(VEC_STRT_FREQ(ML),VEC_VDOP_MIN(ML),
	1             CHIL,LOC_ED,LOC_T,VTURB_VEC,ND,PROF_TYPE(ML),PROF_LIST_LOCATION(ML),
	1             VEC_FREQ(ML),MNL,MNUP,
	1             VEC_SPEC(ML),AT_MASS(SPECIES_LNK(ID)), ATM(ID)%ZXzV,
	1             VEC_ARAD(ML),T2,TDOP,AMASS_DOP,VTURB,                   !T2 is garbage
	1             DOP_PROF_LIMIT,VOIGT_PROF_LIMIT,
	1             V_PROF_LIMIT,MAX_PROF_ED,SET_PROF_LIMS_BY_OPACITY)
!
	          T1=0.01D0*C_KMS/VEC_FREQ(ML)
	          WRITE(134,'(A,T10,A,T25,2I5,2F20.8,2ES15.4)')VEC_SPEC(ML),PROF_TYPE(ML),MNL,MNUP,
	1                  T1,LAMVACAIR(VEC_FREQ(ML)),VEC_VDOP_MIN(ML),
	1                  C_KMS*(VEC_STRT_FREQ(ML)/VEC_FREQ(ML)-1.0D0)
	        END IF
	      END DO
	    END DO
	  END IF
	END DO
!
! 
!
! Get lines and arrange in numerically decreasing frequency. This will
! allow us to consider line overlap, and to include lines with continuum
! frequencies so that the can be handled automatically.
!
	N_LINE_FREQ=ML
!
	CALL INDEXX(N_LINE_FREQ,VEC_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_VDOP_MIN,VEC_INDX,VEC_DP_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_ID,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NL,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NUP,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,PROF_LIST_LOCATION,VEC_INDX,VEC_INT_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,PROF_TYPE,VEC_INDX,VEC_CHAR_WRK)
!
! Output all lines to TRANS_INFO. Usefule for diagnostic purposes.
!
	I=160	!Record length - allow for long names
	CALL GEN_ASCI_OPEN(LUIN,'TRANS_INFO','UNKNOWN',' ','WRITE',I,IOS)
	  WRITE(LUIN,*)'!'
	  WRITE(LUIN,*)' Wavelengths are in air for Lambda > 2000 A'
	  WRITE(LUIN,*)'!'
	  WRITE(LUIN,*)
	1     '     I    NL_F  NUP_F        Nu',
	1     '       Lam(A)    /\V(km/s)    Transition' 
	  WRITE(LUIN,
	1    '(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,16X,A)')
	1         IONE,VEC_MNL_F(1),VEC_MNUP_F(1),
	1         VEC_FREQ(1),LAMVACAIR(VEC_FREQ(1)),
	1         TRIM(VEC_TRANS_NAME(VEC_INDX(1)))
	  DO ML=2,N_LINE_FREQ
	    T1=LAMVACAIR(VEC_FREQ(ML))
	    T2=2.998D+05*(VEC_FREQ(ML-1)-VEC_FREQ(ML))/VEC_FREQ(ML)
	    IF(T2 .GT. 2.998D+05)T2=2.998D+05
	    IF(T1 .LT. 1.0D+04)THEN
	      WRITE(LUIN,
	1      '(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,2X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	    ELSE             
	      WRITE(LUIN,
	1      '(1X,I6,2(1X,I6),2X,F10.6,1X,1P,E11.4,0P,2X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	    END IF
	  END DO
	CLOSE(UNIT=LUIN)
	DEALLOCATE (VEC_TRANS_NAME)
!
! Get lines and arrange in numerically decreasing frequency according to
! the START frequency of the line. This will allow us to consider line overlap,
! and to include lines with continuum frequencies so that the can be handled 
! automatically.
!
	CALL INDEXX(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_ID,VEC_INDX,VEC_INT_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_VDOP_MIN,VEC_INDX,VEC_DP_WRK)
!
	CALL SORTINT(N_LINE_FREQ,VEC_NL,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NUP,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,PROF_LIST_LOCATION,VEC_INDX,VEC_INT_WRK)
!
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,PROF_TYPE,VEC_INDX,VEC_CHAR_WRK)
!
! 
!
! GLOBAL_LINE_SWITCH provides an option to handle all LINE by the same method.
! The local species setting only takes precedence when it is set to NONE.
!
	IF(GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE')THEN
	  DO I=1,N_LINE_FREQ
	    VEC_TRANS_TYPE(I)=GLOBAL_LINE_SWITCH
	  END DO
	ELSE
	  DO I=1,N_LINE_FREQ
	    CALL SET_CASE_UP(VEC_TRANS_TYPE(I),IZERO,IZERO)
	  END DO
	END IF
!
! If desired, we can set transitions with:
!      wavelengths > FLUX_CAL_LAM_END (in A) to the SOBOLEV option.
!      wavelengths < FLUX_CAL_LAM_BEG (in A) to the SOBOLEV option.
!
! The region defined by FLUX_CAL_LAM_BEG < LAM < FLUX_CAL_LAM_END will be computed using
! transition types determined by the earlier species and global options.
!
! Option has 2 uses:
!
! 1. Allows use of SOBOLEV approximation in IR where details of radiative
!    transfer is unimportant. In this case FLUX_CAL_LAM_BEG should be set to zero.
! 2. Allows a full flux calculation to be done in a limited wavelength region
!    as defined by FLUX_CAL_LAM_END and FLUX_CAL_LAM_BEG. 
!
	IF(SET_TRANS_TYPE_BY_LAM)THEN
	  IF(FLUX_CAL_LAM_END .LT. FLUX_CAL_LAM_BEG)THEN
	    WRITE(LUER,*)'Error in CMFGEN'
	    WRITE(LUER,*)'FLUX_CAL_LAM_END must be > FLUX_CAL_LAM_BEG'
	    STOP
	  END IF
	  IF( (.NOT. FLUX_CAL_ONLY) .AND. FLUX_CAL_LAM_BEG .NE. 0)THEN
	    WRITE(LUER,*)'WARNING in CMFGEN'
	    WRITE(LUER,*)'WARNING in CMFGEN'
	    WRITE(LUER,*)'FLUX_CAL_LAM_BEG is normally zero for non-FLUX'
	    WRITE(LUER,*)'calculations:'
	  END IF
	  GLOBAL_LINE_SWITCH='NONE'
	  T1=SPEED_OF_LIGHT()*1.0D-07
	  DO I=1,N_LINE_FREQ
	    IF(T1/VEC_FREQ(I) .GE. FLUX_CAL_LAM_END)THEN
	      VEC_TRANS_TYPE(I)='SOB'
	    END IF
	    IF(T1/VEC_FREQ(I) .LE. FLUX_CAL_LAM_BEG)THEN
	      VEC_TRANS_TYPE(I)='SOB'
	    END IF
	  END DO
	END IF
!
!
! Read in or calculate continuum frequency values.
!
	IF(RD_CONT_FREQ .OR. IMPURITY_CODE)THEN
	  CALL GEN_ASCI_OPEN(LUIN,'CFDAT','OLD',' ','READ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening CFDAT in CMFGEN, IOS=',IOS
	      STOP
	    END IF
	    TEMP_CHAR=' '
	    DO WHILE(INDEX(TEMP_CHAR,'!Number of continuum frequencies') 
	1                                              .EQ. 0)
	      READ(LUIN,'(A)',IOSTAT=IOS)TEMP_CHAR
	      IF(IOS .NE. 0)THEN
                WRITE(LUER,*)'Error reading in number of continuum ',
	1                        'frequencies in CMFGEN'
	        STOP
	      END IF
	    END DO
	    READ(TEMP_CHAR,*)NCF
	    IF(NCF .GT. NCF_MAX)THEN
	      WRITE(LUER,*)'Error - NCF > NCF_MAX in CMFGEN'
	      STOP
	    END IF
	    READ(LUIN,*)(NU(I),I=1,NCF)
	  CLOSE(UNIT=LUIN)
!
! Now need to set continuum frequencies.
!
	ELSE
	  FIRST=.TRUE.		!Check cross section at edge is non-zero.
	  NCF=0 		!Initialize number of continuum frequencies.
!
	  DO ID=1,NUM_IONS-1
	    CALL SET_EDGE_FREQ_V3(ID,OBS,NCF,NCF_MAX,
	1            ATM(ID)%EDGEXzV_F, ATM(ID)%NXzV_F, ATM(ID)%XzV_PRES,
	1            ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV, ATM(ID)%N_XzV_PHOT)
	  END DO
!
	  IF(XRAYS)THEN
	    DO ID=1,NUM_IONS-1
	      CALL SET_X_FREQ_V2(OBS,NCF,NCF_MAX, MAX_CONT_FREQ,
	1             AT_NO(SPECIES_LNK(ID)),ATM(ID)%ZXzV,
	1             ATM(ID)%XzV_PRES, ATM(ID+1)%XzV_PRES)
	    END DO
	  END IF
!
! Now insert addition points into frequency array. WSCI is used as a
! work array - okay since of length NCF_MAX, and zeroed in QUADSE.
! OBS contains the bound-free edges - its contents are zero on
! subroutine exit. J is used as temporary variable for the number of
! frequencies transmitted to SET_CONT_FREQ. NCF is returned as the number 
! of frequency points. FQW is used a an integer array for the sorting ---
! we know it has the correct length since it is the same size as NU.
! LUIN --- Used as temporary LU (opened and closed).
!
	  J=NCF
	  IF(FREQ_GRID_OPTION .EQ. 1)THEN
	    CALL SET_CONT_FREQ(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        J,NCF,NCF_MAX,LUIN)
	  ELSE IF(FREQ_GRID_OPTION .EQ. 2)THEN
	    CALL SET_CONT_FREQ_V3(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_XRAY,NU_XRAY_END,
	1                        J,NCF,NCF_MAX,LUIN)
	  ELSE
!
! Checks if BIG_FREQ_AMP was set to old definition.
!
	    IF(BIG_FREQ_AMP .LT. 0.0D0 .OR. BIG_FREQ_AMP .GT. 1.0D0)THEN
	      BIG_FREQ_AMP=0.5D0
	    END IF
	    CALL SET_CONT_FREQ_V4(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_XRAY,NU_XRAY_END,
	1                        J,NCF,NCF_MAX,LUIN)
	  END IF
!                                             
	END IF              !End set continuum if
!
!                         
!
! We have found all lines. If we are doing a blanketing calculation for this
! line we insert them into the continuum frequency set, otherwise the
! line is not included.
!
	DO ML=1,NCF                                          
	  FQW(ML)=NU(ML)	!FQW has temporary storage of continuum freq.
	END DO                    
	V_DOP=MINVAL(VEC_VDOP_MIN) !12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
	CALL INS_LINE_V6(  NU,LINES_THIS_FREQ,I,NCF_MAX,
	1		  VEC_FREQ,VEC_STRT_FREQ,VEC_VDOP_MIN,VEC_TRANS_TYPE,
	1                 LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,
	1                 N_LINE_FREQ,FQW,NCF,FRAC_DOP,VINF,
	1                 dV_CMF_PROF,dV_CMF_WING,
	1                 ES_WING_EXT,R_CMF_WING_EXT,L_FALSE )
!
	K=NCF		!# of continuum frequencies: Need for DET_MAIN...
	NCF=I		!Revised
	CALL DET_MAIN_CONT_FREQ(NU,NCF,FQW,K,NU_EVAL_CONT,
	1         V_DOP,DELV_CONT,COMPUTE_ALL_CROSS)
!
	CALL GET_PROFILE_STORAGE_LIMITS(NLINES_PROF_STORE,NFREQ_PROF_STORE,
	1         LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,PROF_TYPE,N_LINE_FREQ,NCF)
	CALL INIT_PROF_MODULE(ND,NLINES_PROF_STORE,NFREQ_PROF_STORE)
!
!	DO ML=1,N_LINE_FREQ
!	  T1=2.998D+05*(VEC_STRT_FREQ(ML)/VEC_FREQ(ML)-1)/VEC_VDOP_MIN(ML)
!	  T2=2.988D+05*( VEC_FREQ(ML)/NU(LINE_END_INDX_IN_NU(ML))-1 )/VEC_VDOP_MIN(ML)
!	  WRITE(157,'(I8,2X,A10,4I8,3F7.2,3ES14.6,)')ML,VEC_SPEC(ML),VEC_NL(NL),VEC_NUP(NL),
!	1                  LINE_ST_INDX_IN_NU(ML),LINE_END_INDX_IN_NU(ML),T1,T2,VEC_VDOP_MIN(ML),
!	1                  VEC_FREQ(ML),VEC_STRT_FREQ(ML),NU(LINE_END_INDX_IN_NU(ML))
!	END DO
!
	WRITE(LUER,'(A,T40,I7)')' Number of line frequencies is:',N_LINE_FREQ
	WRITE(LUER,'(A,T40,I7)')' Number of frequencies is:',NCF
	WRITE(LUER,*)' '
!
! Redefine frequency quadrature weights.
!
	IF(FREQ_GRID_OPTION.EQ. 1)THEN
	  CALL SMPTRP(NU,FQW,NCF)
	ELSE
	  CALL TRAPUNEQ(NU,FQW,NCF)
	END IF
	DO ML=1,NCF                                           
	  FQW(ML)=FQW(ML)*1.0D+15
	END DO
!
! Revise number of lines so only those in frequency grid are included.
!
	I=N_LINE_FREQ
	DO ML=I,1,-1
	  N_LINE_FREQ=ML
	  IF(VEC_FREQ(ML) .GT. MIN_CONT_FREQ*(1.0D0+1.1D0*VINF/2.998D+05))EXIT
	END DO
	I=I-N_LINE_FREQ
	IF(I .NE. 0)THEN
	  WRITE(LUER,*)'Warning from SET_FREQUENCY_GRID'
	  WRITE(LUER,'(1X,I5,A,A)')I,' weak lines in ',
	1        'extreme IR will be ignored as outside continuum range.'
	  WRITE(LUER,*)'Min(Nu_CONT)=',NU(NCF)
	  WRITE(LUER,*)'Min(Nu_LINE)=',VEC_FREQ(I+N_LINE_FREQ)
	END IF
!
! Set observers frequencies. The slight fiddling in setting NU_MAX and NU_MIN 
! is done so that the CMF frequencies encompass all observers frame 
! frequencies. This allows computation of all observers fluxes allowing for 
! velocity effects.
!
! We insert lines into the observers frame frequencies if they are being
! treated in blanketed mode, or if SOBN_FREQ_IN_OBS is set to TRUE.
!
	NU_MAX_OBS=NU(3)		!3 Ensures OBS_FREQ contained inside CMF band.
	T1=NU(NCF)*(1.0D0+2.0D0*VINF/2.998D+05)
	T2=MAXVAL(VTURB_VEC)
	NU_MIN_OBS=MAX(NU(NCF-3),T1)
	CALL INS_LINE_OBS_V5(OBS_FREQ,N_OBS,NCF_MAX,
	1               VEC_FREQ,VEC_STRT_FREQ,VEC_VDOP_MIN,VEC_TRANS_TYPE,
	1               N_LINE_FREQ,SOB_FREQ_IN_OBS,
	1		NU_MAX_OBS,NU_MIN_OBS,VINF,
	1               1.0D0,dV_OBS_PROF,dV_OBS_WING,dV_OBS_BIG,
	1               OBS_PRO_EXT_RAT,ES_WING_EXT,T2)
!
	RETURN
	END
