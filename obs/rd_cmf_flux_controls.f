	SUBROUTINE RD_CMF_FLUX_CONTROLS(ND,LUMOD,LUER)
	USE MOD_CMF_OBS
	USE CMF_FLUX_CNTRL_VAR_MOD
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER LUMOD
	INTEGER LUER
!
! Local variables.
!
	CHARACTER(LEN=80)  TMP_STRING
	CHARACTER(LEN=20)  TMP_KEY
	CHARACTER(LEN=132) TEMP_CHAR
!
	INTEGER I
	INTEGER ISPEC
	INTEGER ID
!
!
!
! All options are now read into store in CMF_FLUX, as some options are needed in that
! routine.
!
	  CALL RD_STORE_DBLE(MIN_CONT_FREQ,'MIN_CF',L_TRUE,
	1            'Minimum continuum frequency if calculating NU')
	  CALL RD_STORE_DBLE(MAX_CONT_FREQ,'MAX_CF',L_TRUE,
	1            'Maximum continuum frequency if calculating NU')
	  CALL RD_STORE_DBLE(SMALL_FREQ_RAT,'FRAC_SP',L_TRUE,
	1            'Fractional spacing for small frequencies')
	  CALL RD_STORE_DBLE(BIG_FREQ_AMP,'AMP_FAC',L_TRUE,
	1            'Amplification factor for large frequency ranges')
	  CALL RD_STORE_DBLE(dFREQ_bf_MAX,'MAX_BF',L_TRUE,
	1            'Maximum frequency spacing close to bf edge')
!
	  CALL RD_STORE_LOG(DO_LEV_DISSOLUTION,'DO_DIS',L_TRUE,
	1            'Allow for level dissolution of upper levels?')
	  CALL RD_STORE_DBLE(dV_LEV_DIS,'dV_LEV',L_TRUE,
	1             'Spacing (in km/s) on low side of bf edge for'//
	1             ' level dissolution')
	  CALL RD_STORE_DBLE(AMP_DIS,'AMP_DIS',L_TRUE,
	1            'Amplification factor on low side bf edge')
	  CALL RD_STORE_DBLE(MIN_FREQ_LEV_DIS,'MIN_DIS',L_TRUE,
	1            'Minimum frequency for level dissolution')
!
	  CALL RD_STORE_LOG(COMPUTE_ALL_CROSS,'CROSS',L_TRUE,
	1            'Compute all photoionization cross-sections?')
	  CALL RD_STORE_DBLE(DELV_CONT,'V_CROSS',L_TRUE,
	1            'Max. vel. sep. (km/s) between evaluations of all'//
	1            '  phot. cross-sections?')
!
	  CALL RD_STORE_LOG(DIF,'DIF',L_TRUE,
	1            'Use Diffusion approximation at inner boundary ?')
	  IF(DIF)THEN
	    INNER_BND_METH='DIFFUSION'
	    CALL RD_STORE_CHAR(INNER_BND_METH,'IB_METH',L_FALSE,
	1           'Inner boundary method (DIFFUSION, HOLLOW, or ZERO_FLUX)')
	    IF(INNER_BND_METH .NE. 'DIFFUSION')THEN
	      WRITE(LUER,*)'Error in RD_CONTROL_VARIABLES'
	      WRITE(LUER,*)'Inconsistent inner boundary condition'
	      WRITE(LUER,*)'DIF option and IB_METH are inconsistent'
	      STOP
	    END IF
	  ELSE
	    CALL RD_STORE_CHAR(INNER_BND_METH,'IB_METH',L_TRUE,
	1           'Inner boundary method (DIFUSION, HOLLOW, or ZERO_FLUX)')
	  END IF
	  IB_STAB_FACTOR=0.1D0
	  IF(INNER_BND_METH .EQ. 'DIFFUSION')IB_STAB_FACTOR=0.0D0
	  CALL RD_STORE_DBLE(IB_STAB_FACTOR,'IB_STAB',L_FALSE,'Inner boundary stabilization factor')
	  OUTER_BND_METH='HONJ'
	  CALL RD_STORE_CHAR(OUTER_BND_METH,'OB_METH',L_FALSE,
	1        'Outer boundary method (HONJ or HALF_MOM)')
!
	  CALL RD_STORE_INT(NUM_ES_ITERATIONS,'NUM_ES',L_TRUE,
	1            'Number of electron scattering iterations?')
	  CALL RD_STORE_LOG(RD_COHERENT_ES,'COH_ES',L_TRUE,
	1            'Assume coherent electron scattering? ')
	  CALL RD_STORE_LOG(USE_OLDJ_FOR_ES,'OLD_J',L_TRUE,
	1            'Use old file to provide initial estimate of J_ES?')
	  COHERENT_ES=RD_COHERENT_ES
!
	  CALL RD_STORE_NCHAR(METHOD,'METHOD',ISIX,L_TRUE,
	1           'Which method for continuum tau'//
	1          ' loglog, loglin, linear or zero ?')
	  CALL RD_STORE_NCHAR(N_TYPE,'N_TYPE',ISIX,L_TRUE,
	1           'Method for to handle N for MOM_J_CMF -- '//
	1           'N_ON_J, MIXED, or G_ONLY')
	  CALL RD_STORE_NCHAR(FG_SOL_OPTIONS,'FG_OPT',ITEN,L_TRUE,
	1           'Solution options for FG_J_CMF: DIFF/INS and INT/INS')
	  CALL RD_STORE_LOG(RDTHK_CONT,'THK_CONT',L_TRUE,'Use thick boundary condition for continuum ? ')
	  REXT_FAC=0.0D0
	  CALL RD_STORE_DBLE(REXT_FAC,'REXT_FAC',L_FALSE,'Factor to extend R by for thick continuum solution')
	  CALL RD_STORE_LOG(TRAPFORJ,'TRAP_J',L_TRUE,
	1           'Use trapazoidal weights to compute J? ')
!
	  WRITE(LUMOD,'()')
	  TURB_LAW='LAW_V1'
	  CALL RD_STORE_NCHAR(TURB_LAW,'TURB_LAW',ITEN,L_FALSE,
	1      'Turbulent velocity law: LAW_V1, LAW_TAU1')
	  CALL RD_STORE_DBLE(VTURB_FIX,'VTURB_FIX',L_TRUE,
	1      'Doppler velocity for DOP_FIX Doppler profiles (km/s)')
	  CALL RD_STORE_DBLE(VTURB_MIN,'VTURB_MIN',L_TRUE,
	1      'Minimum turbulent velocity for Doppler profile (km/s)')
	  CALL RD_STORE_DBLE(VTURB_MAX,'VTURB_MAX',L_TRUE,
	1      'Maximum turbulent velocity for Doppler profile (km/s)')
!
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_DBLE(MAX_DOP,'MAX_DOP',L_TRUE,
	1      'Maximum half-width of resonance zone (in Doppler widths)')
	  CALL RD_STORE_DBLE(FRAC_DOP,'FRAC_DOP',L_TRUE,
	1      'Spacing in resonance zone (in Doppler widths)')
	  CALL RD_STORE_DBLE(dV_CMF_PROF,'dV_CMF_PROF',L_TRUE,
	1      'Spacing across cmf profile (in km/s)')
	  CALL RD_STORE_DBLE(dV_CMF_WING,'dV_CMF_WING',L_TRUE,
	1      'Spacing across e.s. wings of cmf profile(in km/s)')
	  CALL RD_STORE_DBLE(ES_WING_EXT,'ES_WING_EXT',L_TRUE,
	1      'Extent of BLUE e.s. wings from resonance core (in km/s)')
	  CALL RD_STORE_DBLE(R_CMF_WING_EXT,'R_CMF_WING_EXT',L_TRUE,
	1      'Extent of RED e.s. wings from RESONANCE core (in Vinf)')
!
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_DBLE(OBS_PRO_EXT_RAT,'OBS_EXT_RAT',L_TRUE,
	1      'Half width of profile in Vinf.')
	  CALL RD_STORE_DBLE(FRAC_DOP_OBS,'FRAC_DOP_OBS',L_TRUE,
	1      'Spacing across intrinsic profile zone (in Doppler widths)')
	  CALL RD_STORE_DBLE(dV_OBS_PROF,'dV_OBS_PROF',L_TRUE,
	1      'Spacing across observed profile (in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_WING,'dV_OBS_WING',L_TRUE,
	1      'Spacing across e.s. wings of observed profile(in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_BIG,'dV_OBS_BIG',L_TRUE,
	1      'Frequency spacing between lines (in km/s)')
!
	  CALL RD_STORE_DBLE(OBS_TAU_MAX,'TAU_MAX',L_TRUE,
	1    'Optical depth at which observers frame integration is terminated')
	  CALL RD_STORE_DBLE(OBS_ES_DTAU,'ES_DTAU',L_TRUE,
	1      'Maximum increments in e.s. optical depth scale')
	  CALL RD_STORE_NCHAR(OBS_INT_METHOD,'INT_METH',ITEN,L_TRUE,
	1            'Integration method for computing I along ray')
	  CALL RD_STORE_INT(N_INS_OBS,'N_INS_OBS',L_TRUE,
	1          'Mininum number of points to be inserted in '//
	1          ' observers frame grid (>= 0)')
!
	  REVISE_P_GRID=.FALSE.
	  CALL RD_STORE_LOG(REVISE_P_GRID,'REVISE_P',L_FALSE,
	1          'Revise p frid for observer''s calculation')
!
	  CALL RD_STORE_LOG(DO_REL_IN_OBSFRAME,'DO_RELO',L_FALSE,
	1        'Use all relativistic terms in Observer''s frame calculation.')
	  CALL RD_STORE_LOG(DO_CMF_REL_OBS,'DO_CMF_RELO',L_FALSE,
	1        'Use all relativistic terms in CMF Observer''s frame calculation.')
!
	  WRITE(LUMOD,'()')
	  USE_FIXED_J=.FALSE.
	  CALL RD_STORE_LOG(USE_FIXED_J,'USE_FIXED_J',L_FALSE,
	1           'Use previously computed J to evaluate ALL rates?')
	  CALL RD_STORE_LOG(FLUX_CAL_ONLY,'FLUX_CAL_ONLY',L_TRUE,
	1           'Compute the observers frame flux only ?')
	  CALL RD_STORE_LOG(EXTEND_FRM_SOL,'EXT_FRM_SOL',L_TRUE,
	1           'Extrapolate the formal solution to larger radii?')
	  CALL RD_STORE_LOG(INSERT_FREQ_FRM_SOL,'INS_F_FRM_SOL',L_TRUE,
	1           'Extrapolate the formal solution to larger radii?')
	  CALL RD_STORE_NCHAR(CMF_FORM_OPTIONS,'FRM_OPT',ITEN,L_TRUE,
	1           'Solution options for CMF_FORM_SOL')
	  CALL RD_STORE_LOG(DO_SOBOLEV_LINES,'DO_SOB_LINES',L_TRUE,
	1        'Compute Sobolev EWs?')
	  CALL RD_STORE_DBLE(EW_CUT_OFF,'EW_CUT',L_TRUE,
	1        'Output EW info only if ABS(EW) > EW_CUT')
	  CALL RD_STORE_LOG(SOB_FREQ_IN_OBS,'SOB_FREQ_IN_OBS',L_TRUE,
	1        ' Allow for SOB & CMF lines in defining observers'//
	1        ' frequencies?')
!
	  COMPUTE_J=.TRUE.
	  CALL RD_STORE_LOG(COMPUTE_J,'COMP_J',L_FALSE,'Compute the radiation field')
	  CALL RD_STORE_LOG(WRITE_ETA_AND_CHI,'WR_ETA',L_TRUE,'Output ETA and CHI? ')
	  CALL RD_STORE_LOG(WRITE_FLUX,'WR_FLUX',L_TRUE, 'Output Flux as a function of depth? ')
	  CALL RD_STORE_LOG(WRITE_CMF_FORCE,'WR_CMF_FORCE',L_TRUE,
	1        'Output CMF line-force multiplier as a function of depth? ')
	  CALL RD_STORE_LOG(WRITE_SOB_FORCE,'WR_SOB_FORCE',L_TRUE,
	1        'Output SOBOLEV line-force multiplier as a function of depth? ')
	  WR_ION_LINE_FORCE=.FALSE.
	  CALL RD_STORE_LOG(WR_ION_LINE_FORCE,'WR_ION_FORCE',L_FALSE,
	1        'Output line-force multiplier as a function of ion & depth? ')
	  CALL RD_STORE_LOG(WRITE_IP,'WR_IP',L_TRUE,
	1        'Output I as a functio of p and frequency?')
	  WRITE_RTAU=.FALSE.
	  CALL RD_STORE_LOG(WRITE_RTAU,'WR_RTAU',L_FALSE,
	1        'Output R(Tau=Tau_ref) as a function of p and frequency?')
	  WRITE_dFR=.FALSE.
	  CALL RD_STORE_LOG(WRITE_dFR,'WR_dFR',L_FALSE,
	1        'Output dFR as a functioin of R and frequency?')
	  CALL RD_STORE_DBLE(TAU_REF,'TAU_REF',WRITE_RTAU,
	1        'Reference tau for WR_TAU')
!
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_NCHAR(GLOBAL_LINE_SWITCH,'GLOBAL_LINE',ITEN,L_TRUE,
	1            'Global switch to indicate handeling of line')
	  CALL SET_CASE_UP(GLOBAL_LINE_SWITCH,IZERO,IZERO)
	  IF( GLOBAL_LINE_SWITCH(1:3) .NE. 'SOB' .AND.
	1       GLOBAL_LINE_SWITCH(1:3) .NE. 'CMF' .AND.
	1       GLOBAL_LINE_SWITCH(1:4) .NE. 'NONE' .AND.
	1       GLOBAL_LINE_SWITCH(1:5) .NE. 'BLANK')THEN
	    WRITE(LUER,*)'Invalid GLOBAL_LINE SWITCH parameter'
	    STOP
	  END IF
	  CALL RD_STORE_LOG(SET_TRANS_TYPE_BY_LAM,'LAM_SET',L_TRUE,
	1         'Set long wavelengths to SOBOLEV approximation')
	  CALL RD_STORE_DBLE(FLUX_CAL_LAM_BEG,'F_LAM_BEG',L_TRUE,
	1         'Inital wavelength (A) for blanketed flux calculation')
	  CALL RD_STORE_DBLE(FLUX_CAL_LAM_END,'F_LAM_END',L_TRUE,
	1         'Final wavelength (A) for blanketed flux calculation')
!
	  CALL RD_STORE_LOG(THK_LINE,'THK_LINE',L_TRUE,
	1           'Use thick boundary condition for lines?')
	  CALL RD_STORE_LOG(CHECK_LINE_OPAC,'CHK_L_POS',L_TRUE,
	1      'Ensure Line opacity is positive ?')
	  CALL RD_STORE_NCHAR(NEG_OPAC_OPTION,'NEG_OPAC_OPT',ITEN,L_TRUE,
	1            'Method for negative opacities in BLANKETING mode')
	  CALL SET_CASE_UP(NEG_OPAC_OPTION,IZERO,IZERO)
	  IF(NEG_OPAC_OPTION .NE. 'SRCE_CHK' .AND.
	1                           NEG_OPAC_OPTION .NE. 'ESEC_CHK')THEN
	    WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	    WRITE(LUER,*)'Invalid NEG_OPAC_OPTION'
	    WRITE(LUER,*)'Valid options are SRCE_CHK and ESEC_CHK'
	    STOP
	  END IF
!
	  CALL RD_STORE_LOG(INCL_TWO_PHOT,'INC_TWO',L_TRUE,'Include two photon transitions?')
	  TWO_PHOTON_METHOD='USE_RAD'
          CALL RD_STORE_CHAR(TWO_PHOTON_METHOD,'TWO_METH',L_FALSE,'USE_RAD, LTE, NOSTIM or OLD_DEFAULT')
	  CALL RD_STORE_LOG(INCL_RAY_SCAT,'INC_RAY',L_TRUE,'Include Rayeligh scattering?')
!
	  CALL RD_STORE_LOG(XRAYS,'INC_XRAYS',L_TRUE,'Include X-ray emission')
	  CALL RD_STORE_LOG(FF_XRAYS,'FF_XRAYS',XRAYS,'Use free-free processes to compute X-ray emission')
	  CALL RD_STORE_LOG(XRAY_SMOOTH_WIND,'X_SM_WIND',XRAYS,'Ignore clumping when computing X-ray emission')
!
	  VSMOOTH_XRAYS=3000.0D0
	  CALL RD_STORE_DBLE(VSMOOTH_XRAYS,'VS_XRAYS',XRAYS,'X-ray smoothing width for SOB/CMF options')
!
	  FILL_FAC_XRAYS_1=0.D0
	  FILL_FAC_XRAYS_2=0.D0
	  CALL RD_STORE_DBLE(FILL_FAC_XRAYS_1,'FIL_FAC_1',XRAYS,
	1           'Filling factor for X-ray emission [1]')
	  CALL RD_STORE_DBLE(T_SHOCK_1,'T_SHOCK_1',XRAYS,
	1           'Shock T for X-ray emission [1]')
	  CALL RD_STORE_DBLE(V_SHOCK_1,'V_SHOCK_1',XRAYS,
	1           'Cut off velocity for X-ray emission [1]')
	  CALL RD_STORE_DBLE(FILL_FAC_XRAYS_2,'FIL_FAC_2',XRAYS,
	1           'Filling factor for X-ray emission [2]')
	  CALL RD_STORE_DBLE(T_SHOCK_2,'T_SHOCK_2',XRAYS,
	1           'Shock T for X-ray emission [2]')
	  CALL RD_STORE_DBLE(V_SHOCK_2,'V_SHOCK_2',XRAYS,
	1           'Cut off velocity for X-ray emission [2]')
!
! If we add _SPEC, SOB or BLANK is the default option. In this case we need only
! specify the species we wish to change from the default.
!
	  IF(GLOBAL_LINE_SWITCH .EQ. 'SOB_SPEC')THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
                ATM(ID)%XzV_TRANS_TYPE='SOB'
	        TEMP_CHAR='TRANS_'//ION_ID(ID)
	        TMP_STRING='Method for treating '//TRIM(ION_ID(ID))//' lines?'
	        CALL RD_STORE_NCHAR( ATM(ID)%XzV_TRANS_TYPE,TEMP_CHAR,ITEN,L_FALSE,TMP_STRING)
	      END IF
	    END DO
	    GLOBAL_LINE_SWITCH='NONE'
	  ELSE IF(GLOBAL_LINE_SWITCH .EQ. 'BLANK_SPEC')THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        ATM(ID)%XzV_TRANS_TYPE='BLANK'
	        TEMP_CHAR='TRANS_'//ION_ID(ID)
	        TMP_STRING='Method for treating '//TRIM(ION_ID(ID))//' lines?'
	        CALL RD_STORE_NCHAR( ATM(ID)%XzV_TRANS_TYPE,TEMP_CHAR,ITEN,L_FALSE,TMP_STRING)
	      END IF
	    END DO
	    GLOBAL_LINE_SWITCH='NONE'
	  ELSE IF(GLOBAL_LINE_SWITCH .EQ. 'NONE')THEN
	    WRITE(LUMOD,'()')
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        TEMP_CHAR='TRANS_'//ION_ID(ID)
	        TMP_STRING='Method for treating '//TRIM(ION_ID(ID))//' lines?'
	        CALL RD_STORE_NCHAR( ATM(ID)%XzV_TRANS_TYPE,TEMP_CHAR,ITEN,L_TRUE,TMP_STRING)
	      END IF
	    END DO
	  END IF
!
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_NCHAR(GLOBAL_LINE_PROF,'GLOBAL_PROF',ITEN,L_TRUE,
	1        'Global switch for intrinsic line absorption profile')
	  CALL SET_CASE_UP(GLOBAL_LINE_PROF,IZERO,IZERO)
	  IF( GLOBAL_LINE_PROF .NE. 'NONE' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOP_FIX' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOPPLER' .AND.
	1       GLOBAL_LINE_PROF .NE. 'LIST' .AND.
	1       GLOBAL_LINE_PROF .NE. 'LIST_VGT' .AND.
	1       GLOBAL_LINE_PROF .NE. 'VOIGT' .AND.
	1       GLOBAL_LINE_PROF .NE. 'HZ_STARK')THEN
	    WRITE(LUER,*)'Invalid GLOBAL_LINE_PROF parameter'
	    STOP
	  END IF
	  IF(GLOBAL_LINE_PROF .EQ. 'DOP_FIX')THEN
            CALL RD_STORE_DBLE(TDOP,'TDOP',L_TRUE,'Temperature to be used in Doppler profile (10^4K)')
            CALL RD_STORE_DBLE(AMASS_DOP,'AMASS_DOP',L_TRUE,'Atomic mass to be used in Doppler profile (amu''s)')
	  ELSE
	    TDOP=2.0D0; AMASS_DOP=1.0D+06
	  END IF
	  CALL RD_STORE_LOG(SET_PROF_LIMS_BY_OPACITY,'OPAC_LIMS',L_TRUE,
	1           'Set prof limits by line to cont. ratio?')
	  CALL RD_STORE_DBLE(DOP_PROF_LIMIT,'DOP_LIM',L_TRUE,
	1           'Edge limits for Doppler line profile')
	  CALL RD_STORE_DBLE(VOIGT_PROF_LIMIT,'VOIGT_LIM',L_TRUE,
	1           'Edge limits for Voigt line profile')
!
	  MAX_PROF_ED=1.0D+16
	  NORM_PROFILE=.FALSE.
	  V_PROF_LIMIT=5000.0D0
	  CALL RD_STORE_DBLE(MAX_PROF_ED,'MAX_PROF_ED',L_FALSE,
	1           'Maximum electron density for Stark profile computation')
	  CALL RD_STORE_DBLE(V_PROF_LIMIT,'V_PROF_LIM',L_FALSE,
	1           'One-sided profile limit for Stark profiles (km/s)')
	  CALL RD_STORE_LOG(NORM_PROFILE,'NORM_PROF',L_FALSE,
	1           'When true, profiles are normalized to have unit area.')
!
! Verify validity of profile option. We also check whether we need to leed
! in the file which links certain types of profiles to individual lines.
!
	  IF(GLOBAL_LINE_PROF .EQ. 'NONE')THEN
	    WRITE(LUMOD,'()')
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        TEMP_CHAR='PROF_'//ION_ID(ID)
	        TMP_STRING='Intrinsic profile for treating '//TRIM(ION_ID(ID))//' lines?'
	        CALL RD_STORE_NCHAR(ATM(ID)%XzV_PROF_TYPE,TEMP_CHAR,ITEN,
	1               L_TRUE,TMP_STRING)
	        IF( ATM(ID)%XzV_PROF_TYPE .NE. 'DOPPLER' .AND.
	1           ATM(ID)%XzV_PROF_TYPE .NE. 'VOIGT' .AND.
	1           ATM(ID)%XzV_PROF_TYPE .NE. 'LIST' .AND.
	1           ATM(ID)%XzV_PROF_TYPE .NE. 'LIST_VGT' .AND.
	1           ATM(ID)%XzV_PROF_TYPE .NE. 'HZ_STARK')THEN
	          WRITE(LUER,*)'Invalid ATM(ID)%XzV_PROF_TYPE SWITCH parameter'
	          STOP
	        END IF
	        IF(ATM(ID)%XzV_PROF_TYPE(1:4) .EQ. 'LIST')RD_STARK_FILE=.TRUE.
	      END IF
	    END DO
	  END IF
!
	  WRITE(LUMOD,'()')
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_LOG(EDD_CONT,'JC_W_EDD',L_TRUE,
	1        'Compute continuum intensity using Eddington factors')
	  CALL RD_STORE_LOG(EDD_LINECONT,'JBAR_W_EDD',L_TRUE,
	1    'Compute line continuum intensity using Eddington factors')
	  INCL_INCID_RAD=.FALSE.
	  PLANE_PARALLEL_NO_V=.FALSE.
	  CALL RD_STORE_LOG(PLANE_PARALLEL_NO_V,'PP_NOV',L_FALSE,
	1    'Plane-paralle geometry WITHOUT velocity field?')
	  PLANE_PARALLEL=.FALSE.
	  CALL RD_STORE_LOG(PLANE_PARALLEL,'PP_MOD',L_FALSE,
	1    'Plane-paralle geometry with velocity field?')
	  USE_J_REL=.FALSE.
	  USE_FORMAL_REL=.FALSE.
	  INCL_REL_TERMS=.FALSE.
	  INCL_ADVEC_TERMS_IN_TRANS_EQ=.FALSE.
	  USE_LAM_ES=.FALSE.
	  CALL RD_STORE_LOG(USE_J_REL,'USE_J_REL',L_FALSE,'Use relativistic moment solver?')
	  IF( USE_J_REL)THEN
	    INCL_REL_TERMS=.TRUE.
	    INCL_ADVEC_TERMS_IN_TRANS_EQ=.TRUE.
	  END IF
	  CALL RD_STORE_LOG(INCL_REL_TERMS,'INCL_REL',L_FALSE,'Include relativistic terms?')
	  CALL RD_STORE_LOG(INCL_ADVEC_TERMS_IN_TRANS_EQ,'INCL_ADV_TRANS',L_FALSE,
	1    'Include advection terms in transfer equation?')
	  INCL_DJDT_TERMS=.FALSE.
	  CALL RD_STORE_LOG(INCL_DJDT_TERMS,'INCL_DJDT',L_FALSE,'DJDt terms in transfer equaton for SN models?')
	  IF(INCL_DJDT_TERMS)THEN
	    USE_DJDT_RTE=.TRUE.
	    USE_Dr4JDT=.TRUE.
	    CALL RD_STORE_LOG(USE_Dr4JDT,'USE_DR4JDT',L_FALSE,'Difference Dr4JDt')
	  ELSE
	    USE_DJDT_RTE=.FALSE.
	    USE_Dr4JDT=.FALSE.
	    CALL RD_STORE_LOG(USE_DJDT_RTE,'USE_DJDT_RTE',L_FALSE,
	1     'Use solver which has DJDt terms in transfer equaton for SN models?')
	  END IF
	  DJDT_RELAX_PARAM=1.0D0
	  CALL RD_STORE_DBLE(DJDT_RELAX_PARAM,'DJDT_RELAX',L_FALSE,
	1          'Factor to scale DJDT terms to assist initial convergence')
	  IF(USE_DJDT_RTE .OR. USE_J_REL)USE_FORMAL_REL=.TRUE.
	  CALL RD_STORE_LOG(USE_FORMAL_REL,'USE_FRM_REL',L_FALSE,'Use CMF_FORMAL_REL to compute F etc?')
	  CALL RD_STORE_LOG(USE_LAM_ES,'USE_LAM_ES',L_FALSE,'Use formal solution for e.s. (done via Lambda iteration)')

	  WRITE(LUMOD,'()')
	  CALL RD_STORE_LOG(ACCURATE,'INC_GRID',L_TRUE,'Increase grid size to improve accuracy? ')
	  CALL RD_STORE_LOG(ALL_FREQ,'ALL_FREQ',L_TRUE,'Increase accuracy for all frequencies?')
	  CALL RD_STORE_DBLE(ACC_FREQ_END,'ACC_END',L_TRUE,'Increase accuracy for all frequencies > ACC_END?')
	  CALL RD_STORE_INT(NPINS,'N_INS',L_TRUE,'Number of points to be inserted in higher'//
	1          ' accuracy grid (1, 2 or 3) ')
	  ST_INTERP_INDX=1; ST_INTERP_INDX=ND; DEEP=10
	  CALL RD_STORE_INT(ST_INTERP_INDX,'ST_INT',L_FALSE,'Interpolate from ? ')
	  CALL RD_STORE_INT(END_INTERP_INDX,'END_INT',L_FALSE,'Interpolate to ? ')
	  CALL RD_STORE_INT(DEEP,'ND_QUAD',L_FALSE,'Quadratic interpolation from ND-? to ND')
	  CALL RD_STORE_NCHAR(INTERP_TYPE,'INTERP_TYPE',10,L_TRUE,
	1         'Perform interpolations in LOG or LIN plane')
!
	  CALL RD_STORE_DBLE(DELV_FRAC_FG,'DELV_FG',L_TRUE,
	1         'Maximum velocity separation (Doppler widths) for FG_J_CMF')
	  CALL RD_STORE_DBLE(DELV_FRAC_MOM,'DELV_MOM',L_TRUE,
	1         'Maximum velocity separation (Doppler widths) for MOM_J_CMF')
!
! Next two variables apply for both ACCURATE and EDDINGTON.
!
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_LOG(COMPUTE_EDDFAC,'COMP_F',L_TRUE,
	1      'Compute new Eddington factors (f)')
	  CALL RD_STORE_DBLE(ACC_EDD_FAC,'ACC_F',L_TRUE,
	1      'Accuracy with which to compute the eddington factor f')
!
	DO ISPEC=1,NUM_SPECIES
	  TMP_KEY='SCL_'//TRIM(SPECIES(ISPEC))//'_ABUND'
	  CALL RD_STORE_DBLE(ABUND_SCALE_FAC(ISPEC),TMP_KEY,L_FALSE,'Factor to scale abundance by')
	END DO
!
	STOP_IF_BAD_PARAM=.TRUE.
	CALL RD_STORE_LOG(STOP_IF_BAD_PARAM,'STOP_IF_BP',L_FALSE,'Undo corrections at last 5 depths')
!
! Memory and options in STORE are no longer required.
!
	CLOSE(UNIT=LUMOD)
	CALL CLEAN_RD_STORE
!
! Check consistency of parameters.
!
	CALL CHECK_CMF_FLUX_PARAM_CONSIS()
!
	RETURN
	END
