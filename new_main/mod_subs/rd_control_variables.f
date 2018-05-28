	SUBROUTINE RD_CONTROL_VARIABLES(LUIN,LUSCR,LUER,NUM_BNDS)
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered : 01-Sep-2016 : Added variabe MINIMUM_ISO_POP (not required)
!                            TIME_SEQ_NO changed from integer to real.
! Altered : 15-Feb-2015 : Added INSTANTANEOUS_ENERGY_DEPOSITION option (12-Jan-2015 on OSPREY[cur_cmf_gam])
! Altered : 02/07-Nov-2011 : Changed location of COMP_GREY_LST_IT so that set for all model.
! Altered : 05-Apr-2011 : Changed the way REVISE_R_GRID is handled. Variable controlling
!                           R grid revision now read in by routine do the R-grid revision.
! Altered : 16-Jul-2010 : Added FIX_ALL_SPECIES variable, and assoicated options.
! Altered : 01-Feb-2010 : GAMRAY_TRANS read inserted
! Altered : 31-Jan-2010 : INNER_BND_METH and OUTER_BND_METH options installed.
! Altered : 23-Nov-2007 : Optional LAM_SCALE_OPT variable included.
! Altered : 29-Jan-2006 : Control variables for relativistic transfer and time
!                          dependent statistical equilibrium equations installed.
! 
	INTEGER LUIN
	INTEGER LUSCR
	INTEGER LUER
	INTEGER NUM_BNDS
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
	REAL*8 ATOMIC_MASS_UNIT
	EXTERNAL ATOMIC_MASS_UNIT
!
	CALL GEN_ASCI_OPEN(LUIN,'VADAT','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening VADAT in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)
C 
C
C Input model parameters and modelling specifications.
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(RP,'RSTAR',L_TRUE,'Stellar radius (in 10^10 cm)')
	  CALL RD_STORE_DBLE(RMAX,'RMAX',L_TRUE,'Maximum radius (in R*)')
	  RMAX=RMAX*RP
	  DO_HYDRO=.FALSE.
	  CALL RD_STORE_LOG(DO_HYDRO,'DO_HYDRO',L_FALSE,'Adjust hydrostatic structure')
!
	  CALL RD_STORE_INT(VELTYPE,'VEL_LAW',L_TRUE,'Velocity Law to be used')
	  IF(VELTYPE .EQ. 1 .OR. VELTYPE .EQ. 2)THEN
	    CALL RD_STORE_DBLE(VRP,'VRP',L_TRUE,'First velocity component')
	    CALL RD_STORE_DBLE(RN,'RN',L_TRUE,' ')
	    CALL RD_STORE_DBLE(VINF,'VINF',L_TRUE,' ')
	    CALL RD_STORE_DBLE(EPPS1,'EPSS1',L_TRUE,' ')
	    CALL RD_STORE_DBLE(GAMMA1,'GAMMA1',L_TRUE,' ')
!
	    CALL RD_STORE_DBLE(RP2,'RP2',L_TRUE,'Second velocity component')
	    CALL RD_STORE_DBLE(VRP2,'VRP2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(RN2,'RN2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(VINF2,'VINF2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(EPPS2,'EPPS2',L_TRUE,' ')
	    CALL RD_STORE_DBLE(GAMMA2,'GAMMA2',L_TRUE,' ')
	    RN=RN*RP
	    RP2=RP*RP2
	    RN2=RN2*RP2
	  ELSE IF(VELTYPE .EQ. 3)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,
	1           'Core velocity (km/s)')
	    CALL RD_STORE_DBLE(VPHOT,'VPHOT',L_TRUE,
	1           'Photospheric velocity (km/s)')
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,
	1           'Terminal velocity (km/s)')
	    CALL RD_STORE_DBLE(SCL_HT,'SCL_HT',L_TRUE,
	1           'Scale Height (in R*) of photosphere')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA',L_TRUE,
	1           'Speed of velocity Law')
	    V_EPPS1=1.0D0
	    VINF=VINF1
	    VINF2=VINF1                !i.e. no 2nd component
	    V_BETA2=1.0D0
	    V_EPPS2=1.0D0
	    N_OB_INS=1                 !Old default
	    CONS_FOR_R_GRID=1.0D0
	    EXP_FOR_R_GRID=0.0D0
	  ELSE IF(VELTYPE .EQ. 4)THEN
!
! No parameters required
!
	  ELSE IF(VELTYPE .EQ. 5)THEN
	    WRITE(LUER,*)'Velocity law 5 not implemented in this version',
	1                 ' of CMFGEN'
	    STOP
	  ELSE IF(VELTYPE .EQ. 6)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,
	1           'Core velocity (km/s)')
	    CALL RD_STORE_DBLE(VPHOT,'VPHOT',L_TRUE,
	1           'Photospheric velocity (km/s)')
	    CALL RD_STORE_DBLE(SCL_HT,'SCL_HT',L_TRUE,
	1           'Scale Height (in R*) of photosphere')
	    CALL RD_STORE_DBLE(VINF1,'VINF1',L_TRUE,
	1            'Terminal velocity (km/s) if no 2nd comp.')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,
	1           'Speed of 1st Beta velocity Law')
	    CALL RD_STORE_DBLE(V_EPPS1,'EPPS1',L_TRUE,
	1           'Scale factor for 1s Beta velocity Law')
	    CALL RD_STORE_DBLE(VINF2,'VINF2',L_TRUE,
	1           'True terminal velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA2,'BETA2',L_TRUE,
	1           'Speed of 2nd Beta velocity Law')
	    CALL RD_STORE_DBLE(V_EPPS2,'EPPS2',L_TRUE,
	1           'Scale factor for 2nd Beta V law')
!
	    N_OB_INS=1          !Old default
	    N_IB_INS=2
	    CONS_FOR_R_GRID=-1.0D0
	    CALL RD_STORE_INT(N_OB_INS,'NBND_INS',L_FALSE,
	1           'Number of additional points to insert in radius grid at boundary')
	    CALL RD_STORE_DBLE(CONS_FOR_R_GRID,'C_R_GRID',L_FALSE,
	1           'Constant to allow improved shoice of R grid')
	    IF(CONS_FOR_R_GRID .GT. 0)THEN
	      CALL RD_STORE_DBLE(EXP_FOR_R_GRID,'E_R_GRID',L_TRUE,
	1           'Constant to allow improved shoice of R grid')
	    ELSE
	      CONS_FOR_R_GRID=1.0D0
	      EXP_FOR_R_GRID=0.0D0
	    END IF
C
C !Required by routines other than STARPCYG
C
	    VINF=VINF2	
	  ELSE IF(VELTYPE .EQ. 7)THEN
	    CALL RD_STORE_NCHAR(VEL_OPTION,'VEL_OPT',ITEN,L_TRUE,
	1                        'Velocity option: RVSIG_COL or deKOTER')
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,'Terminal velocity (km/s)')
	    VCORE=0.0D0		!Not used but initialized
	    VPHOT=0.0D0
	    SCL_HT=0.0D0
	    V_BETA1=0.0D0
	    V_EPPS1=1.0D0
	    VINF=VINF1
	    VINF2=VINF1		!i.e. no 2nd component
	    V_BETA2=1.0D0
	    V_EPPS2=1.0D0
	  ELSE IF(VELTYPE .EQ. 10)THEN
	    SN_MODEL=.TRUE.
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,'Initial velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,'Power of velocity Law')
	    RHO_ZERO=0.0D0
	    CALL RD_STORE_DBLE(RHO_ZERO,'RHO_ZERO',L_FALSE,'Initial density (gm/cm^3)')
	    IF(RHO_ZERO .EQ. 0.0D0)THEN
	      CALL RD_STORE_DBLE(RCUB_RHO_ZERO,'RCUB_RHO',L_TRUE,'r^3 . initial density 10^{-30} gm)')
	      RHO_ZERO=RCUB_RHO_ZERO/(RP**3)
	    END IF
	    RHO_ZERO=RHO_ZERO/ATOMIC_MASS_UNIT()
	    CALL RD_STORE_DBLE(N_RHO,'N_RHO',L_TRUE,'Density exponent (+ve)')
	    VINF=VCORE*(RMAX/RP)**V_BETA1
	  ELSE IF(VELTYPE .EQ. 11)THEN
	    CALL RD_STORE_DBLE(VINF1,'VINF',L_TRUE,'Terminal velocity (km/s)')
	    SN_MODEL=.TRUE.
	    SN_HYDRO_MODEL=.TRUE.
	  ELSE IF(VELTYPE .EQ. 12)THEN
	    CALL RD_STORE_DBLE(VCORE,'VCORE',L_TRUE,'Initial velocity (km/s)')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA1',L_TRUE,'Power of velocity Law')
	  ELSE
	    WRITE(LUER,*)'Velocity law ',VELTYPE, ' not implemented',
	1                ' in this version of CMFGEN'
	    STOP
	  END IF
!
	  VAR_MDOT=.FALSE.
	  IF(VELTYPE .EQ. 7)THEN
	    CALL RD_STORE_LOG(VAR_MDOT,'VAR_MDOT',L_FALSE,'Variable mass-loss rate model?') 
	    CALL RD_STORE_NCHAR(VAR_MDOT_FILE,'VM_FILE',ITEN,VAR_MDOT,
	1                        'File with density and clumping info')
	  END IF
C
	  IF(VAR_MDOT .OR. SN_MODEL)THEN
	    RMDOT=1.0D-20
	  ELSE
	    CALL RD_STORE_DBLE(RMDOT,'MDOT',L_TRUE,'Mass Loss rate (Msun/yr) ')
	  END IF
	  CALL RD_STORE_DBLE(LUM,'LSTAR',L_TRUE,'Stellar luminosity (Lsun)')
!
! TEFF and LOGG only need to be present if DO_HYDRO is TRUE. Values will still be read
! in when available if DO_HYDRO is FALSE. 
!
	  TEFF=0.0D0; LOGG=0.0D0
	  IF(DO_HYDRO)THEN
	    CALL RD_STORE_DBLE(TEFF,'TEFF',DO_HYDRO,'Effective temperature (10^4 K)')
	    CALL RD_STORE_DBLE(LOGG,'LOGG',DO_HYDRO,'Log surface gravity (cgs units)')
	    CALL RD_STORE_DBLE(V_BETA1,'BETA',DO_HYDRO,'Beta exponent for velocity law')
	    PRESSURE_VTURB=0.0D0
	    CALL RD_STORE_DBLE(PRESSURE_VTURB,'P_VTURB',L_FALSE,'Ipressure trubulent velocity')
	  END IF
	  IF(TEFF .NE. 0.0D0 .AND. .NOT. DO_HYDRO)THEN
	    WRITE(LUER,'(A)')' Possible Error in VADAT file.'
	    WRITE(LUER,'(A)')' You have set TEFF in VADAT but DO_HYDRO is false.'
            WRITE(LUER,'(A)')' TEFF will not be VALID unless hydro-iterations are performed.'
            WRITE(LUER,'(A)')' Remove TEFF if you are not going to do a hdyro iteration'
            WRITE(LUER,'(A)')' In this case TEFF is set by L, R and Mdot'
	    STOP
	  END IF

	  CALL RD_STORE_DBLE(STARS_MASS,'MASS',L_TRUE,'Stellar mass (Msun)')
C
C All clumping parameters are read in, even when CLUMPING is switched off.
C
	  CALL RD_STORE_LOG(DO_CLUMP_MODEL,'DO_CL',L_TRUE,
	1            'Calculate a model with clumping?')
	  CALL RD_STORE_NCHAR(CLUMP_LAW,'CL_LAW',ISIX,L_TRUE,
	1      'Which clumping law is being utilized?')
	  CALL SET_CASE_UP(CLUMP_LAW,IZERO,IZERO)
	  CALL RD_STORE_INT(N_CLUMP_PAR,'N_CL_PAR',L_TRUE,
	1          'Number of clumping parameters')
	  IF(N_CLUMP_PAR .GT. N_CLUMP_PAR_MAX)THEN
	    WRITE(LUER,*)'Error in CMFGEN'
	    WRITE(LUER,*)'N_CLUMP_PAR too large: N_CLUMP_PAR=',N_CLUMP_PAR
	    STOP
	  END IF
	  CLUMP_PAR(:)=0.0D0
	  DO I=1,N_CLUMP_PAR			!Should be less than 10
	    TEMP_CHAR='CL_PAR_'
	    WRITE(TEMP_CHAR(8:8),'(I1)')I
	    CALL RD_STORE_DBLE(CLUMP_PAR(I),TEMP_CHAR(1:8),L_TRUE,
	1             'Clumping parameters:')
	  END DO
!
	  REVISE_R_GRID=.FALSE.
	  CALL RD_STORE_LOG(REVISE_R_GRID,'REV_RGRID',L_FALSE,'Automatically revise R grid?')
!
! These calls read in parameters which automatically control the
! revision of the R grid automatically in SN models that have sharp
! ionization fronts.
!
	  JGREY_WITH_V_TERMS=.FALSE.
	  PURE_HUBBLE_FLOW=.FALSE.
	  COMP_GREY_LST_IT=.TRUE.
	  N_RG_PAR=0
	  TIME_SEQ_NO=0.0D0
	  IF(SN_MODEL)THEN
!
	    JGREY_WITH_V_TERMS=.TRUE.
	    CALL RD_STORE_LOG(JGREY_WITH_V_TERMS,'JG_W_V',L_FALSE,
	1            'Include V terms when computing Grey temperature?')
	    DO_CO_MOV_DDT=.FALSE.
	    CALL RD_STORE_LOG(DO_CO_MOV_DDT,'DO_DDT',L_FALSE,
	1            'Include comoving derivative in SE equations?')
	    CALL RD_STORE_DBLE(TIME_SEQ_NO,'TS_NO',DO_CO_MOV_DDT,
	1            'Time sequence number: 1 for inital model')
	    CALL RD_STORE_DBLE(SN_AGE_DAYS,'SN_AGE',DO_CO_MOV_DDT,'Age of SN in days')
	    CALL RD_STORE_LOG(PURE_HUBBLE_FLOW,'PURE_HUB',L_TRUE,
	1            'Force a pre hubble flow using age and radii of SN?')
	    CALL RD_STORE_LOG(INCL_RADIOACTIVE_DECAY,'INC_RAD_DECAYS',L_TRUE,
	1            'Allow for radiactive decays')
!
! Control parameters for handling non-thermal electrons.
!
	    CALL RD_STORE_LOG(TREAT_NON_THERMAL_ELECTRONS,'TRT_NON_TE',L_TRUE,
	1            'Treat non-thermal electrons')
	    SCL_NT_CROSEC=.FALSE.
	    SCL_NT_ION_CROSEC=.FALSE.
	    NT_OMIT_ION_SCALE=1.0D-03
	    NT_OMIT_LEV_SCALE=1.0D-04
	    NT_NKT=1000                 !Number of energy bins
	    NT_EMAX=1000.0D0            !eV
	    NT_EMIN=1.0D0               !ev
	    NON_THERMAL_IT_CNTRL=1
	    NT_SOURCE_TYPE='BELL_SHAPE'
	    IF(TREAT_NON_THERMAL_ELECTRONS)THEN
	      CALL RD_STORE_INT(NT_NKT,'NT_NKT',L_FALSE,'Number of energy bins')
	      CALL RD_STORE_DBLE(NT_EMIN,'NT_EMIN',L_FALSE,'Minimum energy of non-thermal electrons in eV')
	      CALL RD_STORE_DBLE(NT_EMAX,'NT_EMAX',L_FALSE,'Maximum energy of non-thermal electrons in eV')
	      CALL RD_STORE_LOG(SCL_NT_CROSEC,'SCL_NT_CROSEC',L_FALSE,
	1               'Scale the nonthermal excitation cross sections?')
	      CALL RD_STORE_LOG(SCL_NT_ION_CROSEC,'SCL_NT_ION_CROSEC',L_FALSE,
	1               'Scale the nonthermal ionization cross sections?')
	      CALL RD_STORE_INT(NON_THERMAL_IT_CNTRL,'NT_IT_CNTRL',L_FALSE,
	1               'Controls how often we update the non-thermal electron distiution')
	      CALL RD_STORE_DBLE(NT_OMIT_LEV_SCALE,'NT_OMIT_LEV_SCALE',L_FALSE,
	1               'Fractional populations below this level are excluded when computing the non-thermal electron spectrum')
	      CALL RD_STORE_DBLE(NT_OMIT_ION_SCALE,'NT_OMIT_ION_SCALE',L_FALSE,
	1               'Excludes ions with population NT_OMIT_SCALE_FRAC*(species of ion pop)')
	      CALL RD_STORE_CHAR(NT_SOURCE_TYPE,'NT_SOURCE',L_FALSE,'Non-thermal source type - INJECT_DIRAC, CONSTANT or BELL_SHAPE')
	    ENDIF
!
	    ADD_DEC_NRG_SLOWLY=.FALSE.
	    CALL RD_STORE_LOG(ADD_DEC_NRG_SLOWLY,'GAMMA_SLOW',L_FALSE,
	1            'Add radioactivity decay energy slowly?')
	    IF(ADD_DEC_NRG_SLOWLY)THEN
	       CALL RD_STORE_DBLE(DEC_NRG_SCL_FAC_BEG,'DECNRG_SCLFAC_BEG',ADD_DEC_NRG_SLOWLY,
	1            'Initial Scale factor for adding decay energy')
            END IF
	    MINIMUM_ISO_POP=1.0D-20
	    CALL RD_STORE_DBLE(MINIMUM_ISO_POP,'MIN_ISO_POP',L_FALSE,'Minimum population for ant ISOTOPE')
!
	    CALL RD_STORE_LOG(COMP_GREY_LST_IT,'COMP_GREY_LST_IT',L_FALSE,'Compute grey solution on last iteration?')
!
	    CALL RD_STORE_NCHAR(SN_T_OPTION,'SN_T_OPT',ITEN,L_TRUE,
	1           'Method to get T with non-GRID option (USE_T_IN or USE_HYDRO)')
	    CALL RD_STORE_NCHAR(GAMRAY_TRANS,'GAMRAY_TRANS',ITEN,L_TRUE,
	1           'NONLOCAL or LOCAL gamma-ray enegry transport')
	    INSTANTANEOUS_ENERGY_DEPOSITION=.FALSE.
	    CALL RD_STORE_LOG(INSTANTANEOUS_ENERGY_DEPOSITION,'INS_DEP',L_FALSE,'Instantaneous gamma-ray energy deposition?')
!
	    N_IB_INS=2
	    N_OB_INS=3
	    RMAX_ON_RCORE=-1.0D0		!Implies use default.
	    CALL RD_STORE_INT(N_IB_INS,'N_IB_INS',L_FALSE,'# of points for fine grid at inner boundary')
	    CALL RD_STORE_INT(N_OB_INS,'N_OB_INS',L_FALSE,'# of points for fine grid at outer boundary')
	    CALL RD_STORE_DBLE(RMAX_ON_RCORE,'RMAX_ON_RCORE',L_FALSE,'RMAX/RCORE for SN if shrinking radius')
	  END IF	  
	  DO_FULL_REL_OBS=.FALSE.
	  DO_FULL_REL_CMF=.FALSE.
	  CALL RD_STORE_LOG(DO_FULL_REL_OBS,'REL_OBS',L_FALSE,
	1          'Include all reltivistic terms in observer''s frame solution')
	  CALL RD_STORE_LOG(DO_FULL_REL_CMF,'REL_CMF',L_FALSE,
	1        'Include all reltivistic terms in CMF solution for observed intensity')
!
! These two options are installed only to retain consistency with earlier versions.
!
	  USE_OLD_MF_SCALING=.FALSE.
	  USE_OLD_MF_OUTPUT=.FALSE.
	  CALL RD_STORE_LOG(USE_OLD_MF_SCALING,'OLD_MFS',L_FALSE,
	1        'Use old mass scaling when reading SN_HYDRO_DATA')
	  CALL RD_STORE_LOG(USE_OLD_MF_OUTPUT,'OLD_MFO',L_FALSE,
	1        'Use old mass scaling when reading SN_HYDRO_DATA')
!
C
C Read in the un-normalized fractional abundances.
C
	  DO ISPEC=1,NUM_SPECIES
	    TMP_KEY=TRIM(SPECIES(ISPEC))//'/X'
	    TMP_STRING=TRIM(SPECIES(ISPEC))//
	1            '/X fractional abundance by number (un-normalized)'
	    CALL RD_STORE_DBLE(AT_ABUND(ISPEC),TMP_KEY,SPECIES_PRES(ISPEC),
	1            TMP_STRING)
	  END DO
	  WRITE(LUSCR,'()')
C
	  CALL RD_STORE_LOG(RD_CONT_FREQ,'RD_CF_FILE',L_TRUE,
	1            'Read in continuum frequencies from file')
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
! Installed to allow earlier frequency grids to be used. 0 uses the
! the latest default grid.
!
	  FREQ_GRID_OPTION=0
	  CALL RD_STORE_INT(FREQ_GRID_OPTION,'FR_GRID',L_FALSE,
	1            'Which method to compute frequency grid?')
C
	  CALL RD_STORE_LOG(DO_LEV_DISSOLUTION,'DO_DIS',L_TRUE,
	1            'Allow for level dissolution of upper levels?')
	  CALL RD_STORE_DBLE(dV_LEV_DIS,'dV_LEV',L_TRUE,
	1             'Spacing (in km/s) on low side of bf edge for'//
	1             ' level dissolution')
	  CALL RD_STORE_DBLE(AMP_DIS,'AMP_DIS',L_TRUE,
	1            'Amplification factor on low side bf edge')
	  CALL RD_STORE_DBLE(MIN_FREQ_LEV_DIS,'MIN_DIS',L_TRUE,
	1            'Minimum frequency for level dissolution')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(COMPUTE_ALL_CROSS,'CROSS',L_TRUE,
	1            'Compute all photoionization cross-sections?')
	  CALL RD_STORE_DBLE(DELV_CONT,'V_CROSS',L_TRUE,
	1            'Max. vel. sep. (km/s) between evaluations of all'//
	1            '  phot. cross-sections?')
	  CALL RD_STORE_LOG(DIE_AS_LINE,'DIE_AS_LINE',L_TRUE,
	1        'Treat the dielectronic transitions as individual lines?')
	  CALL RD_STORE_DBLE(VSM_DIE_KMS,'VSM_DIE',L_TRUE,
	1        'Velocity (km/s) for smoothing dielectronic transitions')
	  CALL RD_STORE_DBLE(SIG_GAU_KMS,'SIG_GAU_KMS',L_TRUE,
	1        'Sigma of Gaussian used to smooth photoionization data')
	  FRAC_SIG_GAU=0.25D0
	  CALL RD_STORE_DBLE(FRAC_SIG_GAU,'FRAC_SIG',L_FALSE,
	1        'Fractional spacing a across smoothing Gauusian (use 0.25)')
	  CUT_ACCURACY=0.02D0
	  CALL RD_STORE_DBLE(CUT_ACCURACY,'CUT_ACC',L_FALSE,
	1        'Accuracy to retain data when omitting data points to save space (use 0.02)')
	  ABOVE_EDGE=.TRUE.
	  CALL RD_STORE_LOG(ABOVE_EDGE,'ABV_EDGE',L_FALSE,
	1        'Use only data above edge when smoothing (TRUE)')
!
	  CALL RD_STORE_DBLE(EXT_LINE_VAR,'EXT_LINE_VAR',L_TRUE,
	1        'Extent of line variation zone (V/INF) beyond'//
	1        'the resonancze zone')
          IF(EXT_LINE_VAR .LT. 0.0D0 .OR. EXT_LINE_VAR .GT. 2.0D0)THEN
	    WRITE(LUER,*)'Error in CMFGEN --- invalid range for EXT_LINE_VAR'
	    STOP
	  END IF
	  NEW_VAR_STORAGE_METHOD=.TRUE.
	  CALL RD_STORE_LOG(NEW_VAR_STORAGE_METHOD,'NEW_VST_METH',L_FALSE,
	1         'Use new storage method for lines (different storage for each line)?')
	!
! NB: An ideal vale for ZNET_VAR_LIMIT is probably 0.01 or 0.001. If 
! ZNET_VAR_LIMIT is zero, all depths will be included in the linearization, 
! independent of ZNET. A very large value of ZNET (i.e. 10^4), will imply
! an interation on the NET_RATES, with no linearization.
!
	  CALL RD_STORE_DBLE(ZNET_VAR_LIMIT,'ZNET_VAR_LIM',L_TRUE,
	1            'Include lines in full varaition when '//
	1            ' ABS(ZNET-1) > ZNET_VAR_LIM')
	  CALL RD_STORE_LOG(WEAK_WITH_NET,'WNET',L_TRUE,'Use Lambda iteration for weak lines?')
	  CALL RD_STORE_DBLE(WEAK_LINE_LIMIT,'WK_LIM',L_TRUE,'Maximum opacity ratio for weak lines (0.01)?')
!
	  USE_WEAK_TAU_LIM=.FALSE.; WEAK_TAU_LINE_LIMIT=0.01D0
	  CALL RD_STORE_LOG(USE_WEAK_TAU_LIM,'WK_TAU',L_FALSE,'Use TAU(SOB) to decide weak lines?')
	  CALL RD_STORE_DBLE(WEAK_TAU_LINE_LIMIT,'WK_TAU_LIM',USE_WEAK_TAU_LIM,
	1                      'Maximum TAU for weak lines (0.01)?')
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
	  CALL RD_STORE_CHAR(OUTER_BND_METH,'OB_METH',L_FALSE,'Outer boundary method (HONJ or HALF_MOM)')
C
	  CALL RD_STORE_LOG(RD_COHERENT_ES,'COH_ES',L_TRUE,
	1            'Assume coherent electron scattering? ')
	  CALL RD_STORE_LOG(USE_OLDJ_FOR_ES,'OLD_J',L_TRUE,
	1            'Use old file to provide initial estimate of J_ES?')
	  COHERENT_ES=RD_COHERENT_ES
	  CALL RD_STORE_LOG(MIXED_ES_VAR,'MIX_COH',L_TRUE,
	1            'Mix coherent/non-coherent e.s. in linearization?')
	  CALL RD_STORE_DBLE(ES_VAR_FAC,'ES_FAC',L_TRUE,
	1            'Fractional proximity of RJ and RJ_ES for coherent'//
	1            ' variation')
!
	  LTE_MODEL=.FALSE.
	  CALL RD_STORE_LOG(LTE_MODEL,'LTE_MOD',L_FALSE,
	1            'Force populations to LTE and iterate T only?')
!
	  CALL RD_STORE_NCHAR(METHOD,'METHOD',ISIX,L_TRUE,
	1         'Which method for continuum tau'//
	1         ' loglog, loglin, linear or zero ?')
	  LUM_FROM_ETA_METHOD='LINMON'
	  CALL RD_STORE_NCHAR(LUM_FROM_ETA_METHOD,'LUM_METH',ISIX,L_FALSE,
	1         'Which method for computing L from ETA '//
	1         ' loglog, loglin, linear or zero ?')
	  CALL RD_STORE_NCHAR(N_TYPE,'N_TYPE',ISIX,L_TRUE,
	1         'Method for to handle N for MOM_J_CMF -- '//
	1         'N_ON_J, MIXED, or G_ONLY')
	  H_CHK_OPTION='NONE'
	  CALL RD_STORE_CHAR(H_CHK_OPTION,'CHK_H',L_FALSE,'NONE, AV_VAL, or MAX_VAL')
	  CALL RD_STORE_NCHAR(FG_SOL_OPTIONS,'FG_OPT',ITEN,L_TRUE,
	1         'Solution options for FG_J_CMF: DIFF/INS and INT/INS')
	  CALL RD_STORE_DBLE(DELV_FRAC_FG,'VFRAC_FG',L_TRUE,
	1         'Maximum velocity spacing (Doppler units) in FG_J_CMF_V10')
	  CALL RD_STORE_DBLE(DELV_FRAC_MOM,'VFRAC_MOM',L_TRUE,
	1         'Maximum velocity spacing (Doppler units) in MOM_J_CMF_V10')
!
	  CALL RD_STORE_LOG(RDTHK_CONT,'THK_CONT',L_TRUE,'Use thick boundary condition for continuum ? ')
	  RD_OUT_BC_TYPE=1
	  OUT_BC_TYPE=1 
	  OUT_BC_PARAM_ONE=0.299794D0
	  CALL RD_STORE_INT(RD_OUT_BC_TYPE,'OBC_TYPE',L_FALSE,'Outer boundary condition type: 1=def=old')
	  CALL RD_STORE_INT(OUT_BC_PARAM_ONE,'BC_PAR1',L_FALSE,'Frequency to switch to new BC')
	  REXT_FAC=0.0D0
	  CALL RD_STORE_DBLE(REXT_FAC,'REXT_FAC',L_FALSE,'Factor ot extend R grid by for thick continuum')
!
	  INCL_INCID_RAD=.FALSE.
	  CALL RD_STORE_LOG(INCL_INCID_RAD,'INCID_RAD',L_FALSE,
	1           'Include incident radiation for plane-parellel mod with V?')
	  CALL RD_STORE_LOG(TRAPFORJ,'TRAP_J',L_TRUE,
	1           'Use trapazoidal weights to compute J? ')
!
! This section of the code should work with old VADAT  files when using a Doppler profile of
! fixed width (the default).
!
	  WRITE(LUSCR,'()')
	  FIX_DOP=.TRUE.
	  AMASS_DOP=1.0D0
	  CALL RD_STORE_LOG(FIX_DOP,'FIX_DOP',L_FALSE,
	1      'Use the same turbulent velocity for all species?')
	  CALL RD_STORE_DBLE(TDOP,'TDOP',L_TRUE,
	1      'Temperature to be used in Doppler profile (10^4K)')
	  IF(FIX_DOP)THEN
	    CALL RD_STORE_DBLE(AMASS_DOP,'AMASS_DOP',L_TRUE,
	1      'Atomic mass to be used in Doppler profile (amu''s)')
	    CALL RD_STORE_DBLE(VTURB,'VTURB',L_TRUE,
	1      'Turbulent velocity to be used in Doppler profile (km/s)')
	    VTURB_MIN=VTURB; VTURB_MAX=VTURB
	    GLOBAL_LINE_PROF='DOP_FIX'
	  ELSE
	    CALL RD_STORE_DBLE(VTURB_MIN,'VTURB_MIN',L_TRUE,
	1      'Minimum turbulent velocity for Doppler profile (km/s)')
	    CALL RD_STORE_DBLE(VTURB_MAX,'VTURB_MAX',L_TRUE,
	1      'Maximum turbulent velocity for Doppler profile (km/s)')
	    VTURB=VTURB_MIN
!
	    WRITE(LUSCR,'()')
	    CALL RD_STORE_NCHAR(GLOBAL_LINE_PROF,'GLOBAL_PROF',ITEN,L_TRUE,
	1        'Global switch for intrinsic line absorption profile')
	    CALL SET_CASE_UP(GLOBAL_LINE_PROF,IZERO,IZERO)
	    IF( GLOBAL_LINE_PROF .NE. 'NONE' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOP_FIX' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOP_SPEC' .AND.
	1       GLOBAL_LINE_PROF .NE. 'DOPPLER' .AND.
	1       GLOBAL_LINE_PROF .NE. 'LIST' .AND.
	1       GLOBAL_LINE_PROF .NE. 'LIST_VGT' .AND.
	1       GLOBAL_LINE_PROF .NE. 'VOIGT' .AND.
	1       GLOBAL_LINE_PROF .NE. 'HZ_STARK')THEN
	      WRITE(LUER,*)'Invalid GLOBAL_LINE_PROF parameter'
	      STOP
	    END IF
	    CALL RD_STORE_LOG(SET_PROF_LIMS_BY_OPACITY,'OPAC_LIMS',L_TRUE,
	1           'Set prof limits by line to cont. ratio?')
	    CALL RD_STORE_DBLE(DOP_PROF_LIMIT,'DOP_LIM',L_TRUE,
	1           'Edge limits for Doppler line profile')
	    CALL RD_STORE_DBLE(VOIGT_PROF_LIMIT,'VOIGT_LIM',L_TRUE,
	1           'Edge limits for Voigt line profile')
	  END IF
!
! The following options do not apply to a Doppler profile.
!
	  V_PROF_LIMIT=3000.0D0
	  MAX_PROF_ED=1.0D+16
	  NORM_PROFILE=.TRUE.
	  CALL RD_STORE_DBLE(MAX_PROF_ED,'MAX_PROF_ED',L_FALSE,
	1           'Maximum electron density for Stark profile computation')
	  CALL RD_STORE_DBLE(V_PROF_LIMIT,'V_PROF_LIM',L_FALSE,
	1           'One-sided profile limit for Stark profiles (km/s)')
	  CALL RD_STORE_LOG(NORM_PROFILE,'NORM_PROF',L_FALSE,
	1           'When true, profiles are normalized to have unit area.')
!
	  WRITE(LUSCR,'()')
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
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_DBLE(OBS_PRO_EXT_RAT,'OBS_EXT_RAT',L_TRUE,
	1      'Half width of profile in Vinf.')
	  CALL RD_STORE_DBLE(dV_OBS_PROF,'dV_OBS_PROF',L_TRUE,
	1      'Spacing across observed profile (in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_WING,'dV_OBS_WING',L_TRUE,
	1      'Spacing across e.s. wings of observed profile(in km/s)')
	  CALL RD_STORE_DBLE(dV_OBS_BIG,'dV_OBS_BIG',L_TRUE,
	1      'Frequency spacing between lines (in km/s)')
C
	  WRITE(LUSCR,'()')
	  USE_FIXED_J=.FALSE.
	  CALL RD_STORE_LOG(USE_FIXED_J,'USE_FIXED_J',L_FALSE,
	1           'Use previously computed J to evaluate ALL rates?')
	  IF(USE_FIXED_J .AND. .NOT. RD_LAMBDA)THEN
	    WRITE(LUER,'(A)')' Warning: RD_LAMBDA must be TRUE when USE_FIXED_J is TRUE.'
	    WRITE(LUER,'(A)')' Please change the setting in IN_ITS which is read after each iteration'
	    STOP
	  END IF
	  CALL RD_STORE_LOG(FLUX_CAL_ONLY,'FLUX_CAL_ONLY',L_TRUE,
	1           'Compute the observers frame flux only ?')
	  CALL RD_STORE_LOG(EXTEND_FRM_SOL,'EXT_FRM_SOL',L_TRUE,
	1           'Extrapolate the formal solution to larger radii?')
	  CALL RD_STORE_LOG(INSERT_FREQ_FRM_SOL,'INS_F_FRM_SOL',L_TRUE,
	1           'Insert extra frequencies for formal solution?')
	  CALL RD_STORE_NCHAR(CMF_FORM_OPTIONS,'FRM_OPT',ITEN,L_TRUE,
	1           'Solution options for CMF_FORM_SOL')
	  CALL RD_STORE_LOG(DO_SOBOLEV_LINES,'DO_SOB_LINES',L_TRUE,
	1        'Compute Sobolev rates and EWs for flux calculation?')
	  CALL RD_STORE_LOG(SOB_FREQ_IN_OBS,'SOB_FREQ_IN_OBS',L_TRUE,
	1        ' Allow for SOB & CMF lines in defining observers'//
	1        ' frequencies?')
!
	  VERBOSE_OUTPUT=.FALSE.
	  CALL RD_STORE_LOG(VERBOSE_OUTPUT,'VERBOSE_OUT',L_FALSE,
	1        'Switch on enhanced diagnostic output')
	  CALL SET_VERBOSE_INFO(VERBOSE_OUTPUT)
	  WRITE_RATES=.FALSE.
	  CALL RD_STORE_LOG(WRITE_RATES,'WRITE_RATES',L_FALSE,'Write out NETRATE, TOTRATE, etc')
	  WRITE_JH=.FALSE.
	  IF(SN_MODEL)WRITE_JH=.TRUE.
	  CALL RD_STORE_LOG(WRITE_JH,'WRITE_JH',L_FALSE,'Write out JH_AT_CURRENT_TIME')
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_NCHAR(GLOBAL_LINE_SWITCH,'GLOBAL_LINE',ISIX,L_TRUE,
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
	  CALL RD_STORE_DBLE(GF_CUT,'GF_CUT',L_TRUE,
	1          'gf value to omit transitions')
	  CALL RD_STORE_DBLE(AT_NO_GF_CUT,'AT_CUT',L_TRUE,
	1          'Only omit transitions if AT_NO >= AT_CUT')
	  CALL RD_STORE_INT(GF_LEV_CUT,'GF_LEV_CUT',L_TRUE,
	1          'Level above whit transitions omitted if gf < GF_CUT')
	  CALL RD_STORE_INT(MIN_NUM_TRANS,'MIN_TRANS',L_TRUE,
	1          'Minimum number of transitions from each level')
C
	  CALL RD_STORE_LOG(THK_LINE,'THK_LINE',L_TRUE,
	1           'Use thick boundary condition for lines?')
	  CALL RD_STORE_LOG(CHECK_LINE_OPAC,'CHK_L_POS',L_TRUE,
	1      'Ensure Line opacity is positive (SOB & CMF modes only)?')
	  CALL RD_STORE_NCHAR(NEG_OPAC_OPTION,'NEG_OPAC_OPT',ITEN,L_TRUE,
	1            'Method for negative opacities in BLANKETING mode')
	  CALL SET_CASE_UP(NEG_OPAC_OPTION,IZERO,IZERO)
	  IF(NEG_OPAC_OPTION .NE. 'SRCE_CHK' .AND. 
	1                           NEG_OPAC_OPTION .NE. 'ESEC_CHK')THEN
	     WRITE(LUER,*)'Error in CMFGEN_SUB'
	     WRITE(LUER,*)'Invalid NEG_OPAC_OPTION'
	     WRITE(LUER,*)'Valid options are SRCE_CHK and ESEC_CHK'
	     STOP
	  END IF
	  CALL RD_STORE_LOG(SETZERO,'He2_RES=0',L_TRUE,
	1           'Set rates in He2 resonance lines to zero ?')
C
	  CALL RD_STORE_LOG(OVERLAP,'ALLOW_OL',L_TRUE,
	1           'Allow for overlap of close lines (SOB only) ?')
	  CALL RD_STORE_DBLE(OVER_FREQ_DIF,'OL_DIF',L_TRUE,
	1           'Max. difference (in km/s) for overlap')
	  OVER_FREQ_DIF=OVER_FREQ_DIF/2.998D+05
!
	  CALL RD_STORE_LOG(INCL_CHG_EXCH,'INC_CHG',L_TRUE,'Include charge exchange reactions?')
	  CALL RD_STORE_LOG(INCL_TWO_PHOT,'INC_TWO',L_TRUE,'Include two photon transitions?')
	  TWO_PHOTON_METHOD='USE_RAD'
	  CALL RD_STORE_CHAR(TWO_PHOTON_METHOD,'TWO_METH',L_FALSE,'USE_RAD, LTE, NOSTIM or OLD_DEFAULT')
	  INCL_PENNING_ION=.FALSE.
	  CALL RD_STORE_LOG(INCL_PENNING_ION,'INC_PEN',SN_MODEL,'Include Penning ionization?')
!
	  CALL RD_STORE_LOG(INCL_RAY_SCAT,'INC_RAY',L_TRUE,
	1           'Include opacity due to Rayleigh scattering?')
	  CALL RD_STORE_LOG(INCL_ADVECTION,'INC_ADV',L_TRUE,
	1           'Include advection terms in rate equations?')
	  CALL RD_STORE_LOG(INCL_ADIABATIC,'INC_AD',L_TRUE,
	1           'Include adiabatic cooling in energy equation')
	  SCL_SL_LINE_OPAC=.FALSE.
	  CALL RD_STORE_LOG(SCL_SL_LINE_OPAC,'SCL_SL_OPAC',L_FALSE,
	1            'Scale SL line opacities/emissivties for heating/cooling consistency?')
	  CALL RD_STORE_LOG(SCL_LINE_COOL_RATES,'SCL_LN',L_TRUE,
	1            'Scale line cooling rate for Rad. Eq. equation?')
	  IF(SCL_SL_LINE_OPAC .AND. SCL_LINE_COOL_RATES)THEN
	    WRITE(LUER,*)'Inconsitent use of SCL_SL_OPAC and SCL_LN'
	    WRITE(LUER,*)'Only one of the two options can be true'
	    STOP
	  END IF
	  CALL RD_STORE_DBLE(SCL_LINE_HT_FAC,'SCL_LN_FAC',L_TRUE,
	1            'Scale line cooling rate for Rad. Eq. equation?')
	  SCL_LINE_DENSITY_LIMIT=1.0D+30
	  CALL RD_STORE_DBLE(SCL_LINE_DENSITY_LIMIT,'SCL_DEN_LIM',L_FALSE,
	1            'Density beyond which line cooling scaling is switched off')
	  INCLUDE_dSLdT=.FALSE.
	  CALL RD_STORE_LOG(INCLUDE_dSLdT,'INCL_dSLdT',L_FALSE,
	1            'Include variation in distribution of level populations in a SL with T?')
	  NEW_LINE_BA=.FALSE.
	  IF(SN_MODEL)NEW_LINE_BA=.TRUE.
	  INDX_BA_METH_RD=45
	  CALL RD_STORE_LOG(NEW_LINE_BA,'NEW_LINE_BA',L_FALSE,
	1            'Use the new technique for comuting BA_T (i.e., dRE/dT)?')
	  CALL RD_STORE_INT(INDX_BA_METH_RD,'INDX_BA_METH',L_FALSE,'Depth location to begin NEW_LINE_BA technique')
!
	  LINEAR_ADV=.TRUE.
	  CALL RD_STORE_LOG(LINEAR_ADV,'LIN_ADV',L_FALSE,
	1           'Compute advection terms using derivatives in linear plane?')
	  ADVEC_RELAX_PARAM=1.0D0
	  CALL RD_STORE_DBLE(ADVEC_RELAX_PARAM,'ADV_RELAX',L_FALSE,
	1           'Parameter to allow advection terms to be included slowly')
!
! Except for the X-ray switch, the X-ray options are only needed if we 
! are including X-rays.
!
	  VSMOOTH_XRAYS=3000.0D0
	  CALL RD_STORE_LOG(XRAYS,'INC_XRAYS',L_TRUE,
	1           'Include X-ray emission')
	  CALL RD_STORE_LOG(FF_XRAYS,'FF_XRAYS',XRAYS,
	1           'Use free-free processes to compute X-ray emission')
	  CALL RD_STORE_LOG(XRAY_SMOOTH_WIND,'X_SM_WIND',XRAYS,
	1           'Ignore clumping when computing X-ray emission')
	  CALL RD_STORE_DBLE(VSMOOTH_XRAYS,'VS_XRAYS',XRAYS,
	1           'X-ray smoothing width for SOB/CMF options')
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
	  CALL RD_STORE_LOG(ADD_XRAYS_SLOWLY,'XSLOW',XRAYS,
	1           'Add X-rays by slowly increasing filling factors?')
	  CALL RD_STORE_DBLE(FILL_FAC_X1_BEG,'XFI1_BEG',ADD_XRAYS_SLOWLY,
	1           'Initial filling factor for X-ray emission [1]')
	  CALL RD_STORE_DBLE(FILL_FAC_X2_BEG,'XFI2_BEG',ADD_XRAYS_SLOWLY,
	1           'Initial filling factor for X-ray emission [2]')
	  CALL RD_STORE_DBLE(SLOW_XRAY_SCL_FAC,'XSCL_FAC',ADD_XRAYS_SLOWLY,
	1           'Rate to increase X-ray filling factor')
!
	  SCALE_XRAY_LUM=.FALSE.
	  ALLOWED_XRAY_FLUX_ERROR=0.1D0 
	  CALL RD_STORE_LOG(SCALE_XRAY_LUM,'SCL_XLUM',L_FALSE,
	1           'Scale X-ray emissivities to get a specified X-ray luminosiity')
	  CALL RD_STORE_DBLE(DESIRED_XRAY_LUM,'XRAY_LUM',SCALE_XRAY_LUM,
	1           'X-ray luminosity in units of LSTAR')
	  CALL RD_STORE_DBLE(ALLOWED_XRAY_FLUX_ERROR,'XRAY_ERR',L_FALSE,
	1           'Fractional error allowed in observed X-ray luminosity')
!
	  DELV_XRAY=0.5D0*VSMOOTH_XRAYS
	  CALL RD_STORE_DBLE(DELV_XRAY,'V_XRAY',L_FALSE,
	1            'Max. vel. sep. (km/s) between evaluations of '//
	1            '  phot. cross-sections in X-ray region?')
	  NU_XRAY_END=100.0D0
	  CALL RD_STORE_DBLE(NU_XRAY_END,'NU_XRAY',L_FALSE,
	1            'End of X-ray region for continuum definition')
	  
C
	  WRITE(LUSCR,'()')
	  WRITE(LUSCR,'()')
	  INTERP_DC_SPH_TAU=.FALSE.
	  SET_LTE_AS_INIT_ESTIMATES=.FALSE.
	  CALL RD_STORE_LOG( RDINR,'RD_IN_R_GRID',L_TRUE,
	1        'Read in a predetermined R grid ?')
	  CALL RD_STORE_LOG(GRID,'LIN_INT',L_TRUE,
	1        'Use direct linear interpolation if  new model ?')
	  DC_INTERP_METHOD='ED'
	  CALL RD_STORE_NCHAR(DC_INTERP_METHOD,'DC_METH',ITEN,L_FALSE,
	1        'Interpolation method: ED, R, LTE, SPH_TAU')
	  CALL RD_STORE_LOG(INTERP_DC_SPH_TAU,'DC_SPH_TAU',L_FALSE,
	1        'Interpolate d.c''s on the spherical TAU scale?')
	  IF(INTERP_DC_SPH_TAU)DC_INTERP_METHOD='SPH_TAU'
	  CALL RD_STORE_LOG(SET_LTE_AS_INIT_ESTIMATES,'LTE_EST',L_FALSE,
	1        'Use LTE for the initial estimates')
!
	  IF(SET_LTE_AS_INIT_ESTIMATES)DC_INTERP_METHOD='LTE'
	  IF(DC_INTERP_METHOD .EQ. 'LTE')THEN
	    CALL RD_STORE_DBLE(T_EXCITE_MIN,'T_EXC',L_TRUE,'Minimum temp. when using LTE option')
	  END IF
	  CALL RD_STORE_LOG(DO_POP_SCALE,'POP_SCALE',L_TRUE,
	1        'Scale populations so that cons. Eq. satisfied ?')
	  CALL RD_STORE_DBLE(T_INIT_TAU,'T_INIT_TAU',L_TRUE,
	1        'Tau above which T is set exactly to T(spherical)')
	  CALL RD_STORE_LOG(ITERATE_INIT_T,'IT_ON_T',L_TRUE,
	1        'Improve initial T estimate by iteration ?')
	  INTERP_T_ON_R_GRID=.TRUE.
	  IF(DC_INTERP_METHOD .EQ. 'R')INTERP_T_ON_R_GRID=.FALSE.
	  CALL RD_STORE_LOG(INTERP_T_ON_R_GRID,'T_ON_R',L_FALSE,
	1       'Interpolate T on R grid - default is SPH tau grid')
	  CALL RD_STORE_DBLE(GREY_PAR,'GREY_TAU',L_TRUE,
	1        'SpecifysTau above which T is set to TGREY in iterative process')
!
! We only need to specify the TRANSITION type if GLOBAL_LINE_SWITCH is NONE.
!
	  IF(GLOBAL_LINE_SWITCH(1:4) .EQ. 'NONE')THEN
	    WRITE(LUSCR,'()')
	    DO ID=1,NUM_IONS-1
	      TMP_KEY='TRANS_'//TRIM(ION_ID(ID))
	      TMP_STRING='Method for treating '//TRIM(ION_ID(ID))//' lines?'
	      CALL RD_STORE_NCHAR(ATM(ID)%XzV_TRANS_TYPE,TMP_KEY,
	1          ISIX,ATM(ID)%XZV_PRES,TMP_STRING)
	    END DO
	  END IF
C
	  WRITE(LUSCR,'()')
	  DO ID=1,NUM_IONS-1
	    TMP_KEY='DIE_'//TRIM(ION_ID(ID))
	    TMP_STRING='Include (?) LTDR AUOT, WI calc.s for '//TRIM(ION_ID(ID))
	    CALL RD_STORE_2LOG(ATM(ID)%DIE_AUTO_XzV,ATM(ID)%DIE_WI_XzV,
	1         TMP_KEY,ATM(ID)%XZV_PRES,TMP_STRING)
	  END DO
!
! When adding new species / ioization stages it is sometime useful to hold
! some levels fixed. This can be done in two ways. If we set FIX_ALL_SPEC true,
! all levels will be held fixed UNLESS we do an UNFIX_XzV command. This will be the quickest,
! and safest approach when adding a few ioization stages. We still need to set FIX_T and FIX_NE.
!
! Alternatively, we can set FIX_XzV for every ionization stage.
!
! NB: None of  the keywords need be present.
!
	  FIX_ALL_SPECIES=.FALSE.
	  CALL RD_STORE_INT(FIX_ALL_SPECIES,'FIX_ALL_SPEC',L_FALSE,'Fix all species?')
!
! Since we don't use FIX_XzV very much, these have been changed to a hidden variable.
!
	  IF(FIX_ALL_SPECIES)THEN
	    DO ISPEC=1,NUM_SPECIES
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        ATM(ID)%FIX_NXzV=ATM(ID)%NXzV
	      END DO
	      FIX_SPECIES(ISPEC)=1
	    END DO
!
! We now allow the possibility of unfixing some levels.
! The use of UFIX_XzV avoids confusion and abiguity with the FIX_XzV command.
!
	    DO ISPEC=1,NUM_SPECIES
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        TMP_KEY='UNFIX_'//TRIM(ION_ID(ID))
	        TMP_STRING='Fix ? levels of '//TRIM(ION_ID(ID))//' : set to 0 to unfix'
	        I=-1; CALL RD_STORE_INT(I,TMP_KEY,L_FALSE,TMP_STRING)
	        IF(I .NE. -1)ATM(ID)%FIX_NXzV=I
	      END DO
	      TMP_KEY='FIX_'//TRIM(SPECIES(ISPEC))
	      TMP_STRING='Fix (?) highest ionization stage in '//TRIM(SPECIES(ISPEC))
	      I=-1; CALL RD_STORE_INT(I,TMP_KEY,L_FALSE,TMP_STRING)
	      IF(I .NE. -1)FIX_SPECIES(ISPEC)=I
	    END DO
!
	  ELSE
!
! Using the fix command to fix levels.
!
	    DO ISPEC=1,NUM_SPECIES
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	        IF(ID .EQ. SPECIES_BEG_ID(ISPEC))WRITE(LUSCR,'()')
	        TMP_KEY='FIX_'//TRIM(ION_ID(ID))
	        TMP_STRING='Fix ? levels of '//TRIM(ION_ID(ID))
	        ATM(ID)%FIX_NXzV=0
	        CALL RD_STORE_INT(ATM(ID)%FIX_NXzV,TMP_KEY,L_FALSE,TMP_STRING)
	      END DO
	      TMP_KEY='FIX_'//TRIM(SPECIES(ISPEC))
	      TMP_STRING='Fix (?) highest ionization stage in '//TRIM(SPECIES(ISPEC))
	      FIX_SPECIES(ISPEC)=0
	      CALL RD_STORE_INT(FIX_SPECIES(ISPEC),TMP_KEY,L_FALSE,TMP_STRING)
	    END DO
	  END IF
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(FIXED_NE,'FIX_NE',L_TRUE,'Fix the electron density ?')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(RD_FIX_IMP,'FIX_IMP',L_TRUE,'Automatically fix impurity species?')
C
	  WRITE(LUSCR,'()')                             
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(RD_FIX_T,'FIX_T',L_TRUE,
	1            'Keep the Temperature fixed ?')
          FIX_IN_BOUND_T=.FALSE.
          CALL RD_STORE_LOG(FIX_IN_BOUND_T,'FIX_INB_T',L_FALSE,
	1            'Fix the Temperature at the inner boundary ?')
          FIX_LST_X_DPTHS=1
	  CALL RD_STORE_INT(FIX_LST_X_DPTHS,'FIX_X_DPTH',L_FALSE,
	1            'Fix the Temperature at the last ? depths')
	  CALL RD_STORE_LOG(VARFIXT,'FIX_T_AUTO',L_TRUE,
	1            'Fix the Temperature automatically ?')
	  CALL RD_STORE_DBLE(TAU_SCL_T,'TAU_SCL_T',L_TRUE,
	1      'Electron scattering optical depth from which to fix T')
	  IF(TAU_SCL_T .EQ. 0.0D0)THEN
	    CON_SCL_T=0.0D0
	  ELSE
	    CON_SCL_T=1000.0D0
	  END IF
	  T_MIN=0.0D0
	  CALL RD_STORE_DBLE(T_MIN,'T_MIN',L_FALSE,'Minimum electron temperature')
	  DO_SRCE_VAR_ONLY=.FALSE. 
	  CALL RD_STORE_LOG(DO_SRCE_VAR_ONLY,'SRCE_ONLY',L_FALSE,
	1            'Allow only ths source function to vary?')
	  ADD_ADDITIONAL_OPACITY=.FALSE.
	  ADD_OPAC_SCL_FAC=0.0D0
	  CALL RD_STORE_LOG(ADD_ADDITIONAL_OPACITY,'ADD_OPAC',L_FALSE,
	1            'Add additional pacity to help converge model?')
	  CALL RD_STORE_DBLE(ADD_OPAC_SCL_FAC,'OP_SCL_FAC',ADD_ADDITIONAL_OPACITY,
	1            'Scale factor for addition of extra opacity')
C
	  CALL RD_STORE_NCHAR(METH_SOL,'SOL_METH',ISIX,L_TRUE,
	1            'Which Method To solve Matrix Equations'//
	1            ' DIAG, TRI, PEN, GSIT or MIN')
	  IF(NUM_BNDS .EQ. 1 .AND. METH_SOL .NE. 'DIAG')THEN
	    WRITE(LUER,*)'****************************************'
	    WRITE(LUER,*)'******WARNING in CMFGEN*****************'
	    WRITE(LUER,*)'Solution method inconsistent with NUM_BNDS'
	    WRITE(LUER,*)'METH_SOL=',METH_SOL,'NUM_BNDS=',NUM_BNDS
	    WRITE(LUER,'(X,A,/)')'Using diagnoal solution'
	    METH_SOL='DIAG'
	  END IF
	  CALL RD_STORE_NCHAR(SCALE_OPT,'SCALE_OPT',ISIX,L_TRUE,
	1           'Scale option (LOCAL, NONE or GLOBAL) ? ')
	  LAM_SCALE_OPT='LIMIT'
	  CALL RD_STORE_NCHAR(LAM_SCALE_OPT,'LAM_SCALE_OPT',ISIX,L_FALSE,
	1           'Only has an effect if set to LIMIT')
	  CALL RD_STORE_DBLE(EPS,'EPS_TERM',L_TRUE,
	1      'If maximum fractional % change < EPS terminate model ')
	  CALL RD_STORE_DBLE(MAX_LIN_COR,'MAX_LIN',L_TRUE,
	1      'Maximum fractional change for linearization ')
	  CALL RD_STORE_DBLE(MAX_LAM_COR,'MAX_LAM',L_TRUE,
	1      'Maximum fractional change for lambda iteration ')
	  MAX_dT_COR=0.2D0
	  CALL RD_STORE_DBLE(MAX_dT_COR,'MAX_dT',L_FALSE,
	1      'Maximum fractional change in the temperature')
	  IF(MAX_dT_COR .LE. 0.0D0 .OR. MAX_dT_COR .GT. 0.201D0)THEN
	    WRITE(LUER,*)' Error: MAX_dT in VADAT has an invalid value of',MAX_dT_COR
	    WRITE(LUER,*)' Require 0 < MAX_dT_COR < 0.2'
	    STOP
	  END IF
	  CALL RD_STORE_DBLE(MAX_CHNG_LIM,'MAX_CHNG',L_TRUE,
	1      'If maximum % fractional change > MAX_CHNG terminate model ')
!
	  CALL RD_STORE_LOG(COMPUTE_BARDIN,'COMP_BA',L_TRUE,
	1            'Compute BA matrix ?')
	  CALL RD_STORE_LOG(WRBAMAT_RDIN,'STORE_BA',L_TRUE,
	1      'Store the BA matrix for/during each iteration ? ')
	  CALL RD_STORE_LOG(WR_BA_INV,'STORE_BA_INV',L_TRUE,
	1      'Store the INVERSE of the BA matrix on each iteration ? ')
	  CALL RD_STORE_LOG(WR_PART_OF_INV,'WR_PRT_INV',L_TRUE,
	1      'Store part of the INVERSE to reduce storage TRIDIAG only)?')
	  CALL RD_STORE_INT(N_ITS_TO_FIX_BA,'N_FIX_BA',L_TRUE,
	1      'Number of iterations to hold BA fixed')
	  CALL RD_STORE_DBLE(BA_CHK_FAC,'BA_CHK_FAC',L_TRUE,
	1      'If dJ < BA_CHK_FAC*RJ, ignore correction to BA')
	  CALL RD_STORE_DBLE(VAL_FIX_BA,'FIX_BA',L_TRUE,
	1      'Switch off BA computation if MAXCH< VAL_FIX_BA ')
	  CALL RD_STORE_DBLE(VAL_DO_LAM,'LAM_VAL',L_TRUE,
	1      'Do Lambda iterations if MAXCH > VAL_DO_LAM')
	  CALL RD_STORE_INT(RD_CNT_LAM,'NUM_LAM',L_TRUE,
	1      '# of Lambda iterations if MAXCH > VAL_DO_LAM')
	  CNT_LAM=0
	  CALL RD_STORE_LOG(RDINSOL,'RD_SOL',L_TRUE,
	1            'RD in solution vector to update populations')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(EDD_CONT,'JC_W_EDD',L_TRUE,
	1        'Compute continuum intensity using Eddington factors')
	  NO_VEL_FOR_CONTINUUM=.FALSE.
	  CALL RD_STORE_LOG(NO_VEL_FOR_CONTINUUM,'NOV_CONT',L_FALSE,
	1        'In non-blanketed mode, ignore velocity terms when computing continuum')
	  CALL RD_STORE_LOG(EDD_LINECONT,'JBAR_W_EDD',L_TRUE,
	1    'Compute line continuum intensity using Eddington factors')
!
	  PLANE_PARALLEL_NO_V=.FALSE.
	  PLANE_PARALLEL=.FALSE.
	  CALL RD_STORE_LOG(PLANE_PARALLEL_NO_V,'PP_NOV',L_FALSE,'Plane-paralle geometry WITHOUT velocity field?')
	  CALL RD_STORE_LOG(PLANE_PARALLEL,'PP_MOD',L_FALSE,'Plane-paralle geometry with velocity field?')
	  IF(PLANE_PARALLEL_NO_V .AND. PLANE_PARALLEL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'PP_NOV and PP_MOD cannot both be true at the sdame time'
	    STOP
	  END IF
!
	  INCL_DJDT_TERMS=.FALSE.
	  CALL RD_STORE_LOG(INCL_DJDT_TERMS,'INCL_DJDT',L_FALSE,'DJDt terms in transfer equaton for SN models?')
	  IF(INCL_DJDT_TERMS)THEN
	    USE_DJDT_RTE=.TRUE.
	    USE_Dr4JDT=.TRUE.
	    CALL RD_STORE_LOG(USE_Dr4JDT,'USE_DR4JDT',L_FALSE,'Difference Dr4JDt')
	  ELSE
	    JGREY_WITH_V_TERMS=.FALSE.
	    USE_DJDT_RTE=.FALSE.
	    USE_Dr4JDT=.FALSE.
	    CALL RD_STORE_LOG(USE_DJDT_RTE,'USE_DJDT_RTE',L_FALSE,
	1    'Use solver which has DJDt terms in transfer equaton for SN models?')
	  END IF
	  DJDT_RELAX_PARAM=1.0D0
	  CALL RD_STORE_DBLE(DJDT_RELAX_PARAM,'DJDT_RELAX',L_FALSE,
	1          'Factor to scale DJDT terms to assist initial convergence')
!
	  USE_LAM_ES=.FALSE.
	  USE_J_REL=.FALSE.
	  INCL_REL_TERMS=.FALSE.
	  INCL_ADVEC_TERMS_IN_TRANS_EQ=.FALSE.
	  USE_FORMAL_REL=.FALSE.
	  CALL RD_STORE_LOG(USE_LAM_ES,'USE_LAM_ES',L_FALSE,'Use lambda iteration with ray solutions?')
	  CALL RD_STORE_LOG(USE_J_REL,'USE_J_REL',L_FALSE,'Use MOM_J_REL_VN to solve the moment equations?')
	  IF(USE_DJDT_RTE .OR. USE_J_REL)USE_FORMAL_REL=.TRUE.
	  CALL RD_STORE_LOG(USE_FORMAL_REL,'USE_FRM_REL',L_FALSE,'Use CMF_FORMAL_REL to compute F etc?')
	  CALL RD_STORE_LOG(INCL_REL_TERMS,'INCL_REL',USE_J_REL,'Include relativistic terms in the transfer equaton?')
	  IF(INCL_REL_TERMS .AND. .NOT. USE_J_REL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'Can only include relativistic terms if USE_J_REL is set to TRUE'
	    STOP
	  END IF
	  IF(USE_J_REL .AND. INCL_REL_TERMS)INCL_ADVEC_TERMS_IN_TRANS_EQ=.TRUE.
	  CALL RD_STORE_LOG(INCL_ADVEC_TERMS_IN_TRANS_EQ,'INCL_ADV_TRANS',USE_J_REL,
	1    'Include advection terms in the transfer equaton?')
	  IF(INCL_ADVEC_TERMS_IN_TRANS_EQ .AND. .NOT. USE_J_REL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'Can only include advection terms in transfer equation if USE_J_REL is set to TRUE'
	    STOP
	  END IF
	  IF(USE_DJDT_RTE .AND. USE_J_REL)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'USE_DJDT_RTE and USE_J_REL cannot be TRUE at the same time'
	    STOP
	  END IF
	  IF( (PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V) .AND. (USE_DJDT_RTE .OR. USE_J_REL) )THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'USE_DJDT_RET and USE_J_REL cannot be TRUE at the same time'
	    STOP
	  END IF
!
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(ACCURATE,'INC_GRID',L_TRUE,
	1          'Increase grid size to improve accuracy? ')
	  CALL RD_STORE_LOG(ALL_FREQ,'ALL_FREQ',L_TRUE,
	1          'Increase accuracy for all frequencies?')
	  CALL RD_STORE_DBLE(ACC_FREQ_END,'ACC_END',L_TRUE,
	1          'Increase accuracy for all frequencies < ACC_END?')
	  CALL RD_STORE_INT(NPINS,'N_INS',L_TRUE,
	1          'Number of points to be inserted in higher'//
	1          ' accuracy grid (1, 2 or 3) ')
	  CALL RD_STORE_INT(ST_INTERP_INDX,'ST_INT',L_TRUE,
	1          'Interpolate from ? ')
	  CALL RD_STORE_INT(END_INTERP_INDX,'END_INT',L_TRUE,
	1          'Interpolate to ? ')
	  CALL RD_STORE_INT(DEEP,'ND_QUAD',L_TRUE,
	1         'Quadratic interpolation from ND-? to ND')
	  CALL RD_STORE_NCHAR(INTERP_TYPE,'INTERP_TYPE',10,L_TRUE,
	1         'Perform interpolations in LOG or LIN plane')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_INT(N_PAR,'N_PAR',L_TRUE,
	1    'Rate of BA incrementation by BA_PAR in cont. loop (# of freq)')
C
C Next two variables apply for both ACCURATE and EDDINGTON.
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(COMPUTE_EDDFAC,'COMP_F',L_TRUE,
	1      'Compute new Eddington factors (f)')
	  CALL RD_STORE_DBLE(ACC_EDD_FAC,'ACC_F',L_TRUE,
	1      'Accuracy with which to compute the eddington factor f')
C
	  WRITE(LUSCR,'()')
	  CALL RD_STORE_LOG(NG_DO,'DO_NG',L_TRUE,
	1         'Perform NG acceleration when applicable ?')
	  CALL RD_STORE_DBLE(VAL_DO_NG,'BEG_NG',L_TRUE,
	1       'Percentage accuracy at which to begin NG acceleration')
	  CALL RD_STORE_INT(IT_TO_BEG_NG,'IBEG_NG',L_TRUE,
	1       'Iteration at which to begin NG acceleration')
	  CALL RD_STORE_INT(NG_BAND_WIDTH,'BW_NG',L_TRUE,
	1       'Depth band width for NG acceleration')
	  CALL RD_STORE_INT(ITS_PER_NG,'ITS/NG',L_TRUE,
	1         'Number of iterations between NG accelerations (>=4)')
	  IF(ITS_PER_NG .LT. 4)THEN
	     WRITE(LUER,*)'Error in CMFGEN - ITS_PER_NG too small'
	     STOP
	  END IF
!
	  WRITE(LUSCR,'()')
	  AVERAGE_DO=.FALSE.; NUM_OSC_AV=4; ITS_PER_AV=8
	  CALL RD_STORE_LOG(AVERAGE_DO,'DO_AV',L_FALSE,'Perform averaging of oscillating pops')
	  CALL RD_STORE_INT(NUM_OSC_AV,'NOSC_AV',L_FALSE,'# of consecquitive oscillations')
	  CALL RD_STORE_INT(ITS_PER_AV,'ITS/AV',L_FALSE,'# of iterations between averaging')
!
!
	  UNDO_LAST_IT=.FALSE.
	  CALL RD_STORE_LOG(UNDO_LAST_IT,'DO_UNDO',L_FALSE,'Undo corrections at last 5 depths')
!
	  STOP_IF_BAD_PARAM=.TRUE.
	  CALL RD_STORE_LOG(STOP_IF_BAD_PARAM,'STOP_IF_BP',L_FALSE,'Undo corrections at last 5 depths')
!
	  CALL CLEAN_RD_STORE()
!
	CLOSE(UNIT=LUIN)
!
! Check consistency of parameters.
!
	CALL CHECK_PARAM_CONSISTENCY()
!
	END
