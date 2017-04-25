	MODULE CMF_FLUX_CNTRL_VAR_MOD
! 
! Set constants that are regularly used, and passed to subroutines.
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: ITHREE=3
	INTEGER, PARAMETER :: IFOUR=4
	INTEGER, PARAMETER :: IFIVE=5
	INTEGER, PARAMETER :: ISIX=6
	INTEGER, PARAMETER :: ISEV=7
	INTEGER, PARAMETER :: ITEN=10
!
	REAL*8, PARAMETER :: RZERO=0.0D0
	REAL*8, PARAMETER :: RONE=1.0D0
	REAL*8, PARAMETER :: RTWO=2.0D0
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
! Used to tack INPUT/OUPUT erros.
!
	INTEGER IOS
	LOGICAL STOP_IF_BAD_PARAM
!
! Variables used to control the spacing in the comoving frequency grid.
!
	REAL*8 MIN_CONT_FREQ            !Minimum continuum frequency.
	REAL*8 MAX_CONT_FREQ            !Maximum continuum frequency.
	REAL*8 SMALL_FREQ_RAT           !Fractional spacing for small frequencies'
	REAL*8 dFREQ_bf_MAX             !Maximum spacing close to bf edge.
	REAL*8 BIG_FREQ_AMP             !Amplification factor
	REAL*8 dV_LEV_DIS               !dV on low side of bound-free edge.
	REAL*8 AMP_DIS                  !Amplification factor
	REAL*8 MIN_FREQ_LEV_DIS         !Minimum frequency for lev d
!
! Variables for treating lines simultaneously with the continuum.
!
	REAL*8 V_DOP			!Variable
	REAL*8 MAX_DOP			!Maximum half-width of resonance zone (in Doppler widths)
	REAL*8 FRAC_DOP			!Spacing in resonance zone (in Doppler widths)
	REAL*8 dV_CMF_PROF		!Spacing across cmf profile (in km/s)
	REAL*8 dV_CMF_WING		!Spacing across e.s. wings of cmf profile (in km/s)
	REAL*8 ES_WING_EXT		!Extent of BLUE e.s. wings from resonance core (in km/s)
	REAL*8 R_CMF_WING_EXT		!Extent of RED e.s. wings from RESONANCE core (in Vinf)
!
! FRAC_DOP_OBS is used to ensure a fine grid across intrinsic photospheric profile which are
! unaffected by the wind velocity field.
!
        REAL*8 OBS_PRO_EXT_RAT		!Used to set maximum extent of line profile (in units of Vinf).
	REAL*8 FRAC_DOP_OBS		!Spacing to apply across intrinsic profile zone (in Doppler widths)
	REAL*8 dV_OBS_PROF		!Spacing to apply across observed profile (in km/s)
	REAL*8 dV_OBS_WING		!Spacing to apply across e.s. wings of observed profile(in km/s)
	REAL*8 dV_OBS_BIG		!Frequency spacing to apply between lines (in km/s)
!
! Indicates whether lines treated in SOB and CMF mode are allowed for when
! constructing the observers frame grid.
!
	LOGICAL SOB_FREQ_IN_OBS
!
	REAL*8    GF_CUT
	INTEGER GF_LEV_CUT
	INTEGER MIN_NUM_TRANS
	LOGICAL ONLY_OBS_LINES
	LOGICAL ONLY_UNOBS_LINES
	CHARACTER*10 GF_ACTION
!
! Controls for dynamic smoothing of the photoionization cross-sections, and for
! treating dielectronic lines.
!
	REAL*8 SIG_GAU_KMS		!Sigma of Gaussian used to smooth photoionization data
	REAL*8 FRAC_SIG_GAU		!Fractional spacing a across smoothing Gauusian (use 0.25)
	REAL*8 CUT_ACCURACY		!Accuracy to retain data when omitting data points to save space (use 0.02)'
	LOGICAL ABOVE_EDGE		!Use only data above edge when smoothing?
	REAL*8 VSM_DIE_KMS		!Sigma of Gaussian used for dielectronic lines that are read in.
	LOGICAL DIE_AS_LINE		!Treats dielectronic lines as individual lines: not included in spectrum calc. 
!
! DELV_CONT indicates the minimum continuum spacing in km/s between evaluations of
! the continuum opcities. It should be less than the smoothing applied to the 
! photoionization and dielectronic line cross-sections.
!
! When COMPUTE_ALL_CROSS is TRUE, new cros-sections are computed at every wavelength.
!
	REAL*8 DELV_CONT
	LOGICAL COMPUTE_ALL_CROSS
	LOGICAL COMPUTE_NEW_CROSS
!
	REAL*8 TDOP
	REAL*8 AMASS_DOP
	REAL*8 DOP_PROF_LIMIT
	REAL*8 VOIGT_PROF_LIMIT
	REAL*8 MAX_PROF_ED
	REAL*8 V_PROF_LIMIT
	LOGICAL NORM_PROFILE
	LOGICAL SET_PROF_LIMS_BY_OPACITY
	LOGICAL RD_STARK_FILE
	CHARACTER*10 GLOBAL_LINE_PROF
!
! GLOBAL_LINE_SWITCH provides an option to handle all LINE by the same method.
! The local species setting only takes precedence when it is set to NONE.
!
	CHARACTER*10 GLOBAL_LINE_SWITCH
	REAL*8 FLUX_CAL_LAM_BEG
	REAL*8 FLUX_CAL_LAM_END
	LOGICAL SET_TRANS_TYPE_BY_LAM
	LOGICAL DO_SOBOLEV_LINES
!
	REAL*8 VTURB_FIX		!Turbelent velocity used with fixed DOPPLER profile.
	REAL*8 VTURB_MIN		!Minimum turbuleb velocity
	REAL*8 VTURB_MAX		!Maximum turbulent velocity (VTURB = VTURB_MIN + (VTURB_MAX-VTURB_MIN)*V(r)/VInf
	CHARACTER(LEN=10) TURB_LAW	!LAW_V1 or LAW_TAU1
!
	INTEGER NUM_ES_ITERATIONS
	CHARACTER(LEN=10) FG_SOL_OPTIONS
	CHARACTER(LEN=10) CMF_FORM_OPTIONS
	CHARACTER(LEN=10) TWO_PHOTON_METHOD
	CHARACTER(LEN=10) H_CHK_OPTION
	CHARACTER(LEN=6)  N_TYPE
	CHARACTER(LEN=6)  METHOD
!
! FLUX_CAL_ONLY provides a method for computing the continuous spectrum
! only (i.e. no linearization or population corrections):
!    To get a BLANKETED spectrum FLUX_CAL_ONLY should be set to
!       TRUE and GLOBAL_LINE_SWITCH to BLANK
!    To get a pure UNBLANKETED spectrum FLUX_CAL_ONLY should be set to
!       TRUE and GLOBAL_LINE_SWITCH to SOB
!
	LOGICAL FLUX_CAL_ONLY
	LOGICAL USE_FIXED_J
	LOGICAL SN_MODEL
!
! Options to control optional output.
!
	LOGICAL COMPUTE_J                       !Compute J (default is TRUE)
	LOGICAL WRITE_ETA_AND_CHI		!Write Eta and CHI out (EDDFACTOR file format)
	LOGICAL WRITE_dFR			!Write dF as a function of R and frequency.
	LOGICAL WRITE_RTAU			!Write R(Tau=Tau_ref) (i.e., I as a function of impct parameter)
	LOGICAL WRITE_IP			!Write I(p) (i.e., I as a function of impct parameter)
	LOGICAL WRITE_FLUX			!Write flux file
	LOGICAL WRITE_CMF_FORCE
	LOGICAL WRITE_SOB_FORCE
	LOGICAL WR_ION_LINE_FORCE
	REAL*8 EW_CUT_OFF
	REAL*8 TAU_REF
	REAL*8 DJDt_RELAX_PARAM
!
	LOGICAL DO_REL_IN_OBSFRAME
	LOGICAL DO_CMF_REL_OBS
	LOGICAL USE_J_REL
	LOGICAL USE_FORMAL_REL
	LOGICAL USE_LAM_ES
	LOGICAL INCL_REL_TERMS
	LOGICAL INCL_ADVEC_TERMS_IN_TRANS_EQ
	LOGICAL INCL_DJDT_TERMS
	LOGICAL USE_DJDT_RTE
	LOGICAL USE_DR4JDT
!
	LOGICAL PLANE_PARALLEL_NO_V
	LOGICAL PLANE_PARALLEL
	LOGICAL INCL_INCID_RAD
!
	LOGICAL EDDINGTON
	LOGICAL EDD_CONT
	LOGICAL EDD_LINECONT
!
	REAL*8 ACC_EDD_FAC
	LOGICAL COMPUTE_EDDFAC
!
! General boundary condition options.
!
	REAL*8 REXT_FAC
	REAL*8 IB_STAB_FACTOR
	LOGICAL RDTHK_CONT
	LOGICAL DIF
	LOGICAL THK_CONT
	LOGICAL THK_LINE
	CHARACTER(LEN=10) OUTER_BND_METH
	CHARACTER(LEN=10) INNER_BND_METH
!
! Minimum velocity step size (Doppler widths) for FG_J_CMF and MOM_J_CMF.
! These are used to insert extra points alonga  ray, and provide additional
! freedom to using the ACCURATE option.
!
	REAL*8 DELV_FRAC_FG
	REAL*8 DELV_FRAC_MOM
!
! Used to control integration in Observer's frame.
!
	REAL*8 OBS_TAU_MAX			!Cut integration off when TAU > OBS_TAU_MAX
	REAL*8 OBS_ES_DTAU			!Maximum grid spacing, in e.s. optical depth, along ray.
	INTEGER N_INS_OBS               	!# of additional points inserted/per zone for observer's frame calculation.
	LOGICAL REVISE_P_GRID			!Revise p grid for observer's fram calculation?
	CHARACTER(LEN=10) OBS_INT_METHOD	!STAU (integrate over S) or ETAZ (integrate over eta)
!
! Used in the CMF frame formal solution.
! DEsigned to give higher accuracy in CMF spectra calculation.
!
	LOGICAL EXTEND_FRM_SOL
	LOGICAL INSERT_FREQ_FRM_SOL
!
	LOGICAL DO_CLUMP_MODEL
	LOGICAL DO_LEV_DISSOLUTION
	LOGICAL INCL_TWO_PHOT
	LOGICAL INCL_RAY_SCAT
	LOGICAL COHERENT_ES
	LOGICAL RD_COHERENT_ES
!
! Variables etc for computation of continuum in comoving frame.
!
	LOGICAL CONT_VEL
	LOGICAL USE_OLDJ_FOR_ES
!
	LOGICAL TRAPFORJ			!Use trapaoidal weight for moment intgerations
	LOGICAL CHECK_LINE_OPAC			!Ensure line opacity is +ve
	LOGICAL AT_LEAST_ONE_NEG_OPAC		!
	CHARACTER(LEN=10) NEG_OPAC_OPTION	!SRCE_CHK or SECE_CHK
!
! X-ray variables.
! XRAYS indicates whether or not to include X-ray emission
! from shocks in the wind? If FF_XRAYS is true, we assume the X-rays are
! generated by free-free proceses. If XRAY_SMOOTH_WIND is set, we factor out the
! effect of clumping when computing the effect of the X-ray emissivity.
!
	LOGICAL XRAYS
	LOGICAL FF_XRAYS
	LOGICAL XRAY_SMOOTH_WIND
!
! These shocks are characterized by up to 2 temperatures.
!
	REAL*8 FILL_FAC_XRAYS_1,T_SHOCK_1,V_SHOCK_1
	REAL*8 FILL_FAC_XRAYS_2,T_SHOCK_2,V_SHOCK_2
	REAL*8 XRAY_EMISS_1,XRAY_EMISS_2
	REAL*8 VSMOOTH_XRAYS
!
! Variables required to insert additional points into the depth grid. 
! Allows an increase in program accuracy to overcome rapid ionization 
! changes. 
!
	REAL*8 ACC_FREQ_END
	LOGICAL ALL_FREQ
	LOGICAL ACCURATE
	LOGICAL INACCURATE
	LOGICAL THIS_FREQ_EXT           !Frequency specific.
!
	INTEGER NPINS                   !Points inserted for calc (CMF) when ACCURATE=.TRUE..
	INTEGER ST_INTERP_INDX  	!Interp from ST_INT.. to END_INTERP..
	INTEGER END_INTERP_INDX
	CHARACTER*10 INTERP_TYPE
!
! ND-DEEP to DEEP we use a quadratic interpolation scheme so as to try
! and preserve "FLUX" in the diffusion approximation.
!
	INTEGER DEEP
!
	END MODULE CMF_FLUX_CNTRL_VAR_MOD
