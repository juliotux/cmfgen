	MODULE CONTROL_VARIABLE_MOD
!
! Incorporated: 02-Jan-2014: Changed to allow depth dependent profiles.
! Altered : 05-Apr-2011 : Added vriable R_GRD_REVISED (10-Feb-2011).
! Altered : 16-Jul-2010 : Added FIX_ALL_SPECIES variable.
! Altered : 23-Nov-2007 : LAM_SCALE_OPT included.
! Altered : 20-Feb-2006 : ABOVE_EDGE changed to LOGICAL from REAL*8
! Altered : 29-Jan-2006 : Control variable fors relativistic transfer and time
!                          dependent statistical equilibrium equations installed.
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
!
!*******************************************************************************
!
! Star and velocity specifications
!
	REAL*8 RP			!Stellar radius in units of 10^10 cm
	REAL*8 RMAX
	REAL*8 RMDOT
      	REAL*8 VRP,RN,VINF,EPPS1,GAMMA1
	REAL*8 LUM
	REAL*8 TEFF
	REAL*8 LOGG
	REAL*8 PRESSURE_VTURB
	REAL*8 RP2,VRP2,RN2,VINF2,EPPS2,GAMMA2
	REAL*8 CONS_FOR_R_GRID
	REAL*8 EXP_FOR_R_GRID
!
	INTEGER N_IB_INS			!Number of points for fine grid at inner boundary
	INTEGER N_OB_INS			!Number of points for fine grid at outer boundary
	INTEGER, PARAMETER :: NUM_V_OPTS=1
	CHARACTER*10 VEL_OPTION(NUM_V_OPTS)
!
! Parameters for VELTYPE=3 or VELTYPE=6 (STARPCYG_V2)
!
	REAL*8 SCL_HT,VCORE,VPHOT
	REAL*8 VINF1,V_BETA1,V_EPPS1
	REAL*8 V_BETA2,V_EPPS2		!VINF2 above
	INTEGER VELTYPE
!
! Variable Mdot parameters
!
	LOGICAL VAR_MDOT
	CHARACTER*80 VAR_MDOT_FILE
!
! Supernova variables
!
	REAL*8 RHO_ZERO			!Density at core in gm/cm^3
	REAL*8 RCUB_RHO_ZERO            !r^3 . density at core in (in 10^{-30} gm}
	REAL*8 N_RHO			!Exponent for density (+ve)
	REAL*8 DJDT_RELAX_PARAM         !Factor to assist inclusion of DJDT terms.
	REAL*8 SN_AGE_DAYS              !Age of SN in days.
	REAL*8 RMAX_ON_RCORE
!
	INTEGER TIME_SEQ_NO             !Number of model in time sequence.
	LOGICAL SN_MODEL
	LOGICAL SN_HYDRO_MODEL		!Use HYDRO model for SN input
	LOGICAL PURE_HUBBLE_FLOW        !Forces a pure Hubble flow
	LOGICAL INCL_RADIOACTIVE_DECAY	!Allow for radioactive decay.
	LOGICAL JGREY_WITH_V_TERMS 	!Include velocity terms when comuting T(grey).
	LOGICAL INCL_DJDT_TERMS         !Include DJDt terms in radiative transfer equation
	LOGICAL USE_DJDT_RTE            !Use the radiative transfer equation solver which includes DJDt terms
	LOGICAL DO_CO_MOV_DDT		!Include comoving drivative in SE equations.
	LOGICAL DO_FULL_REL_OBS         !Include all relativistic terms in obs. frame computation.
	LOGICAL DO_FULL_REL_CMF         !Include all relativistic terms in CMF obs. frame computation.
!
! Non thermal parameters
!
	LOGICAL TREAT_NON_THERMAL_ELECTRONS     !
	LOGICAL ADD_DEC_NRG_SLOWLY
	REAL*8 DEC_NRG_SCL_FAC_BEG
	REAL*8 DEC_NRG_SCL_FAC
	REAL*8 NT_OMIT_ION_SCALE        !Ion omitted if < NT_OMIT_ION_SCL*(largest ion pop.).
	REAL*8 NT_OMIT_LEV_SCALE        !Fractional pops omitted if < NT_OMIT_SCL
        REAL*8 NT_EMAX			!Maximum energy of non-thermal electrons
        REAL*8 NT_EMIN			!Minimum energy of non-thermal electrons
	INTEGER NT_NKT
	INTEGER NON_THERMAL_IT_CNTRL    !Controls how often we update the nonthermal electron distribution.
	LOGICAL COMP_GREY_LST_IT        !Comput J(GREY) on last iteration [DEFAULT is TRUE].
	LOGICAL SCL_NT_CROSEC
	LOGICAL SCL_NT_ION_CROSEC
	CHARACTER(LEN=12) NT_SOURCE_TYPE
!
	LOGICAL USE_J_REL
	LOGICAL USE_FORMAL_REL
	LOGICAL USE_LAM_ES
	LOGICAL INCL_REL_TERMS
	LOGICAL INCL_ADVEC_TERMS_IN_TRANS_EQ
	LOGICAL USE_OLD_MF_SCALING
	LOGICAL USE_OLD_MF_OUTPUT
!
! Used to revise the R_GRID automatically after each iteration.
!
	INTEGER, PARAMETER :: N_RG_PAR_MAX=4
	INTEGER N_RG_PAR
	REAL*8 RG_PAR(N_RG_PAR_MAX)
	LOGICAL REVISE_R_GRID
	LOGICAL R_GRID_REVISED
	CHARACTER(LEN=10) NEW_RGRID_TYPE
	CHARACTER(LEN=10) SN_T_OPTION
	CHARACTER(LEN=10) GAMRAY_TRANS
!
! Variables for including clumping under the assumption that the clumping
! occurs on scales much smaller than any of the transfer scales. Valid
! for the continuum --- likely to be invalid for lines due to the
! shortness of the SObolev length.
!
! CLUMP_FAC represents the fractional volume occupied by material. Thus its
! value is always less than unity. The remaining volume is assumed to be a 
! vacuum. In practice we are assuming dense spherical shells separated by
! a vacuum.
!
	INTEGER, PARAMETER ::  N_CLUMP_PAR_MAX=4
	REAL*8 CLUMP_PAR(N_CLUMP_PAR_MAX)
	INTEGER N_CLUMP_PAR
	LOGICAL DO_CLUMP_MODEL
	CHARACTER(LEN=6) CLUMP_LAW
!
!************************************************************************
!
! IF GRID is TRUE, variables are interpolated directly on old grid.
!
	LOGICAL GRID
	LOGICAL INTERP_DC_SPH_TAU
	LOGICAL SET_LTE_AS_INIT_ESTIMATES
	LOGICAL DO_HYDRO
	LOGICAL DONE_HYDRO_REVISION
        CHARACTER(LEN=10) DC_INTERP_METHOD
!
! Used when constructing the Temperature distribution on the first
! iteration.
!
	REAL*8 GREY_PAR
	REAL*8 T_INIT_TAU
	LOGICAL ITERATE_INIT_T
	LOGICAL T_MIN_BA_EXTRAP
!
! Used to limit the temperature while iterating on T.
!
	REAL*8 T_MIN
!
! When you use the LTE option to set the departure coefficients, T_EXCITE_MIN
! is used in the outer region to set the escitation coefficient. 

	REAL*8 T_EXCITE_MIN
!
! Indicates that the R grid for the NEWMODEL should be rad in from a file --
! ist is not computed.
!
	LOGICAL RDINR
!
! Reads in a modified from of STEQ_VALS so at to apply corrections to the populations.
! For debuggung purposes only.
!
	LOGICAL RDINSOL
!
! Controls for dynamic smoothing of the photoioization cross-sections, and for
! treating dielectronic lines.
!
	REAL*8 SIG_GAU_KMS
	REAL*8 FRAC_SIG_GAU
	REAL*8 CUT_ACCURACY
	REAL*8 VSM_DIE_KMS
	LOGICAL DIE_AS_LINE
	LOGICAL ABOVE_EDGE
!
! Variables for controlling the frequency grid.
!
	REAL*8 MIN_CONT_FREQ 		!Minimum continuum frequency.
	REAL*8 MAX_CONT_FREQ    	!Maximum continuum frequency.
	REAL*8 SMALL_FREQ_RAT 		!Fractional spacing for small frequencies'
	REAL*8 dFREQ_bf_MAX		!Maximum spacing close to bf edge.
	REAL*8 BIG_FREQ_AMP		!Amplification factor
	REAL*8 dV_LEV_DIS		!dV on low side of bound-free edge.
	REAL*8 AMP_DIS			!Amplification factor
	REAL*8 MIN_FREQ_LEV_DIS		!Minimum frequency for lev dissolution.
	LOGICAL RD_CONT_FREQ		!Read in cont. frequencies from file.
	INTEGER FREQ_GRID_OPTION	!Which frequency grid
!
! Used to describe the Doppler lines profile, which is currently assumed
! to be constant at all depths.
!
	REAL*8 TDOP
	REAL*8 VDOP_FIX
	REAL*8 VTURB
	REAL*8 AMASS_DOP
	REAL*8 VTURB_MIN
	REAL*8 VTURB_MAX
	LOGICAL FIX_DOP
!
	REAL*8 DOP_PROF_LIMIT
	REAL*8 VOIGT_PROF_LIMIT
	REAL*8 V_PROF_LIMIT
	REAL*8 MAX_PROF_ED
	LOGICAL SET_PROF_LIMS_BY_OPACITY
	LOGICAL RD_STARK_FILE
	LOGICAL NORM_PROFILE
	CHARACTER(LEN=10) GLOBAL_LINE_PROF
!
! Variables used to define what transitions are neglected. Cut is
! presently by the gf value, and the lower level of the transition.
!
	REAL*8 GF_CUT
	REAL*8 AT_NO_GF_CUT
	INTEGER GF_LEV_CUT
	INTEGER MIN_NUM_TRANS
!
! Variables for treating lines simultaneously with the continuum.
!
	REAL*8 V_DOP
	REAL*8 MAX_DOP
	REAL*8 FRAC_DOP
	REAL*8 dV_CMF_PROF
	REAL*8 dV_CMF_WING
	REAL*8 ES_WING_EXT
	REAL*8 R_CMF_WING_EXT
	REAL*8 EXT_LINE_VAR
	REAL*8 ZNET_VAR_LIMIT
!
! Variables for controlling the frequency grid in the Observer's frame.
!
	REAL*8 OBS_PRO_EXT_RAT
	REAL*8 dV_OBS_PROF
	REAL*8 dV_OBS_WING
	REAL*8 dV_OBS_BIG
!
! Compute all photoionization cross-sections at all frequencies. If not,
! we assume thei value from an earlier frequncy. Saves considerable computational 
! effort. DELV_CONT gives the maximum velocity separation in km/s between 
! evaluations of the phot. cross-sections.

	REAL*8 DELV_CONT
	LOGICAL COMPUTE_ALL_CROSS
	LOGICAL COMPUTE_NEW_CROSS
	LOGICAL USE_FIXED_J
!
!***************************************************************************************
!
! Method to handle N moment  (N_ON_J, MIXED, or G_ONLY)
!
	CHARACTER(LEN=6)  N_TYPE
!
! These insert extra points into the grid when solving for the radition field.
! These insertions are done internally, and do not directly effect the returned
! results. USed mainly with CMF_FLUX.
!
	REAL*8 DELV_FRAC_FG
	REAL*8 DELV_FRAC_MOM
!
! Indicate solution method used to ocmpute the Eddington factors by FG_J_CMF_V?.
! Options is usually  INT/INS (or DIFF/INS).
!
	CHARACTER(LEN=10) FG_SOL_OPTIONS
!
! Indicates type of BA matrix (DIAG, TRIDIAG etc).
!
	CHARACTER(LEN=6)  METH_SOL
	CHARACTER(LEN=10) CMF_FORM_OPTIONS		!Used for formal solution.
	CHARACTER(LEN=10) NEG_OPAC_OPTION
!
! Specifies method used to compute optical depth.
!
	CHARACTER(LEN=6)  METHOD
	CHARACTER(LEN=6)  LUM_FROM_ETA_METHOD
!
        REAL*8 OVER_FREQ_DIF
        REAL*8 WEAK_LINE_LIMIT
	REAL*8 WEAK_TAU_LINE_LIMIT
!
	LOGICAL WEAK_WITH_NET
        LOGICAL USE_WEAK_TAU_LIM
	LOGICAL OVERLAP
	LOGICAL NEW_VAR_STORAGE_METHOD
	INTEGER NUM_OF_WEAK_LINES    	!Counter in CMFGEN_SUB only
!
! Indicates whether how ALL lines are to be treated.
!
	CHARACTER*6 GLOBAL_LINE_SWITCH
	LOGICAL SET_TRANS_TYPE_BY_LAM
	LOGICAL DO_SOBOLEV_LINES
!
	REAL*8 FLUX_CAL_LAM_BEG
	REAL*8 FLUX_CAL_LAM_END
!
! FLUX_CAL_ONLY provides a method for computing the continuous spectrum
! only (i.e. no linearization or population corrections):
!    To get a BLANKETED spectrum FLUX_CAL_ONLY should be set to
!       TRUE and GLOBAL_LINE_SWITCH to BLANK
!    To get a pure UNBLANKETED spectrum FLUX_CAL_ONLY should be set to
!       TRUE and GLOBAL_LINE_SWITCH to SOB
!
	LOGICAL FLUX_CAL_ONLY
!
! Indicates whether lines treated in SOB and CMF mode are allowed for when
! constructing the observers frame grid.
!
	LOGICAL SOB_FREQ_IN_OBS
!
! Variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	LOGICAL EDD_CONT
	LOGICAL COMPUTE_EDDFAC
	LOGICAL EDD_LINECONT
	LOGICAL NO_VEL_FOR_CONTINUUM
!
! Indicate whether the atmosphere is to be treated in plane-parallel mode
!
	LOGICAL PLANE_PARALLEL_NO_V
	LOGICAL PLANE_PARALLEL
!
! Variables required for adding additional points into the depth 
! grid. Allows an increase in program accuracy to overcome
! rapid ionization changes.
!
	REAL*8 ACC_FREQ_END
	INTEGER NPINS			!Points inserted for error calc.
        INTEGER ST_INTERP_INDX        !Interp from ST_INT.. to END_INTERP..
        INTEGER END_INTERP_INDX
        LOGICAL ACCURATE 
	LOGICAL ALL_FREQ
	CHARACTER*10 INTERP_TYPE
!
! These two variables used only in CMFGEN_SUB.
!
	LOGICAL THIS_FREQ_EXT
        LOGICAL INACCURATE
!
! ND-DEEP to DEEP we use a quadratic interpolation scheme so as to try
! and preserve "FLUX" in the diffusion approximation.
!
	INTEGER DEEP
	REAL*8 ACC_EDD_FAC
!
	LOGICAL EXTEND_FRM_SOL
	LOGICAL INSERT_FREQ_FRM_SOL
!
	REAL*8 ES_VAR_FAC
	LOGICAL MIXED_ES_VAR
!
!***********************************************************************************
! Variables for performing NG acceleration.
!
	INTEGER LAST_NG   		!Indicates iteration on which last
!                                          NG acceleration occurred.
	INTEGER NEXT_NG   		!Indicates iteration on which next
!                                          NG acceleration is to occur.
	INTEGER IT_TO_BEG_NG          !Iteration to beg NG acceleration
	INTEGER ITS_PER_NG  		!Iterations between NG accelerations.
	INTEGER NG_BAND_WIDTH         !Number of depths to acclerate simultaneously
	REAL*8 VAL_DO_NG 		!Begin NG when MAXCH < VAL_DO_NG
	LOGICAL NG_DONE			!iIndicates successfull completion of NG
	LOGICAL NG_DO			!Switch on NG acceleration.
!
! Variables for performing averaging of oscilating populations.
!
	INTEGER LAST_LAMBDA
	INTEGER LAST_AV
	INTEGER NEXT_AV
	INTEGER ITS_PER_AV
	INTEGER NUM_OSC_AV
	LOGICAL AVERAGE_DONE
	LOGICAL AVERAGE_DO
	LOGICAL UNDO_LAST_IT
!
! N_PAR is used to indicate how often the BA matrices should be incremented
! by BA_PAR. After the incrementation, the PAR matrices are zeroed.
!
	INTEGER N_PAR
	REAL*8 MAX_LAM_COR	!Maximum fractional change for Lambda iteration.
	REAL*8 MAX_LIN_COR	!Maximum fractional change for linearization.
	REAL*8 MAX_dT_COR       !Maximum allowed change in the temperature.
	REAL*8 MAX_CHNG_LIM
!
! Indicates how to scale the corections to the populations, so that
! populations don't go negative, and MAX_LIN_COR (or MAX_LAM_COR) is
! satisfied. Options are NONE, LOCAL, MAJOR, and GLOBAL.
! MAJOR is generally the preferred choice.
!
	CHARACTER(LEN=6)  SCALE_OPT
!
! Introduced Nov-2003 to limit LAMBDA changes prior to normal scaling of corrections.
! Option only has an effect if LAM_SCALE_OPT='LIMIT'
!
	CHARACTER(LEN=6)  LAM_SCALE_OPT
!
! Used to determine the accuracy with which the BA and BAION matrices are
! computed. 1.0D-10 < BA_CHK_FAC < 0.1. Smaller number means higher acuracy,
! but slower computation. Too large a value may affect convergence. Most
! models run with 1.0D-04. Basically updates to BA ignored if their effect
! on RJ is less than BA_CHK_FAC*RJ.
!
	REAL*8 BA_CHK_FAC
	LOGICAL INCLUDE_dSLdT
	LOGICAL NEW_LINE_BA
	INTEGER INDX_BA_METH_RD
!
! Variables for scaling the line cooling rates in oder that the radiative
! equilibrium equation is more consistent with the electron heating/cooling 
! equation. The scaling is done when the line frequency is with a fraction
! of SCL_LINE_HT_FAC of the tmean frequency for the super-level under 
! consideration. 0.5 is presently the prefered value.
!
! SCL_SL_LINE_OPAC effects the opacities to achieve consistency.
! SCL_LINE_COOL_RATES only modifies the line cooling rates to achieve consistency.
! Only one of these two options can be true at any one time.
!
	LOGICAL SCL_SL_LINE_OPAC
	LOGICAL SCL_LINE_COOL_RATES
	REAL*8 SCL_LINE_HT_FAC
	REAL*8 SCL_LINE_DENSITY_LIMIT
	REAL*8 EPS                    !If maximum fractional % change < EPS we terminate model
!
! Fix BA variation matrix if % change less than VAL_FIX_BA.
!
	REAL*8 VAL_FIX_BA
	INTEGER N_ITS_TO_FIX_BA
	INTEGER CNT_FIX_BA
	LOGICAL COMPUTE_BARDIN,COMPUTE_BA
	LOGICAL WRBAMAT,WRBAMAT_RDIN
!
! Indicates whether inverse of BA matrix should be output. This takes up
! disk space, but saves computaional effort, especially wen BA matrix is
! being held fixed.

	LOGICAL WR_BA_INV
	LOGICAL WR_PART_OF_INV
!
! Performs a lambda iteration if % change > VAL_DO_LAM
!
	REAL*8 VAL_DO_LAM
	INTEGER CNT_LAM
	INTEGER RD_CNT_LAM
	LOGICAL OLD_RD_LAMBDA
	LOGICAL RD_LAMBDA
	LOGICAL DO_LAMBDA_AUTO
	LOGICAL DO_GREY_T_AUTO
	LOGICAL DO_T_AUTO
	LOGICAL LAMBDA_ITERATION
	LOGICAL STOP_IF_BAD_PARAM
!
! Variables for determining whether some populations are held fixed
! when the new populations are solved for.
!
	REAL*8 CON_SCL_T
	REAL*8 TAU_SCL_T
	LOGICAL RD_FIX_T,RD_FIX_IMP
	LOGICAL FIXED_T,FIXED_NE,FIX_IMPURITY
	LOGICAL VARFIXT
	LOGICAL DO_SRCE_VAR_ONLY
	LOGICAL FIX_ALL_SPECIES
	LOGICAL SET_POPS_D2_EQ_D1
!
	REAL*8 ADD_OPAC_SCL_FAC
 	LOGICAL ADD_ADDITIONAL_OPACITY
!
! 
! Indicates whether level dissolution should be included in the model.
!
	LOGICAL DO_LEV_DISSOLUTION
!
	LOGICAL DIF		  	! Use diffusion approximation?
!
	LOGICAL INCL_ADIABATIC   	!Include adiabatic cooling
	LOGICAL INCL_CHG_EXCH		!Include charge exchange reactions.
	LOGICAL INCL_TWO_PHOT		!Include two-photon transitions
	LOGICAL INCL_RAY_SCAT           !Include Rayleigh scattering.
	LOGICAL INCL_PENNING_ION        !Include H/He Pening ionization
	LOGICAL LINEAR_ADV              !Comput derivatives using linear approximation.
!
	LOGICAL INCL_ADVECTION          !Include advection terms in rate equations.
	REAL*8  ADVEC_RELAX_PARAM       !Allows advection terms to be added slowly
!
! Variables specifying computaion of e.s. souce function.
!
	LOGICAL COHERENT_ES
	LOGICAL RD_COHERENT_ES
	LOGICAL USE_OLDJ_FOR_ES
!
	REAL*8 TSTAR
	REAL*8 IC
!
! X-ray variables. XRAYS indicates whether or not to include X-ray emission 
! from shocks in the wind? If FF_XRAYS is true, we assume the X-rays are 
! generated by free-free proceses. If XRAY_SMOOTH_WIND is set, we factor out the
! effect of cluming when computing the effect of the X-ray emissivity.
!
	LOGICAL XRAYS
	LOGICAL ADD_XRAYS_SLOWLY
	LOGICAL FF_XRAYS
	LOGICAL XRAY_SMOOTH_WIND
!
! The shocks are characterized by up to 2 temperatures.
!
	REAL*8 FILL_FAC_XRAYS_1,T_SHOCK_1,V_SHOCK_1
	REAL*8 FILL_FAC_XRAYS_2,T_SHOCK_2,V_SHOCK_2
	REAL*8 FILL_FAC_X1_BEG,FILL_X1_SAV
	REAL*8 FILL_FAC_X2_BEG,FILL_X2_SAV
	REAL*8 SLOW_XRAY_SCL_FAC
	REAL*8 XRAY_EMISS_1,XRAY_EMISS_2
	REAL*8 VSMOOTH_XRAYS
!
	REAL*8 NU_XRAY_END
	REAL*8 DELV_XRAY
!
! Indicates number of time POPS array is to be written to scratch
! file per iteration.
!
        INTEGER RITE_N_TIMES
        PARAMETER (RITE_N_TIMES=1)
!
	INTEGER NLBEGIN
        LOGICAL NEWMOD
        LOGICAL WRITE_RVSIG
!
	CHARACTER(LEN=10) OUTER_BND_METH
	CHARACTER(LEN=10) INNER_BND_METH
!
	REAL*8 IB_STAB_FACTOR
	REAL*8 OUT_BC_PARAM_ONE
	REAL*8 REXT_FAC
	INTEGER RD_OUT_BC_TYPE
	INTEGER OUT_BC_TYPE
	LOGICAL RDTHK_CONT
	LOGICAL THK_CONT
	LOGICAL THK_LINE
	LOGICAL INCL_INCID_RAD
!
	LOGICAL SETZERO
	lOGICAL DO_POP_SCALE
	LOGICAL TRAPFORJ
	LOGICAL CHECK_LINE_OPAC
        LOGICAL SOBOLEV
        LOGICAL VERBOSE_OUTPUT
	LOGICAL WRITE_RATES
!
	END MODULE CONTROL_VARIABLE_MOD
