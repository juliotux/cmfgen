
! Subroutine to compute the emergent flux in the observer's frame using an
! Observer's frame flux calculation. The program first computes the
! emissivities, opacities, and radiation field in the CMF frame. COHERENT
! or INCOHERENT electron scattering can be assumed. This data is then passed
! to OBS_FRAME_SUB for computation of the OBSERVER's flux.
!
! This routine is a heavily stripped down and modified version of
! CMFGEN.
!
	SUBROUTINE PAR_OPAC_SUB(ND,NC,NP,NDMAX,NPMAX,NT,NLINE_MAX)
	USE MOD_PAR_OPAC
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
!
! Altered 03-Jan-2001 : PROF_TYPE and VEC_CHAR_WRK now of length 12.
! Altered: 10-Mar-2000 : Variable type ATM installed to simplify access to
!                          different model atoms. Installation reduces size
!                          of executable, and allows for faster compilation.
!                          Changed to V4.
! Altered: 16-Jan-2000 : Changed to V2; Call to OBS_FRAME_V2 changed.
! Altered: 04-Jan-2000 : We now check that arrays can be allocated.
! Finalized: 5-Jan-1999
!
! NCF_MAX is the maximum number of continuum points (which will include line
! frequencies in blanketing mode) which can be treated. For small DOPPLER widths
! and large frequency ranges this might need to be increased.
!
	INTEGER*4, PARAMETER :: NCF_MAX=1000000
!
! Maximum number of lines whose profile overlap at a given frequency.
!
	INTEGER*4, PARAMETER :: MAX_SIM=500
!
! Allow space to be set aside for the intrinsic line profiles. This saves 
! computational efort. An error message will be printed if these values
! are too small.
!
	INTEGER*4, PARAMETER :: NLINES_PROF_STORE=60
        INTEGER*4, PARAMETER :: NFREQ_PROF_STORE=5000
C
	INTEGER*4 NCF
	INTEGER*4 ND,NC,NP
	INTEGER*4 NDMAX,NPMAX
	INTEGER*4 NT,NLINE_MAX
C
	REAL*8 POPS(NT,ND)
C
C Constants for opacity etc.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
C
	INTEGER*4, PARAMETER :: IZERO=0
	INTEGER*4, PARAMETER :: IONE=1
	INTEGER*4, PARAMETER :: ITWO=2
	INTEGER*4, PARAMETER :: ITHREE=3
	INTEGER*4, PARAMETER :: IFOUR=4
	INTEGER*4, PARAMETER :: IFIVE=5
	INTEGER*4, PARAMETER :: ISIX=6
	INTEGER*4, PARAMETER :: ITEN=10
C
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
C
	REAL*8, PARAMETER :: RZERO=0.0
	REAL*8, PARAMETER :: RONE=1.0
	REAL*8, PARAMETER :: RTWO=2.0
C
C Internally used variables
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
	REAL*8 DTDR,DBB,DDBBDT
	REAL*8 MEAN_ATOMIC_WEIGHT
	REAL*8 TSTAR,S1,IC,MAXCH
	REAL*8 C_KMS
	REAL*8 T1,T2,T3,T4
	REAL*8 FL,FL_OLD
	REAL*8 FG_COUNT
C
C
C REC_SIZE     is the (maximum) record length in bytes.
C UNIT_SIZE    is the number of bytes per unit that is used to specify
C                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
C WORD_SIZE    is the number of bytes used to represent the number.
C N_PER_REC    is the # of POPS numbers to be output per record.
C
 	INTEGER*4 REC_SIZE
	INTEGER*4 UNIT_SIZE
	INTEGER*4 WORD_SIZE
	INTEGER*4 N_PER_REC
C
C 
C Logical Unit assignments. Those indicated with a # after the ! are open in
C  large sections of the code. Other units generally used temprarily.
C
	INTEGER*4               LUER        !Output/Error file.
	INTEGER*4, PARAMETER :: LUIN=7      !General input unit (closed after accesses).
	INTEGER*4, PARAMETER :: LUMOD=8     !Model Description file.
C
	INTEGER*4, PARAMETER :: LU_FLUX=10   	!Flux/Luminosity Data (OBSFLUX)
	INTEGER*4, PARAMETER :: LU_OPAC=18   	!Rosseland mean opacity etc.
	INTEGER*4, PARAMETER :: LU_EW=20     	!# EW data.
C
	INTEGER*4, PARAMETER :: LU_EDD=35       !Continuum Eddington factors.
	INTEGER*4, PARAMETER :: LU_JCOMP=37     !J_COMP
	INTEGER*4, PARAMETER :: LU_ES=38        !ES_J_CONV
C
C For listing of transitions with TOTAL negative opacity values at some depths.
C
	INTEGER*4, PARAMETER :: LU_NEG=75
C
C 
C
	INTEGER*4 NL,NUP
	INTEGER*4 MNL,MNUP
	INTEGER*4 MNL_F,MNUP_F
	INTEGER*4 I,J,K,L,ML,LS,IOS,LINE_INDX
	INTEGER*4 ID,ISPEC
	INTEGER*4 ES_COUNTER
	INTEGER*4 NUM_ES_ITERATIONS
C
C Functions called
C
	INTEGER*4 ICHRLEN,ERROR_LU
	REAL*8 LAMVACAIR
	REAL*8 ATOMIC_MASS_UNIT
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL JWEIGHT,HWEIGHT,KWEIGHT,NWEIGHT
	EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
	EXTERNAL ICHRLEN,ERROR_LU,SPEED_OF_LIGHT
	EXTERNAL LAMVACAIR,ATOMIC_MASS_UNIT
C
C 
	CHARACTER FG_SOL_OPTIONS*10
	CHARACTER CMF_FORM_OPTIONS*10
	CHARACTER NEG_OPAC_OPTION*10
	CHARACTER METHOD*6
	CHARACTER FMT*120,N_TYPE*6
	CHARACTER SECTION*20
	CHARACTER TMP_KEY*20
	CHARACTER STRING*132
	CHARACTER EW_STRING*132
	CHARACTER TEMP_CHAR*132
C
C Global vectors:
C
	REAL*8 AMASS_ALL(NT)
C
C Arrays for performing LINE frequencies in numerical order
C
	REAL*8 VEC_FREQ(NLINE_MAX)
	REAL*8 VEC_STRT_FREQ(NLINE_MAX)
	REAL*8 VEC_OSCIL(NLINE_MAX)
	REAL*8 VEC_EINA(NLINE_MAX)
	REAL*8 VEC_ARAD(NLINE_MAX)
	REAL*8 VEC_DP_WRK(NLINE_MAX)
	REAL*8 VEC_VDOP_MIN(NLINE_MAX)
C
	INTEGER*4 VEC_INDX(NLINE_MAX)
	INTEGER*4 VEC_NL(NLINE_MAX)
	INTEGER*4 VEC_NUP(NLINE_MAX)
	INTEGER*4 VEC_MNL_F(NLINE_MAX)
	INTEGER*4 VEC_MNUP_F(NLINE_MAX)
	INTEGER*4 VEC_INT_WRK(NLINE_MAX)
	INTEGER*4 PROF_LIST_LOCATION(NLINE_MAX)
	CHARACTER*6 VEC_SPEC(NLINE_MAX)
	CHARACTER*6 VEC_TRANS_TYPE(NLINE_MAX)
	CHARACTER*12 PROF_TYPE(NLINE_MAX)
	CHARACTER*12 VEC_CHAR_WRK(NLINE_MAX)
	CHARACTER*80, ALLOCATABLE :: VEC_TRANS_NAME(:)
	INTEGER*4 N_LINE_FREQ
!
	CHARACTER*10 GLOBAL_LINE_PROF
	REAL*8 DOP_PROF_LIMIT
	REAL*8 VOIGT_PROF_LIMIT
	LOGICAL SET_PROF_LIMS_BY_OPACITY
	LOGICAL RD_STARK_FILE
C
C Arrays and variables for treating lines simultaneously.
C
	REAL*8 EINA(MAX_SIM)
	REAL*8 OSCIL(MAX_SIM)
	REAL*8 GLDGU(MAX_SIM)
	REAL*8 AMASS_SIM(MAX_SIM)
	REAL*8 FL_SIM(MAX_SIM)
	INTEGER*4 SIM_NL(MAX_SIM)
	INTEGER*4 SIM_NUP(MAX_SIM)
C
	REAL*8 CHIL_MAT(ND,MAX_SIM)
	REAL*8 ETAL_MAT(ND,MAX_SIM)
	REAL*8 BB_COR(ND,MAX_SIM)
!
	INTEGER*4 NUM_SIM_LINES
	INTEGER*4 SIM_INDX
	INTEGER*4 TMP_MAX_SIM
	REAL*8 OVER_FREQ_DIF
	REAL*8 EW
	REAL*8 EW_CUT_OFF
	REAL*8 CONT_INT
	LOGICAL OVERLAP
	LOGICAL SOBOLEV
C
C L refers to the lower level, U to the upper level.
C
	REAL*8 L_STAR_RATIO(ND,MAX_SIM)
	REAL*8 U_STAR_RATIO(ND,MAX_SIM)
C
C GLOBAL_LINE_SWITCH provides an option to handle all LINE by the same method.
C The local species setting only takes precedence when it is set to NONE.
C
	CHARACTER*10 GLOBAL_LINE_SWITCH
	REAL*8 FLUX_CAL_LAM_BEG
	REAL*8 FLUX_CAL_LAM_END
	LOGICAL SET_TRANS_TYPE_BY_LAM
	LOGICAL DO_SOBOLEV_LINES
C
C FLUX_CAL_ONLY provides a method for computing the continuous spectrum
C only (i.e. no linearization or population corrections):
C    To get a BLANKETED spectrum FLUX_CAL_ONLY should be set to
C       TRUE and GLOBAL_LINE_SWITCH to BLANK
C    To get a pure UNBLANKETED spectrum FLUX_CAL_ONLY should be set to
C       TRUE and GLOBAL_LINE_SWITCH to SOB
C
	LOGICAL FLUX_CAL_ONLY
C
C Indicates whether lines treated in SOB and CMF mode are allowed for when
C constructing the observers frame grid.
C
	LOGICAL SOB_FREQ_IN_OBS
	LOGICAL WRITE_ETA_AND_CHI
	LOGICAL WRITE_IP
	LOGICAL WRITE_FLUX
	LOGICAL WRITE_CMF_FORCE
	LOGICAL WRITE_SOB_FORCE
C
	REAL*8, ALLOCATABLE :: ETA_CMF_ST(:,:)
	REAL*8, ALLOCATABLE :: CHI_CMF_ST(:,:)
	REAL*8, ALLOCATABLE :: RJ_CMF_ST(:,:)
C
C
C Variables, vectors and rrays for treating lines simultaneously with the
C continuum.
C
	REAL*8 V_DOP 
	REAL*8 MAX_DOP
	REAL*8 FRAC_DOP
	REAL*8 dV_CMF_PROF
	REAL*8 dV_CMF_WING
	REAL*8 ES_WING_EXT
	REAL*8 R_CMF_WING_EXT
	INTEGER*4 LINES_THIS_FREQ(NCF_MAX)
	INTEGER*4 LINE_ST_INDX_IN_NU(NLINE_MAX)
	INTEGER*4 LINE_END_INDX_IN_NU(NLINE_MAX)
	REAL*8 LINE_PROF_SIM(ND,MAX_SIM)
C
	REAL*8 NU_MAX_OBS
	REAL*8 NU_MIN_OBS
	REAL*8 OBS_PRO_EXT_RAT
	REAL*8 FRAC_DOP_OBS
	REAL*8 dV_OBS_PROF
	REAL*8 dV_OBS_WING
	REAL*8 dV_OBS_BIG
C
	REAL*8 CONT_FREQ
	REAL*8 DELV_CONT
	LOGICAL COMPUTE_ALL_CROSS
	LOGICAL COMPUTE_NEW_CROSS
C
	LOGICAL EXTEND_FRM_SOL
	LOGICAL INSERT_FREQ_FRM_SOL
C
	CHARACTER*50 TRANS_NAME_SIM(MAX_SIM)
C
	LOGICAL RESONANCE_ZONE(MAX_SIM)
	LOGICAL END_RES_ZONE(MAX_SIM)
	LOGICAL LINE_STORAGE_USED(MAX_SIM)
C
	INTEGER*4 FREQ_INDX
	INTEGER*4 FIRST_LINE
	INTEGER*4 LAST_LINE
	INTEGER*4 LINE_LOC(NLINE_MAX)
	INTEGER*4 SIM_LINE_POINTER(MAX_SIM)
C
C 
C
C Opacity/emissivity
	REAL*8 CHI(ND)			!Continuum opacity (all sources)
	REAL*8 ETA(ND)			!Continuum emissivity (all sources)
	REAL*8 CHIL(ND)			!Line opacity (without prof.)
	REAL*8 ETAL(ND)			!Line emissivity (without prof.)
	REAL*8 ESEC(ND)			!Continuum electron scattering coef.
	REAL*8 ZETA(ND)			!Source func. (all except elec. scat.)
	REAL*8 THETA(ND)		!Elec. scat. source coef.
	REAL*8 SOURCE(ND)		!Complete source function.
	REAL*8 DTAU(NDMAX)		!Optical depth (used in error calcs)
C		DTAU(I)=0.5*(CHI(I)+CHI(I+1))*(Z(I)-Z(I+1))
	REAL*8 dCHIdR(NDMAX) 		!Derivative of opacity.
C
	REAL*8 P(NP)
C
	REAL*8 CHI_CONT(ND)
	REAL*8 ETA_CONT(ND)
C
C These parameters are used when computing J and the variation of J.
C
	LOGICAL DO_CLUMP_MODEL
	REAL*8 CHI_CLUMP(ND)		!==CHI(I)*CLUMP_FAC(I)
	REAL*8 ETA_CLUMP(ND)		!==ETA(I)*CLUMP_FAC(I)
	REAL*8 ESEC_CLUMP(ND)		!==ESEC(I)*CLUMP_FAC(I)
C
C Variables to limit the computation of the continuum opacities and 
C emissivities. 
C
	REAL*8 EMHNUKT_CONT(ND)
	REAL*8 ETA_C_EVAL(ND)
	REAL*8 CHI_C_EVAL(ND)
C
C 
C
	LOGICAL DO_LEV_DISSOLUTION
	REAL*8 Z_POP(NT)		!Ionic charge for each species
C
C Variables etc for computation of continuum in comoving frame.
C
	LOGICAL CONT_VEL
	LOGICAL FIRST_FREQ
	LOGICAL COHERENT_ES
	LOGICAL RD_COHERENT_ES
	LOGICAL USE_OLDJ_FOR_ES
	LOGICAL NEW_FREQ
	REAL*8 dLOG_NU			!Step in frequency in Log plane
	REAL*8 FEDD_PREV(NDMAX)
	REAL*8 GEDD_PREV(NDMAX)
	REAL*8 N_ON_J(NDMAX)
	REAL*8 N_ON_J_PREV(NDMAX)
	REAL*8 JNU_PREV(NDMAX)
	REAL*8 RSQHNU_PREV(NDMAX)
	REAL*8 GEDD(NDMAX)
	REAL*8 RSQHNU(NDMAX)
	REAL*8 CHI_PREV(ND)
	REAL*8 ETA_PREV(ND)
	REAL*8 HBC_CMF(3),HBC_PREV(3)
	REAL*8 NBC_CMF(3),NBC_PREV(3),INBC_PREV
C
C Quadrature weights.
	REAL*8 FQW(NCF_MAX)		!Frequency weights
	REAL*8 AQW(ND,NP)		!Angular quad. weights. (indep. of v)
	REAL*8 HQW(ND,NP)		!Angular quad. weights. (indep. of v)
                                        !for flux integration.
	REAL*8 KQW(ND,NP)		!Angular quad. weights for K integration.
	REAL*8 HMIDQW(ND,NP)		!Angular quad. weights. (indep. of v)
                                        !for flux integration. Defined at the
                                        !mid points of the radius mesh.
	REAL*8 NMIDQW(ND,NP)		!Angular quad. weights. (indep. of v)
                                        !for N integration. Defined at the
                                        !mid points of the radius mesh.
C
C Continuum matrices
	REAL*8 WM(ND,ND)		!Coef. matrix of J & %J vector
	REAL*8 FB(ND,ND)		!Coef. of J & %J vects in angular equ.
C
C Transfer equation vectors
	REAL*8 TA(NDMAX)
	REAL*8 TB(NDMAX)
	REAL*8 TC(NDMAX)
	REAL*8 XM(NDMAX)		!R.H.S. (SOURCE VECTOR)
C
C 
C
C Arrays and variables for computation of the continuum intensity
C using Eddington factors. This is separate to the "inclusion of
C additional points".
C
	LOGICAL EDDINGTON
	LOGICAL EDD_CONT
	LOGICAL EDD_LINECONT
	REAL*8 FEDD(NDMAX)
	REAL*8 QEDD(NDMAX)
C
	LOGICAL AT_LEAST_ONE_NEG_OPAC
	LOGICAL DIF
	LOGICAL MID
	LOGICAL INCL_TWO_PHOT
	LOGICAL CHECK_LINE_OPAC
	LOGICAL FIRST
	LOGICAL THK_CONT
	LOGICAL THK_LINE
	LOGICAL TRAPFORJ
	LOGICAL RDTHK_CONT
	LOGICAL NEG_OPACITY(ND)
	LOGICAL FIRST_NEG
	LOGICAL LAMBDA_ITERATION
	LOGICAL LST_ITERATION
	LOGICAL LST_DEPTH_ONLY
C 
C
C
C X-ray variables.
C
	LOGICAL XRAYS
	REAL*8 FILL_FAC_XRAYS,T_SHOCK,V_SHOCK
	REAL*8 FILL_VEC_SQ(ND)
	REAL*8 XRAY_LUM(ND)
	REAL*8 GFF,XCROSS
	EXTERNAL GFF,XCROSS
C
C
C
C ACESS_F is the current record we are writing in EDDFACTOR.
C EDD_CONT_REC is the record in EDDFACTOR which points to the first
C record containing the continuum values.
C
	INTEGER*4 ACCESS_F
	INTEGER*4, PARAMETER :: EDD_CONT_REC=3
	LOGICAL COMPUTE_EDDFAC
C
C Arrays and variables required for  additional points into
C the depth grid. Allows an increase in program accuracy to overcome
C rapid ionization changes. The vectors are used throughout,
C and hence should not be equivalenced or put into scratch.
C
	LOGICAL ACCURATE,INACCURATE
	LOGICAL THIS_FREQ_EXT  		!Frequency specific.
	LOGICAL ALL_FREQ
        REAL*8 ACC_FREQ_END
	INTEGER*4 NPINS			!Points inserted for error calc.
	INTEGER*4 ST_INTERP_INDX	!Interp from ST_INT.. to END_INTERP..
	INTEGER*4 END_INTERP_INDX
	CHARACTER*10 INTERP_TYPE
C
C ND-DEEP to DEEP we use a quadratic interpolation scheme so as to try
C and preserve "FLUX" in the diffusion approximation.
C
	INTEGER*4 DEEP
	INTEGER*4 NDEXT,NCEXT,NPEXT
	INTEGER*4 INDX(NDMAX),POS_IN_NEW_GRID(ND)
	REAL*8 COEF(0:3,NDMAX)
	REAL*8 INBC,HBC_J,HBC_S			!Bound. Cond. for JFEAU
	REAL*8 ACC_EDD_FAC
C
	REAL*8 REXT(NDMAX),PEXT(NPMAX),VEXT(NDMAX)
	REAL*8 TEXT(NDMAX),SIGMAEXT(NDMAX)
	REAL*8 CHIEXT(NDMAX),ESECEXT(NDMAX),ETAEXT(NDMAX)
	REAL*8 ZETAEXT(NDMAX),THETAEXT(NDMAX)
	REAL*8 RJEXT(NDMAX),RJEXT_ES(NDMAX)
	REAL*8 FOLD(NDMAX),FEXT(NDMAX),QEXT(NDMAX),SOURCEEXT(NDMAX)
C
C
	REAL*8 F2DAEXT(NDMAX,NDMAX)     !These arrays don't need to be
C
C If required, these arrays shoukd have size NDEXT*NPEXT
C
	REAL*8, ALLOCATABLE :: AQWEXT(:,:)	!Angular quad. weights. (indep. of v)
	REAL*8, ALLOCATABLE :: HQWEXT(:,:)	!Angular quad. weights for flux integration.
	REAL*8, ALLOCATABLE :: KQWEXT(:,:)	!Angular quad. weights for K integration.
	REAL*8, ALLOCATABLE :: HMIDQWEXT(:,:)	!Angular quad. weights for flux integration.
	REAL*8, ALLOCATABLE :: NMIDQWEXT(:,:)	!Angular quad. weights for flux integration.
C
C Arrays for calculating mean opacities.
C
	REAL*8 FLUXMEAN(ND) 		!Flux mean opacity
	REAL*8 LINE_FLUXMEAN(ND) 	!Flux mean opacity due to lines.
	REAL*8 ROSSMEAN(ND)  		!Rosseland mean opacity
	REAL*8 INT_dBdT(ND)  		!Integral of dB/dT over nu 
C                                            (to calculate ROSSMEAN)
	REAL*8 FORCE_MULT(ND)
	REAL*8 NU_FORCE
	REAL*8 NU_FORCE_FAC
	INTEGER*4 N_FORCE
	INTEGER*4 ML_FORCE
	LOGICAL TMP_LOG
C
C Other arrays
	REAL*8 Z(NDMAX)			!Z displacement along a given array
	REAL*8 EMHNUKT(ND)		!EXP(-hv/kT)
	REAL*8 RLUMST(ND)		!Luminosity as a function of depth
	REAL*8 J_INT(ND)		!Frequency integrated J
	REAL*8 K_INT(ND)		!Frequency integrated K
	REAL*8 K_MOM(ND)		!Frequency dependent K moment
	REAL*8 SOB(ND)   	    	!Used in computing continuum flux
	REAL*8 RJ(ND)			!Mean intensity
	REAL*8 RJ_ES(ND)		!Convolution of RJ with e.s. R(v'v')
C
C Line variables.
C
	REAL*8 VAL_DO_NG
	REAL*8 RP
	REAL*8 VINF
	REAL*8 VTURB_FIX
	REAL*8 VTURB_MIN
	REAL*8 VTURB_MAX
	REAL*8 VTURB_VEC(ND)
	REAL*8 MAX_DEL_V_RES_ZONE(ND)
!
	REAL*8 OBS_TAU_MAX
	REAL*8 OBS_ES_DTAU
	CHARACTER*10 OBS_INT_METHOD
C
C Continuum frequency variables and arrays.
C
	REAL*8 NU(NCF_MAX)		!Continuum and line frequencies
	REAL*8 NU_EVAL_CONT(NCF_MAX)	!Frequencies to evaluate continuum
	REAL*8 OBS(NCF_MAX)		!Observers spectrum
	REAL*8 MIN_CONT_FREQ 		!Minimum continuum frequency.
	REAL*8 MAX_CONT_FREQ    	!Maximum continuum frequency.
	REAL*8 SMALL_FREQ_RAT 		!Fractional spacing for small frequencies'
	REAL*8 dFREQ_bf_MAX		!Maximum spacing close to bf edge.
	REAL*8 BIG_FREQ_AMP		!Amplification factor
	REAL*8 dV_LEV_DIS		!dV on low side of bound-free edge.
	REAL*8 AMP_DIS			!Amplification factor
	REAL*8 MIN_FREQ_LEV_DIS		!Minimum frequency for lev dissolution.
C
C Parameters, vectors, and arrays for computing the observed flux.
C
	INTEGER*4, PARAMETER :: NST_CMF=2000
	INTEGER*4 NP_OBS_MAX
	INTEGER*4 NP_OBS
	REAL*8  NU_STORE(NST_CMF)
	REAL*8 V_AT_RMAX		!Used if we extend the atmosphere.
	REAL*8 RMAX_OBS
	REAL*8 H_OUT,H_IN
C
C We allocate memory for the following vectors as we use them for the regular
C flux computation, and when extra depth points are inserted (ACCURATE=.TRUE.)
C
	REAL*8, ALLOCATABLE :: IPLUS_STORE(:,:)
	REAL*8, ALLOCATABLE :: P_OBS(:)
	REAL*8, ALLOCATABLE :: IPLUS(:)
	REAL*8, ALLOCATABLE :: MU_AT_RMAX(:)
	REAL*8, ALLOCATABLE :: HQW_AT_RMAX(:)
C
C Supercedes OBS
C
	INTEGER*4 N_OBS
	REAL*8 OBS_FREQ(NCF_MAX)		!Since N_OBS < NCF =< NCF_MAX
	REAL*8 OBS_FLUX(NCF_MAX)
	LOGICAL FIRST_OBS_COMP
C
C Indicates approximate frequencies for which TAU at outer boundary is written
C to OUTGEN on the last iteration.
C
C They are the He2 ege, NIII/CIII egde, HeI, HI, HI(N=2).
C
	INTEGER*4, PARAMETER :: N_TAU_EDGE=5
	REAL*8 TAU_EDGE(N_TAU_EDGE)
	DATA TAU_EDGE/13.16D0,11.60D0,5.95D0,3.29D0,0.83D0/
C
C
C Open output file for all errors and comments.
C
	LUER=ERROR_LU()
	CALL GEN_ASCI_OPEN(LUER,'OUT_FLUX','UNKNOWN','APPEND',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening OUTGEN in CMFGEN, IOS=',IOS
	  STOP
	END IF
C
C Check whether EQUATION LABELLING is consistent. ' I ' is used as the
C number of the current equation. We also set the variable SPEC_PRES which 
C indicates whether at least one ioization stage of a species is present.
C It is used to determine, foe example,  whether a number conservation 
C equation is required.
C
	I=1
!
!???????????????????????
!
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    SPECIES_PRES(ISPEC)=.FALSE.
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      CALL CHK_EQ_NUM( ATM(ID)%XzV_PRES, ATM(ID)%EQXzV,
	1                      ATM(ID)%NXzV,
	1                      I,SPECIES_PRES(ISPEC),TRIM(ION_ID(ID)))
	    END DO
	    CALL CHK_EQ_NUM(SPECIES_PRES(ISPEC),EQ_SPECIES(ISPEC),IONE,I,
	1                          SPECIES_PRES(ISPEC),TRIM(SPECIES(ISPEC)))
	  END IF
	END DO
!
	IF(EQNE .NE. I)THEN
	  WRITE(LUER,*)'Error - EQNE has wrong value in CMFGEN'
	  STOP
	END IF
	IF(NT .NE. I+1)THEN
	  WRITE(LUER,*)'Error - NT has wrong value in CMFGEN'
	  STOP
	END IF
C 
C
!	LAMBDA_ITERATION=.FALSE.
	LAMBDA_ITERATION=.TRUE.
	LST_ITERATION=.FALSE.
	LST_DEPTH_ONLY=.FALSE.
	ACCESS_F=5
	RP=R(ND)
	VINF=V(1)
	TSTAR=T(ND)
	C_KMS=SPEED_OF_LIGHT()/1.0D+05
C
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
C
C MAXCH and VAL_DO_NG are set so that they defined for the TEST in
C COMP_JCONT_V?.INC whether to do an accurate flux calculation. An 
C accurate flux calculation can be avoided by doing a LAMBDA iteration.
C
	MAXCH=0.0D0
	VAL_DO_NG=5.0D0
C
C Set the vector Z_POP to contain the ionic charge for each species.
C
	Z_POP(1:NT)=0.0D0
C
	DO ID=1,NUM_IONS
	  CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1                       ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
	END DO
C
C Store atomic masses in vector of LENGTH NT for later use by line 
C calculations. G_ALL and LEVEL_ID  are no longer used due to the use
C of super levels.
C
	AMASS_ALL(1:NT)=0.0D0
!
! We also set the mass of the ion corresponding to XzV to AT_MASS for each
! species . Thus the range is EQXzV:EQXzV+NXzV, NOT EQXzV:EQXzV+NXzV-1. Note 
! that HeIII_PRES is always false, even when HeI and HeII are both present.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)AMASS_ALL(ATM(ID)%EQXzV:ATM(ID)%EQXzV+ATM(ID)%NXzV)=
	1              AT_MASS(SPECIES_LNK(ID))
	END DO
C
C
	CALL GEN_ASCI_OPEN(LUIN,'CMF_FLUX_PARAM','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening CMF_FLUX_PARAM in CMFGEN, IOS=',IOS
	  STOP
	END IF
	CALL GEN_ASCI_OPEN(LUMOD,'OUT_PARAMS','UNKNOWN',' ','WRITE',IZERO,IOS)
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUMOD)
C 
C
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
	  CALL RD_STORE_LOG(COMPUTE_ALL_CROSS,'CROSS',L_TRUE,
	1            'Compute all photoionization cross-sections?')
	  CALL RD_STORE_DBLE(DELV_CONT,'V_CROSS',L_TRUE,
	1            'Max. vel. sep. (km/s) between evaluations of all'//
	1            '  phot. cross-sections?')
C
	  CALL RD_STORE_LOG(DIF,'DIF',L_TRUE,
	1            'Use Diffusion approximation at inner boundary ?')
	  CALL RD_STORE_INT(NUM_ES_ITERATIONS,'NUM_ES',L_TRUE,
	1            'Number of electron scattering iterations?')
	  CALL RD_STORE_LOG(RD_COHERENT_ES,'COH_ES',L_TRUE,
	1            'Assume coherent electron scattering? ')
	  CALL RD_STORE_LOG(USE_OLDJ_FOR_ES,'OLD_J',L_TRUE,
	1            'Use old file to provide initial estimate of J_ES?')
	  COHERENT_ES=RD_COHERENT_ES
C
	  CALL RD_STORE_NCHAR(METHOD,'METHOD',ISIX,L_TRUE,
	1           'Which method for continuum tau'//
	1          ' loglog, loglin, linear or zero ?')
	  CALL RD_STORE_NCHAR(N_TYPE,'N_TYPE',ISIX,L_TRUE,
	1           'Method for to handle N for MOM_J_CMF -- '//
	1           'N_ON_J, MIXED, or G_ONLY')
	  CALL RD_STORE_NCHAR(FG_SOL_OPTIONS,'FG_OPT',ITEN,L_TRUE,
	1           'Solution options for FG_J_CMF: DIFF/INS and INT/INS')
	  CALL RD_STORE_LOG(RDTHK_CONT,'THK_CONT',L_TRUE,
	1           'Use thick boundary condition for continuum ? ')
	  CALL RD_STORE_LOG(TRAPFORJ,'TRAP_J',L_TRUE,
	1           'Use trapazoidal weights to compute J? ')
C
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_DBLE(VTURB_FIX,'VTURB_FIX',L_TRUE,
	1      'Doppler velocity for DOP_FIX Doppler profiles (km/s)')
	  CALL RD_STORE_DBLE(VTURB_MIN,'VTURB_MIN',L_TRUE,
	1      'Minimum turbulent velocity for Doppler profile (km/s)')
	  CALL RD_STORE_DBLE(VTURB_MAX,'VTURB_MAX',L_TRUE,
	1      'Maximum turbulent velocity for Doppler profile (km/s)')
C
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
C
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
C
	  CALL RD_STORE_DBLE(OBS_TAU_MAX,'TAU_MAX',L_TRUE,
	1    'Optical depth at which observers frame integration is terminated')
	  CALL RD_STORE_DBLE(OBS_ES_DTAU,'ES_DTAU',L_TRUE,
	1      'Maximum increments in e.s. optical depth scale')
	  CALL RD_STORE_NCHAR(OBS_INT_METHOD,'INT_METH',ITEN,L_TRUE,
	1            'Integration method for computing I along ray')
C
	  WRITE(LUMOD,'()')
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
	  CALL RD_STORE_LOG(WRITE_ETA_AND_CHI,'WR_ETA',L_TRUE,
	1        'Output ETA and CHI? ')
	  CALL RD_STORE_LOG(WRITE_FLUX,'WR_FLUX',L_TRUE,
	1        'Output Flux as a function of depth? ')
	  CALL RD_STORE_LOG(WRITE_CMF_FORCE,'WR_CMF_FORCE',L_TRUE,
	1        'Output CMF line-force multiplier as a function of depth? ')
	  CALL RD_STORE_LOG(WRITE_SOB_FORCE,'WR_SOB_FORCE',L_TRUE,
	1        'Output SOBOLEV line-force multiplier as a function of depth? ')
	  CALL RD_STORE_LOG(WRITE_IP,'WR_IP',L_TRUE,
	1        'Output I as a functio of p and frequency?')
C
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_NCHAR(GLOBAL_LINE_SWITCH,'GLOBAL_LINE',ISIX,L_TRUE,
	1            'Global switch to indicate handeling of line')
	  CALL SET_CASE_UP(GLOBAL_LINE_SWITCH,IZERO,IZERO)
	  IF( GLOBAL_LINE_SWITCH(1:3) .NE. 'SOB' .AND.
	1       GLOBAL_LINE_SWITCH(1:3) .NE. 'CMF' .AND.
	1       GLOBAL_LINE_SWITCH(1:4) .NE. 'LIST' .AND.
	1       GLOBAL_LINE_SWITCH(1:8) .NE. 'LIST_VGT' .AND.
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
C
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
C
	  CALL RD_STORE_LOG(INCL_TWO_PHOT,'INC_TWO',L_TRUE,
	1           'Include two photon transitions?')
	  CALL RD_STORE_LOG(XRAYS,'INC_XRAYS',L_TRUE,
	1           'Include X-ray emission')
	  CALL RD_STORE_DBLE(FILL_FAC_XRAYS,'FIL_FAC',L_TRUE,
	1           'Filling factor for X-ray emission')
	  CALL RD_STORE_DBLE(T_SHOCK,'T_SHOCK',L_TRUE,
	1           'Shock T for X-ray emission')
	  CALL RD_STORE_DBLE(V_SHOCK,'V_SHOCK',L_TRUE,
	1           'Cut off velocity for X-ray emission')
C
	  IF(GLOBAL_LINE_SWITCH .EQ. 'NONE')THEN
	    WRITE(LUMOD,'()')
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        TEMP_CHAR='TRANS_'//ION_ID(ID)
	        ATM(ID)%XzV_TRANS_TYPE='SOB'
	        STRING='Method for treating '//TRIM(ION_ID(ID))//' lines?'
	        CALL RD_STORE_NCHAR( ATM(ID)%XzV_TRANS_TYPE,TEMP_CHAR,
	1                          ITEN,L_FALSE,STRING)
	      END IF
	    END DO
	  ELSE
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        ATM(ID)%XzV_TRANS_TYPE=GLOBAL_LINE_SWITCH
	      END IF
	    END DO
	  END IF
C
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
	  CALL RD_STORE_LOG(SET_PROF_LIMS_BY_OPACITY,'OPAC_LIMS',L_TRUE,
	1           'Set prof limits by line to cont. ratio?')
	  CALL RD_STORE_DBLE(DOP_PROF_LIMIT,'DOP_LIM',L_TRUE,
	1           'Edge limits for Doppler line profile')
	  CALL RD_STORE_DBLE(VOIGT_PROF_LIMIT,'VOIGT_LIM',L_TRUE,
	1           'Edge limits for Voigt line profile')
!
! Verify validity of profile option. We also check whether we need to leed
! in the file which links certain types of profiles to individual lines.
!
	  IF(GLOBAL_LINE_PROF .EQ. 'NONE')THEN
	    WRITE(LUMOD,'()')
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        TEMP_CHAR='PROF_'//ION_ID(ID)
	        STRING='Intrinsic profile for treating '//TRIM(ION_ID(ID))//' lines?'
	        CALL RD_STORE_NCHAR(ATM(ID)%XzV_PROF_TYPE,TEMP_CHAR,ITEN,
	1               L_TRUE,STRING)
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
	
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_LOG(ACCURATE,'INC_GRID',L_TRUE,
	1          'Increase grid size to improve accuracy? ')
	   CALL RD_STORE_LOG(ALL_FREQ,'ALL_FREQ',L_TRUE,
	1          'Increase accuracy for all frequencies?')
	   CALL RD_STORE_DBLE(ACC_FREQ_END,'ACC_END',L_TRUE,
	1          'Increase accuracy for all frequencies > ACC_END?')
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
C Next two variables apply for both ACCURATE and EDDINGTON.
C
	  WRITE(LUMOD,'()')
	  CALL RD_STORE_LOG(COMPUTE_EDDFAC,'COMP_F',L_TRUE,
	1      'Compute new Eddington factors (f)')
	  CALL RD_STORE_DBLE(ACC_EDD_FAC,'ACC_F',L_TRUE,
	1      'Accuracy with which to compute the eddington factor f')
!
	  DO ISPEC=1,NUM_SPECIES
	    TMP_KEY='SCL_'//TRIM(SPECIES(ISPEC))//'_ABUND'
	    CALL RD_STORE_DBLE(ABUND_SCALE_FAC(ISPEC),TMP_KEY,L_FALSE,
	1      'Factor to scale abundance by')
	  END DO
!
C
	CLOSE(UNIT=7)
C 
!
! Scale abunances of species if desired. This should only be done for exploratory
! spectral calculations, and only for IMPURITY species (i.e. not H or He).
!
	  DO ISPEC=1,NUM_SPECIES
	    IF(SPECIES_PRES(ISPEC) .AND. ABUND_SCALE_FAC(ISPEC) .NE. 1.0D0)THEN
	      WRITE(LUER,'(A)')' '
	      WRITE(LUER,'(A)')'******************Warning**********************'
	      WRITE(LUER,'(A)')'Abundance of species ',TRIM(SPECIES(ISPEC)),' scaled'
	      WRITE(LUER,'(A)')'******************Warning**********************'
	      WRITE(LUER,'(A)')' '
	      AT_ABUND(ISPEC)=AT_ABUND(ISPEC)*ABUND_SCALE_FAC(ISPEC)
	      POP_SPECIES(:,ISPEC)=POP_SPECIES(:,ISPEC)*ABUND_SCALE_FAC(ISPEC)
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        ATM(ID)%XzV_F=ATM(ID)%XzV_F*ABUND_SCALE_FAC(ISPEC)
	        ATM(ID)%DXzV_F=ATM(ID)%DXzV_F*ABUND_SCALE_FAC(ISPEC)
	      END DO
	    END IF
	  END DO
!
! Evaluate constants used for level-dissolution calculations.
!
	  CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! Read in file contaning line links to STARK tables.
!
	IF(RD_STARK_FILE .OR. GLOBAL_LINE_PROF(1:4) .EQ. 'LIST')THEN
	  CALL RD_STRK_LIST(LUIN)
	END IF
!                     
! We now need to compute the populations for the model atom with Super-levels.
! We do this in reverse order (i.e. highest ionization stage first) in order
! that we the ion density for the lower ionization stage is available for
! the next call.
!
! For 1st call to FULL_TO_SUP, Last line contains FeX etc as FeXI not installed.
! We only do to NUM_IONS-1 to avoid access arror, and since the ion
! corresponding to NUM_IONS contains only 1 level and is done with NUM_IONS-1
!
	  DO ID=1,NUM_IONS-1
	    CALL FULL_TO_SUP( 
	1        ATM(ID)%XzV,        ATM(ID)%NXzV,   ATM(ID)%DXzV,
	1        ATM(ID)%XzV_PRES,   ATM(ID)%XzV_F,
	1        ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,
	1        ATM(ID+1)%XzV,      ATM(ID+1)%NXzV, ATM(ID+1)%XzV_PRES, ND)
	  END DO
C
C Store all quantities in POPS array. 
C
	  DO ID=1,NUM_IONS
	    CALL IONTOPOP(POPS, ATM(ID)%XzV,        ATM(ID)%DXzV,  ED,T,
	1        ATM(ID)%EQXzV, ATM(ID)%NXzV,NT,ND, ATM(ID)%XzV_PRES)
	  END DO
C
C This routine not only evaluates the LTE populations of both model atoms, but
C it also evaluates the dln(LTE Super level Pop)/dT.
C
	INCLUDE 'EVAL_LTE_INC_V4.INC'
C
C compute the turbulent velocity as a function of depth.
C
	VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*V(1:ND)/V(1)
C
	TA(1:ND)=ABS( CLUMP_FAC(1:ND)-1.0D0 )
	T1=MAXVAL(TA)
	DO_CLUMP_MODEL=.FALSE.
	IF(T1 .GT. 1.0D-05)DO_CLUMP_MODEL=.TRUE.
C
C
C Compute profile frequencies such that for the adopted doppler
C velocity the profile ranges from 5 to -5 doppler widths.
C This section needs to be rewritten if we want the profile to
C vary with depth.
	  FIRST=.TRUE.		!Check cross section at edge is non-zero.
	  NCF=0 		!Initialize number of continuum frequencies.
!
	  DO ID=1,NUM_IONS
	    CALL SET_EDGE_FREQ_V3(ID,OBS,NCF,NCF_MAX,
	1           ATM(ID)%EDGEXzV_F,  ATM(ID)%NXzV_F, ATM(ID)%XzV_PRES,
	1           ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV,
	1           ATM(ID)%N_XzV_PHOT)
	  END DO

	  IF(XRAYS)THEN
	    DO ID=1,NUM_IONS-1
	      CALL SET_X_FREQ(OBS,NCF,NCF_MAX,AT_NO(SPECIES_LNK(ID)),
	1           ATM(ID)%ZXzV, ATM(ID)%XzV_PRES, ATM(ID+1)%XzV_PRES)
	    END DO
	  END IF
C
C Now insert addition points into frequency array. WSCI is used as a
C work array - okay since of length NCF_MAX, and zeroed in QUADSE.
C OBSF contains the bound-free edges - its contents are zero on
C subroutine exit. J is used as temporary variable for the number of
C frequencies transmitted to SET_CONT_FREQ. NCF is returned as the number 
C of frequency points. FQW is used a an integer array for the sorting ---
C we know it has the correct length since it is the same size as NU.
C LUIN --- Used as temporary LU (opened and closed).
C
	  J=NCF
	  CALL SET_CONT_FREQ(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        J,NCF,NCF_MAX,LUIN)
C                                             
C                         
C 
C Set up lines that will be treated with the continuum calculation.
C This section of code is also used by the code treating purely lines
C (either single transition Sobolev or CMF, or overlapping Sobolev).
C
C To define the line transitions we need to operate on the FULL atom models.
C We thus perform separate loops for each species. VEV_TRANS_NAME is
C allocated temporaruly so that we can output the full transitions name
C to TRANS_INFO.
C
	ML=0			!Initialize line counter.
	ALLOCATE (VEC_TRANS_NAME(NLINE_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for VEC_TRANS_NAME'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
C
	ESEC(1:ND)=6.65D-15*ED(1:ND)
!
! The onlye species not present is the ion corresponding
! to the last ioization stage considered.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO MNL=1, ATM(ID)%NXzV_F-1
	      DO MNUP=MNL+1, ATM(ID)%NXzV_F
	        NL=ATM(ID)%F_TO_S_XzV(MNL)+ATM(ID)%EQXzV-1
	        NUP=ATM(ID)%F_TO_S_XzV(MNUP)+ATM(ID)%EQXzV-1
	        IF(ATM(ID)%AXzV_F(MNL,MNUP) .NE. 0)THEN
	          ML=ML+1
	          IF(ML .GT. NLINE_MAX)THEN
	            WRITE(LUER,*)'NLINE_MAX is too small in CMFGEN'
	            STOP
	          END IF
	          VEC_FREQ(ML)=ATM(ID)%EDGEXzV_F(MNL)-ATM(ID)%EDGEXzV_F(MNUP)
	          VEC_SPEC(ML)=ION_ID(ID)
	          VEC_NL(ML)=NL
	          VEC_NUP(ML)=NUP     
	          VEC_MNL_F(ML)=MNL
	          VEC_MNUP_F(ML)=MNUP     
	          VEC_OSCIL(ML)=ATM(ID)%AXzV_F(MNL,MNUP)
	          VEC_EINA(ML)=ATM(ID)%AXzV_F(MNUP,MNL)
	          VEC_ARAD(ML)= ATM(ID)%ARAD(MNL)+ATM(ID)%ARAD(MNUP)
	          VEC_TRANS_TYPE(ML)=ATM(ID)%XzV_TRANS_TYPE
	          VEC_TRANS_NAME(ML)=TRIM(VEC_SPEC(ML))//
	1             '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP))//'-'//
	1             TRIM(ATM(ID)%XzVLEVNAME_F(MNL))//')'
	          T1=VEC_OSCIL(ML)*OPLIN
	          T2=ATM(ID)%GXzV_F(MNL)/ATM(ID)%GXzV_F(MNUP)
	          DO I=1,ND
	            CHIL(I)=ABS(T1*(ATM(ID)%XzV_F(MNL,I)-T2*ATM(ID)%XzV_F(MNUP,I)))
	          END DO
	          PROF_TYPE(ML)=ATM(ID)%XzV_PROF_TYPE
	          IF(GLOBAL_LINE_PROF .NE. 'NONE')PROF_TYPE(ML)=GLOBAL_LINE_PROF
	          T1=0.0D0; T2=0.0D0
	          CALL SET_PROF_LIMITS_V2(VEC_STRT_FREQ(ML),VEC_VDOP_MIN(ML),
	1             CHIL,ED,T,VTURB_VEC,ND,PROF_TYPE(ML),PROF_LIST_LOCATION(ML),
	1             VEC_FREQ(ML),MNL,MNUP,
	1             VEC_SPEC(ML),AT_MASS(SPECIES_LNK(ID)), ATM(ID)%ZXzV,
	1             VEC_ARAD(ML),T2,VTURB_FIX,                !T1,T2: Garbage at presnet
	1             DOP_PROF_LIMIT,VOIGT_PROF_LIMIT,SET_PROF_LIMS_BY_OPACITY)
	        END IF
	      END DO
	    END DO
	  END IF
	END DO
	N_LINE_FREQ=ML
C
C 
C
C GLOBAL_LINE_SWITCH provides an option to handle all LINE by the same method.
C The local species setting only takes precedence when it is set to NONE.
C
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
	CALL INIT_PROF_MODULE(ND,NLINES_PROF_STORE,NFREQ_PROF_STORE)
C
C If desired, we can set transitions with:
C      wavelengths > FLUX_CAL_LAM_END (in A) to the SOBOLEV option.
C      wavelengths < FLUX_CAL_LAM_BEG (in A) to the SOBOLEV option.
C
C The region defined by FLUX_CAL_LAM_BEG < LAM < FLUX_CAL_LAM_END will be computed using
C transition types determined by the earlier species and global options.
C
C Option has 2 uses:
C
C 1. Allows use of SOBOLEV approximation in IR where details of radiative
C    transfer is unimportant. In this case FLUX_CAL_LAM_BEG should be set to zero.
C 2. Allows a full flux calculation to be done in a limited wavelength region
C    as defined by FLUX_CAL_LAM_END and FLUX_CAL_LAM_BEG. 
C
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
C
	DO ML=1,N_LINE_FREQ
	  IF(VEC_TRANS_TYPE(ML) .NE. 'BLANK')VEC_STRT_FREQ(ML)=VEC_FREQ(ML)
	END DO
C
C Sort lines into numerically decreaing frequency. This is used for
C outputing TRANS_INFO file. Need to sort all the VECTORS, as they
C are linked.
C
	CALL INDEXX(N_LINE_FREQ,VEC_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_VDOP_MIN,VEC_INDX,VEC_DP_WRK)
C
	CALL SORTINT(N_LINE_FREQ,VEC_NL,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NUP,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,PROF_LIST_LOCATION,VEC_INDX,VEC_INT_WRK)
!
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,PROF_TYPE,VEC_INDX,VEC_CHAR_WRK)
C
	I=160	!Record length - allow for long names
	CALL GEN_ASCI_OPEN(LUIN,'TRANS_INFO','UNKNOWN',' ','WRITE',I,IOS)
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
	    T2=C_KMS*(VEC_FREQ(ML-1)-VEC_FREQ(ML))/VEC_FREQ(ML)
	    IF(T2 .GT. C_KMS)T2=C_KMS
	    IF(T1 .LT. 1.0E+04)THEN
	      WRITE(LUIN,
	1      '(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,2X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	    ELSE             
	      WRITE(LUIN,
	1      '(1X,I6,2(1X,I6),2X,F10.6,2X,1P,E10.4,0P,2X,F10.2,4X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	    END IF
	  END DO
	CLOSE(UNIT=LUIN)
	DEALLOCATE (VEC_TRANS_NAME)
C
C Get lines and arrange in numerically decreasing frequency according to
C the START frequency of the line. This will allow us to consider line overlap,
C and to include lines with continuum frequencies so that the can be handled 
C automatically.
C
	CALL INDEXX(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_VDOP_MIN,VEC_INDX,VEC_DP_WRK)
C
	CALL SORTINT(N_LINE_FREQ,VEC_NL,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_NUP,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNL_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,VEC_MNUP_F,VEC_INDX,VEC_INT_WRK)
	CALL SORTINT(N_LINE_FREQ,PROF_LIST_LOCATION,VEC_INDX,VEC_INT_WRK)
!
	CALL SORTCHAR(N_LINE_FREQ,VEC_SPEC,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,VEC_TRANS_TYPE,VEC_INDX,VEC_CHAR_WRK)
	CALL SORTCHAR(N_LINE_FREQ,PROF_TYPE,VEC_INDX,VEC_CHAR_WRK)
C
C
C
C We have found all lines. If we are doing a blanketing calculation for this
C line we insert them into the continuum frequency set, otherwise the
C line is not included.
C
	DO ML=1,NCF                                          
	  FQW(ML)=NU(ML)	!FQW has temporary storage of continuum freq.
	END DO                    
	V_DOP=MINVAL(VEC_VDOP_MIN)
	CALL INS_LINE_V5(  NU,LINES_THIS_FREQ,I,NCF_MAX,
	1		  VEC_FREQ,VEC_STRT_FREQ,VEC_VDOP_MIN,VEC_TRANS_TYPE,
	1                 LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,N_LINE_FREQ,
	1                 FQW,NCF,FRAC_DOP,VINF,dV_CMF_PROF,dV_CMF_WING,
	1                 ES_WING_EXT,R_CMF_WING_EXT,L_FALSE )
C
	K=NCF		!# of continuum frequencies: Need for DET_MAIN...
	NCF=I		!Revised
C
	WRITE(LUER,*)' '
	WRITE(LUER,'(A,1X,I7)')' Number of line frequencies is:',N_LINE_FREQ
	WRITE(LUER,'(A,6X,I7)')' Number of continuum frequencies is:',NCF
	WRITE(LUER,*)' '
C
	V_DOP=MAXVAL(VEC_VDOP_MIN)
	CALL DET_MAIN_CONT_FREQ(NU,NCF,FQW,K,NU_EVAL_CONT,
	1             V_DOP,DELV_CONT,COMPUTE_ALL_CROSS)
C
C Redefine frequency quadrature weights.
C
	CALL SMPTRP(NU,FQW,NCF)
	DO ML=1,NCF                                           
	  FQW(ML)=FQW(ML)*1.0D+15
	END DO
C
C 
C Need to calculate impact parameters, and angular quadrature weights here
C as these may be required when setting up the initial temperature
C distribution of the atmosphere (i.e. required by JGREY).
C
C
C Compute impact parameter values P
C
	CALL IMPAR(P,R,RP,NC,ND,NP)
C
C Compute the angular quadrature weights
C
	IF(TRAPFORJ)THEN
	  CALL NORDANGQW(AQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,JTRPWGT)
	  CALL NORDANGQW(HQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,HTRPWGT)
	  CALL NORDANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,KTRPWGT)
	  MID=.TRUE.
	  CALL GENANGQW(HMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,HTRPWGT,MID)
	  CALL GENANGQW(NMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,NTRPWGT,MID)
	ELSE
	  CALL NORDANGQW(AQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,JWEIGHT)
	  CALL NORDANGQW(HQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,HWEIGHT)
	  CALL NORDANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,KWEIGHT)
	  MID=.TRUE.
	  CALL GENANGQW(HMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,HWEIGHT,MID)
	  CALL GENANGQW(NMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,NWEIGHT,MID)
	END IF
C
	IF(ACCURATE)THEN
!
! We first verify that the interpolation range is valid.
!
	  IF(END_INTERP_INDX .GT. ND)END_INTERP_INDX=ND
	  IF(DEEP .GT. ND)DEEP=MIN(5,ND)
	  NDEXT=(END_INTERP_INDX-ST_INTERP_INDX)*NPINS+ND
	  IF(NDEXT .GT. NDMAX)THEN
	    WRITE(LUER,*)' Error - NDEXT larger than NDMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NDMAX in CMFGEN'
	    STOP
	  END IF
	  NCEXT=NC
C
C NB: The following expression guarentees that NPEXT has the same relationship
C to NDEXT and NCEXT as does NP to ND and NC.
C
	  NPEXT=NDEXT+NCEXT+(NP-ND-NC)
	  IF(NPEXT .GT. NPMAX)THEN
	    WRITE(LUER,*)' Error - NPEXT larger than NPMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NPMAX in CMFGEN'
	    STOP
	  END IF
	  I=ND-DEEP
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1         ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1         V,T,SIGMA,ND)
	  CALL IMPAR(PEXT,REXT,RP,NCEXT,NDEXT,NPEXT)
C
	  ALLOCATE (AQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (KQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	    WRITE(LUER,*)'Unable to allocate memory for AQWEXT'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
C
C Note that the F2DAEXT vectors (here used as dummy variables) must be at least
C NPEXT long.
C
	  IF(TRAPFORJ)THEN
	    CALL NORDANGQW(AQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,JTRPWGT)
	    CALL NORDANGQW(HQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,HTRPWGT)
	    CALL NORDANGQW(KQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,KTRPWGT)
	    MID=.TRUE.
	    CALL GENANGQW(HMIDQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,HTRPWGT,MID)
	    CALL GENANGQW(NMIDQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,NTRPWGT,MID)
	  ELSE
	    CALL NORDANGQW(AQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,JWEIGHT)
	    CALL NORDANGQW(HQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,HWEIGHT)
	    CALL NORDANGQW(KQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,KWEIGHT)
	    MID=.TRUE.
	    CALL GENANGQW(HMIDQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,HWEIGHT,MID)
	    CALL GENANGQW(NMIDQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,NWEIGHT,MID)
	  END IF
	ELSE
	  NDEXT=ND ; NCEXT=NC; NPEXT=NP
	  TEXT(1:ND)=T(1:ND)
	END IF
C
C Allocate arrays and vectors for computing observed fluxes.
C
	IF(ACCURATE)THEN
	  NP_OBS_MAX=NPEXT+12
	ELSE
	  NP_OBS_MAX=NP+12
	END IF
	ALLOCATE (IPLUS_STORE(NST_CMF,NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (P_OBS(NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (IPLUS(NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (MU_AT_RMAX(NP_OBS_MAX),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (HQW_AT_RMAX(NP_OBS_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for IPLUS_STORE'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
C
C Used when computing the observed fluxes. These will get overwritten
C if we do an accurate comoving frame soluton using CMF_FORM_SOL.
C
	IF(ACCURATE)THEN
	  DO LS=1,NPEXT
	    MU_AT_RMAX(LS)=SQRT( 1.0D0 -(PEXT(LS)/REXT(1))**2 )
	    HQW_AT_RMAX(LS)=HQWEXT(1,LS)
	  END DO
	ELSE
	  DO LS=1,NP
	    MU_AT_RMAX(LS)=SQRT( 1.0D0 -(P(LS)/R(1))**2 )
	    HQW_AT_RMAX(LS)=HQW(1,LS)
	  END DO
	END IF
C
	ALLOCATE (ETA_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CHI_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RJ_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for IPLUS_STORE'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
C 
C
C Set 2-photon data with current atomic models and populations.
C
	DO ID=1,NUM_IONS
	  CALL SET_TWO_PHOT(TRIM(ION_ID(ID)), 
	1       ATM(ID)%XzVLTE,   ATM(ID)%NXzV,
	1       ATM(ID)%XzVLTE_F, ATM(ID)%XzVLEVNAME_F, ATM(ID)%EDGEXzV_F,
	1       ATM(ID)%GXzV_F,   ATM(ID)%F_TO_S_XzV,   ATM(ID)%NXzV_F, ND,
	1       ATM(ID)%ZXzV,     ATM(ID)%EQXzV,        ATM(ID)%XzV_PRES)
	END DO
C
	DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND))
	COHERENT_ES=.TRUE.
C 
!
!
C***************************************************************************
C***************************************************************************
C
C                         CONTINUUM LOOP
C
C***************************************************************************
C***************************************************************************
C
	FIRST_OBS_COMP=.TRUE.
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the BLANKETING portion of the code.
!
	LAST_LINE=0			!Updated as each line is done
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	        LAST_LINE=LAST_LINE+1
	END DO
	DO SIM_INDX=1,MAX_SIM
	  LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	END DO
C
	CONT_FREQ=0.0D0
!
! Define parameters to allow the Cummulative force multipler to be output at
! a function of frequency. We presently output the force multiplier every 500km/s.
!
	N_FORCE=DLOG(NU(NCF)/NU(1))/DLOG(1.0D0-500.0D0/C_KMS)
	NU_FORCE=NU(1)
	NU_FORCE_FAC=(1.0D0-500.0D0/C_KMS)
	ML_FORCE=1
C                                                                    
C Enter loop for each continuum frequency.
C
	CALL TUNE(IONE,'MLCF')
	DO 10000 ML=1,NCF
	  FREQ_INDX=ML
	  FL=NU(ML)
	  IF(ML .EQ. 1)THEN
	    FIRST_FREQ=.TRUE.
	  ELSE
	    FIRST_FREQ=.FALSE.
	  END IF
	  SECTION='CONTINUUM'
C
	  IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	    COMPUTE_NEW_CROSS=.TRUE.
	    CONT_FREQ=NU_EVAL_CONT(ML)
	  ELSE
	    COMPUTE_NEW_CROSS=.FALSE.
	  END IF
C
C 
C
C Section to include lines automatically with the continuum.
C
C
C  LINES_THIS_FREQ --- Logical vector [NCF] indicating whether this frequency
C                        is part of the resonance zone (i.e. Doppler profile) of 
C                        one (or more) lines.
C
C LINE_ST_INDX_IN_NU --- Integer vector [N_LINES] which specifies the starting
C                          frequency index for this lines resonance zone.
C
C LINE_END_INDX_IN_NU --- Integer vector [N_LINES] which specifies the final
C                          frequency index for this lines resonance zone.
C
C FIRST_LINE   ---- Integer specifying the index of the highest frequency
C                         line which we are taking into account in the
C                         transfer.
C
C LAST_LINE  ---- Integer specifying the index of the lowest frequency
C                         line which we are taking into account in the
C                         transfer.
C                                        
C LINE_LOC   ---- Integer array. Used to locate location of a particular line
C                         in the SIM vectors/arrays.
C
C SIM_LINE_POINTER --- Integer array --- locates the line corresponding to
C                         the indicated storage location in the SIM vectors/
C                         arrays.
C
C Check whether we have to treat another line. We use a DO WHILE, rather
C than an IF statement, to handle lines which begin at the same (upper)
C frequency.
C
C
	DO WHILE( LAST_LINE .LT. N_LINE_FREQ .AND.
	1                ML .EQ. LINE_ST_INDX_IN_NU(LAST_LINE+1) )
C
C Have another line --- need to find its storage location.
C
	  I=1
	  DO WHILE(LINE_STORAGE_USED(I))
	    I=I+1
	    IF(I .GT. MAX_SIM)THEN
	      FIRST_LINE=N_LINE_FREQ
	      DO SIM_INDX=1,MAX_SIM		!Not 0 as used!
	        FIRST_LINE=MIN(FIRST_LINE,SIM_LINE_POINTER(SIM_INDX)) 
	      END DO
	      IF( ML .GT. LINE_END_INDX_IN_NU(FIRST_LINE))THEN
C
C Free up storage location for line.
C
	        I=LINE_LOC(FIRST_LINE)
	        LINE_STORAGE_USED(I)=.FALSE.
	        SIM_LINE_POINTER(I)=0
	      ELSE
	        WRITE(LUER,*)'Too many lines have overlapping '//
	1                      'resonance zones'
	        WRITE(LUER,*)'Current frequency is:',FL
	        STOP
	      END IF
	    END IF
	  END DO
C
	  SIM_INDX=I
	  LAST_LINE=LAST_LINE+1
	  LINE_STORAGE_USED(SIM_INDX)=.TRUE.
	  LINE_LOC(LAST_LINE)=SIM_INDX
	  SIM_LINE_POINTER(SIM_INDX)=LAST_LINE
C
C Have located a storage location. Now must compute all relevant quantities
C necessary to include this line in the transfer calculations.
C
	  SIM_NL(SIM_INDX)=VEC_NL(LAST_LINE)
	  SIM_NUP(SIM_INDX)=VEC_NUP(LAST_LINE)
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
C
	  EINA(SIM_INDX)=VEC_EINA(LAST_LINE)
	  OSCIL(SIM_INDX)=VEC_OSCIL(LAST_LINE)
	  FL_SIM(SIM_INDX)=VEC_FREQ(LAST_LINE)
C
	  TRANS_NAME_SIM(SIM_INDX)=TRIM(VEC_SPEC(LAST_LINE))
	  AMASS_SIM(SIM_INDX)=AMASS_ALL(NL)
C
C 
C
C Compute U_STAR_RATIO and L_STAR_RATIO which are used to switch from
C the opacity/emissivity computed with a FULL_ATOM to an equivalent form
C but written in terms of the SUPER-LEVELS. 
C
C L refers to the lower level of the transition.
C U refers to the upper level of the transition.
C
C At present we must treat each species separately (for those with both FULL
C and SUPER_LEVEL model atoms).
C
C MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
C MNL (MNUP) denotes the lower (upper) level in the super level model atom.
C
	MNL_F=VEC_MNL_F(LAST_LINE)
	MNUP_F=VEC_MNUP_F(LAST_LINE)
	DO K=1,ND
	  L_STAR_RATIO(K,SIM_INDX)=1.0D0
	  U_STAR_RATIO(K,SIM_INDX)=1.0D0
	END DO
C
C T1 is used to represent b(level)/b(super level). If no interpolation of
C the b values in a super level has been performed, this ratio will be unity .
C This ratio is NOT treated in the linearization.
C
	DO ID=1,NUM_IONS
	  IF(VEC_SPEC(LAST_LINE) .EQ. ION_ID(ID))THEN
	    MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	    MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	    DO K=1,ND
	      T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzVLTE_F(MNL_F,K)) /
	1             (ATM(ID)%XzV(MNL,K)/ATM(ID)%XzVLTE(MNL,K))
	      L_STAR_RATIO(K,SIM_INDX)=T1*ATM(ID)%W_XzV_F(MNUP_F,K)*
	1        ATM(ID)%XzVLTE_F(MNL_F,K)/ATM(ID)%XzVLTE(MNL,K)/
	1        ATM(ID)%W_XzV_F(MNL_F,K)
	      T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzVLTE_F(MNUP_F,K)) /
	1             (ATM(ID)%XzV(MNUP,K)/ATM(ID)%XzVLTE(MNUP,K))
	      U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F(MNUP_F,K)/
	1              ATM(ID)%XzVLTE(MNUP,K)
	    END DO
	    GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	    TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1      '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1           TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	    EXIT
	  END IF
	END DO
C 
C
C Compute line opacity and emissivity for this line.
C
	  T1=OSCIL(SIM_INDX)*OPLIN
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=1,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    IF(CHIL_MAT(I,SIM_INDX) .EQ. 0)THEN
	      CHIL_MAT(I,SIM_INDX)=0.01*T1*POPS(NL,I)*L_STAR_RATIO(I,SIM_INDX)
	      WRITE(LUER,*)'Zero line opacity in CMFGEN_SUB'
	      WRITE(LUER,*)'This needs to be fixed'
	      J=ICHRLEN(TRANS_NAME_SIM(SIM_INDX))
	      WRITE(LUER,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)(1:J)
	    END IF
	  END DO
C
C Ensure that LAST_LINE points to the next LINE that is going to be handled 
C in the BLANKETING portion of the code.
C
	  DO WHILE(LAST_LINE .LT. N_LINE_FREQ.AND.
	1            VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	       LAST_LINE=LAST_LINE+1
	  END DO
C	   
	END DO	!Checking whether a  new line is being added.
C
C 
C
C Check whether current frequency is a resonance frequency for each line.
C
	DO SIM_INDX=1,MAX_SIM
	  RESONANCE_ZONE(SIM_INDX)=.FALSE.
	  END_RES_ZONE(SIM_INDX)=.FALSE.
	  IF(LINE_STORAGE_USED(SIM_INDX))THEN
	    L=SIM_LINE_POINTER(SIM_INDX)
	    IF( FREQ_INDX .GE. LINE_ST_INDX_IN_NU(L) .AND.
	1          FREQ_INDX .LT. LINE_END_INDX_IN_NU(L))THEN
	      RESONANCE_ZONE(SIM_INDX)=.TRUE.
	    ELSE IF(FREQ_INDX .EQ. LINE_END_INDX_IN_NU(L))THEN
 	      RESONANCE_ZONE(SIM_INDX)=.TRUE.
	      END_RES_ZONE(SIM_INDX)=.TRUE.
	    END IF
	  END IF
	END DO
!
! Compute profile: Doppler or Stark: T2 and T3 are presently garbage.
!
	  TA(1:ND)=0.0D0; TB(1:ND)=0.0D0
	  T2=0.0D0; T3=0.0D0
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
              J=SIM_LINE_POINTER(SIM_INDX); I=ML
	      T1=Z_POP(VEC_NL(J))+1.0D0
              CALL SET_PROF_V3(TA,NU,I,
	1               LINE_ST_INDX_IN_NU(J),LINE_END_INDX_IN_NU(J),
	1               ED,TA,TB,T,VTURB_VEC,ND,
	1               PROF_TYPE(J),PROF_LIST_LOCATION(J),
	1               VEC_FREQ(J),VEC_MNL_F(J),VEC_MNUP_F(J),
	1               AMASS_SIM(SIM_INDX),T1,VEC_ARAD(J),T3,VTURB_FIX,
	1               END_RES_ZONE(SIM_INDX),L_TRUE,LUIN)
              LINE_PROF_SIM(1:ND,SIM_INDX)=TA(1:ND)
	    ELSE              
	      LINE_PROF_SIM(1:ND,SIM_INDX)=0.0D0
	    END IF
	  END DO                                    
C
C 
C
C Determine which method will be used to compute continuum intensity.
C
	  IF(ACCURATE .AND. ALL_FREQ)THEN       
	    THIS_FREQ_EXT=.TRUE.
	  ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE        
	    THIS_FREQ_EXT=.FALSE.
	  END IF
C
C Compute opacity and emissivity.
C
	  CALL TUNE(IONE,'C_OPAC')
	  INCLUDE 'OPACITIES_V4.INC'
	  CALL TUNE(ITWO,'C_OPAC')
C
C Since resonance zones included, we must add the line opacity and 
C emissivity to the raw continuum values. We first save the pure continuum 
C opacity and emissivity. These are used in carrying the variation of J from 
C one frequency to the next.
C
	  DO I=1,ND
	    CHI_CONT(I)=CHI(I)
	    ETA_CONT(I)=ETA(I)
	  END DO
C
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      DO I=1,ND
	        CHI(I)=CHI(I)+CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	        ETA(I)=ETA(I)+ETAL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	      END DO
	    END IF
	  END DO
C
C CHECK for negative line opacities. We do not distinguish between lines.  
C
	  AT_LEAST_ONE_NEG_OPAC=.FALSE.
	  NEG_OPACITY(1:ND)=.FALSE.
	  IF(NEG_OPAC_OPTION .EQ. 'SRCE_CHK')THEN
	    DO I=1,ND
	      IF(CHI(I) .LT. CHI_CONT(I) .AND.
	1            CHI(I) .LT. 0.1D0*ETA(I)*(CHI_CONT(I)-ESEC(I))/ETA_CONT(I) )THEN
	        CHI(I)=0.1D0*ETA(I)*(CHI_CONT(I)-ESEC(I))/ETA_CONT(I)
	        NEG_OPACITY(I)=.TRUE.
	        AT_LEAST_ONE_NEG_OPAC=.TRUE.
	      ELSE IF(CHI(I) .LT. 0.1D0*ESEC(I))THEN
	        CHI(I)=0.1D0*ESEC(I)
	        NEG_OPACITY(I)=.TRUE.
	        AT_LEAST_ONE_NEG_OPAC=.TRUE.
	      END IF
	    END DO
	  ELSE IF(NEG_OPAC_OPTION .EQ. 'ESEC_CHK')THEN
	    DO I=1,ND
	      IF(CHI(I) .LT. 0.1D0*ESEC(I))THEN
	        T1=CHI(I)
	        CHI(I)=0.1D0*ESEC(I)
	        NEG_OPACITY(I)=.TRUE.
	        AT_LEAST_ONE_NEG_OPAC=.TRUE.
	      END IF
	    END DO
	  END IF
	  IF(LST_ITERATION .AND. AT_LEAST_ONE_NEG_OPAC)THEN
	    WRITE(LU_NEG,'(A,1P,E14.6)')
	1      ' Neg opacity for transition for frequency ',FL
	    DO SIM_INDX=1,MAX_SIM
	      IF(RESONANCE_ZONE(SIM_INDX))THEN
	        WRITE(LU_NEG,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)
	      END IF
	    END DO
	    J=0
	    K=0
	    DO I=1,ND
	     IF(NEG_OPACITY(I) .AND. K .EQ. 0)K=I
	     IF(NEG_OPACITY(I))J=I
	    END DO
	    WRITE(LU_NEG,'(A,2X,I3,5X,A,2XI3)')
	1        ' 1st depth',K,'Last depth',J
	  END IF
C
	  DO I=1,ND
	    ZETA(I)=ETA(I)/CHI(I)
	    THETA(I)=ESEC(I)/CHI(I)
	  END DO
C
	  IF(LST_ITERATION .AND. ML .NE. NCF)THEN
	    DO I=1,N_TAU_EDGE
	      IF(NU(ML) .GE. TAU_EDGE(I) .AND. 
	1                     NU(ML+1) .LT. TAU_EDGE(I))THEN
	        WRITE(LUER,'(A,1P,E10.4,A,E10.3)')' Tau(Nu=',NU(ML),
	1          ') at outer boundary is:',CHI_CONT(1)*R(1)
	      END IF
	    END DO
	  END IF
C
C
C Free up LINE storage locations. As we are only computing the line flux,
C and not its variation, we can free up the memory space as soon as we
C exit the resonance zone. This procedure is much simpler than in CMFGEN,
C
	DO SIM_INDX=1,MAX_SIM
	  IF(END_RES_ZONE(SIM_INDX))THEN
	    LINE_LOC( SIM_LINE_POINTER(SIM_INDX) )=0
	    LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	    SIM_LINE_POINTER(SIM_INDX)=0
	  END IF
	END DO
C
C
!
C mean opacities.
C
	IF(ML .EQ. 1)THEN		!Need to move to main loop imit.
	  DO I=1,ND
	    RLUMST(I)=0.0D0
	    J_INT(I)=0.0D0
	    K_INT(I)=0.0D0
	    FLUXMEAN(I)=0.0D0
	    LINE_FLUXMEAN(I)=0.0D0
	    ROSSMEAN(I)=0.0D0
	    INT_dBdT(I)=0.0d0
	  END DO
	END IF
	T1=TWOHCSQ*HDKT*FQW(ML)*(NU(ML)**4)
	DO I=1,ND		              !(4*PI)**2*Dex(+20)/L(sun)
	  J_INT(I)=J_INT(I)+RJ(I)*FQW(ML)*4.1274D-12
	  K_INT(I)=K_INT(I)+K_MOM(I)*FQW(ML)*4.1274D-12
	  LINE_FLUXMEAN(I)=LINE_FLUXMEAN(I)+T2*(CHI(I)-CHI_CONT(I))
	  T2=T1*EMHNUKT(I)/(  ( (1.0D0-EMHNUKT(I))*T(I) )**2  )
	  ROSSMEAN(I)=ROSSMEAN(I)+T2/CHI(I)
	END DO
C
C The current opacities and emissivities are stored for the variation of the
C radiation field at the next frequency.
C                                                  
	DO I=1,ND
	  CHI_PREV(I)=CHI_CONT(I)
	  ETA_PREV(I)=ETA_CONT(I)
	END DO
!
! Store opacities and emissivities for use in observer's frame 
! calculation.
!
	ETA_CMF_ST(1:ND,ML)=ETA(1:ND)
	CHI_CMF_ST(1:ND,ML)=CHI(1:ND)
	RJ_CMF_ST(1:ND,ML)=RJ(1:ND)
!
10000	CONTINUE
	CALL TUNE(ITWO,'MLCF')
!
!
!
! Compute ROSSELAND and FLUX mean opacities. Compute the respective
! optical depth scales; TA for the FLUX mean optical depth scale, 
! and TB for the ROSSELAND mean optical depth scale.
!
! T1=4 * [STEFAN BOLTZMAN CONS] * 1.0E+16 / pi
!
	T1=7.218771D+11
	DO I=1,ND
	  ROSSMEAN(I)=T1*( T(I)**3 )/ROSSMEAN(I)
	END DO
	CALL DERIVCHI(dCHIdR,ROSSMEAN,R,ND,METHOD)
        CALL NORDTAU(TA,ROSSMEAN,R,R,dCHIdR,ND)
	CALL DERIVCHI(dCHIdR,ESEC,R,ND,METHOD)
        CALL NORDTAU(DTAU,ESEC,R,R,dCHIdR,ND)
!
! 
	IF(WRITE_ETA_AND_CHI)THEN
	  ACCESS_F=5
	  I=WORD_SIZE*(ND+1)/UNIT_SIZE; J=82
	  CALL WRITE_DIRECT_INFO(ND,I,'20-Aug-2000','CHI_DATA',J)
	  CALL WRITE_DIRECT_INFO(ND,I,'20-Aug-2000','ETA_DATA',J)
	  OPEN(UNIT=82,FILE='ETA_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	  OPEN(UNIT=83,FILE='CHI_DATA',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I,IOSTAT=IOS)
	  WRITE(82,REC=EDD_CONT_REC)ACCESS_F,NCF,ND
	  WRITE(83,REC=EDD_CONT_REC)ACCESS_F,NCF,ND
	  DO I=1,ND
	    TA(I)=6.65D-15*ED(I)
	  END DO
	  DO ML=1,NCF
	    WRITE(82,REC=ACCESS_F-1+ML)(ETA_CMF_ST(I,ML),I=1,ND),NU(ML)
	    WRITE(83,REC=ACCESS_F-1+ML)(CHI_CMF_ST(I,ML)/TA(I),I=1,ND),NU(ML)
	  END DO
	  CLOSE(UNIT=82)
	  CLOSE(UNIT=83)
	END IF
	CALL TUNE(3,' ')
!
	STOP
	END
