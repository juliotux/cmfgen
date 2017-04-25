!                                                     
! Main Subroutine to self consistently solve the equation of transfer and the
! equations of statistical equilibrium for a spherically extended atmosphere 
! in the presence of outflows.
!                                 
! At present the atmosphere considered to consist of 
!      H, He, plus other species such as C, N, O and Fe.
! Eithe H or He must be present.
!
! Abundances are un-normalized relative fractional abundances (i.e. specified 
! specified with respect to some arbitrary species X. Generally we have used
! He as X. If negative they are interpreted as mass-fractions.
!
! Several options are available for handeling both the lines and continuum.
! Some options have not been recently tested.
!
! Inclusion of new species is generally straightforward. Only the
! calling routine, CMFGEN, needs to modified. The main work is in
! collating the atomic data. Dynamical memory allocation is used.
!
! 
!
	SUBROUTINE LTE_SUB(ND,NC,NP,NT,
	1                     NUM_BNDS,NION,DIAG_INDX,
	1                     NDMAX,NPMAX,NCF_MAX,NLINE_MAX,
	1                     TX_OFFSET,MAX_SIM,NM,NM_KI,NLF)
	USE MOD_CMFGEN
	USE ANG_QW_MOD
	USE CMF_SOB_MOD
	USE CONTROL_VARIABLE_MOD
	USE OPAC_MOD
	USE STEQ_DATA_MOD
	USE MOD_LEV_DIS_BLK
	USE LINE_VEC_MOD
	USE LINE_MOD
	USE RADIATION_MOD
	USE VAR_RAD_MOD
	IMPLICIT NONE
!
! Altered 04-Apr-2014 : Bug fix -- needed to move computation of V(r) before the
!                         computation of VTURB_VEC.
! Incorporated 02-Jan-2013: Chagend to allow depth dependent line profiles.
! Altered 20-Oct-2010 : Commented out statements accessing DIFFW.
!                         No longer call CALL SET_VAR_RAD_MOD_V2
!                         This reduces memory requirements by not declaring 
!                         CMFGEN variation arrays. DIFFW had to be commented out, as
!                         declared in SET_VAR_RAD_MOD_V2.
!
! Altered 16-Feb-2006 : Changed and modified over 2 month period. Section solving for
!                         populations and performing NG acceleration etc removed to
!                         subroutine (SOLVE_FOR_POPS). Routines added to allow time 
!                         variability of the statistical equilibrium equations. 
!                         Currently these are only for a Hubble law flow. Relativistic 
!                         terms added to COMP_OBS (now COMP_OBS_V2).
! Altered 20-Feb-2005 : Changed to use FLUX_MEAN & ROSS_MEAN which are defined in
!                         MOD_CMFGEN. Previusly used FLUXMEAN & ROSSMEAN defined in
!                         RADIATION_MOD.
! 
	INTEGER ND,NC,NP,NT
	INTEGER NUM_BNDS,NION,DIAG_INDX
	INTEGER NDMAX,NPMAX
	INTEGER NCF_MAX,NLINE_MAX
	INTEGER NM,NM_KI,MAX_SIM,NLF
	INTEGER TX_OFFSET
!
	INTEGER NCF
	LOGICAL, PARAMETER :: IMPURITY_CODE=.FALSE.
!
	CHARACTER*12 PRODATE
	PARAMETER (PRODATE='16-Feb-2006')	!Must be changed after alterations
!
! 
!
	REAL*8 SOL(NT,ND)		!Temp. stor. area for ST. EQ.
	INTEGER DST,DEND
!
! Constants for opacity etc. These are set in CMFGEN.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
!
! Internally used variables
!
	REAL*8 S1,REPA
	REAL*8 MAXCH,MAXCH_SUM
	REAL*8 T1,T2,T3,T4,SRAT
	REAL*8 FL,AMASS,FL_OLD
	REAL*8 FG_COUNT
!
	LOGICAL LST_DEPTH_ONLY
!
! REC_SIZE     is the (maximum) record length in bytes.
! UNIT_SIZE    is the number of bytes per unit that is used to specify
!                 the record length (thus RECL=REC_SIZ_LIM/UNIT_SIZE).
! WORD_SIZE    is the number of bytes used to represent the number.
! N_PER_REC    is the # of POPS numbers to be output per record.
!
 	INTEGER REC_SIZE
	INTEGER UNIT_SIZE
	INTEGER WORD_SIZE
	INTEGER N_PER_REC
!
! 
!
! Logical Unit assignments. Those indicated with a # after the ! are open in
!  large sections of the code. Other units generally used temprarily.
!
	INTEGER               LUER        !Output/Error file.
	INTEGER, PARAMETER :: LUIN=7      !General input unit (closed after accesses).
	INTEGER, PARAMETER :: LUMOD=8     !Model Description file.
!
	INTEGER, PARAMETER :: LU_DC=9      	!Departure coefficient Output.
	INTEGER, PARAMETER :: LU_FLUX=10   	!Flux/Luminosity Data (OBSFLUX)
	INTEGER, PARAMETER :: LU_SE=16     	!Statistical equilibrium and Solution Arrays.
	INTEGER, PARAMETER :: LU_NET=17    	!# Line Netrate data.
	INTEGER, PARAMETER :: LU_OPAC=18   	!Rosseland mean opacity etc.
	INTEGER, PARAMETER :: LU_DR=19     	!# Downward rate (Nu. Z. A).
	INTEGER, PARAMETER :: LU_EW=20     	!# EW data.
	INTEGER, PARAMETER :: LU_REC_CHK=21	!# EW data.
!
! For writing scratch file (SCRTEMP). Also used in reading in MODEL data.
!
	INTEGER, PARAMETER :: LUSCR=26
!
	INTEGER, PARAMETER :: LU_HT=27     !#LINEHEAT (i.e Line heating term in R.E. equation)
!
! Used for RVTJ file and POPCARB, POPNIT etc.
!
	INTEGER, PARAMETER :: LU_POP=30
!
	INTEGER, PARAMETER :: LU_IMP=34       !J and CHI for impurity calculation.
	INTEGER, PARAMETER :: LU_EDD=35       !Continuum Eddington factors.
	INTEGER, PARAMETER :: LU_JEW=36       !J for EW computation (JEW)
	INTEGER, PARAMETER :: LU_JCOMP=37     !J_COMP
	INTEGER, PARAMETER :: LU_ES=38        !ES_J_CONV
!
! Following is used output the BA matrix, and its associated
! pointer file.
!
	INTEGER, PARAMETER :: LU_BA=40 
!
! For listing of transitions with TOTAL negative opacity values at some depths.
!
	INTEGER, PARAMETER :: LU_NEG=75
!
! 
!
	INTEGER NNM				!Include cont. var in line var.
	INTEGER NL,NUP
	INTEGER MNL,MNUP
	INTEGER MNL_F,MNUP_F
	INTEGER PHOT_ID
	INTEGER I,J,K,L,ML,LS,LINE_INDX,NEXT_LOC
	INTEGER IREC,MATELIM
!
	CHARACTER*80 TMP_STRING
	CHARACTER*20 TMP_KEY
	INTEGER ID,LOC_ID,ID_SAV,JJ
        INTEGER  L1,L2,U1,U2
	INTEGER IT,MNT,NIV
	INTEGER ISPEC
!
! Main iteration loop variables.
!
	INTEGER MAIN_COUNTER,NITSF,NUM_ITS_RD,NUM_ITS_TO_DO
	LOGICAL LST_ITERATION
!
! Functions called
!
	INTEGER ICHRLEN,ERROR_LU
	REAL*8 DOP_PRO
	REAL*8 S15ADF
	REAL*8 LAMVACAIR
	REAL*8 ATOMIC_MASS_UNIT
	REAL*8 SPEED_OF_LIGHT
	REAL*8 RAD_SUN
	REAL*8 TEFF_SUN
	LOGICAL EQUAL
	EXTERNAL ICHRLEN,ERROR_LU,SPEED_OF_LIGHT,RAD_SUN,TEFF_SUN
!
! Photoionization cross-section routines.
!
! Collisional routines.
!
	EXTERNAL OMEGA_GEN_V3
!
! Wind variablity arrays.
!
	REAL*8 POPS(NT,ND)		!Population for all species.
	REAL*8 MEAN_ATOMIC_WEIGHT	!Mean atomic weight of atoms  (neutrals
!                          		! and ions) in atomic mass units.
	REAL*8 ABUND_SUM
!
!
! Arrays for improving on the initial T structure --- partition functions. 
! Need one for each atomic species.
!
	REAL*8, ALLOCATABLE :: U_PAR_FN(:,:)
	REAL*8, ALLOCATABLE :: PHI_PAR_FN(:,:)
	REAL*8, ALLOCATABLE :: Z_PAR_FN(:)
!
	REAL*8 TGREY(ND)
	REAL*8 T_SAVE(ND)
!
! Variables for scaling the line cooling rates in oder that the radiative
! equilibrium equation is more consistent with the electron heating/cooling 
! equation. The scaling is done when the line frequency is with a fraction
! of SCL_LINE_HT_FAC of the tmean frequency for the super-level under 
! consideration. 0.5 is presently the prefered value.
!
!	REAL*8 AVE_ENERGY(NT)		!Average energy of each super level
	REAL*8 STEQ_T_SCL(ND)
	REAL*8 STEQ_T_NO_SCL(ND)
! 
!
! Dielectronic recombination variables and arrays.
!
	INTEGER NMAXDIE
	PARAMETER (NMAXDIE=500)
!
	REAL*8 EDGEDIE(NMAXDIE)		!Ionization frequency (negative)
	REAL*8 EINADIE(NMAXDIE)		!Einstein A coefficient
	REAL*8 GUPDIE(NMAXDIE)		!Stat. weight of autoionizing level.
!
	INTEGER LEVDIE(NMAXDIE)  	!Indicates MNL of low state
	INTEGER INDXDIE(NMAXDIE)
!
! Used for species identification as INDXDIE is not unique (specied not
! present can have same index as a species thats present.)
!
	CHARACTER*10 SPECDIE(NMAXDIE)
	CHARACTER*35 DIENAME(NMAXDIE)
!
	INTEGER NDIETOT
!
! Arrays and variables used for both Dielectronic recombination, and
! the implicit recombination.
!
	INTEGER EQION,EQSPEC
	REAL*8 GLOW,GION
	REAL*8 NUST(ND)			!LTE autoionizing population.
	REAL*8 DION(ND)			!Ion population
	REAL*8, ALLOCATABLE :: DIERECOM(:,:)   !Chk for all species.
	REAL*8, ALLOCATABLE :: DIECOOL(:,:)    !Dielec. cooling check for all spec.
	REAL*8, ALLOCATABLE :: ADDRECOM(:,:)     
! 
!
! Opacity/emissivity
!
	REAL*8 CHIL(ND)                 !Line opacity (without prof.)
	REAL*8 ETAL(ND)                 !Line emissivity (without prof.)
!
! Quadrature weights.
!
	REAL*8 FQW(NCF_MAX)		!Frequency weights
!
! Transfer equation vectors
	REAL*8 R_OLD(NDMAX)		!Used to store previous R grid in SN models.
!
! Line vectors
	REAL*8 AV(ND)
	REAL*8 VB(NDMAX)		!Used for error calculations
	REAL*8 VC(NDMAX)		!Used for error calculations
	REAL*8 H(ND)
	REAL*8 Q(ND)			!FREQ DEPENDENT.
	REAL*8 QH(ND)			!  "      "
	REAL*8 GAM(ND)			!FREQ INDEPENDENT
	REAL*8 GAMH(ND)			!  "      "
! 
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	LOGICAL EDDINGTON
!
! Variables for EW's and LINE blanketing.
!
	REAL*8 CONT_INT,EW
	INTEGER ACCESS_JEW
	LOGICAL COMPUTE_EW,COMPUTE_JEW,COMPUTE_LAM,MID,FULL_ES
!
! ACESS_F is the current record we are writing in EDDFACTOR.
! EDD_CONT_REC is the record in EDDFACTOR which points to the first
! record containing the continuum values.
!
	INTEGER ACCESS_F
	INTEGER, PARAMETER :: EDD_CONT_REC=3
!
	INTEGER NDEXT,NCEXT,NPEXT
!
	REAL*8 CNM(NDMAX,NDMAX)		!For collisions cross-section in
	REAL*8 DCNM(NDMAX,NDMAX)	!STEQGEN
!
!
	INTEGER, PARAMETER :: N_FLUX_MEAN_BANDS=12
	REAL*8     LAM_FLUX_MEAN_BAND_END(N_FLUX_MEAN_BANDS)
	REAL*8     BAND_FLUX_MEAN(ND,N_FLUX_MEAN_BANDS)
	REAL*8     BAND_FLUX(ND,N_FLUX_MEAN_BANDS)
	DATA LAM_FLUX_MEAN_BAND_END/100.0D0,150.0D0,200.0D0,227.83D0,258.90D0,300.0D0,504.25D0,911.75D0,
	1                         1200.0D0,1500.0D0,2000.0D0,1.0D+08/
!
!
! Continuum frequency variables and arrays.
!
	REAL*8 NU(NCF_MAX)		!Continuum and line frequencies
	REAL*8 NU_EVAL_CONT(NCF_MAX)	!Frequencies to evaluate continuum
	REAL*8 OBS(NCF_MAX)		!Observers spectrum
!
! Vectors and arrays used for the observed flux.
!
	INTEGER N_OBS
	REAL*8 OBS_FREQ(NCF_MAX)		!Since N_OBS < NCF =< NCF_MAX
	REAL*8 OBS_FLUX(NCF_MAX)
	LOGICAL FIRST_OBS_COMP
!
	CHARACTER TIME*20
	CHARACTER FMT*120
	CHARACTER*20 SECTION,FORMAT_DATE*20
	CHARACTER STRING*132
	CHARACTER EW_STRING*132
	CHARACTER TEMP_CHAR*132
	CHARACTER*2 FORMFEED
!
! Global vectors:
!
	REAL*8 AMASS_ALL(NT)
	INTEGER N_LINE_FREQ
!
	INTEGER LINES_THIS_FREQ(NCF_MAX)
!
	REAL*8 NU_DOP
	REAL*8 NU_MAX_OBS
	REAL*8 NU_MIN_OBS
!
	INTEGER FREQ_INDX
	INTEGER X_INDX
	INTEGER FIRST_LINE
	INTEGER LAST_LINE
!
! Variables to limit the computation of the continuum opacities and
! emissivities.
!
	REAL*8 JREC(ND)
	REAL*8 dJRECdT(ND)
	REAL*8 JPHOT(ND)
	REAL*8 JREC_CR(ND)
	REAL*8 JPHOT_CR(ND)
	REAL*8 BPHOT_CR(ND)
!
	REAL*8 CONT_FREQ
	LOGICAL FINAL_CONSTANT_CROSS
!
! Indicates whether APRXzV, FFXzZ etc should be zeroed.
!
	LOGICAL ZERO_REC_COOL_ARRAYS 
!
! 
!
	REAL*8 Z_POP(NT)		!Ionic charge for each species
!
! Variables etc for computation of continuum in comoving frame.
!
	LOGICAL FIRST_FREQ
	LOGICAL RAT_TOO_BIG
	LOGICAL NEW_FREQ
!
! 
!
! X-ray variables.
! We dimension from 0 so that we can access a Null vector for the 1st included
! ioinization stage of each species.
! 
	REAL*8 X_RECOM(ND,0:NION)			!Next X-ray recombination rate
	REAL*8 X_COOL(ND,0:NION)			!Next X-ray cooling
!
	REAL*8 OBS_XRAY_LUM_0P1
	REAL*8 OBS_XRAY_LUM_1KEV
	REAL*8 GFF,XCROSS_V2
	EXTERNAL GFF,XCROSS_V2
!
	REAL*8 SPEC_DEN(ND,NUM_SPECIES)		!Used by ELEC_PREP
	REAL*8 AT_NO_VEC(ND,NUM_SPECIES)
!
	REAL*8 AD_COOL_V(ND)
	REAL*8 AD_COOL_DT(ND)
	REAL*8 ARTIFICIAL_HEAT_TERM(ND)
!
	LOGICAL FIRST
	LOGICAL CHK,SUCCESS
        LOGICAL VAR_SOB_JC
	LOGICAL NEG_OPACITY(ND),FIRST_NEG
	LOGICAL AT_LEAST_ONE_NEG_OPAC
	LOGICAL FILE_OPEN
!
! Inidicates approximate frequencies for which TAU at outer boundary is written
! to OUTGEN on the last iteration.
!
! They are the He2 ege, NIII/CIII egde, HeI, HI, HI(N=2).
!
	INTEGER, PARAMETER :: N_TAU_EDGE=5
	REAL*8 TAU_EDGE(N_TAU_EDGE)
	DATA TAU_EDGE/13.16D0,11.60D0,5.95D0,3.29D0,0.83D0/
!
! 
!
!****************************************************************************
!
! Initialization section
!
	LUER=ERROR_LU()
	ACCESS_F=5
	COMPUTE_LAM=.FALSE.
	COMPUTE_EW=.TRUE.
	FULL_ES=.TRUE.
	SN_MODEL=.FALSE.
	ZERO_REC_COOL_ARRAYS=.TRUE.
	I=12
	FORMFEED=' '//CHAR(I)
	CNT_FIX_BA=0
	MAXCH_SUM=0.0D0
	LST_ITERATION=.TRUE.		!.FALSE.
!
!
! When TRUE, FIXED_T indicated that T is to be heled fixed (at least at some
! depths) in the linearization. This variable is set automatically by the 
! code depending on the magnitude of the corrections.
!
	FIXED_T=.FALSE.
!
! Set NDIETOT to zero in case no dielectronic lines are include in
! the dielectronic section.
!
	NDIETOT=0
!
! Set so that it is defined for the TEST whether to do an accurate flux
! calculation.
!
	MAXCH=100.D0
!
! A value of 1000 is used to indicate that the last change was greater them
! VAL_DO_NG, or a NEW_MODEL.
!
! A value of 2000 indicates a continuing model. In this case LAST_NG 
! take precedence. NB: NEXT_NG is reset to 1000 for a new model in the
! new model section.
!
	NEXT_NG=2000			!Initial value Indicate model just bega
!
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
! 
!
! Open a scratch file to record model parameters. This file will eventually
! be renamed MODEL.
!
	CALL GEN_ASCI_OPEN(LUSCR,'MODEL_SCR','UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODEL in CMFGEN, IOS=',IOS
	  STOP
	END IF
	WRITE(LUSCR,'()')
!
	REPA=1.2
	MATELIM=1
!
! This section does the following :-
!	1) If new model read/determine the new radius scale  and
! populations.
!	2) If old model, populations are read in from scratch file
! previous radius scale etc are used. Ponit1 and point2
! point to the input data record (Note : Single Record)
!
	CALL RD_CONTROL_VARIABLES(LUIN,LUSCR,LUER,NUM_BNDS)
	TWO_PHOTON_METHOD='LTE'
!
! RMDOT is the density at R=10dex10 cm and V=1km/s (atomic mass units)
!
	RMDOT=RMDOT*3.02286E+23
!
! 
!
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
!
! NB: The following expression guarentees that NPEXT has the same relationship
! to NDEXT and NCEXT as does NP to ND and NC.
!
	  NPEXT=NDEXT+NCEXT+(NP-ND-NC)
	  IF(NPEXT .GT. NPMAX)THEN
	    WRITE(LUER,*)' Error - NPEXT larger than NPMAX in CMFGEN'
	    WRITE(LUER,*)' Need to increase NPMAX in CMFGEN'
	    STOP
	  END IF
	ELSE
	  NDEXT=ND; NCEXT=NC; NPEXT=NP
	END IF
!
	CALL SET_RADIATION_MOD(ND,NDMAX,NPMAX)
	CALL SET_LINE_MOD(ND,NT,MAX_SIM,NM)
	CALL SET_CMF_SOB_MOD(ND,NUM_BNDS,NT,NM_KI,NLF,LUER)
!
!       CALL SET_VAR_RAD_MOD_V2(ND,NDEXT,
!	1        NT,NUM_BNDS,NM,MAX_SIM,NM_KI,ACCURATE,L_FALSE)
!
	T1=10.0/(NLF-1)
	DO ML=1,NLF
	  PF(ML)=5.0D0-T1*(ML-1)
	END DO
! 
!
! Read in bound-free gaunt factors for individual n states of hydrogen,
! and hydrogenic cross-sections for individual l states (n =0 to 30,
! l=0 to n-1)
!
	CALL RD_HYD_BF_DATA(LUIN,LUSCR,LUER)
!
! Read in atomic data for 2-photon transitions.
!
	CALL RD_TWO_PHOT(LUIN,INCL_TWO_PHOT)
!
! Read in data for charge exchange reactions. 
!
	CALL RD_CHG_EXCH_V3(LUIN,INCL_CHG_EXCH)
!
! Read in X-ray photoionization cross-sections.
!
	CALL RD_XRAY_FITS(LUIN)
!
	IF(XRAYS .AND. .NOT. FF_XRAYS)THEN
	  CALL RD_XRAY_SPEC(T_SHOCK_1,T_SHOCK_2,LUIN)
	END IF
!
	IF(XRAYS .AND. ADD_XRAYS_SLOWLY .AND. RD_LAMBDA)THEN
	   FILL_X1_SAV=FILL_FAC_XRAYS_1
	   FILL_X2_SAV=FILL_FAC_XRAYS_2
	   FILL_FAC_XRAYS_1=FILL_FAC_X1_BEG
	   FILL_FAC_XRAYS_2=FILL_FAC_X2_BEG
	END IF
!
! 
!
! Read in oscillator strengths, the photoionization cross section data,
! dielectronic data, and implicit recombination data for carbon. 
! Individual species are grouped together (rather than grouping all the
! oscillator reads) so that the headers of the INPUT files are grouped
! in the MODEL output file.
!
! We do this in reverse order so that GIONXzV can be correctly set.
!
! Note Well - in GENOSICL  - T1 is returned with the ionization energy.
!                            T2 is returned with screened nuclear charge.
!                             I is returned with the number of transitions.
!
! RDGENDIE returns the dielectronic transitions as a line list. These lines
! are treated as individual lines.
!
! RD_XzV_PHOT_DIE associates the dielectronic transitions with the
! photoionization cross-sections. They are then handeled as part of
! the continuum cross-sections.
!
! NB: The passed GF_CUT is set to zero if the Atomic NO. of the species
!      under consideration is less than AT_NO_GF_CUT.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	    IF( ATM(ID)%XzV_PRES)THEN
	      IF( NINT(AT_NO(SPECIES_LNK(ID))) .LT. NINT(AT_NO_GF_CUT) )THEN
	        T2=0.0D0
	      ELSE
	        T2=GF_CUT
	      END IF
	      TMP_STRING=TRIM(ION_ID(ID))//'_F_OSCDAT'
	      CALL GENOSC_V9( ATM(ID)%AXzV_F, ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,ATM(ID)%XzVLEVNAME_F,
	1                 ATM(ID)%ARAD,ATM(ID)%GAM2,ATM(ID)%GAM4,ATM(ID)%OBSERVED_LEVEL,
	1                 T1, ATM(ID)%ZXzV,
	1                 ATM(ID)%XzV_OSCDATE, ATM(ID)%NXzV_F,I,
	1                 'SET_ZERO',T2,GF_LEV_CUT,MIN_NUM_TRANS,L_FALSE,L_FALSE,
	1                 LUIN,LUSCR,TMP_STRING)
	      TMP_STRING=TRIM(ION_ID(ID))//'_F_TO_S'
	      CALL RD_F_TO_S_IDS( ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1           ATM(ID)%XzVLEVNAME_F, ATM(ID)%NXzV_F, ATM(ID)%NXzV,
	1           LUIN,TMP_STRING)
	      CALL RDPHOT_GEN_V2( ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1           ATM(ID)%GIONXzV_F,AT_NO(SPECIES_LNK(ID)),
	1           ATM(ID)%ZXzV, ATM(ID)%NXzV_F,
	1           ATM(ID)%XzV_ION_LEV_ID, ATM(ID)%N_XzV_PHOT,  NPHOT_MAX,
	1           ATM(ID+1)%XzV_PRES,     ATM(ID+1)%EDGEXzV_F, ATM(ID+1)%GXzV_F,
	1           ATM(ID+1)%F_TO_S_XzV,   ATM(ID+1)%XzVLEVNAME_F, ATM(ID+1)%NXzV_F,
	1           SIG_GAU_KMS,FRAC_SIG_GAU,CUT_ACCURACY,ABOVE_EDGE,
	1           XRAYS,ID,ION_ID(ID),LUIN,LUSCR)   
              IF(ATM(ID+1)%XzV_PRES) ATM(ID)%GIONXzV_F= ATM(ID+1)%GXzV_F(1)
 	      IF(DIE_AS_LINE .AND. (ATM(ID)%DIE_AUTO_XzV .OR.  ATM(ID)%DIE_WI_XzV) )THEN
	        TMP_STRING='DIE'//TRIM(ION_ID(ID))
	        CALL RDGENDIE_V4( ATM(ID)%XzVLEVNAME_F, ATM(ID)%INDX_XzV,
	1             ATM(ID)%NXzV_F,
	1             EDGEDIE,EINADIE,GUPDIE,
	1             LEVDIE,INDXDIE,SPECDIE,DIENAME, ATM(ID)%GIONXzV_F,
	1             ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LUSCR,L_TRUE,TMP_STRING,NMAXDIE,NDIETOT)
	      ELSE IF( ATM(ID)%DIE_AUTO_XzV .OR.  ATM(ID)%DIE_WI_XzV)THEN
	        TMP_STRING='DIE'//TRIM(ION_ID(ID))
	        CALL RD_PHOT_DIE_V1(ID,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1             ATM(ID)%NXzV_F,    ATM(ID)%GIONXzV_F,
	1             VSM_DIE_KMS, ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LUSCR,TMP_STRING)
	      END IF
	    END IF
	  END DO
	END DO
!
! 
!
! We open a new MODEL file so that the information is at the head of the
! file.
!
	CALL GEN_ASCI_OPEN(LUMOD,'MODEL','UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error opening MODEL in CMFGEN, IOS=',IOS
	  STOP
	END IF
!
! Output description of model. This is done after the reading of most data
! since some of the model information is read in.
!
	  CALL DATE_TIME(TIME)
	  WRITE(LUMOD,'(//,'' Model Started on:'',15X,(A))')TIME
	  WRITE(LUMOD,
	1       '('' Main program last changed on:'',3X,(A))')PRODATE
	  WRITE(LUMOD,'()')
	  FMT='(5X,I8,5X,''!Number of depth points'')'
	  WRITE(LUMOD,FMT)ND
	  FMT='(5X,I8,5X,''!Number of core rays'')'
	  WRITE(LUMOD,FMT)NC
	  FMT='(5X,I8,5X,''!Total number of rays'')'
	  WRITE(LUMOD,FMT)NP
	  FMT='(5X,I8,5X,''!Total number of variables'')'
	  WRITE(LUMOD,FMT)NT
	  FMT='(5X,I8,5X,''!Maximum number of frequencies'')'
	  WRITE(LUMOD,FMT)NCF_MAX
	  FMT='(5X,I8,5X,''!Number of bands'')'
	  WRITE(LUMOD,FMT)NUM_BNDS
!	  
	  WRITE(LUMOD,'()')
	  CALL RITE_ATMHD_V2(LUMOD)
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      ISPEC=SPECIES_LNK(ID)
	      CALL RITE_ATMDES_V2( ATM(ID)%XzV_PRES, ATM(ID)%NXzV, 
	1          ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1          ATM(ID)%NXzV_F, ATM(ID)%GIONXzV_F, ATM(ID)%N_XzV_PHOT,
	1          AT_NO(ISPEC),AT_MASS(ISPEC),LUMOD,ION_ID(ID))
	    END IF
	  END DO
!
! Append VADAT information and atomic data headers to model file.
! Thus only have 1 model file output.
!
	  REWIND(LUSCR)
	  IOS=0
	  TEMP_CHAR=FORMFEED
	  DO WHILE(IOS .EQ. 0)
	    I=ICHRLEN(TEMP_CHAR)
	    IF(I .GT. 0)THEN
	      WRITE(LUMOD,'(A)')TEMP_CHAR(1:I)
	    ELSE
	      WRITE(LUMOD,'()')
	    END IF
	    READ(LUSCR,'(A)',IOSTAT=IOS)TEMP_CHAR
	 END DO
!
! Finished all data read, so can close LUMOD (output descriptor).
!
	CLOSE(UNIT=LUSCR,STATUS='DELETE')
	CLOSE(UNIT=LUMOD)
! 
!
! Set the vector Z_POP to contain the ionic charge for each species.
!
	DO I=1,NT
	  Z_POP(I)=0.0D0
	END DO
!
	DO ID=1,NUM_IONS-1
	  CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1              ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
	END DO
	DO I=1,NT
	  WRITE(6,*)Z_POP(I)
	END DO
!
! Store atomic masses in vector of LENGTH NT for later use by line 
! calculations. G_ALL and LEVEL_ID  are no longer used due to the use
! of super levels.
!
	AMASS_ALL(1:NT)=0.0D0
	DO ID=1,NUM_IONS-1
	  IF(ATM(ID)%XzV_PRES)AMASS_ALL( ATM(ID)%EQXzV: ATM(ID)%EQXzV+ATM(ID)%NXzV-1)=
	1         AT_MASS(SPECIES_LNK(ID))
	END DO
!
! These don't effect the LTE opacities but are needed in some places
! (for consitenct with CMFGEN).
!
	DO I=1,ND
	  R(ND-I+1)=1.0D0+(I-1)*0.01D0
	  V(I)=1.0D0
	  SIGMA(I)=0.0D0
	END DO
!
! 
!
! Compute profile frequencies such that for the adopted doppler
! velocity the profile ranges from 5 to -5 doppler widths.
! This section needs to be rewritten if we want the profile to
! vary with depth.
!
! ERF is used in computing the Sobolev incident intensity at the
! outer boundary. ERF = int from "x" to "inf" of -e(-x^2)/sqrt(pi).
! Note that ERF is not the error function. ERF is related to the
! complementary error function by ERF =-0.5D0 . erfc(X).
! S15ADF is a NAG routine which returns erfc(x).
!
! The incident Sobolev intensity is S[ 1.0-exp(tau(sob)*ERF) ]
! NB -from the definition, -1<erf<0 .
!
	T1=4.286299D-05*SQRT( TDOP/AMASS_DOP + (VTURB/12.85)**2 )
	J=0
	DO I=1,NLF
	  ERF(I)=-0.5D0*S15ADF(PF(I),J)
	  PF(I)=PF(I)*T1
	END DO
        VDOP_VEC(1:ND)=12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
        VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*V(1:ND)/V(1)
!
        IF(GLOBAL_LINE_PROF(1:4) .EQ. 'LIST')THEN
          CALL RD_STRK_LIST(LUIN)
        END IF
!
! Compute the frequency grid for CMFGEN. Routine also allocates the vectors 
! needed for the line data, sets the line data, and puts the line data into 
! numerical order.
!
	CALL SET_FREQUENCY_GRID_V2(NU,FQW,LINES_THIS_FREQ,NU_EVAL_CONT,
	1               NCF,NCF_MAX,N_LINE_FREQ,ND,
	1               OBS_FREQ,OBS,N_OBS,LUIN,IMPURITY_CODE)
!
! Define the average energy of each super level. At present this is
! depth independent, which should be adequate for most models.
! This average energy is used to scale the line cooling rates in
! the radiative equilibrium equation so that is more consistent
! with the electron cooling rate. The need for this scaling
! arises when levels within a super level have a 'relatively large'
! energy separation, and the dominat rates are scattering.
!
	AVE_ENERGY(:)=0.0D0
	DO ID=1,NUM_IONS-1
	   CALL AVE_LEVEL_ENERGY(AVE_ENERGY, ATM(ID)%EDGEXzV_F,
	1         ATM(ID)%GXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%EQXzV, 
	1         ATM(ID)%NXzV,   ATM(ID)%NXzV_F, NT, ATM(ID)%XzV_PRES)
	END DO
!
! 
!
	NLBEGIN=0		! Initialize for lines.
!
! Compute CLUMP_FAC(1:ND) which allow for the possibility that the wind is
! clumped. At the sime time, we compute the vectors which give the density,
! the atom density, and the species density at each depth.
!
	CALL SET_ABUND_CLUMP(MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,ND)
!
! 
!
	IF(ACCURATE)THEN
	  I=ND-DEEP
	  IF(INTERP_TYPE .NE. 'LOG')THEN
	    WRITE(LUER,*)'Error in CMFGEN_SUB'
	    WRITE(LUER,*)'The INTERP_TYPE currently implemented is LOG'
	    STOP
	  END IF
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1         ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	  TA(1:ND)=1.0D0	!TEXT not required, T currently zero
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,TA,SIGMA,ND)
!
          VDOP_VEC_EXT(1:NDEXT)=12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
!
	END IF
!
! Need to calculate impact parameters, and angular quadrature weights here
! as these may be required when setting up the initial temperature
! distribution of the atmosphere (i.e. required by JGREY).
!
	CALL SET_ANG_QW(R,NC,ND,NP,REXT,NCEXT,NDEXT,NPEXT,TRAPFORJ,ACCURATE)
!
! Allocate memory for opacities.
!
        CALL INIT_OPAC_MOD(ND,NT,L_TRUE)
! 
!
!		'NEW MODEL'
!
! Read in old estimates for the departure coefficents for the new model.
! T and ED are also estimated.
!
	NEWMOD=.TRUE.
	IF(NEWMOD)THEN
	  NITSF=0
	  IREC=0
	  LAST_NG=-1000  			!Must be -1000
	  NEXT_NG=1000				!Must be initialized to 1000
	  CALL SET_LTE_EST(POPS,Z_POP,MEAN_ATOMIC_WEIGHT,ABUND_SUM,LUER,LUIN,ND,NT)
	END IF
!
! VEXT and SIGMAEXT have already been computed. We need TEXT for
! convolving J with the electron scattering redistribution function.
!
	IF(ACCURATE)THEN
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1        V,T,SIGMA,ND)
	ELSE
	  TEXT(1:ND)=T(1:ND)
	END IF
!
! 
!
! Set TSTAR which is required by TOTOPA_JILA. Its precise value is
! irrelevant (provided >0) as we always adopt the diffusion
! approximation. 
!
	TSTAR=T(ND)
!
! 
!
! Associate charge exchange reactions with levels in the model atoms.
!
	CALL SET_CHG_LEV_ID_V4(ND,LUMOD)
	CALL VERIFY_CHG_EXCH_V3()
!
! Determine the number of important variables for each species, and
! set the links.
!
	CALL DETERMINE_NSE(NION,XRAYS)
        CALL CREATE_IV_LINKS_V2(NT,NION)
!
! 
!
!**************************************************************************
!**************************************************************************
!
!                    MAIN ITERATION LOOP
!
!**************************************************************************
!**************************************************************************
!
! MAIN_COUNTER is an integer variable which keeps track of the TOTAL number
! of iterations performed. NB - A NG acceleration is counted as a single
! acceleration.
!
! NUM_ITS_TO_DO indicates the number of iterations left to do. For the
! last iteration this will be zero in the "2000" LOOP.
!
! LST_ITERATION is a logical variable which indicate that the current
! iteration is the last one, and hence DEBUGING and INTERPRETATION data
! should be written out. Its equivalent to NUM_ITS_TO_DO=0 in the loop.
!
	FIRST_OBS_COMP=.TRUE.
!
	IF(LST_ITERATION)THEN
	  CALL GEN_ASCI_OPEN(LU_NEG,'NEG_OPAC','UNKNOWN',' ',' ',IZERO,IOS)
	END IF
!
! 
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	 DO J=1,ND
	    POPION(J)=0.0D0
	    DO I=1,NT
	      IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	    END DO
	 END DO
	WRITE(26,*)POPION(ND)
!
! 
!
! This routine not only evaluates the LTE populations of both model atoms, but
! it also evaluates the dln(LTE Super level Pop)/dT.
!
	CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
!
! 
!
! Set 2-photon data with current atomic models and populations.
!
	DO ID=1,NUM_IONS-1
	  ID_SAV=ID
	  CALL SET_TWO_PHOT_V3(ION_ID(ID), ID_SAV,
	1       ATM(ID)%XzVLTE,          ATM(ID)%NXzV,
	1       ATM(ID)%XzVLTE_F_ON_S,   ATM(ID)%XzVLEVNAME_F,
	1       ATM(ID)%EDGEXzV_F,       ATM(ID)%GXzV_F,
	1       ATM(ID)%F_TO_S_XzV,      ATM(ID)%NXzV_F, ND,
	1       ATM(ID)%ZXzV,            ATM(ID)%EQXzV,  ATM(ID)%XzV_PRES)
	END DO
!
! 
!
! Section to compute DT/DR and D(DT/DR)/D? for use in the
! diffusion approximation. DIFFW is used to store the second
! derivatives.
!
! This was based on the subroutine DTSUB. No longer a subroutine as
! too many variables need to be included.
!
! Zero variation of DTDR vector
!
!	DO I=1,NT
!	  DIFFW(I)=0.0D0
!	END DO
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the portion of the code that computes dTdR.
!
	LAST_LINE=0			!Updated as each line is done
	DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	        LAST_LINE=LAST_LINE+1
	END DO
	DO SIM_INDX=1,MAX_SIM
	  LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	END DO
!
	CALL TUNE(1,'FULL_OPAC')
	DTDR=0.0D0
	SECTION='DTDR'
	IF(IONE .EQ. IONE)THEN
!
! We only need to compute the opacity at the innermost depth, but to save
! programing we will compute it at all depths. As this is only done once
! per iteration, not much time will be wasted.
!
! Setting LST_DEPTH_ONLY to true limits the computation of CHI, ETA, and
! dCHI and dETA to the inner boundary only (in some cases).
!
	  LST_DEPTH_ONLY=.FALSE.            !.TRUE.
!
! RJ is used in VARCONT to compute the varaition of ETA. In this section
! we only want the variation of CHI, so we initialize its value to zero.
! This prevents a floating point exception.
!
	  RJ(1:ND)=0.0D0
	  INT_dBdT=0.0D0
	  ROSS_MEAN=0.0D0
!
	  CONT_FREQ=0.0D0
	  DO ML=1,NCF
	    FREQ_INDX=ML
! 
	    FL=NU(ML)
	    IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	      COMPUTE_NEW_CROSS=.TRUE.
	      CONT_FREQ=NU_EVAL_CONT(ML)
	    ELSE
	      COMPUTE_NEW_CROSS=.FALSE.
	    END IF
!
	    CALL TUNE(1,'DTDR_OPAC')
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
!	    INCLUDE 'OPACITIES_V4.INC'
	    CALL TUNE(2,'DTDR_OPAC')
!
! 
!
! Compute variation of opacity/emissivity. Store in VCHI and VETA.
!
!	    CALL TUNE(1,'DTDR_VOPAC')
!	    INCLUDE 'VAROPAC_V4.INC'
!	    CALL COMP_VAR_OPAC(POPS,RJ,FL,CONT_FREQ,FREQ_INDX,
!	1                SECTION,ND,NT,LST_DEPTH_ONLY)
!	    CALL TUNE(2,'DTDR_VOPAC')
! 
!
! Compute contribution to CHI and VCHI by lines.
!
! Section to include lines automatically with the continuum.
! Only computes line opacity at final depth point. This is used in the
! computation of dTdR.
!
! NB: Care must taken to ensure that this section remains consistent
!      with that in continuum calculation section.
!
	    CALL SET_LINE_OPAC(POPS,NU,FREQ_INDX,LAST_LINE,N_LINE_FREQ,
	1          LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)
!
! Add in line opacity.
!
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      DO L=1,ND
	         CHI(L)=CHI(L)+CHIL_MAT(L,SIM_INDX)*LINE_PROF_SIM(L,SIM_INDX)
	      END DO
	    END IF
	  END DO
!
! Now do the line variation. This presently ignores the effect of a 
! temperature variation.
!
!	  DO SIM_INDX=1,MAX_SIM
!	    IF(RESONANCE_ZONE(SIM_INDX))THEN
!	      NL=SIM_NL(SIM_INDX)
!	      NUP=SIM_NUP(SIM_INDX)
!	      VCHI(NL,ND)=VCHI(NL,ND)+LINE_PROF_SIM(SIM_INDX)*
!	1        LINE_OPAC_CON(SIM_INDX)*L_STAR_RATIO(ND,SIM_INDX)
!	      VCHI(NUP,ND)=VCHI(NUP,ND)-LINE_PROF_SIM(SIM_INDX)*
!	1        LINE_OPAC_CON(SIM_INDX)*U_STAR_RATIO(ND,SIM_INDX)*
!	1        GLDGU(SIM_INDX)
!	    END IF
!	  END DO
!
! 
!
! Set TA = to the variation vector at the inner boundary.
!
!	    CALL TUNE(1,'DTDR_VEC')
!	    DO I=1,NT
!	      TA(I)=VCHI(I,ND)
!	    END DO
!	    T1=HDKT*NU(ML)/T(ND)
!
! Increment Parameters
!
!	    T3=FQW(ML)*TWOHCSQ*( NU(ML)**3 )*T1*EMHNUKT(ND)/
!	1         CHI(ND)/T(ND)/(1.0D0-EMHNUKT(ND))**2
!	    DTDR=DTDR+T3
!	    DO I=1,NT-1
!	      DIFFW(I)=DIFFW(I)+T3*TA(I)/CHI(ND)
!	    END DO
!	    DIFFW(NT)=DIFFW(NT)+T3*(TA(NT)/CHI(ND)-(T1*(1.0D0+EMHNUKT(ND))
!	1           /(1.0D0-EMHNUKT(ND))-2.0D0)/T(ND))
!	    CALL TUNE(2,'DTDR_VEC')
!
	    T1=TWOHCSQ*HDKT*FQW(ML)*(NU(ML)**4)
	    DO L=1,ND
	      T2=T1*EMHNUKT(L)/(  ( (1.0D0-EMHNUKT(L))*T(L) )**2  )
	      INT_dBdT(L)=INT_dBdT(L)+T2
	      ROSS_MEAN(L)=ROSS_MEAN(L)+T2/CHI(L)
	    END DO
!
	    IF(MOD(ML,1000) .EQ. 0)WRITE(6,*)'ML=',ML
	  END DO
!
! Output the Rosseland mean data, first determining the number
! of uniques temperatures.
!
	  T1=1.0D-12
	  I=1
	  DO L=2,ND
	    IF(EQUAL(T(1),T(L),T1))THEN
	      I=L-1
	      EXIT
	    END IF
	  END DO
	  OPEN(UNIT=150,FILE='ROSSELAND_LTE_TAB',STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(150,'(A)')'!'
	    WRITE(150,'(A,F8.5)')'! Mean atomic weight (atoms/ions only)',MEAN_ATOMIC_WEIGHT
	    WRITE(150,'(A)')'!'
	    WRITE(150,'(I5,20X,A)')I,'!Number of temperatures'
	    WRITE(150,'(I5,20X,A)')ND/I,'!Number of densities'
	    WRITE(150,'(A)')'!'
	    WRITE(150,'(A,8(5X,A))')'!','       T','  Density','Atom Pop.',
	1    '       Ne','Chi(Ross)','  Chi(es)','Kap(Ross)','  Kap(es)'
	    DO L=1,ND
	      ROSS_MEAN(L)=INT_dBdT(L)/ROSS_MEAN(L)
	      WRITE(150,'(8ES14.4)')T(L),DENSITY(L),POP_ATOM(L),ED(L),
	1                   ROSS_MEAN(L),6.65D-15*ED(L),1.0D-10*ROSS_MEAN(L)/DENSITY(L),
	1                   6.65D-25*ED(L)/DENSITY(L)
	    END DO
	  CLOSE(UNIT=150)
!
! The luminosity of the Sun is 3.826E+33 ergs/sec. For convenience
! DTDR will have the units  (E+04K)/(E+10cm) .
!
!	  T1=LUM*7.2685D+11/R(ND)/R(ND)
!	  DTDR=T1/DTDR
!	  T1=( DTDR**2 )/T1
!	  DO I=1,NT
!	    DIFFW(I)=DIFFW(I)*T1
!	  END DO
	END IF
	CALL TUNE(2,'FULL_OPAC')
!	CALL TUNE(3,'  ')
!
!	LST_DEPTH_ONLY=.FALSE.
!	WRITE(LUER,*)'The value of DTDR is :',DTDR
! 
! We first, however, need to compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	DO J=1,ND
	  POPION(J)=0.0D0
	  DO I=1,NT
	    IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	  END DO
	END DO
!
! Revise vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD.
!
!	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
!	CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!
! 
! Output brief summary of the model. This is to facilate the creation
! of compact model logs.
!
	IF(LST_ITERATION .AND. IONE .EQ. ITWO)THEN
!
	  CALL GEN_ASCI_OPEN(LUMOD,'MOD_SUM','UNKNOWN',' ',' ',IZERO,IOS)
!
	  WRITE(LUMOD,'(/,''Model Started on:'',15X,(A))')TIME
	  CALL DATE_TIME(TIME)
	  WRITE(LUMOD,'(''Model Finalized on:'',13X,(A))')TIME
	  WRITE(LUMOD,
	1       '(''Main program last changed on:'',3X,(A))')PRODATE
	  WRITE(LUMOD,'()')
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'ND',ND)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NC',NC)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NP',NP)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NT',NT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NUM_BNDS',NUM_BNDS)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NCF',NCF)
	  CALL WR_INT_INFO(STRING,NEXT_LOC,'NLINES',N_LINE_FREQ)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  WRITE(LUMOD,'(A)')' '
!
! Output brief summary of atomic models.
!
	  DO ISPEC=1,NUM_SPECIES
	    STRING=' '
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      IF(ATM(ID)%XzV_PRES)THEN
	       CALL WR_SL_INFO(STRING,ATM(ID)%NXzV,ATM(ID)%NXzV_F,
	1                        ATM(ID)%ZXzV,ION_ID(ID),LUMOD)
	      END IF
	    END DO
	    IF(STRING .NE. ' ')WRITE(LUMOD,'(A)')TRIM(STRING)
	  END DO
	  WRITE(LUMOD,'(A)')' '
!
! Output stellar parameters.
!
	  STRING=' '
	  NEXT_LOC=1
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'L*',LUM)
	  T1=RMDOT/3.02286D+23
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Mdot',T1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*  ',RP)
	  T1=RMAX/RP   ; CALL WR_VAL_INFO(STRING,NEXT_LOC,'RMAX/R*',T1)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
! Compute the Radius and Velocity at Tau=10, and at Tau=2/3.
!
	  TCHI(1:ND)=ROSS_MEAN(1:ND)*CLUMP_FAC(1:ND)
	  CALL DERIVCHI(dCHIdR,TCHI,R,ND,METHOD)
          CALL NORDTAU(DTAU,TCHI,R,R,dCHIdR,ND)
	  TA(1:ND)=0.0D0 ; DO I=2,ND ; TA(I) = TA(I-1)+DTAU(I-1) ; END DO
	  TB(1)=MIN(2.0D0/3.0D0,TA(ND))  ; TB(2)=MIN(10.0D0,TA(ND))
	  TB(3)=MIN(20.0D0,TA(ND))
	  CALL MON_INTERP(TC,ITHREE,IONE,TB,ITHREE,R,ND,TA,ND)
	  CALL MON_INTERP(AV,ITHREE,IONE,TB,ITHREE,V,ND,TA,ND)
!
! Output summary of Teff, R, and V at RSTAR, Tau=10, and TAU=2/3.
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TA(ND))
	  T1=1.0D+10*RP/RAD_SUN(); CALL WR_VAL_INFO(STRING,NEXT_LOC,'R*/Rsun',T1)
	  T1=TEFF_SUN()*(LUM/T1**2)**0.25
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'T*  ',T1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',V(ND))
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  IF(TA(ND) .GT. 20.0D0)THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(3))		!20.0D0
	    T1=1.0D+10*TC(3)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(LUM/T1**2)**0.25
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(3))
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  IF(TA(ND) .GT. 10.0D0)THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(2))		!10.0D0
	    T1=1.0D+10*TC(2)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(LUM/T1**2)**0.25
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(2))
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  IF(TA(ND) .GT. 0.67D0)THEN
	    NEXT_LOC=1  ;   STRING=' '
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Tau',TB(1))		!0.67D0
	    T1=1.0D+10*TC(1)/RAD_SUN()
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'R /Rsun',T1)
	    T1=TEFF_SUN()*(LUM/T1**2)**0.25
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Teff',T1)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'V(km/s)',AV(1))
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	  END IF
!
	  STRING=' '
	  NEXT_LOC=1
	  T1=4.9376D+07*(RMDOT/3.02286D+23)*V(1)/LUM
	  T2=8.235D+03*(RMDOT/3.02286D+23)*V(1)*V(1)/LUM
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Eta',T1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Ek/L(%)',T2)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  WRITE(LUMOD,'(A)')' '
!
! Velocity law information.
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Vinf1',VINF1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'Beta1',V_BETA1)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'SCL_HT/RP',SCL_HT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
!
	  NEXT_LOC=1  ;   STRING=' '
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'VCORE',VCORE)
	  CALL WR_VAL_INFO(STRING,NEXT_LOC,'VPHOT',VPHOT)
	  WRITE(LUMOD,'(A)')TRIM(STRING)
	  IF(VELTYPE .EQ. 6)THEN
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Vinf2',VINF2)
	    CALL WR_VAL_INFO(STRING,NEXT_LOC,'Beta2',V_BETA2)
	    WRITE(LUMOD,'(A)')TRIM(STRING)
	    WRITE(LUMOD,'(A)')' '
	  END IF
!
! Output abundance information.
!
	  DO ISPEC=1,NUM_SPECIES
	    CALL WR_ABUND_INFO_V2(SPECIES(ISPEC),AT_MASS(ISPEC),
	1           AT_ABUND(ISPEC),ABUND_SUM,
	1           MEAN_ATOMIC_WEIGHT,SOL_MASS_FRAC(ISPEC),LUMOD)
	  END DO
!
	  WRITE(LUMOD,'(A)')' '
	  IF(DO_CLUMP_MODEL)THEN
	    WRITE(LUMOD,'(A,A)')
	1            'Running clumped model: ',TRIM(CLUMP_LAW)
	    WRITE(LUMOD,'(A,1PE10.3)')
	1            'Filling factor at boundary is: ',CLUMP_FAC(1)
	    STRING=' '
	    NEXT_LOC=1
	    DO I=1,N_CLUMP_PAR
	      TEMP_CHAR='CL_P_'
	      WRITE(TEMP_CHAR(6:6),'(I1)')I
	      CALL WR_VAL_INFO(STRING,NEXT_LOC,TEMP_CHAR,CLUMP_PAR(I))
	      IF(NEXT_LOC .GT. 80)THEN
	        WRITE(LUMOD,'(A)')TRIM(STRING)
	        STRING=' '
	        NEXT_LOC=1
	      END IF
	    END DO
	    IF(STRING .NE. ' ')WRITE(LUMOD,'(A)')TRIM(STRING)
	    WRITE(LUMOD,'(A)')' '
	  END IF
!
	  WRITE(LUMOD,'(A,1PE10.3)')
	1            'Maximum correcion (%) on last iteration: ',MAXCH
	  CLOSE(LUMOD)
!
	END IF
!
! 
9999	CONTINUE
!
! NB - Need GE because of NG acceleration.
!
	IF(LST_ITERATION)THEN
!
!	  CALL TUNE(ITWO,'GIT')
	  CALL TUNE(ITHREE,' ')
!
! This file is a direct access file and contains the models
! output (i.e. T,density,population levels etc). No longer
! assumes that RVTJ exists (28-Sep-1990). Altered (1-Dec-1991) to be
! an asci file for CRAY/VAX compatibility. NB. Obs may be ZERO if program
! was inadvertently started during the last iteration.
!
	  CALL DATE_TIME(TIME)
!
	  CALL GEN_ASCI_OPEN(LU_POP,'RVTJ','UNKNOWN',' ',' ',IZERO,IOS)
	    FORMAT_DATE='15-Jun-2000'
	    WRITE(LU_POP,'(1X,A,T30,A)')'Output format date:',FORMAT_DATE
	    WRITE(LU_POP,'(1X,A,T30,A)')'Completion of Model:',TIME
	    WRITE(LU_POP,'(1X,A,T30,A)')'Program Date:',PRODATE
!
	    WRITE(LU_POP,'(1X,A,T30,I5)')'ND:',ND
	    WRITE(LU_POP,'(1X,A,T30,I5)')'NC:',NC
	    WRITE(LU_POP,'(1X,A,T30,I5)')'NP:',NP
	    WRITE(LU_POP,'(1X,A,T30,I6)')'NCF:',N_OBS
!
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'Mdot(Msun/yr):',
	1                                   RMDOT/3.02286D+23
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'L(Lsun):',LUM
	    WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')'H/He abundance:',AT_ABUND(1)
	    WRITE(LU_POP,'(1X,A,T30,L1)')'Was T fixed?:',RD_FIX_T
	    WRITE(LU_POP,'(1X,A,T30,A)')'Species naming convention:',
	1                                        NAME_CONVENTION
!
	    WRITE(LU_POP,'(A)')' Radius (10^10 cm)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')R
	    WRITE(LU_POP,'(A)')' Velocity (km/s)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')V
	    WRITE(LU_POP,'(A)')' dlnV/dlnr-1'
	    WRITE(LU_POP,'(1X,1P8E16.7)')SIGMA
	    WRITE(LU_POP,'(A)')' Electron density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')ED
	    WRITE(LU_POP,'(A)')' Temperature (10^4K)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')T
!
! These are written to OBSFLUX, and hence do not need to be output to
! RVTJ. They are not accessed by DISPGEN.
!
!	    WRITE(LU_POP,'(A)')' Continuum Frequencies (10^15 Hz)'
!	    WRITE(LU_POP,'(1X,1P6E20.12)')(OBS_FREQ(I),I=1,N_OBS)
!	    WRITE(LU_POP,'(A)')' Continuum Fluxes (Jy kpc^2)'
!	    WRITE(LU_POP,'(1X,1P8E16.7)')(OBS_FLUX(I),I=1,N_OBS)
!
	    WRITE(LU_POP,'(A)')' Rosseland Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(ROSS_MEAN(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Flux Mean Opacity'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(FLUX_MEAN(I),I=1,ND)
!
! Compute the ion population at each depth.                          
! These are required when evaluation the occupation probabilities.
!
	    DO J=1,ND
	      POPION(J)=0.0D0
	      DO I=1,NT
	        IF(Z_POP(I) .GT. 0)POPION(J)=POPION(J)+POPS(I,J)
	      END DO
	    END DO
	    WRITE(LU_POP,'(A)')' Atom Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(POP_ATOM(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Ion Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(POPION(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Mass Density (gm/cm^3)'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(DENSITY(I),I=1,ND)
	    WRITE(LU_POP,'(A)')' Clumping Factor'
	    WRITE(LU_POP,'(1X,1P8E16.7)')(CLUMP_FAC(I),I=1,ND)
!
	    WRITE(LU_POP,'(A)')' Hydrogen Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,1)
	    WRITE(LU_POP,'(A)')' Helium Density'
	    WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,2)
!
	    IF(ATM(1)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')TRIM(ION_ID(1))//' populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Hydrogen levels:',ATM(1)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(1))//
	1                ' oscillator date:',ATM(1)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(1)%XzV_F
	      WRITE(LU_POP,'(A)')' DHYD population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(1)%DXzV_F
	    END IF
!
	    IF(ATM(3)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')' HeI populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Helium I levels:',ATM(3)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(3))//
	1                ' oscillator date:',ATM(3)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(3)%XzV_F
	      WRITE(LU_POP,'(A)')' DHeI population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(3)%DXzV_F
	    END IF
!
	    IF(ATM(4)%XzV_PRES)THEN
	      WRITE(LU_POP,'(A)')' He2 populations'
	      WRITE(LU_POP,'(A,T30,I3)')' Number of Helium II levels:',ATM(4)%NXzV_F
	      WRITE(LU_POP,'(A,T30,A)')TRIM(ION_ID(4))//
	1                ' oscillator date:',ATM(1)%XzV_OSCDATE
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(4)%XzV_F
	      WRITE(LU_POP,'(A)')' DHe2 population'
	      WRITE(LU_POP,'(1X,1P8E16.7)')ATM(4)%DXzV_F
	    END IF
!
	  CLOSE(UNIT=LU_POP)
! 
!
	  FORMAT_DATE='27-JAN-1992'
	  DO ISPEC=1,NUM_SPECIES
	    IF(POP_SPECIES(ND,ISPEC) .NE. 0)THEN
	      TMP_STRING='POP'//TRIM(SPECIES(ISPEC))
	      CALL GEN_ASCI_OPEN(LU_POP,TMP_STRING,'UNKNOWN',' ',' ',IZERO,IOS)
	      WRITE(LU_POP,'(1X,A,T30,A)')'Output format date:',FORMAT_DATE
	      WRITE(LU_POP,'(1X,A,T30,A)')'Completion of Model:',TIME
	      WRITE(LU_POP,'(1X,A,T30,I5)')'ND:',ND
	      WRITE(LU_POP,'(1X,A,T30,1P,E12.5)')
	1              TRIM(SPECIES(ISPEC))//'/He abundance:',AT_ABUND(ISPEC)
     	      WRITE(LU_POP,'(1X,1P8E16.7)')POP_SPECIES(1:ND,ISPEC)
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        CALL RITE_ASC( ATM(ID)%XzV_PRES, ATM(ID)%XzV_F,ATM(ID)%DXzV_F,
	1              ATM(ID)%NXzV_F, ND,
	1              ATM(ID)%XzV_OSCDATE,TRIM(ION_ID(ID)),LU_POP)
	      END DO
	    END IF
	  END DO
! 
!
! Write out departure coefficients to ASCI file. 
! NB - 1 refers to dimension of DHYD (i.e. DHYD(1,nd)
!      1 refers to format for output.
!      1,NHY - For use with HeI.
!
	CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
!      
! GAM_SPECIES refers to the number of electrons arising from each species (eg
! carbon).
!
	DO ISPEC=1,NUM_SPECIES
	  GAM_SPECIES(1:ND,ISPEC)=0.0D0
	  FIRST=.TRUE.
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	    IF( ATM(ID)%XzV_PRES)THEN
	      CALL UPDATE_GAM( GAM_SPECIES(1,ISPEC),
	1          ATM(ID)%XzV_F, ATM(ID)%DXzV_F, ATM(ID)%ZXzV,
	1          ATM(ID)%NXzV_F,ND,
	1          ATM(ID+1)%XzV_PRES,FIRST)
	      TMP_STRING=TRIM(ION_ID(ID))//'OUT'
	      CALL WRITEDC_V2( ATM(ID)%XzV_F, ATM(ID)%XzVLTE_F,
	1          ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,IONE,
	1          R,T,ED,V,CLUMP_FAC,LUM,ND,
	1          TRIM(TMP_STRING),'DC',IONE)
	    END IF
	  END DO
	END DO
!
! We only output GAM when it is not zero.
!
	  CALL RITE_GAM_HEAD(R,ED,T,ND,LUIN,'GAMMAS')
	  DO ISPEC=1,NUM_SPECIES
	    CALL RITE_GAM_V2(POP_SPECIES(1,ISPEC),GAM_SPECIES(1,ISPEC),
	1                       AT_NO(ISPEC),SPECIES(ISPEC),ND,LUIN)
	  END DO
	  CLOSE(LUIN)
!
! If we have a time variability model, we output the model so it can be used by the
! next model in the time sequence. IREC is the ouput record, and is thus simply
! the current TIME_SEQ_NO.
!
	  IF(SN_MODEL .AND. TIME_SEQ_NO .NE. 0)THEN
	    IREC=TIME_SEQ_NO
	    CALL RITE_TIME_MODEL(R,V,SIGMA,POPS,IREC,L_FALSE,L_TRUE,
	1               NT,ND,LUSCR,CHK)
	  END IF
!
	  CLOSE(UNIT=LUER)
	  CLOSE(UNIT=LU_SE)
	  STOP
	END IF
!
	END 
