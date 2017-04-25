!
! Subroutine to compute the emergent flux in the observer's frame using an
! Observer's frame flux calculation. The program first computes the
! emissivities, opacities, and radiation field in the CMF frame. COHERENT
! or INCOHERENT electron scattering can be assumed. This data is then passed
! to OBS_FRAME_SUB for computation of the OBSERVER's flux.
!
! This routine is a heavily stripped down and modified version of CMFGEN.
!
	SUBROUTINE CMF_FLUX_SUB_V5(ND,NC,NP,NDMAX,NPMAX,NT,NLINE_MAX)
	USE MOD_CMF_OBS
	USE MOD_FREQ_OBS
	USE CMF_FLUX_CNTRL_VAR_MOD
	USE MOD_LEV_DIS_BLK
	IMPLICIT NONE
!
! Altered: 21-Sep-2016 : Error corrected -- I was writing out KAPPA from RVTJ rather than the from the CMF_FLUX
!                           calculaton. Note that the ROSSELAND mean is only correct (at depth) if we compute
!                           the spectrum over the full frequency range.
! Altered: 09-Sep-2015 : Changed to C4(line)= ABS(C4[upper]) + ABS(C4[lower) [I was just summing values]
! Altered: 18-May-2015 : Changed GAM2, GAM4 to C4 and C6 (quadratic and Van der Waals interacton constants).
!                           C4 is now utilized (read into VEC_C4). C6 is still not used (09-Jun-2015).
! Altered: 04-Apr-2015 : Changed SET_TWO_PHOT_V2 to _V3.
! Altered: 21-Jan-2013 : Change to vectors passed to SET_PRO_V3. Error probably affects IR HeI lines.
!                           No change top UV/optical spectrum was seen.
!                        Placed large vectors (dimension with NCF_MAX and NLINE_MAX) in module MOD_FREQ_OBS.
!                           These were being placed on the stack, and limiting the available stacks size.
!                           Error in STACK found using -g debug option with -Mchkstk.
!
! Altered:  1-Oct-2012 : Fixed bug: NC_OBS was not being set for the normal P-grid which caused
!                           a problem with the observer's frame calculation.
! Altered:  3-Aug-2012 : Changed to allow the computation of a revised grid for shell models.
! Altered: 17-Dec-2011 : Now call OPACITIES_V5.INC and EVAL_LTE_INC_V5.INC 
!                          L_STAR_RATIO and U_STAR_RATIO now computed using XzVLTE_F_ON_S 
!                          Done to allow lower wind temperatures.
! Altered: 20-Oct-2009 : Changed to call OBS_FRAME_SUB_V8 (WRITE_RTAU, TAU_REF added to call).
! Altered: 07-Jul-2008 : Substantial changes over a period of time to include
!                          relativistic radiative transfer. Lots of changes to
!                          COMP_JCONT.INC. NCEXT is now increased according to 
!                          NPINS when ACCURATE option is used.
! Altered: 31-Jul-2007 : Call to MOM_JREL_V3 included in INCLUDE file.
!                          Variabled and control variables for rel. trans. installed.
! Altered: 31-Jul-2007 : Call to MOM_JREL_V3 included in INCLUDE file.
!                          Variabled and control variables for rel. trans. installed.
! Altered: 01-Feb-2006 : COMP_OBS_V2 installed. Frequency transformation between
!                          comoving and observing franes is now accuratle for all speeds.
!                          Frequency scaling of intensities currently switched off.
! Altered: 03-May-2004 : Options are now read into store in CMF_FLUX.
!                          STORE cleaned after use.
! Altered: 05-Jun-2002 : Min bug fix: Changed value of VDOP passed to
!                          DET_MAIN_CONT_FREQ. Reduces number of continuum frequencies.
! Altered: 16-Aug-2002 : Bug fix, NU_OBS_MAX was not always chosen correctly.
! Altered: 03-Jan-2001 : PROF_TYPE and VEC_CHAR_WRK now of length 12.
! Altered: 10-Mar-2000 : Variable type ATM installed to simplify access to
!                          different model atoms. Installation reduces size
!                          of executable, and allows for faster compilation.
!                          Changed to V4.
! Altered: 16-Jan-2000 : Changed to V2; Call to OBS_FRAME_V2 changed.
! Altered: 04-Jan-2000 : We now check that arrays can be allocated.
! Finalized: 5-Jan-1999
!
! Maximum number of lines whose profile overlap at a given frequency.
!
	INTEGER, PARAMETER :: MAX_SIM=3000
!
! Define vinteger variables to set aside for storage for the intrinsic line profiles.
! These are now computed by a subroutine (as of Nov-2013).
!
	INTEGER NLINES_PROF_STORE
	INTEGER NFREQ_PROF_STORE
!
	INTEGER NCF
	INTEGER ND,NC,NP
	INTEGER NDMAX,NPMAX
	INTEGER NT,NLINE_MAX
!
	REAL*8 POPS(NT,ND)
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
!
! Internally used variables
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
	REAL*8 DTDR,BNUE,DBB,DDBBDT
	REAL*8 MEAN_ATOMIC_WEIGHT
	REAL*8 TSTAR,S1,IC,MAXCH
	REAL*8 C_KMS
	REAL*8 T1,T2,T3,T4
	REAL*8 FL,FL_OLD
	REAL*8 FG_COUNT
!
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
! Logical Unit assignments. Those indicated with a # after the ! are open in
!  large sections of the code. Other units generally used temprarily.
!
	INTEGER               LUER        !Output/Error file.
	INTEGER, PARAMETER :: LUIN=7      !General input unit (closed after accesses).
	INTEGER, PARAMETER :: LUMOD=8     !Model Description file.
!
	INTEGER, PARAMETER :: LU_FLUX=10   	!Flux/Luminosity Data (OBSFLUX)
	INTEGER, PARAMETER :: LU_OPAC=18   	!Rosseland mean opacity etc.
	INTEGER, PARAMETER :: LU_EW=20     	!# EW data.
!
	INTEGER, PARAMETER :: LU_EDD=35       !Continuum Eddington factors.
	INTEGER, PARAMETER :: LU_JCOMP=37     !J_COMP
	INTEGER, PARAMETER :: LU_ES=38        !ES_J_CONV
!
! For listing of transitions with TOTAL negative opacity values at some depths.
!
	INTEGER, PARAMETER :: LU_NEG=75
!
! 
!
	INTEGER NL,NUP
	INTEGER MNL,MNUP
	INTEGER MNL_F,MNUP_F
	INTEGER I,J,K,L,ML,LS,LINE_INDX
	INTEGER ID,ID_SAV,ISPEC
	INTEGER ES_COUNTER
!
! Functions called
!
	INTEGER ICHRLEN,ERROR_LU
	REAL*8 LAMVACAIR
	REAL*8 ATOMIC_MASS_UNIT
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL JWEIGHT,HWEIGHT,KWEIGHT,NWEIGHT
	EXTERNAL JTRPWGT,HTRPWGT,KTRPWGT,NTRPWGT
	EXTERNAL ICHRLEN,ERROR_LU,SPEED_OF_LIGHT
	EXTERNAL LAMVACAIR,ATOMIC_MASS_UNIT
!
! 
	CHARACTER FMT*120
	CHARACTER SECTION*20
	CHARACTER TMP_KEY*20
	CHARACTER STRING*132
	CHARACTER EW_STRING*132
	CHARACTER TEMP_CHAR*132
!
! Global vectors:
!
	REAL*8 AMASS_ALL(NT)
!
! Arrays and variables for treating lines simultaneously.
!
	REAL*8 EINA(MAX_SIM)
	REAL*8 OSCIL(MAX_SIM)
	REAL*8 GLDGU(MAX_SIM)
	REAL*8 AMASS_SIM(MAX_SIM)
	REAL*8 FL_SIM(MAX_SIM)
	INTEGER SIM_NL(MAX_SIM)
	INTEGER SIM_NUP(MAX_SIM)
!
	REAL*8 CHIL_MAT(ND,MAX_SIM)
	REAL*8 ETAL_MAT(ND,MAX_SIM)
	REAL*8 BB_COR(ND,MAX_SIM)
!
	INTEGER NUM_SIM_LINES
	INTEGER SIM_INDX
	INTEGER TMP_MAX_SIM
	REAL*8 OVER_FREQ_DIF
	REAL*8 EW
	REAL*8 CONT_INT
	LOGICAL OVERLAP
	LOGICAL SOBOLEV
!
! L refers to the lower level, U to the upper level.
!
	REAL*8 L_STAR_RATIO(ND,MAX_SIM)
	REAL*8 U_STAR_RATIO(ND,MAX_SIM)
!
	REAL*8, ALLOCATABLE :: ETA_CMF_ST(:,:)
	REAL*8, ALLOCATABLE :: CHI_CMF_ST(:,:)
	REAL*8, ALLOCATABLE :: RJ_CMF_ST(:,:)
	REAL*8, ALLOCATABLE :: ION_LINE_FORCE(:,:)
!
! Vectors for treating lines simultaneously with the continuum.
!
	REAL*8 LINE_PROF_SIM(ND,MAX_SIM)
!
	REAL*8 NU_MAX_OBS
	REAL*8 NU_MIN_OBS
!
	REAL*8 CONT_FREQ
!
	CHARACTER*50 TRANS_NAME_SIM(MAX_SIM)
!
	LOGICAL RESONANCE_ZONE(MAX_SIM)
	LOGICAL END_RES_ZONE(MAX_SIM)
	LOGICAL LINE_STORAGE_USED(MAX_SIM)
!
	INTEGER FREQ_INDX
	INTEGER FIRST_LINE
	INTEGER LAST_LINE
	INTEGER SIM_LINE_POINTER(MAX_SIM)
!
! 
!
! Opacity/emissivity
	REAL*8 CHI(ND)			!Continuum opacity (all sources)
	REAL*8 CHI_RAY(ND)
	REAL*8 CHI_SCAT(ND)
	REAL*8 ETA(ND)			!Continuum emissivity (all sources)
	REAL*8 CHIL(ND)			!Line opacity (without prof.)
	REAL*8 ETAL(ND)			!Line emissivity (without prof.)
	REAL*8 ESEC(ND)			!Continuum electron scattering coef.
	REAL*8 ZETA(ND)			!Source func. (all except elec. scat.)
	REAL*8 THETA(ND)		!Elec. scat. source coef.
	REAL*8 SOURCE(ND)		!Complete source function.
	REAL*8 DTAU(NDMAX)		!Optical depth (used in error calcs)
!		DTAU(I)=0.5*(CHI(I)+CHI(I+1))*(Z(I)-Z(I+1))
	REAL*8 dCHIdR(NDMAX) 		!Derivative of opacity.
!
	REAL*8 P(NP)
!
	REAL*8 CHI_CONT(ND)
	REAL*8 ETA_CONT(ND)
!
! These parameters are used when computing J and the variation of J.
!
	REAL*8 CHI_SCAT_CLUMP(ND)
	REAL*8 CHI_RAY_CLUMP(ND)
	REAL*8 CHI_CLUMP(ND)		!==CHI(I)*CLUMP_FAC(I)
	REAL*8 ETA_CLUMP(ND)		!==ETA(I)*CLUMP_FAC(I)
	REAL*8 ESEC_CLUMP(ND)		!==ESEC(I)*CLUMP_FAC(I)
!
! Variables to limit the computation of the continuum opacities and 
! emissivities. 
!
	REAL*8 EMHNUKT_CONT(ND)
	REAL*8 ETA_C_EVAL(ND)
	REAL*8 CHI_C_EVAL(ND)
!
! 
!
	REAL*8 Z_POP(NT)		!Ionic charge for each species
!
! Variables etc for computation of continuum in comoving frame.
!
	LOGICAL FIRST_FREQ
	LOGICAL NEW_FREQ
	REAL*8 dLOG_NU			!Step in frequency in Log plane
	REAL*8 FEDD_PREV(NDMAX)
	REAL*8 GEDD_PREV(NDMAX)
	REAL*8 H_ON_J(NDMAX)
	REAL*8 N_ON_J_NODE(NDMAX)
	REAL*8 dlnJdlnR(NDMAX)
	REAL*8 KMID_ON_J(NDMAX)
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
	REAL*8 HFLUX_AT_IB,HFLUX_AT_OB
!
! Quadrature weights.
	REAL*8 AQW(ND,NP)		!Angular quad. weights. (indep. of v)
	REAL*8 HQW(ND,NP)		!Angular quad. weights. (indep. of v) for flux integration.
	REAL*8 KQW(ND,NP)		!Angular quad. weights for K integration.
	REAL*8 NQW(ND,NP)		!Angular quad. weights. (indep. of v) for N integration.
	REAL*8 HMIDQW(ND,NP)		!Angular quad. weights. (indep. of v)
                                        !for flux integration. Defined at the
                                        !mid points of the radius mesh.
	REAL*8 NMIDQW(ND,NP)		!Angular quad. weights. (indep. of v)
                                        !for N integration. Defined at the
                                        !mid points of the radius mesh.
!
! Continuum matrices
	REAL*8 WM(ND,ND)		!Coef. matrix of J & %J vector
	REAL*8 FB(ND,ND)		!Coef. of J & %J vects in angular equ.
!
! Transfer equation vectors
	REAL*8 TA(NDMAX)
	REAL*8 TB(NDMAX)
	REAL*8 TC(NDMAX)
	REAL*8 XM(NDMAX)		!R.H.S. (SOURCE VECTOR)
!
! 
!
! Arrays and variables for computation of the continuum intensity
! using Eddington factors. This is separate to the "inclusion of
! additional points".
!
	REAL*8 FEDD(NDMAX)
	REAL*8 QEDD(NDMAX)
!
	LOGICAL MID
	LOGICAL FIRST
	LOGICAL NEG_OPACITY(ND)
	LOGICAL FIRST_NEG
	LOGICAL LAMBDA_ITERATION
	LOGICAL LST_ITERATION
	LOGICAL LST_DEPTH_ONLY
! 
!
!
! X-ray variables.
!
	REAL*8 FILL_VEC_SQ(ND)
	REAL*8 XRAY_LUM(ND)
	REAL*8 GFF,XCROSS_V2
	EXTERNAL GFF,XCROSS_V2
!
!
!
! ACESS_F is the current record we are writing in EDDFACTOR.
! EDD_CONT_REC is the record in EDDFACTOR which points to the first
! record containing the continuum values.
!
	INTEGER ACCESS_F
	INTEGER, PARAMETER :: EDD_CONT_REC=3
	CHARACTER(LEN=20) DA_FILE_DATE
!
	INTEGER NDEXT,NCEXT,NPEXT
	INTEGER INDX(NDMAX),POS_IN_NEW_GRID(ND)
	REAL*8 COEF(0:3,NDMAX)
	REAL*8 INBC,HBC_J,HBC_S			!Bound. Cond. for JFEAU
!
	REAL*8 REXT(NDMAX),PEXT(NPMAX),VEXT(NDMAX)
	REAL*8 TEXT(NDMAX),SIGMAEXT(NDMAX)
	REAL*8 CHIEXT(NDMAX),ESECEXT(NDMAX),ETAEXT(NDMAX)
	REAL*8 CHI_RAY_EXT(NDMAX),CHI_SCAT_EXT(NDMAX)
	REAL*8 ZETAEXT(NDMAX),THETAEXT(NDMAX)
	REAL*8 RJEXT(NDMAX),RJEXT_ES(NDMAX)
	REAL*8 FOLD(NDMAX),FEXT(NDMAX),QEXT(NDMAX),SOURCEEXT(NDMAX)
	REAL*8 VDOP_VEC_EXT(NDMAX)
!
!
	REAL*8 F2DAEXT(NDMAX,NDMAX)     !These arrays don't need to be
!
! If required, these arrays shoukd have size NDEXT*NPEXT
!
	REAL*8, ALLOCATABLE :: AQWEXT(:,:)	!Angular quad. weights. (indep. of v)
	REAL*8, ALLOCATABLE :: HQWEXT(:,:)	!Angular quad. weights for flux integration.
	REAL*8, ALLOCATABLE :: KQWEXT(:,:)	!Angular quad. weights for K integration.
	REAL*8, ALLOCATABLE :: NQWEXT(:,:)	!Angular quad. weights for N integration.
	REAL*8, ALLOCATABLE :: HMIDQWEXT(:,:)	!Angular quad. weights for flux integration.
	REAL*8, ALLOCATABLE :: NMIDQWEXT(:,:)	!Angular quad. weights for flux integration.
!
! Arrays for calculating mean opacities.
!
	REAL*8 FLUXMEAN(ND) 		!Flux mean opacity
	REAL*8 LINE_FLUXMEAN(ND) 	!Flux mean opacity due to lines.
	REAL*8 ROSSMEAN(ND)  		!Rosseland mean opacity
	REAL*8 INT_dBdT(ND)  		!Integral of dB/dT over nu 
!                                            (to calculate ROSSMEAN)
	REAL*8 FORCE_MULT(ND)
	REAL*8 NU_FORCE
	REAL*8 NU_FORCE_FAC
	INTEGER N_FORCE
	INTEGER ML_FORCE
	LOGICAL TMP_LOG
!
! Other arrays
	REAL*8 Z(NDMAX)			!Z displacement along a given array
	REAL*8 EMHNUKT(ND)		!EXP(-hv/kT)
	REAL*8 RLUMST(ND)		!Luminosity as a function of depth
	REAL*8 J_INT(ND)		!Frequency integrated J
	REAL*8 K_INT(ND)		!Frequency integrated K
	REAL*8 K_MOM(ND)		!Frequency dependent K moment
	REAL*8 SOB(ND)   	    	!Used in computing continuum flux
	REAL*8 RJ(ND)			!Mean intensity
	REAL*8 RJ_ES(ND)		!Convolution of RJ with e.s. R(v'v')
!
! Line variables.
!
	REAL*8 VAL_DO_NG
	REAL*8 RP
	REAL*8 VINF
	REAL*8 VTURB_VEC(ND)
	REAL*8 VDOP_VEC(ND)
	REAL*8 MAX_DEL_V_RES_ZONE(ND)
!
! Parameters, vectors, and arrays for computing the observed flux.
!
	INTEGER, PARAMETER :: NST_CMF=20000
	INTEGER NP_OBS_MAX
	INTEGER NP_OBS
	INTEGER NC_OBS
	REAL*8  NU_STORE(NST_CMF)
	REAL*8 V_AT_RMAX		!Used if we extend the atmosphere.
	REAL*8 RMAX_OBS
	REAL*8 H_OUT,H_IN
!
! We allocate memory for the following vectors as we use them for the regular
! flux computation, and when extra depth points are inserted (ACCURATE=.TRUE.)
!
	REAL*8, ALLOCATABLE :: IPLUS_STORE(:,:)
	REAL*8, ALLOCATABLE :: P_OBS(:)
	REAL*8, ALLOCATABLE :: IPLUS(:)
	REAL*8, ALLOCATABLE :: MU_AT_RMAX(:)
	REAL*8, ALLOCATABLE :: HQW_AT_RMAX(:)
!
! Supercedes OBS
!
	INTEGER N_OBS
	LOGICAL FIRST_OBS_COMP
!
! Indicates approximate frequencies for which TAU at outer boundary is written
! to OUTGEN on the last iteration.
!
! They are the He2 ege, NIII/CIII egde, HeI, HI, HI(N=2).
!
	INTEGER, PARAMETER :: N_TAU_EDGE=5
	REAL*8 TAU_EDGE(N_TAU_EDGE)
	DATA TAU_EDGE/13.16D0,11.60D0,5.95D0,3.29D0,0.83D0/
!
!
! Check whether EQUATION LABELLING is consistent. ' I ' is used as the
! number of the current equation. We also set the variable SPEC_PRES which 
! indicates whether at least one ioization stage of a species is present.
! It is used to determine, foe example,  whether a number conservation 
! equation is required.
!
!
!???????????????????????
!
	I=1
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
	LUER=ERROR_LU()
	IF(EQNE .NE. I)THEN
	  WRITE(LUER,*)'Error - EQNE has wrong value in CMF_FLUX_SUB_V5'
	  STOP
	END IF
	IF(NT .NE. I+1)THEN
	  WRITE(LUER,*)'Error - NT has wrong value in CMF_FLUX_SUB_V5'
	  STOP
	END IF
!
	CALL INIT_MOD_FREQ_OBS(NLINE_MAX)
! 
!
!	LAMBDA_ITERATION=.FALSE.
	LAMBDA_ITERATION=.TRUE.
	LST_ITERATION=.FALSE.
	LST_DEPTH_ONLY=.FALSE.
	RD_STARK_FILE=.FALSE.
	ACCESS_F=5
	RP=R(ND)
	VINF=V(1)
	TSTAR=T(ND)
	C_KMS=SPEED_OF_LIGHT()/1.0D+05
	DO_REL_IN_OBSFRAME=.FALSE.
	DA_FILE_DATE='20-Aug-2003' 
!
	CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
! MAXCH and VAL_DO_NG are set so that they defined for the TEST in
! COMP_JCONT_V?.INC whether to do an accurate flux calculation. An 
! accurate flux calculation can be avoided by doing a LAMBDA iteration.
!
	MAXCH=0.0D0
	VAL_DO_NG=5.0D0
!
! Set the vector Z_POP to contain the ionic charge for each species.
!
	Z_POP(1:NT)=0.0D0
!
	DO ID=1,NUM_IONS
	  CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,
	1                       ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
	END DO
!
! Store atomic masses in vector of LENGTH NT for later use by line 
! calculations. G_ALL and LEVEL_ID  are no longer used due to the use
! of super levels.
!
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
!
!
! Read in all options controlling the spectrum computation.
! As subroutine exits, it also cleans and releases the STORE.
!
	  CALL RD_CMF_FLUX_CONTROLS(ND,LUMOD,LUER)
!
! 
!
! Scale abunances of species if desired. This should only be done for exploratory
! spectral calculations, and only for IMPURITY species (i.e. not H or He).
!
	  K=0
	  DO ISPEC=1,NUM_SPECIES
	    IF(SPECIES_PRES(ISPEC) .AND. ABUND_SCALE_FAC(ISPEC) .NE. 1.0D0)THEN
	      IF(K .EQ. 0)THEN
	        WRITE(LUER,'(A)')' '
	        WRITE(LUER,'(A)')'******************Warning**********************'
	      END IF
	      K=1
	      WRITE(LUER,'(A,A,A)')'Abundance of species ',TRIM(SPECIES(ISPEC)),' scaled'
	      AT_ABUND(ISPEC)=AT_ABUND(ISPEC)*ABUND_SCALE_FAC(ISPEC)
	      POP_SPECIES(:,ISPEC)=POP_SPECIES(:,ISPEC)*ABUND_SCALE_FAC(ISPEC)
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        ATM(ID)%XzV_F=ATM(ID)%XzV_F*ABUND_SCALE_FAC(ISPEC)
	        ATM(ID)%DXzV_F=ATM(ID)%DXzV_F*ABUND_SCALE_FAC(ISPEC)
	      END DO
	    END IF
	  END DO
	  IF(K .NE. 0)THEN
	    WRITE(LUER,'(A)')'******************Warning**********************'
	    WRITE(LUER,'(A)')' '
	  END IF
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
! Read in atomic data for 2-photon transitions.
!
	CALL RD_TWO_PHOT(LUIN,INCL_TWO_PHOT)
!
! Read in X-ray photoionization cross-sections.
!
	CALL RD_XRAY_FITS(LUIN)
!
	IF(XRAYS .AND. .NOT. FF_XRAYS)THEN
	  CALL RD_XRAY_SPEC(T_SHOCK_1,T_SHOCK_2,LUIN)
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
	  DO ID=NUM_IONS-1,1,-1
	    CALL FULL_TO_SUP( 
	1        ATM(ID)%XzV,        ATM(ID)%NXzV,   ATM(ID)%DXzV,
	1        ATM(ID)%XzV_PRES,   ATM(ID)%XzV_F,
	1        ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,
	1        ATM(ID+1)%XzV,      ATM(ID+1)%NXzV, ATM(ID+1)%XzV_PRES, ND)
	  END DO
!
! Store all quantities in POPS array. 
!
	  DO ID=1,NUM_IONS
	    CALL IONTOPOP(POPS, ATM(ID)%XzV,        ATM(ID)%DXzV,  ED,T,
	1        ATM(ID)%EQXzV, ATM(ID)%NXzV,NT,ND, ATM(ID)%XzV_PRES)
	  END DO
!
! This routine not only evaluates the LTE populations of both model atoms, but
! it also evaluates the dln(LTE Super level Pop)/dT.
!
	INCLUDE 'EVAL_LTE_INC_V5.INC'
!
! Compute the turbulent velocity and MINIMUM Doppler velocity as a function of depth.
! For the later, the iron mas of 55.8amu is assumed.
!
	IF(TURB_LAW .EQ. 'LAW_V1')THEN
	  VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*V(1:ND)/V(1)
	ELSE IF(TURB_LAW .EQ. 'LAW_TAU1')THEN
	  ESEC(1:ND)=6.65D-15*ED(1:ND)*CLUMP_FAC(1:ND)
	  TA(1)=1.0D0
	  DO I=1,ND
	    TA(I)=TA(1)+0.5D0*(R(I)-R(I+1))*(ESEC(I)+ESEC(I+1))
	  END DO
	  VTURB_VEC(1:ND)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)/(1.0D0+TA(1:ND))
	ELSE
	  WRITE(LUER,'(A)')' Error --- TURBULENT velocity law not recognized'
	  WRITE(LUER,'(A)')TRIM(TURB_LAW)
	  STOP
	END IF
	VDOP_VEC(1:ND)=SQRT( VTURB_VEC(1:ND)**2 + 2.96*T(1:ND) )
!
	TA(1:ND)=ABS( CLUMP_FAC(1:ND)-1.0D0 )
	T1=MAXVAL(TA)
	DO_CLUMP_MODEL=.FALSE.
	IF(T1 .GT. 1.0D-05)DO_CLUMP_MODEL=.TRUE.
!
!
! Compute profile frequencies such that for the adopted doppler
! velocity the profile ranges from 5 to -5 doppler widths.
! This section needs to be rewritten if we want the profile to
! vary with depth.
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
!
! Now insert addition points into frequency array. WSCI is used as a
! work array - okay since of length NCF_MAX, and zeroed in QUADSE.
! OBSF contains the bound-free edges - its contents are zero on
! subroutine exit. J is used as temporary variable for the number of
! frequencies transmitted to SET_CONT_FREQ. NCF is returned as the number 
! of frequency points. FQW is used a an integer array for the sorting ---
! we know it has the correct length since it is the same size as NU.
! LUIN --- Used as temporary LU (opened and closed).
!
	  J=NCF
!	  CALL SET_CONT_FREQ(NU,OBS,FQW,
!	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
!	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
!	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
!	1                        J,NCF,NCF_MAX,LUIN)
	T1=5000.0D0
	  CALL SET_CONT_FREQ_V4(NU,OBS,FQW,
	1                        SMALL_FREQ_RAT,BIG_FREQ_AMP,dFREQ_BF_MAX,
	1                        MAX_CONT_FREQ,MIN_CONT_FREQ,
	1                        dV_LEV_DIS,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_CONT,T1,
	1                        J,NCF,NCF_MAX,LUIN)
!                                             
!                         
! 
! Set up lines that will be treated with the continuum calculation.
! This section of code is also used by the code treating purely lines
! (either single transition Sobolev or CMF, or overlapping Sobolev).
!
! To define the line transitions we need to operate on the FULL atom models.
! We thus perform separate loops for each species. VEV_TRANS_NAME is
! allocated temporaruly so that we can output the full transitions name
! to TRANS_INFO.
!
	ML=0			!Initialize line counter.
	ALLOCATE (VEC_TRANS_NAME(NLINE_MAX),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for VEC_TRANS_NAME'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
!
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
	            WRITE(LUER,*)'NLINE_MAX is too small in CMF_FLUX_SUB_V5'
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
	          VEC_C4(ML)= ABS(ATM(ID)%C4(MNL)) + ABS(ATM(ID)%C4(MNUP))
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
	          T2=0.0D0
	          CALL SET_PROF_LIMITS_V4(VEC_STRT_FREQ(ML),VEC_VDOP_MIN(ML),
	1             CHIL,ED,T,VTURB_VEC,ND,PROF_TYPE(ML),PROF_LIST_LOCATION(ML),
	1             VEC_FREQ(ML),MNL,MNUP,
	1             VEC_SPEC(ML),AT_MASS(SPECIES_LNK(ID)), ATM(ID)%ZXzV,
	1             VEC_ARAD(ML),VEC_C4(ML),TDOP,AMASS_DOP,VTURB_FIX,                !2: Garbage at present
	1             DOP_PROF_LIMIT,VOIGT_PROF_LIMIT,V_PROF_LIMIT,MAX_PROF_ED,
	1             SET_PROF_LIMS_BY_OPACITY)
	          IF(ID .EQ. 1)THEN
	            IF(ATM(ID)%XzVLEVNAME_F(2) .NE. '2___' .AND. INDEX(PROF_TYPE(ML),'DOP') .EQ. 0)THEN
	              WRITE(LUER,*)'Error in CMF_FLUX -- only Doppler profiles are implemented for split H levels'
	              STOP
	            END IF
	          END IF
	        END IF
	      END DO
	    END DO
	  END IF
	END DO
	N_LINE_FREQ=ML
	CALL ADJUST_LINE_FREQ_V2(VEC_FREQ,VEC_STRT_FREQ,VEC_TRANS_NAME,N_LINE_FREQ,LUIN)
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
	    WRITE(LUER,*)'Error in CMF_FLUX_SUB_V5'
	    WRITE(LUER,*)'FLUX_CAL_LAM_END must be > FLUX_CAL_LAM_BEG'
	    STOP
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
	DO ML=1,N_LINE_FREQ
	  IF(VEC_TRANS_TYPE(ML) .NE. 'BLANK')VEC_STRT_FREQ(ML)=VEC_FREQ(ML)
	END DO
!
! Sort lines into numerically decreaing frequency. This is used for
! outputing TRANS_INFO file. Need to sort all the VECTORS, as they
! are linked.
!
	CALL INDEXX(N_LINE_FREQ,VEC_FREQ,VEC_INDX,L_FALSE)
	CALL SORTDP(N_LINE_FREQ,VEC_STRT_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_FREQ,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_OSCIL,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_C4,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
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
	    IF(T1 .LT. 1.0D+04)THEN
	      WRITE(LUIN,'(1X,I6,2(1X,I6),2X,F10.6,2X,F10.3,2X,F10.2,4X,I7,2X,A12,3X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,PROF_LIST_LOCATION(ML),PROF_TYPE(ML),TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
	    ELSE             
	      WRITE(LUIN,'(1X,I6,2(1X,I6),2X,F10.6,2X,1P,E10.4,0P,2X,F10.2,4X,I7,2X,A12,3X,A)')
	1         ML,VEC_MNL_F(ML),VEC_MNUP_F(ML),
	1         VEC_FREQ(ML),T1,T2,PROF_LIST_LOCATION(ML),PROF_TYPE(ML),TRIM(VEC_TRANS_NAME(VEC_INDX(ML)))
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
	CALL SORTDP(N_LINE_FREQ,VEC_ARAD,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_C4,VEC_INDX,VEC_DP_WRK)
	CALL SORTDP(N_LINE_FREQ,VEC_EINA,VEC_INDX,VEC_DP_WRK)
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
! We have found all lines. If we are doing a blanketing calculation for this
! line we insert them into the continuum frequency set, otherwise the
! line is not included.
!
	DO ML=1,NCF                                          
	  FQW(ML)=NU(ML)	!FQW has temporary storage of continuum freq.
	END DO                    
	V_DOP=MINVAL(VEC_VDOP_MIN)
	CALL INS_LINE_V5(  NU,LINES_THIS_FREQ,I,NCF_MAX,
	1		  VEC_FREQ,VEC_STRT_FREQ,VEC_VDOP_MIN,VEC_TRANS_TYPE,
	1                 LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,N_LINE_FREQ,
	1                 FQW,NCF,FRAC_DOP,VINF,dV_CMF_PROF,dV_CMF_WING,
	1                 ES_WING_EXT,R_CMF_WING_EXT,L_FALSE )
!
	K=NCF		!# of continuum frequencies: Need for DET_MAIN...
	NCF=I		!Revised
!
	WRITE(LUER,*)' '
	WRITE(LUER,'(A,1X,I7)')' Number of line frequencies is:',N_LINE_FREQ
	WRITE(LUER,'(A,6X,I7)')' Number of continuum frequencies is:',NCF
!
! Used inthis context, every edge must be with V_DOP km/s of a frequency in the NU
! array.
!
	V_DOP=FRAC_DOP*MINVAL(VEC_VDOP_MIN)
	CALL DET_MAIN_CONT_FREQ(NU,NCF,FQW,K,NU_EVAL_CONT,
	1             V_DOP,DELV_CONT,COMPUTE_ALL_CROSS)
!
! Redefine frequency quadrature weights.
!
	CALL SMPTRP(NU,FQW,NCF)
	DO ML=1,NCF                                           
	  FQW(ML)=FQW(ML)*1.0D+15
	END DO
!
! Set observers frequencies. The slight fiddling in setting NU_MAX and NU_MIN 
! is done so that the CMF frequencies encompass all observers frame 
! frequencies. This allows computation of all observers fluxes allowing for 
! velocity effects.
!
! We always insert lines into the observers frame frequencies. This allows
! a blanketed model spectrum to be divided by an unblanketed model
! spectrum.
!
	T1=NU(1)/(1.0D0+2.0D0*VINF/C_KMS)
	NU_MAX_OBS=MIN(T1,NU(3))
	T1=NU(NCF)*(1.0D0+2.0D0*VINF/C_KMS)
	NU_MIN_OBS=MAX(NU(NCF-3),T1)
	CALL INS_LINE_OBS_V4(OBS_FREQ,N_OBS,NCF_MAX,
	1               VEC_FREQ,VEC_STRT_FREQ,VEC_VDOP_MIN,VEC_TRANS_TYPE,
	1               N_LINE_FREQ,SOB_FREQ_IN_OBS,
	1		NU_MAX_OBS,NU_MIN_OBS,VINF,
	1               FRAC_DOP_OBS,dV_OBS_PROF,dV_OBS_WING,dV_OBS_BIG,
	1               OBS_PRO_EXT_RAT,ES_WING_EXT,VTURB_MAX)
!
	CALL GET_PROFILE_STORAGE_LIMITS(NLINES_PROF_STORE,NFREQ_PROF_STORE,
	1         LINE_ST_INDX_IN_NU,LINE_END_INDX_IN_NU,PROF_TYPE,N_LINE_FREQ,NCF)
	CALL INIT_PROF_MODULE(ND,NLINES_PROF_STORE,NFREQ_PROF_STORE)
!
! 
! Need to calculate impact parameters, and angular quadrature weights here
! as these may be required when setting up the initial temperature
! distribution of the atmosphere (i.e. required by JGREY).
!
!
! Compute impact parameter values P
!
	CALL IMPAR(P,R,RP,NC,ND,NP)
!
! Compute the angular quadrature weights
!
	IF(TRAPFORJ)THEN
	  CALL NORDANGQW(AQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,JTRPWGT)
	  CALL NORDANGQW(HQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,HTRPWGT)
	  CALL NORDANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,KTRPWGT)
	  CALL NORDANGQW(NQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,NTRPWGT)
	  MID=.TRUE.
	  CALL GENANGQW(HMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,HTRPWGT,MID)
	  CALL GENANGQW(NMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,NTRPWGT,MID)
	ELSE
	  CALL NORDANGQW(AQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,JWEIGHT)
	  CALL NORDANGQW(HQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,HWEIGHT)
	  CALL NORDANGQW(KQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,KWEIGHT)
	  CALL NORDANGQW(NQW,R,P,WM(1,1),WM(1,3),WM(1,5),NC,ND,NP,NWEIGHT)
	  MID=.TRUE.
	  CALL GENANGQW(HMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,HWEIGHT,MID)
	  CALL GENANGQW(NMIDQW,R,P,WM(1,1),WM(1,3),WM(1,5),
	1               NC,ND,NP,NWEIGHT,MID)
	END IF
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
	  NCEXT=NC*(NPINS+1)
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
	  I=ND-DEEP
	  CALL REXT_COEF_V2(REXT,COEF,INDX,NDEXT,R,POS_IN_NEW_GRID,
	1         ND,NPINS,L_TRUE,I,ST_INTERP_INDX,END_INTERP_INDX)
	  CALL EXTEND_VTSIGMA(VEXT,TEXT,SIGMAEXT,COEF,INDX,NDEXT,
	1         V,T,SIGMA,ND)
	  CALL IMPAR(PEXT,REXT,RP,NCEXT,NDEXT,NPEXT)
!
	  ALLOCATE (AQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (KQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (NMIDQWEXT(NDEXT,NPEXT),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	    WRITE(LUER,*)'Unable to allocate memory for AQWEXT'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
!	  
! Compute the turbulent velocity and MINIMUM Doppler velocity as a function of depth.
! For the later, the iron mass of 55.8amu is assumed.
!
	  TA(1:NDEXT)=VTURB_MIN+(VTURB_MAX-VTURB_MIN)*VEXT(1:NDEXT)/VEXT(1)
	  VDOP_VEC_EXT(1:NDEXT)=SQRT( TA(1:NDEXT)**2 + 2.96*TEXT(1:NDEXT) )
!
! Note that the F2DAEXT vectors (here used as dummy variables) must be at least
! NPEXT long.
!
	  IF(TRAPFORJ)THEN
	    CALL NORDANGQW(AQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,JTRPWGT)
	    CALL NORDANGQW(HQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,HTRPWGT)
	    CALL NORDANGQW(KQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,KTRPWGT)
	    CALL NORDANGQW(NQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,NTRPWGT)
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
	    CALL NORDANGQW(NQWEXT,REXT,PEXT,F2DAEXT(1,1),F2DAEXT(1,4),
	1               F2DAEXT(1,7),NCEXT,NDEXT,NPEXT,NWEIGHT)
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
        CALL SET_POP_FOR_TWOJ(POS_IN_NEW_GRID,EDD_CONT_REC,LU_EDD,NDEXT)
!
! Allocate arrays and vectors for computing observed fluxes.
!
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
!
! Used when computing the observed fluxes. These will get overwritten
! if we do an accurate comoving frame soluton using CMF_FORM_SOL.
!
	IF(ACCURATE)THEN
	  DO LS=1,NPEXT
	    MU_AT_RMAX(LS)=SQRT( 1.0D0 -(PEXT(LS)/REXT(1))**2 )
	    HQW_AT_RMAX(LS)=HQWEXT(1,LS)
	  END DO
	  IF(PEXT(NPEXT) .EQ. REXT(1))MU_AT_RMAX(NPEXT)=0.0D0
	ELSE
	  DO LS=1,NP
	    MU_AT_RMAX(LS)=SQRT( 1.0D0 -(P(LS)/R(1))**2 )
	    HQW_AT_RMAX(LS)=HQW(1,LS)
	  END DO
	  IF(P(NP) .EQ. R(1))MU_AT_RMAX(NP)=0.0D0
	END IF
!
	ALLOCATE (ETA_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CHI_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (RJ_CMF_ST(ND,NCF),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB'
	  WRITE(LUER,*)'Unable to allocate memory for IPLUS_STORE'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
! 
!
! Set 2-photon data with current atomic models and populations.
!
	DO ID=1,NUM_IONS
	  ID_SAV=ID
	  CALL SET_TWO_PHOT_V3(TRIM(ION_ID(ID)),ID_SAV,
	1       ATM(ID)%XzVLTE,         ATM(ID)%NXzV,
	1       ATM(ID)%XzVLTE_F_ON_S,  ATM(ID)%XzVLEVNAME_F, 
	1       ATM(ID)%EDGEXzV_F,      ATM(ID)%GXzV_F,
	1       ATM(ID)%F_TO_S_XzV,     ATM(ID)%NXzV_F,     ND,
	1       ATM(ID)%ZXzV,           ATM(ID)%EQXzV,      ATM(ID)%XzV_PRES)
	END DO
!
	DTDR=(T(ND)-T(ND-1))/(R(ND-1)-R(ND))
! 
! 
!
! Electron scattering iteration loop. We perform this loop to
! take into account incoherent electron scattering.
!
	IF(RD_COHERENT_ES)NUM_ES_ITERATIONS=1
	DO ES_COUNTER=1,NUM_ES_ITERATIONS
!
! This forces CMF_FORM_SOL_V? to be used for OBSFLUX computation.
!
	  IF(ES_COUNTER .EQ. NUM_ES_ITERATIONS)LST_ITERATION=.TRUE.
!
! NB: NDEXT contains ND if ACCURATE is false. The +1 arrises since
! we also write NU(ML) out.
!
	IF(ACCURATE .OR. EDD_CONT .OR. EDD_LINECONT)THEN
	  I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
	  CALL WRITE_DIRECT_INFO_V3(NDEXT,I,DA_FILE_DATE,'EDDFACTOR',LU_EDD)
	  IF(.NOT. COMPUTE_EDDFAC)THEN
	    OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	    IF(IOS .EQ. 0)THEN
	      READ(LU_EDD,REC=5,IOSTAT=IOS)T1
	      IF(T1 .EQ. 0 .OR. IOS .NE. 0)THEN
	        WRITE(LUER,*)'Error --- All Eddfactors not'//
	1                      ' computed - will compute new F'
	        COMPUTE_EDDFAC=.TRUE.
	      END IF
	    ELSE
	      WRITE(LUER,'(/,1X,A)')'Error opening EDDFACTOR - will compute new F'
	        COMPUTE_EDDFAC=.TRUE.
	    END IF
	  END IF
	  IF(COMPUTE_EDDFAC)THEN
	      OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='NEW',RECL=I)
	    WRITE(LU_EDD,REC=1)0
	    WRITE(LU_EDD,REC=2)0
	    WRITE(LU_EDD,REC=3)0
	    WRITE(LU_EDD,REC=4)0
!
! We set record 5 to zero, to signify that the eddington factors are
! currently being computed. A non zero value signifies that all values
! have successfully been computed. (Consistent with old Eddfactor
! format a EDD_FAC can never be zero : Reason write a real number).
!
	    T1=0.0
	    WRITE(LU_EDD,REC=5)T1
	  END IF
	END IF
!
	IF(.NOT. RD_COHERENT_ES)THEN
	  I=WORD_SIZE*(ND+1)/UNIT_SIZE
	  IF(ACCURATE)I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
	  OPEN(UNIT=LU_ES,FILE='ES_J_CONV',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I,IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening ES_J_CONV'//
	1                    ' - will compute new J'
	      COHERENT_ES=.TRUE.
	    END IF
	END IF
!
! Decide whether to use an file with old J values to provide an initial 
! estimate of J with incoherent electron scattering. The options
!
!             .NOT. RD_COHERENT_ES .AND. COHERENT_ES
!
! indicate that no ES_J_CONV file exists already.
! Note that TEXT and NDEXT contain T and ND when ACCURATE is FALSE.
!
	IF(USE_OLDJ_FOR_ES .AND. .NOT. RD_COHERENT_ES .AND. COHERENT_ES)THEN
	  COHERENT_ES=RD_COHERENT_ES
	  I=ND*NCF
	  CALL COMP_J_CONV_V2(ETA_CMF_ST,I,NU,TEXT,NDEXT,NCF,LUIN,
	1        'OLD_J_FILE',
	1        EDD_CONT_REC,L_FALSE,L_TRUE,LU_ES,'ES_J_CONV')
!
! Now open the file so it can be read in the CONTINUUM loop (read in 
! COMP_JCONT).
!
	  OPEN(UNIT=LU_ES,FILE='ES_J_CONV',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',IOSTAT=IOS)
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error opening ES_J_CONV'//
	1                    ' - will compute new J'
	      COHERENT_ES=.TRUE.
	    END IF
	END IF
!
!***************************************************************************
!***************************************************************************
!
!                         CONTINUUM LOOP
!
!***************************************************************************
!***************************************************************************
!
	EDDINGTON=EDD_CONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  IF(COMPUTE_EDDFAC)THEN
	    WRITE(LU_EDD,REC=EDD_CONT_REC)ACCESS_F,NCF,NDEXT
	  ELSE
	    READ(LU_EDD,REC=EDD_CONT_REC)ACCESS_F
	  END IF
	END IF
!
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
!
	CONT_FREQ=0.0D0
!
! Define parameters to allow the Cummulative force multipler to be output at
! a function of frequency. We presently output the force multiplier every 500km/s.
!
	N_FORCE=DLOG(NU(NCF)/NU(1))/DLOG(1.0D0-500.0D0/C_KMS)
	NU_FORCE=NU(1)
	NU_FORCE_FAC=(1.0D0-500.0D0/C_KMS)
	ML_FORCE=1
!                                                                    
! Enter loop for each continuum frequency.
!
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
!
	  IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	    COMPUTE_NEW_CROSS=.TRUE.
	    CONT_FREQ=NU_EVAL_CONT(ML)
	  ELSE
	    COMPUTE_NEW_CROSS=.FALSE.
	  END IF
!
! 
!
! Section to include lines automatically with the continuum.
!
!
!  LINES_THIS_FREQ --- Logical vector [NCF] indicating whether this frequency
!                        is part of the resonance zone (i.e. Doppler profile) of 
!                        one (or more) lines.
!
! LINE_ST_INDX_IN_NU --- Integer vector [N_LINES] which specifies the starting
!                          frequency index for this lines resonance zone.
!
! LINE_END_INDX_IN_NU --- Integer vector [N_LINES] which specifies the final
!                          frequency index for this lines resonance zone.
!
! FIRST_LINE   ---- Integer specifying the index of the highest frequency
!                         line which we are taking into account in the
!                         transfer.
!
! LAST_LINE  ---- Integer specifying the index of the lowest frequency
!                         line which we are taking into account in the
!                         transfer.
!                                        
! LINE_LOC   ---- Integer array. Used to locate location of a particular line
!                         in the SIM vectors/arrays.
!
! SIM_LINE_POINTER --- Integer array --- locates the line corresponding to
!                         the indicated storage location in the SIM vectors/
!                         arrays.
!
! Check whether we have to treat another line. We use a DO WHILE, rather
! than an IF statement, to handle lines which begin at the same (upper)
! frequency.
!
!
	CALL TUNE(IONE,'ADD_LINE')
	DO WHILE( LAST_LINE .LT. N_LINE_FREQ .AND.
	1                ML .EQ. LINE_ST_INDX_IN_NU(MIN(LAST_LINE+1,N_LINE_FREQ)) )
!
! Have another line --- need to find its storage location.
!
	  I=1
	  DO WHILE(LINE_STORAGE_USED(I))
	    I=I+1
	    IF(I .GT. MAX_SIM)THEN
	      FIRST_LINE=N_LINE_FREQ
	      DO SIM_INDX=1,MAX_SIM		!Not 0 as used!
	        FIRST_LINE=MIN(FIRST_LINE,SIM_LINE_POINTER(SIM_INDX)) 
	      END DO
	      IF( ML .GT. LINE_END_INDX_IN_NU(FIRST_LINE))THEN
!
! Free up storage location for line.
!
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
!
	  SIM_INDX=I
	  LAST_LINE=LAST_LINE+1
	  LINE_STORAGE_USED(SIM_INDX)=.TRUE.
	  LINE_LOC(LAST_LINE)=SIM_INDX
	  SIM_LINE_POINTER(SIM_INDX)=LAST_LINE
!
! Have located a storage location. Now must compute all relevant quantities
! necessary to include this line in the transfer calculations.
!
	  SIM_NL(SIM_INDX)=VEC_NL(LAST_LINE)
	  SIM_NUP(SIM_INDX)=VEC_NUP(LAST_LINE)
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
!
	  EINA(SIM_INDX)=VEC_EINA(LAST_LINE)
	  OSCIL(SIM_INDX)=VEC_OSCIL(LAST_LINE)
	  FL_SIM(SIM_INDX)=VEC_FREQ(LAST_LINE)
!
	  TRANS_NAME_SIM(SIM_INDX)=TRIM(VEC_SPEC(LAST_LINE))
	  AMASS_SIM(SIM_INDX)=AMASS_ALL(NL)
!
! 
!
! Compute U_STAR_RATIO and L_STAR_RATIO which are used to switch from
! the opacity/emissivity computed with a FULL_ATOM to an equivalent form
! but written in terms of the SUPER-LEVELS. 
!
! L refers to the lower level of the transition.
! U refers to the upper level of the transition.
!
! At present we must treat each species separately (for those with both FULL
! and SUPER_LEVEL model atoms).
!
! MNL_F (MNUP_F) denotes the lower (upper) level in the full atom.
! MNL (MNUP) denotes the lower (upper) level in the super level model atom.
!
	MNL_F=VEC_MNL_F(LAST_LINE)
	MNUP_F=VEC_MNUP_F(LAST_LINE)
	DO K=1,ND
	  L_STAR_RATIO(K,SIM_INDX)=1.0D0
	  U_STAR_RATIO(K,SIM_INDX)=1.0D0
	END DO
!
! T1 is used to represent b(level)/b(super level). If no interpolation of
! the b values in a super level has been performed, this ratio will be unity .
! This ratio is NOT treated in the linearization.
!
	DO ID=1,NUM_IONS
	  IF(VEC_SPEC(LAST_LINE) .EQ. ION_ID(ID))THEN
	    MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	    MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	    DO K=1,ND
              T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzV(MNL,K))/ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
              L_STAR_RATIO(K,SIM_INDX)=T1*(ATM(ID)%W_XzV_F(MNUP_F,K)/ATM(ID)%W_XzV_F(MNL_F,K))*ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
              T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzV(MNUP,K))/ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	      U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	    END DO
	    GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	    TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1      '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1           TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	    EXIT
	  END IF
	END DO
! 
!
! Compute line opacity and emissivity for this line.
!
	  T1=OSCIL(SIM_INDX)*OPLIN
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=1,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    IF(CHIL_MAT(I,SIM_INDX) .EQ. 0.0D0)THEN
	      CHIL_MAT(I,SIM_INDX)=0.01D0*T1*POPS(NL,I)*L_STAR_RATIO(I,SIM_INDX)
	      WRITE(LUER,*)'Zero line opacity in CMFGEN_SUB'
	      WRITE(LUER,*)'This needs to be fixed'
	      J=ICHRLEN(TRANS_NAME_SIM(SIM_INDX))
	      WRITE(LUER,'(1X,A)')TRANS_NAME_SIM(SIM_INDX)(1:J)
	    END IF
	  END DO
!
! Ensure that LAST_LINE points to the next LINE that is going to be handled 
! in the BLANKETING portion of the code.
!
	  DO WHILE(LAST_LINE .LT. N_LINE_FREQ.AND.
	1            VEC_TRANS_TYPE(MIN(LAST_LINE+1,N_LINE_FREQ))(1:4) .NE. 'BLAN')
	       LAST_LINE=LAST_LINE+1
	  END DO
!	   
	END DO	!Checking whether a  new line is being added.
	CALL TUNE(ITWO,'ADD_LINE')
!
! 
!
! Check whether current frequency is a resonance frequency for each line.
!
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
! TB is the proton density.
! TC is the He+ density.
!
	  CALL TUNE(IONE,'SET_PROF')
	  TB(1:ND)=0.0D0; TC(1:ND)=0.0D0
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. ION_ID(ID) .EQ. 'HI')THEN
	      TB(1:ND)=ATM(ID)%DxzV(1:ND)
	    ELSE IF(ATM(ID)%XzV_PRES .AND. ION_ID(ID) .EQ. 'HeI')THEN
	      TC(1:ND)=ATM(ID)%DxzV(1:ND)
	    END IF
	  END DO
!
! Paralleizing this loop seems to cost CPU time.
!!$OMP PARALLEL DO SCHEDULE(DYNAMIC,20) PRIVATE (TA,I,J,T1,T3)
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
              J=SIM_LINE_POINTER(SIM_INDX); I=ML
	      T1=Z_POP(VEC_NL(J))+1.0D0; T3=0.0D0
              CALL SET_PROF_V5(TA,NU,I,
	1               LINE_ST_INDX_IN_NU(J),LINE_END_INDX_IN_NU(J),
	1               ED,TB,TC,T,VTURB_VEC,ND,
	1               PROF_TYPE(J),PROF_LIST_LOCATION(J),
	1               VEC_FREQ(J),VEC_MNL_F(J),VEC_MNUP_F(J),
	1               AMASS_SIM(SIM_INDX),T1,VEC_ARAD(J),VEC_C4(J),
	1               TDOP,AMASS_DOP,VTURB_FIX,MAX_PROF_ED,
	1               END_RES_ZONE(SIM_INDX),NORM_PROFILE,LUIN)
              LINE_PROF_SIM(1:ND,SIM_INDX)=TA(1:ND)
	    ELSE              
	      LINE_PROF_SIM(1:ND,SIM_INDX)=0.0D0
	    END IF
	  END DO
!!$OMP END PARALLEL DO                                    
	  CALL TUNE(ITWO,'SET_PROF')
!
! 
!
! Determine which method will be used to compute continuum intensity.
!
	  IF(ACCURATE .AND. ALL_FREQ)THEN       
	    THIS_FREQ_EXT=.TRUE.
	  ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	    THIS_FREQ_EXT=.TRUE.
	  ELSE        
	    THIS_FREQ_EXT=.FALSE.
	  END IF
!
! Compute opacity and emissivity.
!
	  CALL TUNE(IONE,'C_OPAC')
	  INCLUDE 'OPACITIES_V5.INC'
	  CALL TUNE(ITWO,'C_OPAC')
!
! Since resonance zones included, we must add the line opacity and 
! emissivity to the raw continuum values. We first save the pure continuum 
! opacity and emissivity. These are used in carrying the variation of J from 
! one frequency to the next.
!
	  DO I=1,ND
	    CHI_CONT(I)=CHI(I)
	    ETA_CONT(I)=ETA(I)
	  END DO
!
	  CALL TUNE(IONE,'OP_TOT')
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      DO I=1,ND
	        CHI(I)=CHI(I)+CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	        ETA(I)=ETA(I)+ETAL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	      END DO
	    END IF
	  END DO
	  CALL TUNE(ITWO,'OP_TOT')
!
! CHECK for negative line opacities. We do not distinguish between lines.  
!
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
!
	  DO I=1,ND
	    ZETA(I)=ETA(I)/CHI(I)
	    THETA(I)=CHI_SCAT(I)/CHI(I)
	  END DO
!
	  IF(LST_ITERATION .AND. ML .NE. NCF)THEN
	    T1=MAX( LOG(DENSITY(5)/DENSITY(1))/LOG(R(1)/R(5))-1.0D0,1.0D0 )
	    DO I=1,N_TAU_EDGE
	      IF(NU(ML) .GE. TAU_EDGE(I) .AND. NU(ML+1) .LT. TAU_EDGE(I))THEN
	        IF(I .EQ. 1)WRITE(LUER,'(A)')' '
	        WRITE(LUER,'(A,1P,E10.4,A,E10.3)')' Tau(Nu=',NU(ML),
	1          ') at outer boundary is:',CHI_CONT(1)*R(1)/T1
	      END IF
	    END DO
	  END IF
!
! Compute continuum intensity.
!
	  IF(COMPUTE_J)THEN
	    CALL TUNE(IONE,'COMP_JCONT')
	    INCLUDE 'COMP_JCONT_V4.INC'	
	    CALL TUNE(ITWO,'COMP_JCONT')
	  END IF
!
!
! Free up LINE storage locations. As we are only computing the line flux,
! and not its variation, we can free up the memory space as soon as we
! exit the resonance zone. This procedure is much simpler than in CMFGEN,
!
	DO SIM_INDX=1,MAX_SIM
	  IF(END_RES_ZONE(SIM_INDX))THEN
	    LINE_LOC( SIM_LINE_POINTER(SIM_INDX) )=0
	    LINE_STORAGE_USED(SIM_INDX)=.FALSE.
	    SIM_LINE_POINTER(SIM_INDX)=0
	  END IF
	END DO
!
!
! Compute flux distribution and luminosity (in L(sun)) of star. NB: For
! NORDFLUX we always assume coherent scattering.
!
!
	CALL TUNE(IONE,'FLUX_DIST')
	IF(.NOT. COMPUTE_J)THEN
	ELSE IF(THIS_FREQ_EXT .AND. .NOT. CONT_VEL)THEN
!
! Since ETAEXT is not required any more, it will be used
! flux.
!
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCEEXT,ZETAEXT,THETAEXT,RJEXT,NDEXT)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,REXT,Z,PEXT,
	1               SOURCEEXT,CHIEXT,dCHIdR,HQWEXT,ETAEXT,
	1               S1,THK_CONT,DIF,DBB,IC,NCEXT,NDEXT,NPEXT,METHOD)
	  CALL UNGRID(SOB,ND,ETAEXT,NDEXT,POS_IN_NEW_GRID)
	  SOB(2)=ETAEXT(2)				!Special case
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  OBS_FLUX(ML)=6.599341D0*SOB(1)*2.0D0		!2 DUE TO 0.5U
	ELSE IF(CONT_VEL)THEN
!
! TA is a work vector. TB initially used for extended SOB.
!
	   IF(ACCURATE)THEN
	     CALL REGRID_H(TB,REXT,RSQHNU,HFLUX_AT_OB,HFLUX_AT_IB,NDEXT,TA)
	     DO I=1,ND
	       SOB(I)=TB(POS_IN_NEW_GRID(I))
	     END DO
	   ELSE
	     CALL REGRID_H(SOB,R,RSQHNU,HFLUX_AT_OB,HFLUX_AT_IB,ND,TA)
	   END IF
	   IF(PLANE_PARALLEL .OR. PLANE_PARALLEL_NO_V)THEN
	     SOB(1)=HFLUX_AT_OB; SOB(ND)=HFLUX_AT_IB
	     SOB(1:ND)=SOB(1:ND)*R(ND)*R(ND)
	   END IF
!
	   IF(WRITE_FLUX .AND. ES_COUNTER .EQ. NUM_ES_ITERATIONS)THEN
	     IF(ML .EQ. 1)THEN
	       I=WORD_SIZE*(ND+1)/UNIT_SIZE
	       CALL WRITE_DIRECT_INFO_V3(ND,I,DA_FILE_DATE,'FLUX_FILE',LU_FLUX)
	       OPEN(UNIT=LU_FLUX,FILE='FLUX_FILE',FORM='UNFORMATTED',
	1                 ACCESS='DIRECT',STATUS='NEW',RECL=I)
	       WRITE(LU_FLUX,REC=1)0
	       WRITE(LU_FLUX,REC=2)0
	       WRITE(LU_FLUX,REC=3)0
	       WRITE(LU_FLUX,REC=4)0
	       WRITE(LU_FLUX,REC=EDD_CONT_REC)EDD_CONT_REC+1,NCF,ND
	     END IF
	     TA(1:ND)=13.1986D0*SOB(1:ND)
	     WRITE(LU_FLUX,REC=EDD_CONT_REC+ML)(SOB(I),I=1,ND),FL
	   END IF
!
	   CALL COMP_OBS_V2(IPLUS,FL,
	1           IPLUS_STORE,NU_STORE,NST_CMF,
	1           MU_AT_RMAX,HQW_AT_RMAX,OBS_FREQ,OBS_FLUX,N_OBS,
	1           V_AT_RMAX,RMAX_OBS,'IPLUS','LIN_INT',DO_CMF_REL_OBS,
	1           FIRST_OBS_COMP,NP_OBS)
!
	ELSE                          
	  S1=(ETA(1)+RJ(1)*CHI_SCAT(1))/CHI(1)
	  CALL MULTVEC(SOURCE,ZETA,THETA,RJ,ND)
	  CALL NORDFLUX(TA,TB,TC,XM,DTAU,R,Z,P,SOURCE,CHI,THETA,HQW,SOB,
	1               S1,THK_CONT,DIF,DBB,IC,NC,ND,NP,METHOD)
!
! Compute observed flux in Janskys for an object at 1 kpc .
!	(const=dex(23)*2*pi*dex(20)/(3.0856dex(21))**2 )
!
	  OBS_FLUX(ML)=6.599341D0*SOB(1)*2.0D0		!2 DUE TO 0.5U
	END IF
!
! The quantity output is the LINE FORCE multiplier M(t).
!
	IF(FL .LT. NU_FORCE .AND. LST_ITERATION .AND. WRITE_CMF_FORCE)THEN
	  INQUIRE(UNIT=82,OPENED=TMP_LOG)
	  K=5   					!ACCESS_F
	  IF(.NOT. TMP_LOG)THEN
	    I=WORD_SIZE*(ND+1)/UNIT_SIZE; J=82
	    CALL OPEN_DIR_ACC_V1(ND,I,DA_FILE_DATE,'CMF_FORCE_DATA',J)
	    WRITE(82,REC=EDD_CONT_REC)K,N_FORCE,ND
	  END IF
	  WRITE(82,REC=K-1+ML_FORCE)( LINE_FLUXMEAN(I)/STARS_LUM/ESEC(I) ,I=1,ND),NU_FORCE
	  ML_FORCE=ML_FORCE+1
	  NU_FORCE=NU_FORCE*NU_FORCE_FAC
	END IF
!
! Compute the luminosity, the FLUX mean opacity, and the ROSSELAND
! mean opacities.
!
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
	IF(COMPUTE_J)THEN
	  T1=TWOHCSQ*HDKT*FQW(ML)*(NU(ML)**4)
	  DO I=1,ND		              !(4*PI)**2*Dex(+20)/L(sun)
	    T2=SOB(I)*FQW(ML)*4.1274D-12
	    RLUMST(I)=RLUMST(I)+T2
	    J_INT(I)=J_INT(I)+RJ(I)*FQW(ML)*4.1274D-12
	    K_INT(I)=K_INT(I)+K_MOM(I)*FQW(ML)*4.1274D-12
	    FLUXMEAN(I)=FLUXMEAN(I)+T2*CHI(I)
	    LINE_FLUXMEAN(I)=LINE_FLUXMEAN(I)+T2*(CHI(I)-CHI_CONT(I))
	    T2=T1*EMHNUKT(I)/(  ( (1.0D0-EMHNUKT(I))*T(I) )**2  )
	    INT_dBdT(I)=INT_dBdT(I)+T2
	    ROSSMEAN(I)=ROSSMEAN(I)+T2/CHI(I)
	  END DO
	END IF
	CALL TUNE(ITWO,'FLUX_DIST')
!
! Compute and output line force contributed by each species.
!
	IF(WR_ION_LINE_FORCE)THEN
	  IF(.NOT. ALLOCATED(ION_LINE_FORCE))ALLOCATE (ION_LINE_FORCE(ND,NUM_IONS))
	  IF(ML .EQ. 1)ION_LINE_FORCE=0.0D0
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX))THEN
	      K=SIM_LINE_POINTER(SIM_INDX)
	      DO ID=1,NUM_IONS
	        IF(VEC_SPEC(K) .EQ. ION_ID(ID))THEN
	          DO I=1,ND
	            ION_LINE_FORCE(I,ID)=ION_LINE_FORCE(I,ID)+4.1274D-12*FQW(ML)*SOB(I)*
	1                CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	          END DO
	          EXIT
	        END IF
	      END DO
	    END IF
	  END DO
	  IF(ML .EQ. NCF)THEN
	    DO I=1,ND
	       ION_LINE_FORCE(I,:)=ION_LINE_FORCE(I,:)/ESEC(I)/STARS_LUM
	       TA(I)=FLUXMEAN(I)/ESEC(I)/STARS_LUM
	    END DO
	    OPEN(UNIT=LUIN,FILE='ION_LINE_FORCE',STATUS='UNKNOWN',ACTION='WRITE')
	      WRITE(LUIN,'(A)')' '
	      WRITE(LUIN,'(A)')' Summary of line force contributions by individual ions.'
	      WRITE(LUIN,'(A)')' Ion contributions are expressed as a % of total radiation force.'
	      WRITE(LUIN,'(A)')' At depth, continuum opacities will also be important.'
	      WRITE(LUIN,'(A)')' '
	      WRITE(LUIN,'(3X,A,500(A8))')'d',' V(km/s)','    M(t)',(TRIM(ION_ID(ID)),ID=1,NUM_IONS)
	      DO I=1,ND
	        WRITE(LUIN,'(I4,1X,F8.3,500(F8.2))')I,V(I),TA(I),
	1                    (100.0D0*ION_LINE_FORCE(I,ID)/TA(I),ID=1,NUM_IONS)
	      END DO
	      WRITE(LUIN,'(3X,A,500(A8))')'d',' V(km/s)','    M(t)',(TRIM(ION_ID(ID)),ID=1,NUM_IONS)
	    CLOSE(LUIN)
	  END IF
	END IF
!
! The current opacities and emissivities are stored for the variation of the
! radiation field at the next frequency.
!                                                  
	DO I=1,ND
	  CHI_PREV(I)=CHI_CONT(I)
	  ETA_PREV(I)=ETA_CONT(I)
	END DO
!
! Store opacities and emissivities for use in observer's frame 
! calculation. Rayleigh scattering is assumed to be coherent.
!
	IF(INCL_RAY_SCAT)THEN
	  ETA_CMF_ST(1:ND,ML)=ETA(1:ND)+CHI_RAY(1:ND)*RJ(1:ND)
	ELSE
	  ETA_CMF_ST(1:ND,ML)=ETA(1:ND)
	END IF
	CHI_CMF_ST(1:ND,ML)=CHI(1:ND)
	RJ_CMF_ST(1:ND,ML)=RJ(1:ND)
!
10000	CONTINUE
	CALL TUNE(ITWO,'MLCF')
!
! NB: We use K here, rather than ACCESS_F, so that we don't corrupt EDDFACTOR if
!      evaluate EW is set to TRUE.
! 
	IF(WRITE_ETA_AND_CHI .AND. (ES_COUNTER .EQ. NUM_ES_ITERATIONS .OR. .NOT. COMPUTE_J))THEN
	  K=5						!Use for ACCESS_F
	  I=WORD_SIZE*(ND+1)/UNIT_SIZE
	  J=82; CALL OPEN_DIR_ACC_V1(ND,I,DA_FILE_DATE,'ETA_DATA',J)
	  J=83; CALL OPEN_DIR_ACC_V1(ND,I,DA_FILE_DATE,'CHI_DATA',J)
	  WRITE(82,REC=EDD_CONT_REC)K,NCF,ND
	  WRITE(83,REC=EDD_CONT_REC)K,NCF,ND
	  DO ML=1,NCF
	    WRITE(82,REC=K-1+ML)(ETA_CMF_ST(I,ML),I=1,ND),NU(ML)
	    WRITE(83,REC=K-1+ML)(CHI_CMF_ST(I,ML),I=1,ND),NU(ML)
	  END DO
	  CLOSE(UNIT=82)
	  CLOSE(UNIT=83)
	END IF
	IF(.NOT. COMPUTE_J)STOP
!
	COMPUTE_EDDFAC=.FALSE.
!
! Compute the convolution of J with the electron redistribution function.
! If USE_FIXED_J is true, the frequencies in the EDDFACTOR file will not
! match those of this model. In addition, since we are not iterating, there
! is no need to create an ES_J_CONV file.
!
	COHERENT_ES=RD_COHERENT_ES
	IF(.NOT. COHERENT_ES .AND. .NOT. USE_FIXED_J)THEN
	   I=ND*ND
	   CALL COMP_J_CONV_V2(WM,I,NU,TEXT,NDEXT,NCF,LU_EDD,
	1         'EDDFACTOR',
	1          EDD_CONT_REC,L_FALSE,L_FALSE,LU_ES,'ES_J_CONV')
	END IF
!
	END DO		!ES_COUNTER
!
!
!
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFLUX','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(STRING,'(I10)')N_OBS
	  DO WHILE(STRING(1:1) .EQ. ' ')
	       STRING(1:)=STRING(2:)
	  END DO
	  STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	  CALL WRITV_V2(OBS_FREQ,N_OBS,ISEV,TRIM(STRING),LU_FLUX)
	  CALL WRITV_V2(OBS_FLUX,N_OBS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	  CALL WRITV(RLUMST,ND,'Luminosity',LU_FLUX)
	CLOSE(UNIT=LU_FLUX)
!
! Make sure ALL frequencies have been output to CMF_FORCE_MULT
!
	IF(WRITE_CMF_FORCE)THEN
	  K=5   					!ACCESS_F
	  DO ML=ML_FORCE,N_FORCE
	    WRITE(82,REC=K-1+ML)( LINE_FLUXMEAN(I)/STARS_LUM/ESEC(I) ,I=1,ND),NU_FORCE
	    NU_FORCE=NU_FORCE*NU_FORCE_FAC
	  END DO
	  CLOSE(UNIT=82)
	END IF
!
! Output errors that have occurred in MOM_J_CMF
!
        I=1; CALL WRITE_J_CMF_ERR(I)
!
! Compute ROSSELAND and FLUX mean opacities. Compute the respective
! optical depth scales; TA for the FLUX mean optical depth scale, 
! and TB for the ROSSELAND mean optical depth scale.
!
! T1=4 * [STEFAN BOLTZMAN CONS] * 1.0D+16 / pi
!
	T1=7.218771D+11
	DO I=1,ND
	  IF(ABS(RLUMST(I)) .LE. 1.0D-20)RLUMST(I)=1.0D-20
	  FLUXMEAN(I)=FLUXMEAN(I)/RLUMST(I)
	  INT_dBdT(I)=INT_dBdT(I)/ROSSMEAN(I)		!Program rosseland opac.
	  ROSSMEAN(I)=T1*( T(I)**3 )/ROSSMEAN(I)
	END DO
!
! Due to instabilities, the FLUX mean opacity can sometimes be
! negative. If so we can continue, but we note that the
! results may be in error, and need to be checked.
!
	IF(MINVAL(ROSSMEAN) .GT. 0.0D0)THEN
	  CALL DERIVCHI(dCHIdR,ROSSMEAN,R,ND,METHOD)
        ELSE
	  dCHIDR(1:ND)=0.0D0
	  WRITE(LUER,*)'Error in CMF_FLUX_SUB: Check MEANOPAC'
	  WRITE(LUER,*)'Rosseland mean opacity has negative values'
	END IF
	CALL NORDTAU(TA,ROSSMEAN,R,R,dCHIdR,ND)
!
	IF(MINVAL(FLUXMEAN) .GT. 0.0D0)THEN
	  CALL DERIVCHI(dCHIdR,FLUXMEAN,R,ND,METHOD)
        ELSE
	  dCHIDR(1:ND)=0.0D0
	  WRITE(LUER,*)'Warning from CMF_FLUX_SUB: Check MEANOPAC'
	  WRITE(LUER,*)'Flux mean opacity has zero or negative values'
	END IF
        CALL NORDTAU(TB,FLUXMEAN,R,R,dCHIdR,ND)
!
	CALL DERIVCHI(dCHIdR,ESEC,R,ND,METHOD)
        CALL NORDTAU(DTAU,ESEC,R,R,dCHIdR,ND)
!
	TA(ND)=0.0D0; TB(ND)=0.0D0; TC(ND)=0.0D0; DTAU(ND)=0.0D0
	CALL GEN_ASCI_OPEN(LU_OPAC,'MEANOPAC','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(LU_OPAC,
	1  '( ''     R        I    Tau(Ross)   /\Tau   Rat(Ross)'//
	1  '  Chi(Ross)  Chi(ross)  Chi(Flux)   Chi(es) '//
	1  '  Tau(Flux)  Tau(es)  Rat(Flux)  Rat(es)     Kappa   V(km/s)'' )' )
	IF(R(1) .GE. 1.0D+05)THEN
	  FMT='(ES17.10,I4,2ES10.3,ES10.2,4ES11.3,4ES10.2,2ES11.3)'
	ELSE
	  FMT='( F17.10,I4,2ES10.3,ES10.2,4ES11.3,4ES10.2,2ES11.3)'
	END IF
	  DO I=1,ND
	    IF(I .EQ. 1)THEN
	      T1=0.0D0		!Rosseland optical depth scale
	      T2=0.0D0		!Flux optical depth scale
	      T3=0.0D0		!Electron scattering optical depth scale.
	      TC(1:3)=0.0D0
	    ELSE
	      T1=T1+TA(I-1)
	      T2=T2+TB(I-1)
	      T3=T3+DTAU(I-1)
	      TC(1)=TA(I)/TA(I-1)
	      TC(2)=TB(I)/TB(I-1)
	      TC(3)=DTAU(I)/DTAU(I-1)
	    END IF
	    WRITE(LU_OPAC,FMT)R(I),I,T1,TA(I),TC(1),
	1      ROSSMEAN(I),INT_dBdT(I),FLUXMEAN(I),ESEC(I),
	1      T2,T3,TC(2),TC(3),1.0D-10*ROSSMEAN(I)/DENSITY(I),V(I)
	  END DO
	CLOSE(UNIT=LU_OPAC)
!
! Output hydrodynamical terms to allow check on radiation driving of the wind.
!
	I=18
	MEAN_ATOMIC_WEIGHT=DENSITY(ND)/POP_ATOM(ND)/ATOMIC_MASS_UNIT()
	CALL HYDRO_TERMS(POP_ATOM,R,V,T,SIGMA,ED,RLUMST,
	1                 STARS_MASS,MEAN_ATOMIC_WEIGHT,
	1		  FLUXMEAN,ESEC,I,ND)
!
! 
! If requested, convolve J with the electron-scattering redistribution
! funtion. K is used for LUIN, and LUOUIT but is not accessed.
!
	IF(.NOT. COHERENT_ES)THEN
	  I=ND*NCF
	  CALL COMP_J_CONV_V2(RJ_CMF_ST,I,NU,T,ND,NCF,
	1         K,'J PASSED VIA CALL',EDD_CONT_REC,L_FALSE,L_FALSE,
	1         K,'RETURN J VIA CALL')
	END IF
!
	ESEC(1:ND)=6.65D-15*ED(1:ND)
	DO ML=1,NCF
	  ETA_CMF_ST(1:ND,ML)=ETA_CMF_ST(1:ND,ML) +
	1                        RJ_CMF_ST(1:ND,ML)*ESEC(1:ND)
	END DO
	DEALLOCATE (RJ_CMF_ST)
!
! 
! ***************************************************************************
! ***************************************************************************
!
! Determine the number of points inserted along each ray via 
! MAX_DEL_V_RES_ZONE.
!
	DO I=1,ND
	  MAX_DEL_V_RES_ZONE(I)=VTURB_VEC(I)*FRAC_DOP*0.5D0
	END DO
!
	IF(DO_CLUMP_MODEL)THEN
	  DO ML=1,NCF
	    ETA_CMF_ST(:,ML)=ETA_CMF_ST(:,ML)*CLUMP_FAC(:)
	    CHI_CMF_ST(:,ML)=CHI_CMF_ST(:,ML)*CLUMP_FAC(:)
	  END DO
	END IF
!
!
! Ensure that MU and the quadrature weights are correctly defined.
! By definition, p * dp equals R**2 * mu * dmu. Integration over mu is
! more stable, and is to be preferred.
!
! For a plane-parallel model, we only have NC integration angles.
!
	IF(PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL)THEN
	  HQW_AT_RMAX(:)=0.0D0
	  CALL GAULEG(RZERO,RONE,MU_AT_RMAX,HQW_AT_RMAX,NC)
	  HQW_AT_RMAX(1:NC)=HQW_AT_RMAX(1:NC)*MU_AT_RMAX(1:NC)
	  HQW_AT_RMAX(1:NC)=HQW_AT_RMAX(1:NC)*R(ND)*R(ND)/R(1)/R(1)
	  P(1:NC)=R(ND)*SQRT( (1.0D0-MU_AT_RMAX(1:NC))*(1.0D0+MU_AT_RMAX(1:NC)) )
	  NC_OBS=NC; NP_OBS=NP
	  P_OBS(1:NP_OBS)=P(1:NP_OBS)
	ELSE IF(REVISE_P_GRID)THEN
	  DEALLOCATE (P_OBS)
	  I=NP*10; ALLOCATE (P_OBS(I))
	  CALL REVISE_OBS_P(P_OBS,NP_OBS,I,NC_OBS,NC,R,ND,LUIN,LUMOD)
!
	  DEALLOCATE (MU_AT_RMAX,HQW_AT_RMAX)
	  ALLOCATE (MU_AT_RMAX(NP_OBS),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (HQW_AT_RMAX(NP_OBS),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error allocate MU_AT_RMAX in refine p grid section'
	    STOP
	  END IF
	  DO LS=1,NP_OBS
	    MU_AT_RMAX(LS)=SQRT((R(1)-P_OBS(LS))*(R(1)+P_OBS(LS)))/R(1)
	  END DO
	  CALL HWEIGHT(MU_AT_RMAX,HQW_AT_RMAX,NP_OBS)
!
	ELSE
	  DO LS=1,NP
	    MU_AT_RMAX(LS)=SQRT((R(1)-P(LS))*(R(1)+P(LS)))/R(1)
	  END DO
	  CALL HWEIGHT(MU_AT_RMAX,HQW_AT_RMAX,NP)
	  NC_OBS=NC; NP_OBS=NP
	  P_OBS(1:NP_OBS)=P(1:NP_OBS)
	END IF
!
! ***************************************************************************
! ***************************************************************************
!
! Compute the observer's frame fluxes. The fluxes are returned in Janskies.
! V6 of the observer's frame code can now handle a plane-parallel model atmosphere,
! with, or without, a velocity field.
!
! We use TA for V in the call to OBS_FRAME_SUB_V6.
!
	TMP_LOG=.FALSE.
	IF(PLANE_PARALLEL_NO_V .OR. PLANE_PARALLEL)TMP_LOG=.TRUE.
	IF(PLANE_PARALLEL_NO_V)THEN
	   TA(1:ND)=0.0D0
	ELSE
	   TA(1:ND)=V(1:ND)
	END IF
	CALL OBS_FRAME_SUB_V9(ETA_CMF_ST,CHI_CMF_ST,NU,
	1            R,TA,T,ED,ND,NCF,
	1            P_OBS,MU_AT_RMAX,HQW_AT_RMAX,NC_OBS,NP_OBS,
	1            OBS_FREQ,OBS_FLUX,N_OBS,
	1            MAX_DEL_V_RES_ZONE,OBS_TAU_MAX,OBS_ES_DTAU,
	1            N_INS_OBS,OBS_INT_METHOD,
	1            TAU_REF,WRITE_RTAU,WRITE_IP,WRITE_dFR,
	1            DO_REL_IN_OBSFRAME,TMP_LOG)
!
	CALL GEN_ASCI_OPEN(LU_FLUX,'OBSFRAME','UNKNOWN',' ',' ',IZERO,IOS)
	  WRITE(STRING,'(I10)')N_OBS
	  DO WHILE(STRING(1:1) .EQ. ' ')
	       STRING(1:)=STRING(2:)
	  END DO
	  STRING='Continuum Frequencies ( '//TRIM(STRING)//' )'
	  CALL WRITV_V2(OBS_FREQ,N_OBS,ISEV,TRIM(STRING),LU_FLUX)
	  CALL WRITV_V2(OBS_FLUX,N_OBS,IFOUR,'Observed intensity (Janskys)',LU_FLUX)
	CLOSE(UNIT=LU_FLUX)
!
! At present the line EW's are automatically computed when the
! CONTINUM Flux is evaluated.
!
	IF(.NOT. DO_SOBOLEV_LINES)RETURN
!
! This section of the program compute the SOBOLEV line EW's. These
! are only approximate, and can be used to determine which lines are
! imortant contributors to INDIVIDUAL emission features. NB: If a line
! has a P Cygni profile it may be an important contributor but give an
! EW of zero.
!
! Set access counter for Continuum Eddington file, ensuring that file
! is OPENED. FIRST is used as a temporary LOGICAL variable.
!
	COMPUTE_EDDFAC=.TRUE.
	EDDINGTON=EDD_LINECONT
	IF(ACCURATE .OR. EDDINGTON)THEN
	  INQUIRE(UNIT=LU_EDD,OPENED=FIRST)
	  IF(.NOT. FIRST)THEN
	    I=WORD_SIZE*(NDEXT+1)/UNIT_SIZE
	    OPEN(UNIT=LU_EDD,FILE='EDDFACTOR',FORM='UNFORMATTED',
	1       ACCESS='DIRECT',STATUS='OLD',RECL=I)
	  END IF
	  WRITE(LU_EDD,REC=4)ACCESS_F,N_LINE_FREQ,NDEXT
	END IF
!
! Open file to store EW data
!
	CALL GEN_ASCI_OPEN(LU_EW,'EWDATA','UNKNOWN',' ',' ',IZERO,IOS)
	WRITE(LU_EW,'(4X,A,3X,A,6X,A,2X,A,3X,A)')'Lam(Ang)','C. Flux',
	1                  'EW(Ang)','Sob','Trans. Name'
!
! Zero the vector that will be used to store the force-multiplier computed
! for all lines using the SOBOLEV approximation.
!
	FORCE_MULT(1:ND)=0.0D0
	N_FORCE=DLOG(NU(NCF)/NU(1))/DLOG(1.0D0- 500.0D0/C_KMS)
	NU_FORCE=NU(1)
	NU_FORCE_FAC=(1.0D0-500.0D0/C_KMS)
	ML_FORCE=1
!
! Enter line loop.
!
	LINE_INDX=1
	OVERLAP=.FALSE.
	CALL TUNE(1,'SOB_EW')
	DO WHILE (LINE_INDX .LE. N_LINE_FREQ)
!
	IF(OVERLAP)THEN
	  TMP_MAX_SIM=MAX_SIM
	ELSE
	  TMP_MAX_SIM=1
	  OVER_FREQ_DIF=0.0D0
	END IF
!
! Determine next line (or lines), and store line parameters (eg frequency,
! Einstein A coefficient etc) in appropriate locations for the subsequent 
! treatment of line.
!
! The temporary variable J [=MAX(N_LINE_FREQ,LINE_INDX+SIM_INDX) ] is used so that
! VEC_FREQ(N_LINE_FREQ+1) is not accessed (remember that the arguments of
! a "DO WHILE" can be done in any order.
!
	SIM_INDX=0
	J=LINE_INDX
	SOBOLEV=.TRUE.
	DO WHILE( LINE_INDX+SIM_INDX .LE. N_LINE_FREQ .AND.
	1      (VEC_FREQ(LINE_INDX)-VEC_FREQ(J))/VEC_FREQ(LINE_INDX) 
	1          .LE. OVER_FREQ_DIF 
	1          .AND. SIM_INDX .LT. TMP_MAX_SIM )
	  SIM_INDX=SIM_INDX+1
	  SIM_NL(SIM_INDX)=VEC_NL(LINE_INDX+SIM_INDX-1)
	  SIM_NUP(SIM_INDX)=VEC_NUP(LINE_INDX+SIM_INDX-1)
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  EINA(SIM_INDX)=VEC_EINA(LINE_INDX)
	  OSCIL(SIM_INDX)=VEC_OSCIL(LINE_INDX)
	  FL_SIM(SIM_INDX)=VEC_FREQ(LINE_INDX)
	  SIM_LINE_POINTER(SIM_INDX)=LINE_INDX
	  TRANS_NAME_SIM(SIM_INDX)=TRIM(VEC_SPEC(LINE_INDX))
	  J=MIN(LINE_INDX+SIM_INDX,N_LINE_FREQ)
	END DO
	NUM_SIM_LINES=SIM_INDX
!
! Determine center of blend. Used as continuum wavelength.
!
	FL=0.0D0
	DO SIM_INDX=1,NUM_SIM_LINES
	  FL=FL+FL_SIM(SIM_INDX)
	END DO
	FL=FL/NUM_SIM_LINES
!
	CONT_FREQ=FL
	COMPUTE_NEW_CROSS=.TRUE.
!
! Determine method to compute continuum intensity.
!
	IF(ACCURATE .AND. ALL_FREQ)THEN
	  THIS_FREQ_EXT=.TRUE.
	ELSE IF( ACCURATE .AND. FL .GT. ACC_FREQ_END )THEN
	  THIS_FREQ_EXT=.TRUE.
	ELSE
	  THIS_FREQ_EXT=.FALSE.
	END IF
! 
! Compute U_STAR_RATIO and L_STAR_RATIO which are used to switch from
! the opacity/emissivity computed with a FULL_ATOM to an equivalent form
! but written in terms of the SUPER-LEVELS. 
!
! L refers to the lower level of the transition.
! U refers to the upper level of the transition.
!
! At present we must treat each species separately (for those with both FULL
! and SUPER_LEVEL model atoms).
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  I=SIM_LINE_POINTER(SIM_INDX)
	  MNL_F=VEC_MNL_F(I)
	  MNUP_F=VEC_MNUP_F(I)
	  L_STAR_RATIO(1:ND,SIM_INDX)=0.0D0
	  U_STAR_RATIO(1:ND,SIM_INDX)=0.0D0
!
	  DO ID=1,NUM_IONS
	    IF(VEC_SPEC(I) .EQ. ION_ID(ID))THEN
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
! T1 is used to represent b(level)/b(super level). If no interpolation of
! the b values in a super level has been performed, this ratio will be unity .
! This ratio is NOT treated in the linearization.
!
	      DO K=1,ND
	        T1=(ATM(ID)%XzV_F(MNL_F,K)/ATM(ID)%XzV(MNL,K))/ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
	        L_STAR_RATIO(K,SIM_INDX)=T1*(ATM(ID)%W_XzV_F(MNUP_F,K)/ATM(ID)%W_XzV_F(MNL_F,K))*ATM(ID)%XzVLTE_F_ON_S(MNL_F,K)
	        T2=(ATM(ID)%XzV_F(MNUP_F,K)/ATM(ID)%XzV(MNUP,K))/ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
	        U_STAR_RATIO(K,SIM_INDX)=T2*ATM(ID)%XzVLTE_F_ON_S(MNUP_F,K)
    	      END DO
	      GLDGU(SIM_INDX)=ATM(ID)%GXzV_F(MNL_F)/ATM(ID)%GXzV_F(MNUP_F)
	      TRANS_NAME_SIM(SIM_INDX)=TRIM(TRANS_NAME_SIM(SIM_INDX))//
	1             '('//TRIM(ATM(ID)%XzVLEVNAME_F(MNUP_F))//'-'//
	1                  TRIM(ATM(ID)%XzVLEVNAME_F(MNL_F))//')'
	      EXIT
	    END IF
	  END DO
	END DO
!
! Compute line opacity and emissivity. 
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  T1=OSCIL(SIM_INDX)*OPLIN
	  T2=FL_SIM(SIM_INDX)*EINA(SIM_INDX)*EMLIN
	  NL=SIM_NL(SIM_INDX)
	  NUP=SIM_NUP(SIM_INDX)
	  DO I=1,ND
	    CHIL_MAT(I,SIM_INDX)=T1*(L_STAR_RATIO(I,SIM_INDX)*POPS(NL,I)-
	1            GLDGU(SIM_INDX)*U_STAR_RATIO(I,SIM_INDX)*POPS(NUP,I))
	    ETAL_MAT(I,SIM_INDX)=T2*POPS(NUP,I)*U_STAR_RATIO(I,SIM_INDX)
	    NEG_OPACITY(I)=.FALSE.
	  END DO
	END DO
! 
!
	SECTION='LINE'
!
! Compute continuum opacity and emissivity at the line frequency.
!
	INCLUDE 'OPACITIES_V5.INC'
!
! Compute continuum intensity.
!
	INCLUDE 'COMP_JCONT_V4.INC'	
!
! SOURCE is used by SOBEW
!
	DO I=1,ND
	  SOURCE(I)=ZETA(I)+RJ(I)*THETA(I)
	END DO
! 
!
! Scale emissivity because of different frequencies.
!
	DO SIM_INDX=1,NUM_SIM_LINES
	  DO I=1,ND
	    T1=FL**3/(EXP(HDKT*FL/T(I))-1.0D0)
	    T2=FL_SIM(SIM_INDX)**3/(EXP(HDKT*FL_SIM(SIM_INDX)/T(I))-1.0D0)
	    BB_COR(I,SIM_INDX)=T1/T2
	  END DO
	END DO
!
! AT present overlapping lines in the CMF frame is not installed. We still
! use CHIL and ETAL --- not CHIL_MAT etc. 
!
	IF(SOBOLEV .OR. NUM_SIM_LINES .EQ. 1)THEN
	  DO I=1,ND
	    CHIL(I)=0.0D0
	    ETAL(I)=0.0D0
	    DO SIM_INDX=1,NUM_SIM_LINES
	      CHIL(I)=CHIL(I)+CHIL_MAT(I,SIM_INDX)
	      ETAL(I)=ETAL(I)+ETAL_MAT(I,SIM_INDX)*BB_COR(I,SIM_INDX)
	    END DO
	  END DO
	END IF
!
        IF(CHECK_LINE_OPAC .AND. SOBOLEV)THEN
	  T1=3.0D-10/FL
	  FIRST_NEG=.TRUE.
	  DO I=1,ND
	    T2=CHIL(I)*T1*R(I)/V(I)		!Tau(Sobolev) for dlnV/dlnr=1
	    IF(T2 .LT. -0.5)THEN
	      NEG_OPACITY(I)=.TRUE.
	      IF(LST_ITERATION)THEN
	        IF(FIRST_NEG)THEN
	          J=ICHRLEN(TRANS_NAME_SIM(1))
	          WRITE(LU_NEG,'(1X,A)')TRANS_NAME_SIM(1)(1:J)
	          FIRST_NEG=.FALSE.
	        END IF
	        WRITE(LU_NEG,2290)I,T2,ED(I)
2290	        FORMAT(1X,'I= ',I3,'  : TAU(Sob)= ',1P,E9.2,'  : Ne=',E8.2)
	        WRITE(LU_NEG,*)CHIL(I),POPS(NL,I),R(I),V(I),T1
	      END IF
	      CHIL(I)=1.0D0	!Reset after we output its value.
	    END IF
	  END DO
!
! Note that both CHIL and CHIL_MAT are used by SOBJBAR_SIM to compute ZNET.
! The following ensures consistency. Alternative is to make positive
! only those opacities that are negative --- then need  NEG_OPACITY indicator
! for each line?
!
	  DO SIM_INDX=1,NUM_SIM_LINES
	    DO I=1,ND
              IF(NEG_OPACITY(I))CHIL_MAT(I,SIM_INDX)=1.0D0/NUM_SIM_LINES
	    END DO
	  END DO
!
	END IF
!
! Estimate the line EW using a Modified Sobolev approximation.
!
	CHI_CLUMP(1:ND)=CHI(1:ND)*CLUMP_FAC(1:ND)
	ETA_CLUMP(1:ND)=ETA(1:ND)*CLUMP_FAC(1:ND)
	ESEC_CLUMP(1:ND)=ESEC(1:ND)*CLUMP_FAC(1:ND)
	CHI_SCAT_CLUMP(1:ND)=CHI_SCAT(1:ND)*CLUMP_FAC(1:ND)
!
	CHIL(1:ND)=CHIL(1:ND)*CLUMP_FAC(1:ND)
	ETAL(1:ND)=ETAL(1:ND)*CLUMP_FAC(1:ND)
!
	IF(FL .LT. NU_FORCE .AND. WRITE_SOB_FORCE)THEN
	  INQUIRE(UNIT=82,OPENED=TMP_LOG)
	  K=5   					!ACCESS_F
	  IF(.NOT. TMP_LOG)THEN
	    I=WORD_SIZE*(ND+1)/UNIT_SIZE; J=82
	    CALL OPEN_DIR_ACC_V1(ND,I,DA_FILE_DATE,'SOB_FORCE_DATA',J)
	    WRITE(82,REC=EDD_CONT_REC)K,N_FORCE,ND
	  END IF
	  WRITE(82,REC=K-1+ML_FORCE)(FORCE_MULT(I),I=1,ND),NU_FORCE
	  ML_FORCE=ML_FORCE+1
	  NU_FORCE=NU_FORCE*NU_FORCE_FAC
	END IF
!
! We use TA as a temporary vector which indicates the origin
! of the line emission. Not required in this code as used only
! for display purposes.
!
	CALL SOBEW_GRAD_V2(SOURCE,CHI_CLUMP,CHI_SCAT_CLUMP,CHIL,ETAL,
	1            V,SIGMA,R,P,FORCE_MULT,STARS_LUM,AQW,HQW,TA,EW,CONT_INT,
	1            FL,INNER_BND_METH,DBB,IC,THK_CONT,L_FALSE,NC,NP,ND,METHOD)
!
	IF(ABS(EW) .GE. EW_CUT_OFF)THEN
	  T1=LAMVACAIR(FL_SIM(1)) 		!Wavelength(Angstroms)
	  DO SIM_INDX=1,NUM_SIM_LINES
	    NUP=SIM_NUP(SIM_INDX)
	    CALL EW_FORMAT(EW_STRING,TRANS_NAME_SIM(SIM_INDX),T1,
	1                     CONT_INT,EW,SOBOLEV)
	    IF(SIM_INDX .NE. 1)THEN
	      I=INDEX(EW_STRING,'  ')
	      EW_STRING(3:I+1)=EW_STRING(1:I-1)
	      EW_STRING(2:3)='##'
	    END IF
	    L=ICHRLEN(EW_STRING)
	    WRITE(LU_EW,'(A)')EW_STRING(1:L)
	  END DO
	END IF
!
! Update line counter so that we move onto next line
!
	LINE_INDX=LINE_INDX+NUM_SIM_LINES
!
	END DO		!Line loop
!
! Outout the Sobolev Line-Force multiplier
!
	IF(WRITE_SOB_FORCE)THEN
	  K=5   					!ACCESS_F
	  DO ML=ML_FORCE,N_FORCE
	    WRITE(82,REC=K-1+ML)(FORCE_MULT(I),I=1,ND),NU_FORCE
	    NU_FORCE=NU_FORCE*NU_FORCE_FAC
	  END DO
	  CLOSE(UNIT=82)
	  OPEN(UNIT=82,FILE='SOB_FORCE_MULT',STATUS='UNKNOWN')
	   DO I=1,ND
	      WRITE(82,'(1X,I3,2X,3ES14.6)')I,R(I),V(I),FORCE_MULT(I)
	    END DO
	  CLOSE(UNIT=82)
	END IF
!	
	CALL TUNE(2,'SOB_EW')
!
	RETURN
!
	END
