!
! Subroutine to provide population, T, and ED estimates for a new model.
! Routine replaces main IF(NEWMOD) ... END IF section in CMFGEN_SUB.
! This routine will be easier to revise than CMFGEN.
!
! Options:
!
!       GRID:
!            Departure coefficients, T, Ne etc read from old model.
!            Model assumed to be same as current model.
!
!       GRID=.FALSE.
!            Used for a completely new model. Ne is estimated and the 
!            departure coefficients are assumed to be the same on the
!            electron density scale. T is also estimated on the electron
!            density scale. For TAU < TAU_INIT_TAU, T is not changed
!            from that read in (smooth variation around T_INIT_TAU).
!
!            When
!               ITERATE_INIT_T=.TRUE.
!            we compute the grey temperature structure and iterate to
!            improve T. Only done when GRID=.FALSE. If GREY_SCK_FAC_IN
!            is present, we modify TGREY according to
!                        TGREY = TGREY . (T/TGREY)_old
!
	SUBROUTINE SET_NEW_MODEL_ESTIMATES(POPS,Z_POP,NU,NU_EVAL_CONT,FQW,
	1            LUER,LUIN,NC,ND,NP,NT,NCF,N_LINE_FREQ,MAX_SIM)
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE CONTROL_VARIABLE_MOD
	USE LINE_VEC_MOD
	USE LINE_MOD
	IMPLICIT NONE
!
! Altered 31-Jan-2016 : Added SAVED_TWO_PHOT_METHOD.
! Altered 07-Apr-2015 : Changed to SET_TWO_PHOT_V3.
! Altered 05-Apr-2011 : Many changes done in order to facilitate the USE of LTE populations
!                         over a wider dynamic range (18-Dec-2010).
!                         Single routine (REGRID_LOG_DC_V1) now used to read all departure coefficient files.
!                         CNVT_FR_DC_V2 and PAR_FUN_V4 now called.
! Altered 15-Jan-2009 : Changed REGRID_T_ED to V2
!                          'R' option in with non 'GRID' option no longer effects T & ED.
! Created 17-Dec-2004
! Altered 06-Jun-2005 : Call to SUP_TO FULL inserted to get better consistency.
!                          Only done when GRID=.FALSE.
!
	INTEGER NC
	INTEGER ND
	INTEGER NP
	INTEGER NT
!
	INTEGER NCF				!Number of continuum frequencies
	INTEGER N_LINE_FREQ			!Number of lines
	INTEGER LUER				!Unit for error messages
	INTEGER LUIN				!Unit for input
	INTEGER MAX_SIM				!Maximum number of lines that can be treated simultaneously.
!
	REAL*8 POPS(NT,ND)
	REAL*8 Z_POP(NT)			!Vector containing Z of atom/ion (not core)
!
	REAL*8 FQW(NCF)
	REAL*8 NU_EVAL_CONT(NCF)
	REAL*8 NU(NCF)
!
! These are set in CMFGEN.
!
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 OPLIN,EMLIN
!
! Arrays for improving on the initial T structure --- partition functions.
! Need one for each atomic species.
!
        REAL*8, ALLOCATABLE :: U_PAR_FN(:,:)
        REAL*8, ALLOCATABLE :: PHI_PAR_FN(:,:)
        REAL*8, ALLOCATABLE :: Z_PAR_FN(:)
        REAL*8 SPEC_DEN(ND,NUM_SPECIES)         !Used by ELEC_PREP
	REAL*8 AT_NO_VEC(ND,NUM_SPECIES)
!
	REAL*8 TGREY(ND)			!Grey temperature structure
	REAL*8 T_SAVE(ND)
	REAL*8 ROSSMEAN(ND)			!Rosseland mean opacity
	REAL*8 PLANCKMEAN(ND)			!Planck mean opacity
!
! These are all work vectors.
!
	REAL*8 RJ(ND)
	REAL*8 DTAU(ND)
	REAL*8 Z(ND)
	REAL*8 dCHIdR(ND)
!
	REAL*8 TA(ND)
	REAL*8 TB(ND)
	REAL*8 TC(ND)
	REAL*8 QH(ND)
	REAL*8 Q(ND)
	REAL*8 GAM(ND)
	REAL*8 GAMH(ND)
	REAL*8 H(ND)
	REAL*8 SOB(ND)
	REAL*8 XM(ND)
	REAL*8 FEDD(ND)
!
	REAL*8 T1,T2,T3
	REAL*8 HBC_J
	REAL*8 NU_DOP
	REAL*8 FL		!Current frequency
	REAL*8 CONT_FREQ	!Frequency at which current ETA/CHI was evaluated
!
	INTEGER FREQ_INDX 	!Index of current frequency in NU
	INTEGER ML		!Same as FREQ_INDX
	INTEGER LAST_LINE	!Next line to be accessed
	INTEGER GREY_IOS	!Used to return error if GREY_SCL_FAC_IN can't be read.
!
	INTEGER I,J,K,L
	INTEGER ISPEC
	INTEGER ID
	INTEGER ID_SAV
	INTEGER NL,NUP
	INTEGER MNL_F,MNUP_F
	INTEGER MNL,MNUP
	INTEGER MAIN_COUNTER
!
	LOGICAL LST_DEPTH_ONLY
	LOGICAL FIRST
	LOGICAL COMPUTED
!
	CHARACTER*80 TMP_STRING
	CHARACTER*20 SECTION
	CHARACTER(LEN=12) SAVED_TWO_PHOTON_METHOD
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
	LST_DEPTH_ONLY=.FALSE.
	SECTION='CONTINUUM'
	GREY_IOS=0
	SAVED_TWO_PHOTON_METHOD=TWO_PHOTON_METHOD
	TWO_PHOTON_METHOD='LTE'
!
! Compute temperature distribution and populations. The call ' 'WSC are
! too allow NDOLD in the input files to be larger than ND.
! The first call to REGRIDWS is effectively used to compute DHeI only.
!
	IF(GRID) THEN
	  WRITE(LUER,'(/,A,/)')' Using direct interpolation option (i.e. GRID) for new model.'
!
	  DO ID=1,NUM_IONS-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      TMP_STRING=TRIM(ION_ID(ID))//'_IN'
	      ISPEC=SPECIES_LNK(ID)
	      CALL REGRID_LOG_DC_V1( ATM(ID)%XzV_F,R,ED,T, ATM(ID)%DXzV_F,CLUMP_FAC,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1             POP_SPECIES(1,ISPEC),ATM(ID)%NXzV_F,ND,LUIN,'R',TMP_STRING)
	    END IF
	  END DO
!
! Regrid the temperature and the electron density. By using this call the
! last species can be taken from a different model to the H, He populations 
! etc. Normally T_IN can be the same as He2_IN (i.e. any input departure
! coefficient file). 
!
	  CALL REGRID_T_ED(R,ED,T,POP_ATOM,ND,'T_IN')
!
!
! SPEC_DEN will contain the density of each species, while
! AT_NO_VEC will contain the the atomic number. These are set at all depths.
!
!	  IF(VAR_MDOT)THEN
!	    DO I=1,ND
!	      ED(I)=POP_SPECIES(I,1)+2.0D0*POP_SPECIES(I,2)
!	    END DO
!	  END IF
!
	  IF(VAR_MDOT)THEN
	    I=0
	    DO ISPEC=1,NUM_SPECIES
	      CALL ELEC_PREP(SPEC_DEN,AT_NO_VEC,I,NUM_SPECIES,
	1                POP_SPECIES(1,ISPEC),AT_NO(ISPEC),SPECIES_PRES(ISPEC),ND)
	    END DO
	    CALL GETELEC_V2(SPEC_DEN,AT_NO_VEC,I,ED,ND,LUIN,'GAMMAS_IN')
	  END IF
!
! 
	ELSE
	  WRITE(LUER,*)'Using NON-GRID option for new model.'
	  IF(.NOT. DO_POP_SCALE)THEN
	     DO_POP_SCALE=.TRUE.
	     WRITE(LUER,*)'Warning - setting DO_POP_SCALE=.TRUE. in SET_NEW_MODEL_ESTIMATES'
	     WRITE(LUER,*)'DO_POP_SCALE adjusted as non-GRID option.'
	  END IF
!
	  IF(SN_HYDRO_MODEL)THEN
!
! Set T and ED
!
	    IF(SN_T_OPTION .EQ. 'USE_T_IN')THEN
	      CALL REGRID_T_ED_V2(R,ED,T,POP_ATOM,VOL_EXP_FAC,ND,'T_IN')
	    ELSE IF(SN_T_OPTION .EQ. 'USE_HYDRO') THEN
!
! For this options, we need to do nothing as T and ED were set by SET_ABUND_CLUMP.
!
	    ELSE
	      WRITE(LUER,'(A)')'Unrecognized SN option for determin initial T'
	      WRITE(LUER,'(A)')'Availabvle options are USE_T_IN and USE_HYDRO'
	      STOP
	    END IF
	  ELSE
!
	    IF(SET_LTE_AS_INIT_ESTIMATES)THEN
	      OPEN(UNIT=LUIN,STATUS='OLD',ACTION='READ',FILE='FIN_CAL_GRID')
	        READ(LUIN,'(A)')TMP_STRING
		DO I=1,ND
	          READ(LUIN,*)T1,T1,T1,T1,T(I),T2,ED(I)
	          ED(I)=ED(I)*T2
		END DO
	      CLOSE(LUIN)
	    ELSE
!
! SPEC_DEN will contain the density of each species, while
! AT_NO_VEC will contain the the atomic number. These are set at all depths.
!
	      I=0
	      DO ISPEC=1,NUM_SPECIES
	        CALL ELEC_PREP(SPEC_DEN,AT_NO_VEC,I,NUM_SPECIES,
	1                POP_SPECIES(1,ISPEC),AT_NO(ISPEC),SPECIES_PRES(ISPEC),ND)
	      END DO
	      CALL GETELEC_V2(SPEC_DEN,AT_NO_VEC,I,ED,ND,LUIN,'GAMMAS_IN')
!
	      IF(INTERP_T_ON_R_GRID)THEN
!
! We do a simple linera interpolation on R, and assume that the Grey iteration procdure
! will correct T in the inner region. Usefull when large changes in T.
!
	        CALL REGRID_T_ED(R,TA,T,POP_ATOM,ND,'T_IN')
	      ELSE
!                         
! The INIT_TEMP routine assumes that the T can interpolated using 
! a Spherical TAU scale computed using the electron scattering opacity.
!
	        CALL INIT_TEMP_V2(R,ED,CLUMP_FAC,T,LUM,T_INIT_TAU,ND,LUIN,'T_IN')
	      END IF
	    END IF
	  END IF
!
! If clumping is present we interpret on the departure coefficients using the
! electron density as the independent variable. At present we correct the
! electron density for clumping (i.e. ED(I)*CLUMP_FAC(I)), which we store in TC.
! It may be better to regrid on the actual electron densities.
!
! NB: ATM(ID)$XzV_F is to LOG(dep. coef.)
!
	  IF(DC_INTERP_METHOD .EQ. 'LTE')THEN
	    WRITE(LUER,*)'LTE assumed for departure coefficients.'
	    T1=T_EXCITE_MIN
	    CALL DET_LTE_ED(T1,ND,DO_LEV_DISSOLUTION)
	    DO ID=1,NUM_IONS-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        CALL SET_DC_LTE_V2(ATM(ID)%XzV_F,ATM(ID)%DXzV_F,ATM(ID)%EDGEXzV_F,ATM(ID)%NXzV_F,T,T1,ND)
	        ATM(ID)%DXzV_F=1.0D-200
	      END IF
	    END DO
!
	  ELSE IF(DC_INTERP_METHOD .EQ. 'SPH_TAU')THEN
	    WRITE(LUER,*)'Departure coefficients assumed to be function of Tau(spherical).'
	    DO ID=1,NUM_IONS-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        TMP_STRING=TRIM(ION_ID(ID))//'_IN'
	        ISPEC=SPECIES_LNK(ID)
	        CALL REGRID_LOG_DC_V1( ATM(ID)%XzV_F,R,ED,T, ATM(ID)%DXzV_F,CLUMP_FAC,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1             POP_SPECIES(1,ISPEC),ATM(ID)%NXzV_F,ND,LUIN,'SPH_TAU',TMP_STRING)
	      END IF
	    END DO
!
	  ELSE IF(DC_INTERP_METHOD .EQ. 'ED')THEN
	    WRITE(LUER,*)'Departure coefficients assumed to be function of Ne.'
	    TC(1:ND)=ED(1:ND)*CLUMP_FAC(1:ND)
	    DO ID=1,NUM_IONS-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        TMP_STRING=TRIM(ION_ID(ID))//'_IN'
	        ISPEC=SPECIES_LNK(ID)
	        CALL REGRID_LOG_DC_V1( ATM(ID)%XzV_F,R,TC,T, ATM(ID)%DXzV_F,CLUMP_FAC,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1             POP_SPECIES(1,ISPEC),ATM(ID)%NXzV_F,ND,LUIN,'ED',TMP_STRING)
	      END IF
	    END DO
!
	  ELSE IF(DC_INTERP_METHOD .EQ. 'R')THEN
!
! We use TA for ED and TB for T since ED and T have already been set.
!
	    WRITE(LUER,*)'Departure coefficients assumed to be function of R.'
	    DO ID=1,NUM_IONS-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        TMP_STRING=TRIM(ION_ID(ID))//'_IN'
	        ISPEC=SPECIES_LNK(ID)
	        CALL REGRID_LOG_DC_V1( ATM(ID)%XzV_F,R,ED,T, ATM(ID)%DXzV_F,CLUMP_FAC,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1             POP_SPECIES(1,ISPEC),ATM(ID)%NXzV_F,ND,LUIN,'R',TMP_STRING)
	      END IF
	    END DO
!
	  ELSE IF(DC_INTERP_METHOD .EQ. 'RTX')THEN
!
! NB: T must have been previoulsy computed.
!
	    WRITE(LUER,*)'Excitation temperatures assumed to be function of R.'
	    DO ID=1,NUM_IONS-1
	      IF(ATM(ID)%XzV_PRES)THEN
	        TMP_STRING=TRIM(ION_ID(ID))//'_IN'
	        ISPEC=SPECIES_LNK(ID)
	        CALL REGRID_LOG_DC_V1( ATM(ID)%XzV_F,R,ED,T, ATM(ID)%DXzV_F,CLUMP_FAC,
	1             ATM(ID)%EDGEXzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%INT_SEQ_XzV,
	1             POP_SPECIES(1,ISPEC),ATM(ID)%NXzV_F,ND,LUIN,'RTX',TMP_STRING)
	      END IF
	    END DO
	  END IF
!
	END IF				!if(grid)
!
! 
!
! Compute vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD.
!
! As a first estimate of POPION, we assume all species are fully ionized.
! The POPION is just the number of atoms.
!
	DO I=1,ND
	  POPION(I)=POP_ATOM(I)
	END DO
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! Compute the LTE populations and convert from departure coefficients to
! populations. The populations are scaled to ensure number conservation -
! during the scaling the ionization remains fixed. On each call to CNVT_FR_DC,
! TA is incremented by the population. IF FIRST is .TRUE., TA is zeroed first.
! The second flag indicates whether to add in the ION contribution, which
! will also be added as the ground state population of the next species.
! If it is FALSE, the higher ionization stage is assumed not to be present,
! and DION is added in.
!
! TB is used as a dummy vector when we are dealing with the lowest ionization
! stage. It is returned with the ground state population.
!
	DO ISPEC=1,NUM_SPECIES
	  FIRST=.TRUE.
	  DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	    IF(ATM(ID)%XzV_PRES)THEN
	      CALL LTEPOP_WLD_V2(ATM(ID)%XzVLTE_F, ATM(ID)%LOG_XzVLTE_F, ATM(ID)%W_XzV_F,
	1              ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,
	1              ATM(ID)%ZXzV,      ATM(ID)%GIONXzV_F,
	1              ATM(ID)%NXzV_F,    ATM(ID)%DXzV_F,     ED,T,ND)
	      CALL CNVT_FR_DC_V2(ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
	1              ATM(ID)%DXzV_F,    ATM(ID)%NXzV_F,
	1              TB,                TA,ND,
	1              FIRST,             ATM(ID+1)%XzV_PRES)
	      IF(ID .NE. SPECIES_BEG_ID(ISPEC))ATM(ID-1)%DXzV_F(1:ND)=TB(1:ND)
  	    END IF
	  END DO
!
! Now scale the population for EACH species to ensure that the species
! conservation equation is satisfied.
!
! This option should always be set if new model with T iteration 
! and correction. With the GRID=T option it may provide a convenient
! method for reducing the populations temporarily (in conjunction with
! DISPGEN) in the outer layers to overcome a large jump in optical
! depth.
!
	  IF(DO_POP_SCALE)THEN
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      CALL SCALE_POPS(ATM(ID)%XzV_F,ATM(ID)%DXzV_F,
	1              POP_SPECIES(1,ISPEC),TA,ATM(ID)%NXzV_F,ND)
	    END DO
	  END IF
	END DO			!ISPEC
!                     
! We now need to compute the populations for the model atom with Super-levels.
! We do this in reverse order (i.e. highest ionization stage first) in order
! that we the ion density for the lower ionization stage is available for
! the next call.
!
! For 1st call to FULL_TO_SUP, Last line contains FeX etc as FeXI not installed.
!
	DO ID=NUM_IONS-1,1,-1
	   CALL FULL_TO_SUP(
	1      ATM(ID)%XzV,   ATM(ID)%NXzV,       ATM(ID)%DXzV,   ATM(ID)%XzV_PRES,
	1      ATM(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F, ATM(ID)%DXzV_F,
	1      ATM(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES,  ND)
	END DO
!
! Store all quantities in POPS array. This is done here as it enables POPION 
! to be readily computed. It also ensures that POS is correct if we don't
! iterate on T.
!
	DO ID=1,NUM_IONS-1
	    CALL IONTOPOP(POPS,  ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1          ATM(ID)%EQXzV, ATM(ID)%NXzV, NT,ND, ATM(ID)%XzV_PRES)
	END DO
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	DO J=1,ND
	  POPION(J)=0.0D0
	  DO I=1,NT
	     IF(Z_POP(I) .GT. 0.01D0)POPION(J)=POPION(J)+POPS(I,J)
	  END DO
	END DO
!
! Evaluates LTE populations for both the FULL atom, and super levels.
!
	CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
!
! 
! 
! Iterate on the initial temperature distribution so that the
! temperature distribution at depth corresponds to the GREY solution.
! We use the Rosseland mean opacities to evaluate the GREY temperature
! distribution. We then compute non-LTE partition functions which
!
! This page computes the Rosseland mean opacity from the temperature
! distribution and the population levels. TA is a working vector. The
! Rosseland opacity is given in ROSSMEAN. 
!
	CALL TUNE(1,'T_ITERATE')
	MAIN_COUNTER=1
	DO WHILE (ITERATE_INIT_T .AND. .NOT. GRID .AND.
	1                                 MAIN_COUNTER .LE. 5)
!
	    IF(.NOT. ALLOCATED(U_PAR_FN))THEN
	      ALLOCATE (U_PAR_FN(ND,NUM_IONS),STAT=IOS)
	      IF(IOS .EQ. 0)ALLOCATE (PHI_PAR_FN(ND,NUM_IONS),STAT=IOS)
	      IF(IOS .EQ. 0)ALLOCATE (Z_PAR_FN(NUM_IONS),STAT=IOS)
	      IF(IOS .NE. 0)THEN
	        WRITE(LUER,*)'Unable to allocate PHI_PAR_FN in SET_NEW_MODEL_ESTIMATES'
	        STOP
	      END IF
	    END IF
! 
!
! Set 2-photon data with current atomic models and populations.
!
	    DO ID=1,NUM_IONS-1
	       ID_SAV=ID
	       CALL SET_TWO_PHOT_V3(ION_ID(ID), ID_SAV, 
	1          ATM(ID)%XzVLTE,          ATM(ID)%NXzV,
	1          ATM(ID)%XzVLTE_F_ON_S,   ATM(ID)%XzVLEVNAME_F,
	1          ATM(ID)%EDGEXzV_F,       ATM(ID)%GXzV_F,
	1          ATM(ID)%F_TO_S_XzV,      ATM(ID)%NXzV_F, ND,
	1          ATM(ID)%ZXzV,            ATM(ID)%EQXzV,  ATM(ID)%XzV_PRES)
	    END DO
!
! 
!
! We ensure that LAST_LINE points to the first LINE that is going to
! be handled in the BLANKETING portion of the code.
!
	    LAST_LINE=0	    		!Updated as each line is done
	    DO WHILE(LAST_LINE .LT. N_LINE_FREQ .AND.
	1             VEC_TRANS_TYPE(LAST_LINE+1)(1:4) .NE. 'BLAN')
	            LAST_LINE=LAST_LINE+1
	    END DO
!
! ROSSMEAN is initially used to accumulate the integral of 1/chi (weighted 
! by dB/DT). After the frequency loop it is corrected so that it contains
! Rosseland mean opacity.
!
	    CALL DP_ZERO(ROSSMEAN,ND)
	    CALL DP_ZERO(PLANCKMEAN,ND)
	    TSTAR=T(ND)			!Required for IC in OPACITIES
	    CONT_FREQ=0.0D0
	    DO ML=1,NCF
	      FREQ_INDX=ML
	      FL=NU(ML)
!
	      IF(NU_EVAL_CONT(ML) .NE. CONT_FREQ)THEN
	        COMPUTE_NEW_CROSS=.TRUE.
	        CONT_FREQ=NU_EVAL_CONT(ML)
	      ELSE
	        COMPUTE_NEW_CROSS=.FALSE.
	      END IF
!
	      CALL COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
!
! Compute the line opacity. This initializes the storage locations on
! entry when FREQ_INDX=ML=1.
!
              CALL SET_LINE_OPAC(POPS,NU,FREQ_INDX,LAST_LINE,N_LINE_FREQ,
	1            LST_DEPTH_ONLY,LUER,ND,NT,NCF,MAX_SIM)

!
! Now add in line opacity to continuum opacity.
!
	      DO SIM_INDX=1,MAX_SIM
	        IF(RESONANCE_ZONE(SIM_INDX))THEN
	          DO I=1,ND
	            CHI(I)=CHI(I) +
	1             CHIL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	            ETA(I)=ETA(I) +
	1             ETAL_MAT(I,SIM_INDX)*LINE_PROF_SIM(I,SIM_INDX)
	          END DO
	        END IF
	      END DO
!
! CHECK for negative line opacities. 
!
	      DO I=1,ND
	        CHI_NOSCAT(I)=MAX(0.0D0,CHI(I)-ESEC(I))
	        IF(CHI(I) .LT. 0.1D0*ESEC(I))CHI(I)=0.1D0*ESEC(I)
	      END DO
!
! Note division by T**2 is included with Stefan-Boltzman constant.
!
	      T1=-HDKT*NU(ML)
	      T2=FQW(ML)*TWOHCSQ*(NU(ML)**3)
	      T3=-T1*FQW(ML)*TWOHCSQ*(NU(ML)**3)
	      DO I=1,ND
	        PLANCKMEAN(I)=PLANCKMEAN(I) + T2*CHI_NOSCAT(I)*EMHNUKT(I)/(1.0D0-EMHNUKT(I))
	        ROSSMEAN(I)=ROSSMEAN(I) + T3*EMHNUKT(I)/CHI(I)/(1.0D0-EMHNUKT(I))**2
	      END DO
	    END DO
!
! Compute CHI, and then optical depth scale.
! Stefan-Boltzman constant *1D-15*1D+16/PI (T**4/PI). NB --- T1 is a factor
! of 10^15 larger than in MAINGEN as FQW has already been multiplied by
! 10^15 for dv integrations.
!
! If clumping is important, we need to correct the Rosseland mean opacity 
! for clumping. Since it is a simple scale factor at each depth, we can do
! it here, rather than adjust CHI for each frequency.
!
	    T1=1.8047D+11
	    DO I=1,ND
	      ROSSMEAN(I)=4.0D0*CLUMP_FAC(I)*T1*(T(I)**5)/ROSSMEAN(I)
	      PLANCKMEAN(I)=CLUMP_FAC(I)*PLANCKMEAN(I)/T1/(T(I)**4)
	    END DO
!
	    CALL WRITV(ROSSMEAN,ND,'Rosseland Mean Opacity',88)
	    CALL WRITV(PLANCKMEAN,ND,'Planck Mean Opacity',88)
	    TA(1:ND)=1.0D-10*ROSSMEAN(1:ND)/DENSITY(1:ND)
	    TB(1:ND)=1.0D-10*PLANCKMEAN(1:ND)/DENSITY(1:ND)
	    CALL WRITV(TA,ND,'Rosseland mean mass absorption coefficient',88)
	    CALL WRITV(TB,ND,'Planck mean mass absorption coefficient',88)
! 
!
! Check that inner boundary is deep enough so that LTE can be fully recovered. SOURCE and
! TC are used as temporary vectors.
!
	    IF(MAIN_COUNTER .EQ. 1)THEN
	      CALL TORSCL(TA,ROSSMEAN,R,TB,TC,ND,METHOD,' ')
	      CALL ESOPAC(ESEC,ED,ND)
	      CALL TORSCL(TB,ESEC,R,SOURCE,TC,ND,METHOD,' ')
	      WRITE(LUER,*)' '
	      WRITE(LUER,'(A,ES10.3)')' Thompson scattering optical depth at inner boundary is:',TB(ND)
	      WRITE(LUER,'(A,ES10.3)')' Rosseland optical depth at inner boundary is:          ',TA(ND)
	      WRITE(LUER,'(A,ES10.3)')' Rosseland optical depth at outer boundary is:          ',TA(1)
	      WRITE(LUER,*)' '
	      IF(TA(ND) .LT. 10.0D0)THEN
	        WRITE(LUER,*)('*',I=1,70)
	        WRITE(LUER,*)('*',I=1,70)
	        WRITE(LUER,*)' '
	        WRITE(LUER,*)'Warning --- your core optical depth is probably too low'
	        WRITE(LUER,*)'You should use a value in excess of 10'
	        WRITE(LUER,*)' '
	        WRITE(LUER,*)('*',I=1,70)
	        WRITE(LUER,*)('*',I=1,70)
	      END IF
	    END IF
!
! Compute the grey temperature structure and the Rosseland optical depth scale.
! ROSSMEAN already includes the effect of clumping.
!
	    CHI(1:ND)=ROSSMEAN(1:ND)
	    WRITE(LUER,*)'Callng COMP_GREY_V4 in SET_NEW'
	    CALL COMP_GREY_V4(POPS,TGREY,TA,ROSSMEAN,PLANCKMEAN,COMPUTED,LUER,NC,ND,NP,NT)
	    IF(.NOT. COMPUTED)THEN
	      WRITE(LUER,*)'Unable to compute grey temperature structure'
	      WRITE(LUER,*)'As this is needed to provide initial T estimate, stopping code'
	      STOP
	    END IF
!
! SCALE_GREY modifies the computed grey temperature distribution according
! to that computed in a previous model.
!
! i.e. TGREY = TGREY . (T/TGREY)_old
!
	   IF(GREY_IOS .EQ. 0)THEN
	     CALL SCALE_GREY(TGREY,TA,GREY_IOS,LUIN,ND)
	   END IF
!
! Now correct T distribution towards grey value. As we don't require the 
! old T, we can overwrite it straight away. If we multiplied T1 by a number
! less than  unity, this would be equivalent to only a partial correction 
! of T towards TGREY:
!          GREY_PAR=0 set T=TGREY
!          GREY_PAR=INFINITY leaves T=T.
!
! T3 and T2 are used to determine the current largest correction.
!
	    TC(1:ND)=T(1:ND)
	    IF( MAIN_COUNTER .EQ. 1)THEN
	      DO I=1,ND
	        T_SAVE(I)=T(I)		!Save original T for use when
	      END DO                    !correcting T towards TGREY.
	    END IF
 	    T2=0.0D0
	    DO I=1,ND
	      IF(GREY_PAR .LE. 0)then
	        T1=1.0D0
	      ELSE
	        T1=1.0D0-EXP(-TA(I)/GREY_PAR)
	      END IF
	      IF(TA(I) .LT. 0.1D0*GREY_PAR)T1=0.0D0     !Changed T1 to TA(I) [14-Jan-2009]
	      T3=ABS( T1*(TGREY(I)-T(I)) )
	      T(I)=T1*TGREY(I)+(1.0D0-T1)*T_SAVE(I)
	      T(I)=MAX(T(I),0.95D0*T_MIN)
	      T2=MAX(T3/T(I),T2)
	    END DO
	    WRITE(LUER,'('' Largest correction to T in GREY initialization loop '//
	1              'is '',1P,E9.2,'' %'')')100.0*T2
!
! Now compute non-LTE partition functions. These assume that the
! departure coefficients are independent of Temperature. This
! is a good assumption at depth where b is approximately unity.
!
! GAM_SPECIES is used as a storage location for the population of the
! highest ionization stage. Must be done in forward direction.
!
! NB: After calling PAR_FUN_V2, ATM(ID)%XzV_F will contain DCs -
!       NOT populations.
!
! TMP_STRING is used to indicate whether we interpolate in depature
! coefficients (DC) or excitation temperatures (TX).
!
	    TMP_STRING='DC'
	    IF(DC_INTERP_METHOD .EQ. 'RTX')TMP_STRING='TX'
	    DO ID=1,NUM_IONS
	      J=ID-1			!1 is added in PAR_FUN_V2
	      ISPEC=SPECIES_LNK(ID)
	      CALL PAR_FUN_V4(U_PAR_FN, PHI_PAR_FN, Z_PAR_FN,
	1          GAM_SPECIES(1,ISPEC),
	1          ATM(ID)%XzV_F,     ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%W_XzV_F,
	1          ATM(ID)%DXzV_F,    ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,
	1          ATM(ID)%GIONXzV_F, ATM(ID)%ZXzV,T, TC, ED,
	1          ATM(ID)%NXzV_F, ND,J,NUM_IONS, 
	1          ATM(ID)%XzV_PRES,ION_ID(ID),TMP_STRING)
	   END DO
!
! The non-LTE partition functions are density independent, provided
! we assume the departure coefficients remain fixed.
!
! We now evaluate the contribution to the electron density by each
! species, using the non-LTE partition functions.
!
! We use H for ED(est)
! We use QH for dED(est)/dT.
!
	    T1=1.0D0
	    J=0
	    DO WHILE (T1 .GT. 1.0D-04)
	      FIRST=.TRUE.
!
! Recall GAM_SPECIES is set to be the population of the highest ionization
! stage.
!
	      DO ISPEC=1,NUM_SPECIES
	        ID=SPECIES_BEG_ID(ISPEC)
	        J=SPECIES_END_ID(ISPEC)-SPECIES_BEG_ID(ISPEC)+1
	        IF(SPECIES_PRES(ISPEC))THEN
	          CALL EVAL_ED(H,QH,U_PAR_FN(1,ID),PHI_PAR_FN(1,ID),
	1                  Z_PAR_FN(ID),ED,POP_SPECIES(1,ISPEC),
	1                  GAM_SPECIES(1,ISPEC),XM,TB,TC,J,ND,FIRST)
	        END IF
	      END DO
!
	      T1=0.0
	      DO I=1,ND
	        TA(I)=-(H(I)-ED(I))/(QH(I)-1.0D0)/ED(I)
	        T1=MAX(T1,ABS(TA(I)))
	        IF(TA(I) .LT. -0.9D0)TA(I)=-0.9D0
	        IF(TA(I) .GT. 9.0D0)TA(I)=9.0D0
	        ED(I)=ED(I)*(1.0D0+TA(I))
	      END DO
	      J=J+1
	      IF(J .GT. 20)THEN
	        WRITE(LUER,*)'Error --- Computation of ED in EVAL_ED section'//
	1                 ' has taken more than 20 iterations'
	        STOP
	      END IF
	    END DO
! 
!
! Now need to compute LTE populations, and populations.
! Since T and Ne have altered, we revise the vectors for evaluating the 
! level dissolution. These constants are the same for all species. These are 
! stored in a common block, and are required by SUP_TO_FULL and LTE_POP_WLD.
!
! NB: POPION will also alter but in W-R and LBV's all species will be ionized,
! and hence POPION will not change from iteration to iteration. In any event,
! it has a smaller effect than changes in Ne.
!
	    CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! We do low ionization species second, as first need DION.
!
	    DO ISPEC=1,NUM_SPECIES
	      FIRST=.TRUE.   
	      DO ID=SPECIES_END_ID(ISPEC),SPECIES_BEG_ID(ISPEC),-1
	        IF(ATM(ID)%XzV_PRES)THEN
	          CALL LTEPOP_WLD_V2(ATM(ID)%XzVLTE_F, ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%W_XzV_F,
	1               ATM(ID)%EDGEXzV_F,  ATM(ID)%GXzV_F,  ATM(ID)%ZXzV,
	1               ATM(ID)%GIONXzV_F,  ATM(ID)%NXzV_F,  ATM(ID)%DXzV_F,
	1               ED,T,ND)
	          CALL CNVT_FR_DC_V2(ATM(ID)%XzV_F, ATM(ID)%LOG_XzVLTE_F,
	1               ATM(ID)%DXzV_F,   ATM(ID)%NXzV_F,
	1               TB,               TA,ND,FIRST,      ATM(ID+1)%XzV_PRES)
	          IF(ID .NE. SPECIES_BEG_ID(ISPEC))ATM(ID-1)%DXzV_F(1:ND)=TB(1:ND)
	        END IF
	      END DO
!
! We need to scale the populations to ensure that the change in temperature
!   has not causes some population to blow up. We always do this --- the
! DO_POP_SCALE option has no effect.
!
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        CALL SCALE_POPS(ATM(ID)%XzV_F, ATM(ID)%DXzV_F,
	1           POP_SPECIES(1,SPECIES_LNK(ID)),TA, ATM(ID)%NXzV_F,ND)
	      END DO
	    END DO
!
	    MAIN_COUNTER=MAIN_COUNTER+1
!
! We now need to compute the populations for the model atom with Super-levels.
! We do this in reverse order (i.e. highest ionization stage first) in order
! that we the ion density for the lower ionization stage is available for
! the next call.
!
! For 1st call to FULL_TO_SUP, Last line contains FeX etc as FeXI not installed.
!
	    DO ID=NUM_IONS-1,1,-1
	      CALL FULL_TO_SUP(
	1      ATM(ID)%XzV,   ATM(ID)%NXzV,       ATM(ID)%DXzV,      ATM(ID)%XzV_PRES,
	1      ATM(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F,    ATM(ID)%DXzV_F,
	1      ATM(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, ND)
	    END DO
!
! Store all quantities in POPS array. This is done here (rather than
! after final iteration) as it enable POPION to be readily computed.
!
	    DO ID=1,NUM_IONS-1
	      CALL IONTOPOP(POPS, ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1         ATM(ID)%EQXzV, ATM(ID)%NXzV, NT,ND,
	1         ATM(ID)%XzV_PRES)
	    END DO
!
! Compute the ion population at each depth.
! These are required when evaluation the occupation probabilities.
!
	    DO J=1,ND
	      POPION(J)=0.0D0
	      DO I=1,NT        
	        IF(Z_POP(I) .GT. 0.01D0)POPION(J)=POPION(J)+POPS(I,J)
	      END DO
	    END DO
!
! While the following may seem superfolous, it ensures absolute consistency. Its possible
! that, for H, He etc that the upper levels with interpolating sequences may not be fully
! consistent (due to rounding errors, from another model, etc).
!
	    CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!
! Revise ALL LTE populations.
!
	    CALL EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
!
	END DO		!ITERATE_INIT_T
	CALL TUNE(2,'T_ITERATE')
!
	IF(ALLOCATED(U_PAR_FN))THEN
	  DEALLOCATE (U_PAR_FN,STAT=IOS)
	  DEALLOCATE (PHI_PAR_FN,STAT=IOS)
	  DEALLOCATE (Z_PAR_FN,STAT=IOS)
	END IF 
!
! Restore two photon method option.
!
	TWO_PHOTON_METHOD=SAVED_TWO_PHOTON_METHOD
!
	RETURN
	END
