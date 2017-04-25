!
! Program to compute the emergent flux in the observer's frame using an
! Observer's frame flux calculation. The program first computes the
! emissivities, opacities, and radiation field in the CMF frame. COHERENT
! or INCOHERENT electron scattering can be assumed. This data is then passed
! to OBS_FRAME_SUB for computation of the OBSERVER's flux.
!
! CMF_FLUX calling program is partially based on DISPGEN.
!
	PROGRAM LTE
	USE MOD_CMFGEN
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Altered:  1-Nov-2010: ND and NP are now computed internally from the data
!                         in the GRID_PARAMS file. This no longer any need to
!                         alter MODEL_SPEC.
! Altered: 20-Mar-2005: FLUX_MEAN & ROSS_MEAN now zeroed. These vectors now
!                         used in CMFGEN_SUB.
! Altered: 03-Mar-2000: Variable type ATM installed to simplify handling
!	                   of multiple species.
!
! Created:  5-Jan-1998=9 (Progran began late Dec, 1998)
!
	INTEGER ND		!Actual number of depth points in atmosphere
	INTEGER NC		!Actual number of core rays 
	INTEGER NP		!Total number of rays (ND+NC)
	INTEGER NUM_BNDS	!Number of bans in linearization matrix
	INTEGER MAX_SIM	!Maximum number of lines that can be treated sim.
	INTEGER NCF_MAX	!Maximum number of frequencies that can be treated.
!
	INTEGER ND_MAX,NP_MAX
	INTEGER N_LINE_MAX
!
	INTEGER NM
	INTEGER NLF
	INTEGER NM_KI
	INTEGER TX_OFFSET
	INTEGER NION
	INTEGER DIAG_INDX
!
	CHARACTEr*20 TEMP_KEY
!
	REAL*8 T1		!Temporary variable
	INTEGER I,J,IOS,NT
	INTEGER EQ_TEMP
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: LU_IN=7
	INTEGER, PARAMETER :: LU_OUT=8
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	INTEGER NF
	INTEGER NS
	INTEGER NV
	INTEGER ID
	INTEGER ISPEC
	INTEGER NUM_IONS_RD
!
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	LOGICAL AT_LEAST_ONE_ION_PRES
	LOGICAL FND_END_OF_IONS
	LOGICAL DO_TERM_OUT
!
! Set constants.
!
	CHIBF=2.815D-06
	CHIFF=3.69D-29
	HDKT=4.7994145D0
	TWOHCSQ=0.0147452575D0
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
	LUER=ERROR_LU()
!
! Open output file for all errors and comments. Change DO_TERM_OUT to
! have the output go to the terminal/batch log file.
!
	DO_TERM_OUT=.FALSE.
	IF(.NOT. DO_TERM_OUT)THEN
	  CALL GEN_ASCI_OPEN(LUER,'OUTLTE','UNKNOWN','APPEND',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error opening OUTLTE in LTE, IOS=',IOS
	    STOP
	  END IF
	  CALL SET_LINE_BUFFERING(LUER)
	END IF
!
! Set all atomic data. New species can be simple added by insertion.
! Try to add species in order of atomic number. Hydrogen should ALWAYS
! be species 1, Helium should ALWAYS be species 2. 
!
! While this tabulation is Verbose, it is simple to change.
! Note that the Solar abundances are only used for reference in
! the output table MOD_SUM.
!
	ID=1
	AT_NO(1)=1.0D0;		    AT_MASS(ID)=1.0D0
	SPECIES(ID)='HYD';	    SPECIES_ABR(ID)='H'
	SOL_ABUND_HSCL(ID)=12.0D0
!
	ID=ID+1
	AT_NO(ID)=2.0D0; 	    AT_MASS(ID)=4.0D0		!Helium
	SPECIES(ID)='HE';	    SPECIES_ABR(ID)='He'
	SOL_ABUND_HSCL(ID)=11.0D0
!
	ID=ID+1
	AT_NO(ID)=6.0D0;	    AT_MASS(ID)=12.0D0		!Carbon
	SPECIES(ID)='CARB';	    SPECIES_ABR(ID)='C'
	SOL_ABUND_HSCL(ID)=8.56D0
!
	ID=ID+1
	AT_NO(ID)=7.0D0;	    AT_MASS(ID)=14.0D0		!Nitrogen
	SPECIES(ID)='NIT';	    SPECIES_ABR(ID)='N'
	SOL_ABUND_HSCL(ID)=8.05D0
!
	ID=ID+1
	AT_NO(ID)=8.0D0; 	    AT_MASS(ID)=16.0D0		!Oxygen
	SPECIES(ID)='OXY';	    SPECIES_ABR(ID)='O'
	SOL_ABUND_HSCL(ID)=8.93D0
!
	ID=ID+1
	AT_NO(ID)=9.0D0;            AT_MASS(ID)=19.00D0         !Fluorine
	SPECIES(ID)='FLU';          SPECIES_ABR(ID)='F'
	SOL_ABUND_HSCL(ID)=4.56D0
!
	ID=ID+1
	AT_NO(ID)=10.0D0;	    AT_MASS(ID)=20.2D0		!Neon
	SPECIES(ID)='NEON';	    SPECIES_ABR(ID)='Ne'
	SOL_ABUND_HSCL(ID)=8.09D0
!
	ID=ID+1
	AT_NO(ID)=11.0D0;	    AT_MASS(ID)=23.0D0		!Sodium
	SPECIES(ID)='SOD';	    SPECIES_ABR(ID)='Na'
	SOL_ABUND_HSCL(ID)=6.33D0
!
	ID=ID+1
	AT_NO(ID)=12.0D0;	    AT_MASS(ID)=24.3D0		!Magnesium
	SPECIES(ID)='MAG';	    SPECIES_ABR(ID)='Mg'
	SOL_ABUND_HSCL(ID)=7.58D0
!
	ID=ID+1
	AT_NO(ID)=13.0D0;	    AT_MASS(ID)=27.0D0		!Aluminium
	SPECIES(ID)='ALUM';	    SPECIES_ABR(ID)='Al'
	SOL_ABUND_HSCL(ID)=6.47D0
!
	ID=ID+1
	AT_NO(ID)=14.0D0;	    AT_MASS(ID)=28.1D0		!Silicon
	SPECIES(ID)='SIL';	    SPECIES_ABR(ID)='Sk'
	SOL_ABUND_HSCL(ID)=7.55D0
!
	ID=ID+1
	AT_NO(ID)=15.0D0;	    AT_MASS(ID)=31.0D0		!Phosphorous
	SPECIES(ID)='PHOS';	    SPECIES_ABR(ID)='P'
	SOL_ABUND_HSCL(ID)=5.45D0
!
	ID=ID+1
	AT_NO(ID)=16.0D0;	    AT_MASS(ID)=32.1D0		!Sulpher
	SPECIES(ID)='SUL';	    SPECIES_ABR(ID)='S'
	SOL_ABUND_HSCL(ID)=7.21D0
!
	ID=ID+1
	AT_NO(ID)=17.0D0;	    AT_MASS(ID)=35.5D0		!Chlorine
	SPECIES(ID)='CHL';	    SPECIES_ABR(ID)='Cl'
	SOL_ABUND_HSCL(ID)=5.5D0
!
	ID=ID+1
	AT_NO(ID)=18.0D0;	    AT_MASS(ID)=39.9D0		!Argon
	SPECIES(ID)='ARG';	    SPECIES_ABR(ID)='Ar'
	SOL_ABUND_HSCL(ID)=6.56D0
!
	ID=ID+1
	AT_NO(ID)=19.0D0;	    AT_MASS(ID)=39.1D0		!Potassium
	SPECIES(ID)='POT';	    SPECIES_ABR(ID)='K'
	SOL_ABUND_HSCL(ID)=5.12D0
!
	ID=ID+1
	AT_NO(ID)=20.0D0;	    AT_MASS(ID)=40.1D0		!Calcium
	SPECIES(ID)='CAL';	    SPECIES_ABR(ID)='Ca'
	SOL_ABUND_HSCL(ID)=6.36D0
!
	ID=ID+1
	AT_NO(ID)=21.0D0;           AT_MASS(ID)=44.96D0         !Scandium
	SPECIES(ID)='SCAN';         SPECIES_ABR(ID)='Sc'
	SOL_ABUND_HSCL(ID)=3.10D0
!
	ID=ID+1
	AT_NO(ID)=22.0D0;	    AT_MASS(ID)=47.88D0		!Titanium
	SPECIES(ID)='TIT';	    SPECIES_ABR(ID)='Tk'	!Actual symbol is Ti
	SOL_ABUND_HSCL(ID)=4.99D0
!
	ID=ID+1
	AT_NO(ID)=23.0D0;           AT_MASS(ID)=50.94D0         !Vandium
	SPECIES(ID)='VAN';          SPECIES_ABR(ID)='V'         !Actual symbol is V 
	SOL_ABUND_HSCL(ID)=4.00D0
!
	ID=ID+1
	AT_NO(ID)=24.0D0;	    AT_MASS(ID)=52.0D0		!Chromium
	SPECIES(ID)='CHRO';	    SPECIES_ABR(ID)='Cr'
	SOL_ABUND_HSCL(ID)=5.67D0
!
	ID=ID+1
	AT_NO(ID)=25.0D0;	    AT_MASS(ID)=54.9D0		!Maganese
	SPECIES(ID)='MAN';	    SPECIES_ABR(ID)='Mn'
	SOL_ABUND_HSCL(ID)=5.39D0
!
	ID=ID+1
	AT_NO(ID)=26.0D0;	    AT_MASS(ID)=55.8D0		!Iron
	SPECIES(ID)='IRON';	    SPECIES_ABR(ID)='Fe'
	SOL_ABUND_HSCL(ID)=7.54D0        
!
	ID=ID+1
	AT_NO(ID)=27.0D0;	    AT_MASS(ID)=58.9D0		!Cobalt
	SPECIES(ID)='COB';	    SPECIES_ABR(ID)='Co'
	SOL_ABUND_HSCL(ID)=4.92D0
!
	ID=ID+1
	AT_NO(ID)=28.0D0;	    AT_MASS(ID)=58.7D0		!Nickel
	SPECIES(ID)='NICK';	    SPECIES_ABR(ID)='Nk'
	SOL_ABUND_HSCL(ID)=6.25D0
!
	ID=ID+1
	AT_NO(ID)=56.0D0;           AT_MASS(ID)=137.33D0        !Barium
	SPECIES(ID)='BAR';          SPECIES_ABR(ID)='Ba'
	SOL_ABUND_HSCL(ID)=2.13D0
!
	IF(ID .NE. NUM_SPECIES)THEN
	  WRITE(LUER,*)'Error in CMFGEN: Invalid species setup'
	  WRITE(LUER,*)'This likely means a new species has not been added to necesseary files'
	  STOP
	END IF
!
! Convert from the abundance on a logarithmic scale with H=12.0 dex, to
! mass-fractions.
!
	T1=0.0D0
	DO ID=1,NUM_SPECIES
	  SOL_MASS_FRAC(ID)=10.0D0**(SOL_ABUND_HSCL(ID)-12.0D0)
	  SOL_MASS_FRAC(ID)=AT_MASS(ID)*SOL_MASS_FRAC(ID)
	  T1=T1+SOL_MASS_FRAC(ID)
	END DO
	SOL_MASS_FRAC(:)=SOL_MASS_FRAC(:)/T1
!
! Initilaization: These parameters will remain as they are when a species
! is not present.
!
	EQ_SPECIES(:)=0
	SPECIES_PRES(:)=.FALSE.
	AT_ABUND(:)=0.0D0
	SPECIES_BEG_ID(:)=0
	SPECIES_END_ID(:)=-1
!
! Although these will be set later, we set them here to allow us to check
! that the included ionization stages are sequential (for each species).
!
	DO ID=1,MAX_NUM_IONS
	  ATM(ID)%NXzV_F=1; 	ATM(ID)%NXzV=1
	  ATM(ID)%XzV_PRES=.FALSE.
	END DO
!
! Get size of the electron and temperature grid. This will be used to set ND and
! NP. We arbitrarily set NC to 10.
!
	OPEN(UNIT=LU_IN,FILE='GRID_PARAMS',STATUS='OLD',ACTION='READ')
          READ(LU_IN,*)NC,ND
        CLOSE(LU_IN)
	ND=NC*ND
	NC=10
	NP=ND+NC
	NUM_BNDS=1
!
! Get data describing number of depth points, number of atomic levels
! etc.
!
! There is no need to read NC, ND, and NP, as these are ast from GRID_PARAMS.
!
	OPEN(UNIT=LU_IN,FILE='MODEL_SPEC',STATUS='OLD',ACTION='READ')
	CALL RD_OPTIONS_INTO_STORE(LU_IN,LU_OUT)
!
!	CALL RD_STORE_INT(ND,'ND',L_TRUE,'Number of depth points')
!	CALL RD_STORE_INT(NC,'NC',L_TRUE,'Number of core rays')
!	CALL RD_STORE_INT(NP,'NP',L_TRUE,'Number of impact parameters')
!	CALL RD_STORE_INT(NUM_BNDS,'NUM_BNDS',L_TRUE,'Number of bands in linearization matrix (BA)')
!
	CALL RD_STORE_INT(MAX_SIM,'MAX_SIM',L_TRUE,
	1        'Maximum # of lines that can treated simultaneously')
	CALL RD_STORE_INT(N_LINE_MAX,'NLINE_MAX',L_TRUE,
	1        'Maximum # of lines that can be treated')
	CALL RD_STORE_INT(NCF_MAX,'NCF_MAX',L_TRUE,
	1        'Maximum # of frequencies that can be treated')
	CALL RD_STORE_INT(NLF,'NLF',L_TRUE,
	1        'Number of frequencies per Doppler profile in CMF mode (21)')
!
! We now get the number of atomic levels. Old MODEL_SPEC files, with NSF in the
! keyword can be read. ISF take's precident over NSF, and there is no check that 
! there is not both an effectively identical NSF and ISF keyowrd. In this
! case, the number of important variables is assumed to be the same as NS.
!
	ID=0
	NUM_IONS_RD=0
	EQ_TEMP=1
	DO ISPEC=1,NUM_SPECIES
	  AT_LEAST_ONE_ION_PRES=.FALSE.
	  FND_END_OF_IONS=.FALSE.
	  DO J=1,MIN(MAX_IONS_PER_SPECIES,NINT(AT_NO(ISPEC))+1)
	    TEMP_KEY=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))//'_ISF'
	    NS=0
	    CALL RD_STORE_3INT(NV,NS,NF,TEMP_KEY,L_FALSE,' ')
	    IF(NS .EQ .0)THEN
	      TEMP_KEY=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))//'_NSF'
	      CALL RD_STORE_2INT(NS,NF,TEMP_KEY,L_FALSE,' ')
	      NV=NS
	    END IF
	    IF(NS .NE. 0)THEN
	      ID=ID+1
	      ION_ID(ID)=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))
	      ATM(ID)%XzV_PRES=.TRUE.
	      ATM(ID)%NXzV_F=NF
	      ATM(ID)%NXzV=NS
	      ATM(ID)%NXzV_IV=NV
	      AT_LEAST_ONE_ION_PRES=.TRUE.
	      IF(SPECIES_BEG_ID(ISPEC) .EQ. 0)SPECIES_BEG_ID(ISPEC)=ID
	      ATM(ID)%EQXzV=EQ_TEMP
	      EQ_TEMP=EQ_TEMP+ATM(ID)%NXzV
	      SPECIES_LNK(ID)=ISPEC
	      IF(FND_END_OF_IONS)THEN
	        WRITE(LUER,*)'Error in CMFGEN'
	        WRITE(LUER,*)'Ionization stages for species ',
	1                      TRIM(SPECIES_ABR(ISPEC)),' are not consecutive'
	        WRITE(LUER,*)'Check your MODEL_SPEC file for typo''s'
	        STOP
	      END IF
	    ELSE IF(AT_LEAST_ONE_ION_PRES)THEN
	      FND_END_OF_IONS=.TRUE.
	    END IF
	  END DO
!
! We now set the variables for the SINGLE level of the HIGHEST ionization stage.
! XzV_PRES must be false.
!
	  IF(AT_LEAST_ONE_ION_PRES)THEN
	    ID=ID+1
            SPECIES_PRES(ISPEC)=.TRUE.
	    SPECIES_END_ID(ISPEC)=ID
	    SPECIES_LNK(ID)=ISPEC
	    ATM(ID)%XzV_PRES=.FALSE.
	    ATM(ID)%EQXzV=EQ_TEMP
	    EQ_SPECIES(ISPEC)=EQ_TEMP
	    EQ_TEMP=EQ_TEMP+1
	    NUM_IONS_RD=NUM_IONS_RD+1		!i.e. number of species
	  END IF
	END DO
	NUM_IONS=ID
	NUM_IONS_RD=NUM_IONS-NUM_IONS_RD	!
	EQNE=EQ_TEMP
	NT=EQNE+1
!
! Check whether any ION type'os in MODEL_SPEC file --- i.e. check that all ions
! specified by [XzV_NSF] are valid.
!
	CALL CNT_NUM_KEYS(I,'_ISF]')
	IF(I .NE. NUM_IONS_RD)THEN
	  WRITE(LUER,*)'Error in CMFGEN'
	  WRITE(LUER,*)'You have an invalid [XzV_ISF] key in MODEL_SPEC'
	  WRITE(LUER,*)'NUM_IONS=',NUM_IONS
	  WRITE(LUER,*)'NUM_IONS_RD=',NUM_IONS_RD
	  WRITE(LUER,*)'No IONS in MODEL_SPEC',I
	  STOP
	END IF

!
	CALL CLEAN_RD_STORE()
	CLOSE(LU_IN)
!
!
	IF(NP .NE. ND+NC .AND. NP .NE. ND+NC-2)THEN
	  WRITE(LUER,*)'Error in CMFGEN'
	  WRITE(LUER,*)'Invalid NP: should be ND+NC or ND+NC-2'
	  STOP
	END IF
!
! Check whether EQUATION LABELLING is consistent. ' I ' is used as the
! number of the current equation. We also set the variable SPEC_PRES which 
! indicates whether at least one ionization stage of a species is present.
! It is used to determine, for example,  whether a number conservation 
! equation is required.
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
	IF(EQNE .NE. I)THEN
	  WRITE(LUER,*)'Error - EQNE has wrong value in CMFGEN'
	  STOP
	END IF
	IF(NT .NE. I+1)THEN
	  WRITE(LUER,*)'Error - NT has wrong value in CMFGEN'
	  STOP
	END IF
!
!
	ALLOCATE (R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SIGMA(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (T(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ED(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ROSS_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (FLUX_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POP_ATOM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (DENSITY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POPION(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CLUMP_FAC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (VTURB_VEC(ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX'
	  WRITE(LUER,*)'Unable to allocate Atmosphere arrays'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
	FLUX_MEAN(:)=0.0D0
	ROSS_MEAN(:)=0.0D0
!
! Allocate population vectors.
!
	ALLOCATE (POP_SPECIES(ND,NUM_SPECIES)); POP_SPECIES=0.0D0
	ALLOCATE (GAM_SPECIES(ND,NUM_SPECIES)); GAM_SPECIES=0.0D0
!
! Allocate ATM memory
!
	DO ID=1,NUM_IONS
!
	  NF=ATM(ID)%NXzV_F
	  NS=ATM(ID)%NXzV
	                ALLOCATE (ATM(ID)%XzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE(NS,ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%LOG_XzVLTE(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%dlnXzVLTE_dlnT(NS,ND),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzV_F(NF,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F(NF,ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%LOG_XzVLTE_F(NF,ND),STAT=IOS)
          IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F_ON_S(NF,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%W_XzV_F(NF,ND),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%AXzV_F(NF,NF),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%EDGEXzV_F(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%F_TO_S_XzV(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%INT_SEQ_XzV(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLEVNAME_F(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GXzV_F(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%ARAD(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM2(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM4(NF),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%OBSERVED_LEVEL(NF),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV_F(ND),STAT=IOS)       ; ATM(ID)%DXzV_F(:)=0.0D0
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV(ND),STAT=IOS)         ; ATM(ID)%DXzV(:)=0.0D0
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WSXzV(NS,ND,NPHOT_MAX),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WCRXzV(NS,ND,NPHOT_MAX),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%dWSXzVdT(NS,ND,NPHOT_MAX),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WSE_X_XzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%WCR_X_XzV(NS,ND),STAT=IOS)
!
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%APRXzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%ARRXzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%BFCRXzV(NS,ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%FFXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CPRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CRRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CHG_PRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%CHG_RRXzV(ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%COOLXzV(ND),STAT=IOS)
!
	  IF(IOS .NE. 0)THEN
	    WRITE(LUER,*)'Error in CMF_FLUX'
	    WRITE(LUER,*)'Unable to allocate arrays for species XzV'
	    WRITE(LUER,*)'STATUS=',IOS
	    STOP
	  END IF
	END DO
!
	I=0
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    I=I+1
	    ATM(ID)%INDX_XzV=I
	    ATM(ID)%EQXzVBAL=I
	  ELSE
	    ATM(ID)%INDX_XzV=0
	    ATM(ID)%EQXzVBAL=0
	  END IF
	END DO
!
! Create a vector, IMP_VAR, to inidicate which variables are considered as
! important for the atmospheric structure.
!
	ALLOCATE (IMP_VAR(NT))
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	     CALL SET_IMP_VEC(IMP_VAR,ATM(ID)%NXzV,ATM(ID)%NXzV_IV,
	1                        ATM(ID)%EQXzV,NT,ATM(ID)%XzV_PRES)
	  END IF
	END DO
!
! Used when computing dCHI and dETA for CMF mode. Required dimension is 4.
!
	NM_KI=4
!
! Indicates what storage locations in dzNET can be used for lines.
! The first 2 locations are for the current total opacity and emissivity.
! The next 3 locations are for the continuum opacity, emissivity, and
! electron scattering opacity.
! levels can be used for 
!
	TX_OFFSET=5
!
! Total number of storage locations to be set aside. We multiply MAX_SIM by
! 2 to account for both the upper and lower levels.
!
	NM=TX_OFFSET+2*MAX_SIM
!
! NION and NUM_IONS can be used interchangably. NION is passed to CMFGEN_SUB
! for dynamic array allocation.
!
	NION=NUM_IONS
	DIAG_INDX=NUM_BNDS/2+1
!
! Set arrays asside for temoray storage, and in an ACCURATE J calculation
! is to be performed.
!
	ND_MAX=MAX(NT,2*ND)
	NP_MAX=ND_MAX+2*NC
!
	CALL LTE_SUB(ND,NC,NP,NT,
	1            NUM_BNDS,NION,DIAG_INDX,
	1            ND_MAX,NP_MAX,NCF_MAX,N_LINE_MAX,
	1            TX_OFFSET,MAX_SIM,NM,NM_KI,NLF)
!
	CALL TUNE(3,' ')
	STOP
!
	END
