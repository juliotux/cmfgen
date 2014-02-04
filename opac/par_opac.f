!
! Program to compute the emergent flux in the observer's frame using an
! Observer's frame flux calculation. The program first computes the
! emissivities, opacities, and radiation field in the CMF frame. COHERENT
! or INCOHERENT electron scattering can be assumed. This data is then passed
! to OBS_FRAME_SUB for computation of the OBSERVER's flux.
!
! CMF_FLUX calling program is partially based on DISPGEN.
!
	PROGRAM PAR_OPAC
	USE MOD_PAR_OPAC
	USE MOD_USR_OPTION
	IMPLICIT NONE
!
! Altered: 03-Mar-2000: Variable type ATM installed to simplify handling
!	                   of multiple species.
!
! Created:  5-Jan-1998=9 (Progran began late Dec, 1998)
!
	INTEGER*4 ND		!Actula number of depth points in atmosphere
	INTEGER*4 NC		!Actual number of core rays 
	INTEGER*4 NP		!Total number of rays (ND+NC)
!
	INTEGER*4 ND_MAX,NP_MAX,NC_MAX
	INTEGER*4 N_MAX
	INTEGER*4 N_LINE_MAX
!
	REAL*8    GF_CUT
	INTEGER*4 GF_LEV_CUT
	INTEGER*4 MIN_NUM_TRANS
	LOGICAL   ONLY_OBS_LINES
	CHARACTER*10 GF_ACTION
!
	CHARACTER*20 TIME
	CHARACTER*80 FILENAME,BLANK,DIR_NAME,STRING
	CHARACTER*20 FILE_EXTENT
	CHARACTER*11 FORMAT_DATE
	CHARACTER*10 NAME_CONVENTION
!
	LOGICAL ASK  		!Ask of filenames or uset defaults.
	LOGICAL FILE_PRES
	LOGICAL SCRAT
	LOGICAL XRAYS
	LOGICAL DIE_AS_LINE
	INTEGER*4 I,J,IOS,NT,SCRATREC
	INTEGER*4 LEN_DIR
	INTEGER*4 EQ_TEMP
	REAL*8 T1,T2
	REAL*8 RMDOT
	REAL*8 VSM_DIE_KMS
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER*4, PARAMETER :: IZERO=0
	INTEGER*4, PARAMETER :: T_IN=5		!Terminal IO
	INTEGER*4, PARAMETER :: T_OUT=6		!Terminal IO
	INTEGER*4, PARAMETER :: LUIN=8
	INTEGER*4, PARAMETER :: LUMOD=9
	INTEGER*4, PARAMETER :: LU_TMP=10
	INTEGER*4, PARAMETER :: LU=30			!Used for pop file io
	INTEGER*4, PARAMETER :: LUSCRAT=33
!
	INTEGER*4 NF
	INTEGER*4 NS
	INTEGER*4 ID
	INTEGER*4 ISPEC
	LOGICAL PREV_STAGE_PRES
!
	INTEGER*4 LUER
	INTEGER*4 ERROR_LU
	EXTERNAL ERROR_LU
!
!	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
!
	DATA BLANK/' '/
!
! Set constants.
!
	CHIBF=2.815E-06
	CHIFF=3.69E-29
	HDKT=4.7994145
	TWOHCSQ=0.0147452575
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
	OPLIN=2.6540081E+08
	EMLIN=5.27296E-03
	LUER=ERROR_LU()
	FILE_EXTENT=' '
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
	AT_NO(ID)=22.0D0;	    AT_MASS(ID)=47.88D0		!Titanium
	SPECIES(ID)='TIT';	    SPECIES_ABR(ID)='Tk'	!Actual symbol is Ti
	SOL_ABUND_HSCL(ID)=4.99D0
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
	IF(ID .NE. NUM_SPECIES)THEN
	  WRITE(LUER,*)'Error in CMFGEN: Invalid species setup'
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
	ABUND_SCALE_FAC(:)=1.0D0
	DO ID=1,MAX_NUM_IONS
	  ATM(ID)%NXzV_F=1; 	ATM(ID)%NXzV=0
	  ATM(ID)%XzV_PRES=.FALSE.
	END DO
!
! Read in the gaunt factors for individual l states of hydrogen.
!
	CALL RD_HYD_BF_DATA(LUIN,LUMOD,T_OUT)
!
! *************************************************************************
!
! Read in basic model. This includes scaler quantities which describe the
! model (ie MDOT, ND, NC) as well as important vectors [i.e. Density
! structure, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
10	FILENAME=' '
	DIR_NAME=' '
	ASK=.FALSE.
	IOS=2			!Filename has to exits, blank not allowed.
	WRITE(T_OUT,*)'Append (ask) to so that subsequent FILE names '//
	1         'are not defaulted'
	WRITE(T_OUT,*)'Append (scrat) to get scratch output.'
	WRITE(T_OUT,*)' '
	I=80
	CALL RD_NCHAR(FILENAME,'RVTJ',I,T_IN,LUER,'Input RVTJ filename')
	CALL SET_CASE_UP(FILENAME,0,0)
	STRING=FILENAME
	IF( INDEX(STRING,'(SCRAT)') .NE. 0)SCRAT=.TRUE.
	IF( INDEX(STRING,'(ASK)') .NE. 0)ASK=.TRUE.
	I= INDEX(STRING,'(')
	IF(I .NE. 0)FILENAME=FILENAME(1:I-1)//' '
	INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	IF(.NOT. FILE_PRES)GOTO 10
!
! Get default directory.
!
	DIR_NAME=' '		!Valid DIR_NAME if not present.
	LEN_DIR=0
	J=LEN(FILENAME)
	DO WHILE(J .GT. 0)
	  IF( FILENAME(J:J) .EQ. ']' .OR.
	1     FILENAME(J:J) .EQ. ':' .OR.
	1     FILENAME(J:J) .EQ. '/'        )THEN
	    DIR_NAME=FILENAME(1:J)
	    LEN_DIR=J
	    J=0
	  END IF
	  J=J-1
	END DO
!
! Get default extension.
!
	FILE_EXTENT=' '
	J=LEN(FILENAME)
	DO WHILE(J .GT. LEN_DIR)
	  IF( FILENAME(J:J) .EQ. '.' )THEN
	    FILE_EXTENT=FILENAME(J:)
	    J=0
	  END IF
	  J=J-1
	END DO
!
	CALL RD_RVTJ_PARAMS_V2(RMDOT,STARS_LUM,AT_ABUND(1),TIME,NAME_CONVENTION,
	1                    ND,NC,NP,FILENAME,LUIN)
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
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in CMF_FLUX'
	  WRITE(LUER,*)'Unable to allocate Atmosphere arrays'
	  WRITE(LUER,*)'STATUS=',IOS
	  STOP
	END IF
	CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,DENSITY,CLUMP_FAC,ND,LUIN)
!
	CLOSE(UNIT=LUIN)
	WRITE(T_OUT,3)TIME
3	FORMAT(1X,'Model completed on ',A20)
!
! Convert to old naming convention if necessary.
!
	CALL SET_CASE_UP(NAME_CONVENTION,IZERO,IZERO)
	IF(NAME_CONVENTION .EQ. 'K_FOR_I')THEN
	ELSE IF(NAME_CONVENTION .EQ. 'X_FOR_I')THEN
	  DO I=1,NUM_SPECIES
	    IF(SPECIES_ABR(I)(2:2) .EQ. 'k')SPECIES_ABR(I)(2:2)='x'
	  END DO
	ELSE
	  WRITE(T_OUT,*)'Don''t recognize naming convention in DISPGEN'
	  WRITE(T_OUT,*)'NAME_CONVENTION= ',NAME_CONVENTION
	  STOP
	END IF
!
! 
! Open MODEL data file to get N_S, DIE_AUTO, and DIE_WI for each species.
!
	IOS=1
	DO WHILE(IOS .NE. 0)
	  IOS=0
	  IF(ASK .AND. IOS .EQ. 0)THEN
	    IOS=1  		!File has to exist if filename input.
	    CALL GET_FILENAME(FILENAME,1,2,'MODEL','!',IOS)
	  ELSE
	    FILENAME=' '
	    FILENAME=DIR_NAME(1:LEN_DIR)//'MODEL'//FILE_EXTENT
	    INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	    IF(.NOT. FILE_PRES)FILENAME=BLANK
	  END IF
	  CALL RD_MODEL_FILE(FILENAME,LU,IOS)
	END DO
	IF(IOS .NE. 0)STOP
	CALL RD_DBLE(STARS_MASS,'MASS',T_IN,LUER,'Stars mass')
	IF(STARS_MASS .EQ. 0)CALL RD_MODEL_DBLE(STARS_MASS,'MASS')
	CALL RD_MODEL_LOG(DIE_AS_LINE,'DIE_AS_LINE')
	CALL RD_MODEL_DBLE(VSM_DIE_KMS,'VSM_DIE')
!
! Read in other parameters from batch file. These must be in order.
!
	CALL RD_LOG(ONLY_OBS_LINES,'ONLY_OBS_LINES',T_IN,LUER,
	1            'Observed lines only?')
! 
!
!
! Allocate population vectors.
!
	ALLOCATE (POP_SPECIES(ND,NUM_SPECIES)); POP_SPECIES=0.0D0
!
	EQ_TEMP=1
	NT=2
	ID=1
	DO ISPEC=1,NUM_SPECIES
!
! Open input file for each species: This file is an ASCI (new format) of a direct
! access file (old format) and contains the ouput for the entrie atom.
! The file need not exist. We do this for all species including H and He.
!
	  IOS=1
	  DO WHILE(IOS .NE. 0)
	    IOS=0
	    IF(ASK .AND. IOS .EQ. 0)THEN
	      IOS=1  		!File has to exist if filename input.
	      CALL GET_FILENAME(FILENAME,1,2,TRIM(SPECIES(ISPEC)),'!',IOS)
	    ELSE
	      FILENAME=DIR_NAME(1:LEN_DIR)//'POP'//TRIM(SPECIES(ISPEC))//FILE_EXTENT
	      INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	      IF(.NOT. FILE_PRES)FILENAME=BLANK
	    END IF
!
! Read in basic model data (TIMECHK, and ABUNDC)
!
	    IF(FILENAME .NE. BLANK)THEN
	      CALL OP_SPEC_FILE_V2(FILENAME,LU,AT_ABUND(ISPEC),POP_SPECIES(1,ISPEC),
	1          ND,FORMAT_DATE,IOS,TIME,TRIM(SPECIES(ISPEC)))
	    END IF
	  END DO
!
	  PREV_STAGE_PRES=.FALSE.
	  IF(FILENAME .NE. BLANK)THEN
!
	    DO J=1,MAX_IONS_PER_SPECIES
	      ION_ID(ID)=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))
	      CALL RD_POP_DIM(ATM(ID)%NXzV_F, ATM(ID)%XzV_PRES,
	1                     TRIM(ION_ID(ID)),FORMAT_DATE,LU)
!
! Only allocate memory if ion is available. We need, however, to allocate
! for the ion also.
!
	      IF(ATM(ID)%XzV_PRES .OR. PREV_STAGE_PRES)THEN
	        CALL RD_MODEL_SPEC_INFO( TRIM(ION_ID(ID)),
	1          ATM(ID)%NXzV_F,       ATM(ID)%XzV_PRES,  ATM(ID)%NXzV,    
	1          ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV)
                SPECIES_PRES(ISPEC)=.TRUE.
	        IF(SPECIES_BEG_ID(ISPEC) .EQ. 0)SPECIES_BEG_ID(ISPEC)=ID
	        NF=ATM(ID)%NXzV_F
	        NS=ATM(ID)%NXzV
!
	                      ALLOCATE (ATM(ID)%XzV(NS,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE(NS,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%dlnXzVLTE_dlnT(NS,ND),STAT=IOS)
!
	  	IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzV_F(NF,ND),STAT=IOS)
	  	IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F(NF,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%W_XzV_F(NF,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%AXzV_F(NF,NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%EDGEXzV_F(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%F_TO_S_XzV(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%INT_SEQ_XzV(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLEVNAME_F(NF),STAT=IOS)
!
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV_F(ND),STAT=IOS)       ; ATM(ID)%DXzV_F(:)=0.0D0
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV(ND),STAT=IOS)         ; ATM(ID)%DXzV(:)=0.0D0
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GXzV_F(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%ARAD(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM2(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM4(NF),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%OBSERVED_LEVEL(NF),STAT=IOS)
	        IF(IOS .NE. 0)THEN
	          WRITE(LUER,*)'Error in CMF_FLUX'
	          WRITE(LUER,*)'Unable to allocate arrays for species XzV'
	          WRITE(LUER,*)'STATUS=',IOS
	          STOP
	        END IF
!
	        DO I=1,ATM(ID)%NXzV_F ; ATM(ID)%F_TO_S_XzV(I)=I ; END DO
	        CALL RD_ION_POP_V3(ATM(ID)%XzV_F, ATM(ID)%DXzV_F,
	1             ATM(ID)%XzV_PRES,    ATM(ID)%NXzV_F,
	1             ATM(ID)%XzV_OSCDATE, TRIM(ION_ID(ID)),
	1             FORMAT_DATE,LU,ND,SCRAT,LUSCRAT,SCRATREC)
	        SPECIES_LNK(ID)=ISPEC
	        PREV_STAGE_PRES=ATM(ID)%XzV_PRES
	        IF(.NOT. ATM(ID)%XZV_PRES)SPECIES_END_ID(ISPEC)=ID
	        ATM(ID)%EQXzV=EQ_TEMP
	        IF(.NOT. ATM(ID)%XzV_PRES)THEN
	          NT=NT+1
	          EQ_SPECIES(ISPEC)=EQ_TEMP
	          EQ_TEMP=EQ_TEMP+1
	        ELSE
	          NT=NT+ATM(ID)%NXzV
	          EQ_TEMP=EQ_TEMP+ATM(ID)%NXzV
	        END IF
	        ID=ID+1
	      END IF
	    END DO
	  END IF
	END DO			!Over species
!
	NUM_IONS=ID-1
!
	CLOSE(UNIT=LU)
	EQNE=NT-1
!
	CALL CLEAN_MODEL_STORE
!
! 
!
! Read in oscilator and photoionization data.
!
	GF_CUT=0.0D0
	GF_LEV_CUT=1000
	GF_ACTION=' '
	MIN_NUM_TRANS=1000
	XRAYS=.FALSE.
	T2=0.0D0
!
	DO ID=NUM_IONS,1,-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    FILENAME=TRIM(ION_ID(ID))//'_F_OSCDAT'
	    CALL GENOSC_V8(ATM(ID)%AXzV_F,  ATM(ID)%EDGEXzV_F,
	1          ATM(ID)%GXzV_F, ATM(ID)%XzVLEVNAME_F, 
	1          ATM(ID)%ARAD,ATM(ID)%GAM2,ATM(ID)%GAM4,
	1          ATM(ID)%OBSERVED_LEVEL,T1,ATM(ID)%ZXzV,
	1          ATM(ID)%NEW_XzV_OSCDATE, ATM(ID)%NXzV_F, I,
	1          GF_ACTION,GF_CUT,GF_LEV_CUT,MIN_NUM_TRANS,ONLY_OBS_LINES,
	1          LUIN,LU_TMP,FILENAME)
	    IF(ATM(ID)%XzV_OSCDATE .NE. ATM(ID)%NEW_XzV_OSCDATE)THEN
	       WRITE(T_OUT,*)'Warning --- invalid date for ',ION_ID(ID)
	       WRITE(T_OUT,*)'Old oscilator date:',ATM(ID)%XzV_OSCDATE
	       WRITE(T_OUT,*)'New oscilator date:',ATM(ID)%NEW_XzV_OSCDATE
	    END IF
	    FILENAME=TRIM(ION_ID(ID))//'_F_TO_S'
	    CALL RD_F_TO_S_IDS(ATM(ID)%F_TO_S_XzV,ATM(ID)%INT_SEQ_XzV,
	1           ATM(ID)%XzVLEVNAME_F,
	1           ATM(ID)%NXzV_F,ATM(ID)%NXzV,LUIN,FILENAME)
	    CALL RDPHOT_GEN_V1(ATM(ID)%EDGEXzV_F, ATM(ID)%XzVLEVNAME_F,
	1            ATM(ID)%GIONXzV_F,      AT_NO(SPECIES_LNK(ID)),
	1            ATM(ID)%ZXzV,           ATM(ID)%NXzV_F,
	1            ATM(ID)%XzV_ION_LEV_ID, ATM(ID)%N_XzV_PHOT, NPHOT_MAX,
	1            ATM(ID+1)%XzV_PRES,     ATM(ID+1)%EDGEXzV_F,
	1            ATM(ID+1)%GXzV_F,       ATM(ID+1)%F_TO_S_XzV,
	1            ATM(ID+1)%XzVLEVNAME_F, ATM(ID)%NXzV_F, 
	1            XRAYS,ID,ION_ID(ID),LUIN,LU_TMP)
            IF(ATM(ID+1)%XzV_PRES)ATM(ID)%GIONXzV_F=ATM(ID+1)%GXzV_F(1)
 	    IF(.NOT. DIE_AS_LINE .AND. (ATM(ID)%DIE_AUTO_XzV .OR. ATM(ID)%DIE_WI_XzV) )THEN
	      FILENAME='DIE'//TRIM(ION_ID(ID))
	      CALL RD_PHOT_DIE_V1(ID,  ATM(ID)%EDGEXzV_F,    ATM(ID)%XzVLEVNAME_F,
	1             ATM(ID)%NXzV_F,  ATM(ID)%GIONXzV_F,
	1             VSM_DIE_KMS,     ATM(ID)%DIE_AUTO_XzV, ATM(ID)%DIE_WI_XzV,
	1             ION_ID(ID),LUIN,LU_TMP,FILENAME)
 	    ELSE IF(DIE_AS_LINE .AND. (ATM(ID)%DIE_AUTO_XzV .OR. ATM(ID)%DIE_WI_XzV) )THEN
	      WRITE(LUER,*)'Warning in CMF_FLUX_V5'
	      WRITE(LUER,*)'Dielectronic lines from files not included'//
	1                    ' in spectrum calculation.'
	    END IF
	  END IF
	END DO		!Over NUM_SPECIES
!
! 
! Determine the total number of bound-bound transitions in the model.
!
	N_MAX=0
	N_LINE_MAX=0
	DO ID=1,NUM_IONS
	  N_MAX=MAX(N_MAX,ATM(ID)%NXzV_F)
	  IF(ATM(ID)%XzV_PRES)THEN
	    DO I=1,ATM(ID)%NXzV_F-1
	      DO J=I+1,ATM(ID)%NXzV_F
	        IF(ATM(ID)%AXzV_F(I,J) .NE. 0)N_LINE_MAX=N_LINE_MAX+1
	      END DO
	    END DO
	  END IF
	END DO
!
	ND_MAX=4*ND-1
	NC_MAX=4*NC
	NP_MAX=ND_MAX+NC_MAX
!
!*****************************************************************************
!*****************************************************************************
!
! Call program to compute emissivities, opacities, and the mean intensities
! in the CMF. The routine then passes (calls) OBS_FRAM_SUb to compute the
! model spectrum
!
	CALL PAR_OPAC_SUB(ND,NC,NP,ND_MAX,NP_MAX,NT,N_LINE_MAX)
!
	CALL TUNE(3,' ')
	STOP
!
	END
