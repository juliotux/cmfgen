C
	PROGRAM DISPGEN
	USE MOD_DISP
	USE MOD_USR_OPTION
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 19-AUg-2015 : Changed to GENOSC_V9 (from _V8).
! Altered 13-May-2015 : Changed to GENOSC_V8 (from _V5). Updated MOD_DISP.
! Altered 07-Sep-2005 : XRAYS option set to TRUE. Lithium cross-sections will be set.
! Altered 20-Apr-2004 : Use RDHOTGEN_V2 will allows for dynamic smoothing of 
!                         photoioization cross-sections.
!
	INTEGER ND,NP,NC
	INTEGER N_MAX,ND_MAX,NC_MAX,NP_MAX
	INTEGER N_LINE_MAX,N_PLT_MAX
	INTEGER NF_SUM
C
	CHARACTER*20 TIME
	CHARACTER*80 FILENAME,BLANK,DIR_NAME,STRING
	CHARACTER*20 FILE_EXTENT
	CHARACTER*11 FORMAT_DATE
	CHARACTER*11 RVTJ_FORMAT_DATE
	CHARACTER*10 NAME_CONVENTION
	CHARACTER(LEN=10) GF_ACTION
C
	LOGICAL ASK  		!Ask of filenames or uset defaults.
	LOGICAL FILE_PRES
	LOGICAL SCRAT
	LOGICAL XRAYS
	INTEGER I,J,IOS,SCRATREC
	INTEGER FILE_OPT
	INTEGER LEN_DIR
	INTEGER GF_LEV_CUT
	INTEGER MIN_NUM_TRANS
	REAL*8 GF_CUT
	REAL*8 T1,T2
	REAL*8 RMDOT,RLUM
!
! Variables for reading in photoioization data.
!
	REAL*8 SIG_GAU_KMS
	REAL*8 FRAC_SIG_GAU
	REAL*8 CUT_ACCURACY
	LOGICAL ABOVE_EDGE
!
	CHARACTER*80 TMP_STRING
	INTEGER ID,ISPEC
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: T_IN=5		!Terminal Input
	INTEGER, PARAMETER :: T_OUT=6		!Terminal Output
	INTEGER, PARAMETER :: LUIN=8
	INTEGER, PARAMETER :: LUMOD=9
	INTEGER, PARAMETER :: LU_TMP=10
	INTEGER, PARAMETER :: LU=30			!Used for pop file io
	INTEGER, PARAMETER :: LUSCRAT=33
C
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL PREV_STAGE_PRES
C
	DATA BLANK/' '/
C
C Set constants.
C
	CHIBF=2.815D-06
	CHIFF=3.69D-29
	HDKT=4.7994145D0
	TWOHCSQ=0.0147452575D0
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
	OPLIN=2.6540081D+08
	EMLIN=5.27296D-03
!
! 
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
	IF(ID .NE. NSPEC)THEN
	  WRITE(T_OUT,*)'Error in CMFGEN: Invalid species setup'
	  STOP
	END IF
C
C Initilaization: These parameters will remain as they are when a species
C is not present.
C
	ABUND(:)=0.0D0
	DO I=1,NSPEC*NION_MAX
	  ATM(I)%NXzV_F=1
	  ATM(I)%XzV_PRES=.FALSE.
	END DO
	SPECIES_BEG_ID(1:NSPEC)=0
	SPECIES_END_ID(1:NSPEC)=0
!
! Read in the gaunt factors for individual l states of hydrogen. 
!
	CALL RD_HYD_BF_DATA(LUIN,LUMOD,T_OUT)
!
! Read in atomic data for 2-photon transitions.
!
        CALL RD_TWO_PHOT(LUIN,L_TRUE)
!
! Read in X-ray photoionization cross-sections.
!
        CALL RD_XRAY_FITS(LUIN)
!
! *************************************************************************
!
! Read in basic model [i.e. H, HeI, He2, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
!
	ASK=.FALSE.
	SCRAT=.FALSE.
!
	FILENAME='RVTJ'
10	DIR_NAME=' '
	IOS=2			!Filename has to exits, blank not allowed.
	WRITE(T_OUT,*)'Append (ask) to so that subsequent FILE names '//
	1         'are not defaulted'
	WRITE(T_OUT,*)'Append (scrat) to get scratch output.'
	WRITE(T_OUT,*)' '
	CALL GEN_IN(FILENAME,'Structure file')
	STRING=FILENAME
	CALL SET_CASE_UP(STRING,0,0)
	IF( INDEX(STRING,'(SCRAT)') .NE. 0)SCRAT=.TRUE.
	IF( INDEX(STRING,'(ASK)') .NE. 0)ASK=.TRUE.
	I= INDEX(STRING,'(')
	IF(I .NE. 0)FILENAME=FILENAME(1:I-1)//' '
	INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	IF(.NOT. FILE_PRES)THEN
	  FILENAME='../RVTJ'
	  GOTO 10
	END IF
C
C Get default directory.
C
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
C
C Get default extension.
C
	FILE_EXTENT=' '
	J=LEN(FILENAME)
	DO WHILE(J .GT. LEN_DIR)
	  IF( FILENAME(J:J) .EQ. '.' )THEN
	    FILE_EXTENT=FILENAME(J:)
	    J=0
	  END IF
	  J=J-1
	END DO
C
	CALL RD_RVTJ_PARAMS_V3(RMDOT,RLUM,ABUND(1),TIME,NAME_CONVENTION,
	1         ND,NC,NP,RVTJ_FORMAT_DATE,FILENAME,LUIN)
	ALLOCATE (R(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (V(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (SIGMA(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (T(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ED(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (ROSS_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (FLUX_MEAN(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POP_ATOM(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (MASS_DENSITY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (POPION(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CLUMP_FAC(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (CMFGEN_TGREY(ND),STAT=IOS)
	IF(IOS .EQ. 0)ALLOCATE (dE_RAD_DECAY(ND),STAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error in DISPGEN -- error allocating atmospheric vectors'
	  WRITE(T_OUT,*)'STATUS=',IOS
	  STOP
	END IF
	CALL RD_RVTJ_VEC_V2(R,V,SIGMA,ED,T,CMFGEN_TGREY,dE_RAD_DECAY,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,RVTJ_FORMAT_DATE,ND,LUIN)
!
	CALL SET_CASE_UP(NAME_CONVENTION,IZERO,IZERO)
	IF(NAME_CONVENTION .EQ. 'K_FOR_I')THEN
	  Write(T_OUT,*)'Using K convention: i.e. Silicon abrevaition is Sk'
	ELSE IF(NAME_CONVENTION .EQ. 'X_FOR_I')THEN
	  Write(T_OUT,*)'Using X convention: i.e. Silicon abrevaition is Sx'
	  DO I=1,NSPEC
	    IF(SPECIES_ABR(I)(2:2) .EQ. 'k')SPECIES_ABR(I)(2:2)='x'
	  END DO
	ELSE
	  WRITE(T_OUT,*)'Don''t recognize naming convention in DISPGEN'
	  WRITE(T_OUT,*)'NAME_CONVENTION= ',NAME_CONVENTION
	  STOP
	END IF 
!
	WRITE(T_OUT,3)TIME
3	FORMAT(1X,'Model completed on ',A20)
C
	IF(SCRAT)THEN
	  INQUIRE(IOLENGTH=I)R(1)
	  I=ND*I
	  OPEN(UNIT=LUSCRAT,FORM='UNFORMATTED',
	1      ACCESS='DIRECT',STATUS='UNKNOWN',RECL=I)
	  I=0
	  WRITE(LUSCRAT,REC=2)I,ND
	  WRITE(LUSCRAT,REC=3)(R(J),J=1,ND)
	  WRITE(LUSCRAT,REC=4)(V(J),J=1,ND)
	  WRITE(LUSCRAT,REC=5)(SIGMA(J),J=1,ND)
	  SCRATREC=6
	END IF
!
! Allocate population vectors.
!
	ALLOCATE (POPDUM(ND,NSPEC),STAT=IOS) 
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error in DISPGEN -- error allocating POPDUM'
	  WRITE(T_OUT,*)'STATUS=',IOS
	  STOP
	END IF
! 
!
! Open DUM input file: This file is an ASCI (new format) of a direct
! access file (old format) and contains the ouput for the DUM atom.
! The file need not exist. We do this for all species includeing H and He.
!
	ID=0
	NF_SUM=2			!T and Ne
	DO ISPEC=1,NSPEC
	  FILENAME=DIR_NAME(1:LEN_DIR)//'POP'//TRIM(SPECIES(ISPEC))//FILE_EXTENT
	  IF(ASK)THEN
	    WRITE(6,*)' '
	    CALL GEN_IN(FILENAME,'Filename')
	    INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	    IF(.NOT. FILE_PRES)FILENAME=BLANK
	  ELSE
	    INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	    IF(.NOT. FILE_PRES)FILENAME=BLANK
	  END IF
!
! Read in basic model data (TIMECHK, and ABUNDC)
!
	  IF(FILENAME .NE. BLANK)THEN
	    CALL OP_SPEC_FILE_V2(FILENAME,LU,ABUND(ISPEC),
	1          POPDUM(1,ISPEC),ND,
	1          FORMAT_DATE,IOS,TIME,TRIM(SPECIES(ISPEC)))
	  END IF
C
	  IF(FILENAME .NE. BLANK)THEN
	    PREV_STAGE_PRES=.FALSE.
	    DO J=1,NION_MAX
	      ID=ID+1
	      ION_ID(ID)=TRIM(SPECIES_ABR(ISPEC))//TRIM(GEN_ION_ID(J))
	      CALL RD_POP_DIM(ATM(ID)%NXzV_F,ATM(ID)%XzV_PRES,
	1                    TRIM(ION_ID(ID)),FORMAT_DATE,LU)
!
! Only allocate memory if ion is available. We need, however, to allocate
! for the ion also.
!
	      IF(ATM(ID)%XzV_PRES .OR. PREV_STAGE_PRES)THEN
	        IF(SPECIES_BEG_ID(ISPEC) .EQ. 0)SPECIES_BEG_ID(ISPEC)=ID
	        ALLOCATE (ATM(ID)%XzV_F(ATM(ID)%NXzV_F,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLTE_F(ATM(ID)%NXzV_F,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%LOG_XzVLTE_F(ATM(ID)%NXzV_F,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%W_XzV_F(ATM(ID)%NXzV_F,ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV_F(ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%DXzV(ND),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%AXzV_F(ATM(ID)%NXzV_F,ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%EDGEXzV_F(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GXzV_F(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%F_TO_S_XzV(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%ARAD(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM2(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%GAM4(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%OBSERVED_LEVEL(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .EQ. 0)ALLOCATE (ATM(ID)%XzVLEVNAME_F(ATM(ID)%NXzV_F),STAT=IOS)
	        IF(IOS .NE. 0)THEN
	          WRITE(T_OUT,*)'Error in DISPGEN -- error allocating atomic data arrays'
	          WRITE(T_OUT,*)'STATUS=',IOS
	          STOP
	        END IF
!
	        ATM(ID)%DXzV_F(:)=0.0D0
	        ATM(ID)%DXzV(:)=0.0D0
	        DO I=1,ATM(ID)%NXzV_F ; ATM(ID)%F_TO_S_XzV(I)=I ; END DO
!
	        CALL RD_ION_POP_V3(ATM(ID)%XzV_F,ATM(ID)%DXzV_F,
	1       ATM(ID)%XzV_PRES,ATM(ID)%NXzV_F,
	1       ATM(ID)%XzV_OSCDATE,TRIM(ION_ID(ID)),FORMAT_DATE,
	1       LU,ND,SCRAT,LUSCRAT,SCRATREC)
	        SPECIES_LNK(ID)=ISPEC
	        IF(ATM(ID)%XzV_PRES)THEN
	          NF_SUM=NF_SUM+ATM(ID)%NXzV_F
	        ELSE
	          NF_SUM=NF_SUM+1
	        END IF
	        PREV_STAGE_PRES=ATM(ID)%XzV_PRES
	        IF(.NOT. ATM(ID)%XZV_PRES)SPECIES_END_ID(ISPEC)=ID
	      ELSE
	        ID=ID-1
	      END IF
	    END DO
	  END IF
C
	  CLOSE(UNIT=LU)
!
	END DO
	NUM_IONS=ID
! 
!
! Read in oscilator and photoionization data.
!
	GF_LEV_CUT=5
	XRAYS=.TRUE.              !.FALSE.
	T2=0.0D0
!
! Default values for reading in photoioization cross-sections.
!
	SIG_GAU_KMS=1000.0D0
	FRAC_SIG_GAU=0.25
	CUT_ACCURACY=0.02
	ABOVE_EDGE=.TRUE.
	CALL GEN_IN(SIG_GAU_KMS,'Default smoothing for photoionization cross-sections (km/s)')
!
	GF_CUT=0.0D0
	GF_LEV_CUT=5000
	GF_ACTION=' '
	MIN_NUM_TRANS=1000
!
	WRITE(T_OUT,*)' '
	DO ID=NUM_IONS,1,-1
	  IF(ATM(ID)%XzV_PRES)THEN
	    TMP_STRING=TRIM(ION_ID(ID))//'_F_OSCDAT'
	    CALL GENOSC_V9(ATM(ID)%AXzV_F,ATM(ID)%EDGEXzV_F,ATM(ID)%GXzV_F,ATM(ID)%XzVLEVNAME_F,
	1          ATM(ID)%ARAD,ATM(ID)%GAM2,ATM(ID)%GAM4,
	1          ATM(ID)%OBSERVED_LEVEL,T1,ATM(ID)%ZXzV,
	1          ATM(ID)%NEW_XzV_OSCDATE,ATM(ID)%NXzV_F,I,
	1          GF_ACTION,GF_CUT,GF_LEV_CUT,MIN_NUM_TRANS,L_FALSE,L_FALSE,
	1          LUIN,LU_TMP,TRIM(TMP_STRING))
	    IF(ATM(ID)%XzV_OSCDATE .NE. ATM(ID)%NEW_XzV_OSCDATE)THEN
	       WRITE(T_OUT,*)'Warning --- invalid date for ',ION_ID(ID)
	       WRITE(T_OUT,*)'Old oscilator date:',ATM(ID)%XzV_OSCDATE
	       WRITE(T_OUT,*)'New oscilator date:',ATM(ID)%NEW_XzV_OSCDATE
	    END IF
	    CALL RDPHOT_GEN_V2(ATM(ID)%EDGEXzV_F,ATM(ID)%XzVLEVNAME_F,
	1        ATM(ID)%GIONXzV_F,AT_NO(SPECIES_LNK(ID)),
	1        ATM(ID)%ZXzV,ATM(ID)%NXzV_F,
	1        ATM(ID)%XzV_ION_LEV_ID,ATM(ID)%N_XzV_PHOT,NPHOT_MAX,
	1        ATM(ID+1)%XzV_PRES,ATM(ID+1)%EDGEXzV_F,
	1        ATM(ID+1)%GXzV_F,ATM(ID+1)%F_TO_S_XzV,
	1        ATM(ID+1)%XzVLEVNAME_F,ATM(ID)%NXzV_F,
	1        SIG_GAU_KMS,FRAC_SIG_GAU,CUT_ACCURACY,ABOVE_EDGE,
	1        XRAYS,ID,ION_ID(ID),LUIN,LU_TMP)
            IF(ATM(ID+1)%XzV_PRES)ATM(ID)%GIONXzV_F=ATM(ID+1)%GXzV_F(1)
	    WRITE(T_OUT,*)'Successfully read atomic data for species ',ION_ID(ID)
	  END IF
	END DO
	WRITE(T_OUT,*)' '
!
! 
!
	IF(SCRAT)THEN
	  WRITE(LUSCRAT,REC=SCRATREC)0,'Ne and T     '
	  WRITE(LUSCRAT,REC=SCRATREC+1)(ED(I),I=1,ND)
	  WRITE(LUSCRAT,REC=SCRATREC+2)(T(I),I=1,ND)
	  READ(LUSCRAT,REC=2)I
	  WRITE(LUSCRAT,REC=2)I+2,ND
	END IF
!
!
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
	NC_MAX=25
	ND_MAX=3*ND
	NP_MAX=ND_MAX+NC_MAX	
	N_PLT_MAX=MAX(NP_MAX,2*N_LINE_MAX,N_MAX)
!
	WRITE(T_OUT,*)'Total number of levels (including Ne and T) is',NF_SUM
!
!*****************************************************************************
!*****************************************************************************
!
	CALL MAINGEN(RMDOT,RLUM,
	1                ND,NP,NC,
	1                N_MAX,ND_MAX,NC_MAX,NP_MAX,
	1                N_LINE_MAX,N_PLT_MAX)
!
	STOP
!
	END
