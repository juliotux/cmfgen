!
! General routine to guess the departure coefficients for a new species or
! ionization stage. The following files are required:
!
!     EDDFACTOR
!     RVTJ
!     XzV_F_OSCDAT
!     XzIVOUT   (for high ionization species only)
!
! These can come from a similar model but the frequency range must be sufficient.
! For high ionization species, the model should be identical, since the departure
! coefficients are very sensitive to the electron temperature.  Program estimates 
! the ground state departure coefficient, and use the same excitation temperature
! for all other levels. If adding a whole new species, start with lowest ioization
! species.
!
	PROGRAM GUESS_DC
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered 03-Feb-2010 : Simply set DC=1 if RR/PR is with 10% of unity. There
!                         is potentially a problem with the simple iteration
!                         for low iozaton stages with high temperatures.
! Altered 19-Sep-2007 : Species length increased from 5 to 6 (handle ArVIII).
! Altered 14-May-2006 : For a high ionization species, ground state population is 
!                       read in from file. This allows estimate of ion population 
!                       to be made.
!
	INTEGER NCF
	INTEGER ND
	CHARACTER*10 DATA_TYPE
	CHARACTER*40 FILE_DATE
	CHARACTER*80 FILENAME
	CHARACTER*132 STRING
!
	INTEGER NUM_FILES
	INTEGER ID
!
	REAL*8, ALLOCATABLE :: TA(:)
	REAL*8, ALLOCATABLE :: TB(:)
	REAL*8, ALLOCATABLE :: TC(:)
	REAL*8, ALLOCATABLE :: TAU_ROSS(:)
	REAL*8, ALLOCATABLE :: TAU_ES(:)
!
	REAL*8, ALLOCATABLE :: PHOT_SUM(:)
	REAL*8, ALLOCATABLE :: RECOM_SUM(:)
	REAL*8, ALLOCATABLE :: GS_DC(:)
	REAL*8, ALLOCATABLE :: DC(:,:)
	REAL*8, ALLOCATABLE :: T_EXC(:)
!
	REAL*8, ALLOCATABLE :: GS_ION_POP(:)
	REAL*8, ALLOCATABLE :: DC_RUB(:)
!
! Needed when reading EDDFACTOR
!
	REAL*8, POINTER :: RJ(:,:)
	REAL*8, POINTER :: NU(:)
!
! Needed when reading RVTJ.
!
	REAL*8 RMDOT
	REAL*8 RLUM
	REAL*8 ABUND_HYD
	INTEGER ND_ATM,NC_ATM,NP_ATM
	CHARACTER*21 TIME
	REAL*8, ALLOCATABLE :: R(:)
	REAL*8, ALLOCATABLE :: V(:)
	REAL*8, ALLOCATABLE :: SIGMA(:)
	REAL*8, ALLOCATABLE :: T(:)
	REAL*8, ALLOCATABLE :: ION_POP(:)
	REAL*8, ALLOCATABLE :: ED(:)
	REAL*8, ALLOCATABLE :: ROSS_MEAN(:)
	REAL*8, ALLOCATABLE :: FLUX_MEAN(:)
	REAL*8, ALLOCATABLE :: POP_ATOM(:)
	REAL*8, ALLOCATABLE :: MASS_DENSITY(:)
	REAL*8, ALLOCATABLE :: POPION(:)
	REAL*8, ALLOCATABLE :: CLUMP_FAC(:)
!
! Used when reading file containing energy levels, oscillator strengths, etc
! (e.g. CV_F_OSCDAT)
!
	INTEGER, PARAMETER :: N_MAX=5000
	CHARACTER*30 NAME(N_MAX)
	REAL*8 FEDGE(N_MAX)
	REAL*8 ENERGY(N_MAX)
	REAL*8 G(N_MAX)
	REAL*8 ION_EN
	REAL*8 ZION
	CHARACTER*30 EN_DATE
	INTEGER NLEV
!
	CHARACTER*6 METHOD,TYPE_ATM
	CHARACTER*10 NAME_CONVENTION
!
        INTEGER ACCESS_F
        INTEGER, PARAMETER :: EDD_CONT_REC=3
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
! Miscellaneous variables.
!
	INTEGER IOS			!Used for Input/Output errors.
	INTEGER I,J,ML
	INTEGER EDGE_ML
	INTEGER ST_REC
	INTEGER REC_LENGTH
	INTEGER ND_RD,NLEV_RD
	INTEGER GION_LOW,GION_UP
	REAL*8 NU1,NU2
	REAL*8 RJ1,RJ2
	REAL*8 T1,T2,T3,T4
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: T_IN=5		!For terminal input 
	INTEGER, PARAMETER :: T_OUT=6           !For terminal output
	INTEGER, PARAMETER :: LU_IN=10		!For file I/O
	INTEGER, PARAMETER :: LU_OUT=11
	INTEGER, PARAMETER :: LU_HEAD=12
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN               
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
!
	INTEGER GET_INDX_DP
!
	CHARACTER*80 RVTJ_FILE_NAME
	CHARACTER*6 SPECIES
!
	REAL*8 SPEED_OF_LIGHT
	CHARACTER*30 UC
	EXTERNAL SPEED_OF_LIGHT,UC
!
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
!
	METHOD='LOGMON'
	TYPE_ATM=' '
!
        CALL DIR_ACC_PARS(REC_SIZE,UNIT_SIZE,WORD_SIZE,N_PER_REC)
!
!  Read in EDDFACTOR file. This is used to compute the photoionization and
!         recombination rates.
!
	FILENAME='EDDFACTOR'
5	CALL GEN_IN(FILENAME,'First data file')
	CALL READ_DIRECT_INFO_V3(I,REC_LENGTH,FILE_DATE,FILENAME,LU_IN,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error opening/reading INFO file: check format'
	  WRITE(T_OUT,*)'Also check error file or fort.2'
	  GOTO 5
	END IF
	OPEN(UNIT=LU_IN,FILE=FILENAME,STATUS='OLD',ACTION='READ',
	1                 RECL=REC_LENGTH,ACCESS='DIRECT',FORM='UNFORMATTED',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	     WRITE(T_OUT,*)'Error opening ',TRIM(FILENAME)
	     WRITE(T_OUT,*)'IOS=',IOS
	     GOTO 5
	  END IF
	  READ(LU_IN,REC=3)ST_REC,NCF,ND
	  ND=ND; NCF=NCF
	  ALLOCATE (RJ(ND,NCF))
	  ALLOCATE (NU(NCF))
	  DO ML=1,NCF
	    READ(LU_IN,REC=ST_REC+ML-1,IOSTAT=IOS)(RJ(I,ML),I=1,ND),NU(ML)
	    IF(IOS .NE. 0)THEN
	      WRITE(T_OUT,*)'Error reading all frequencies'
	      NCF=ML-1
	      EXIT
	    END IF
	  END DO
	CLOSE(LU_IN)
	WRITE(T_OUT,*)'Successfully read in ',TRIM(FILENAME),' file as MODEL A (default)'
	WRITE(T_OUT,*)'Number of depth points is',ND
	WRITE(T_OUT,*)'Number of frequencies is ',NCF
	WRITE(T_OUT,*)' '
!
!
!
! *************************************************************************
!
! Read in basic model [i.e. R, V, T, SIGMA etc ] from RVTJ file.
!
! The file is a SEQUENTIAL (new version) or DIRECT (old version) ACCESS
! file.
!
! *************************************************************************
!
10	RVTJ_FILE_NAME='RVTJ'
	CALL GEN_IN(RVTJ_FILE_NAME,'File with R, V, T etc (RVTJ)')
	OPEN(UNIT=LU_IN,FILE=RVTJ_FILE_NAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Unable to open RVTJ: IOS=',IOS
	    GOTO 10
	  END IF
	CLOSE(LU_IN)
	CALL RD_RVTJ_PARAMS_V2(RMDOT,RLUM,ABUND_HYD,TIME,NAME_CONVENTION,
	1             ND_ATM,NC_ATM,NP_ATM,RVTJ_FILE_NAME,LU_IN)
	ALLOCATE (R(ND_ATM))
	ALLOCATE (V(ND_ATM))
	ALLOCATE (SIGMA(ND_ATM))
	ALLOCATE (T(ND_ATM))
	ALLOCATE (ION_POP(ND_ATM))
	ALLOCATE (ED(ND_ATM))
	ALLOCATE (ROSS_MEAN(ND_ATM))
	ALLOCATE (FLUX_MEAN(ND_ATM))
	ALLOCATE (POP_ATOM(ND_ATM))
	ALLOCATE (MASS_DENSITY(ND_ATM))
	ALLOCATE (POPION(ND_ATM))
	ALLOCATE (CLUMP_FAC(ND_ATM))
	CALL RD_RVTJ_VEC(R,V,SIGMA,ED,T,ROSS_MEAN,FLUX_MEAN,
	1       POP_ATOM,POPION,MASS_DENSITY,CLUMP_FAC,ND_ATM,LU_IN)
	CLOSE(LU_IN)
!
! Now compute the important optical depth scales.
!
	 ALLOCATE (TA(ND_ATM))
	 ALLOCATE (TB(ND_ATM))
	 ALLOCATE (TC(ND_ATM))
	 ALLOCATE (TAU_ROSS(ND_ATM))
	 ALLOCATE (TAU_ES(ND_ATM))
	 IF(ROSS_MEAN(ND_ATM) .NE. 0)THEN
	   CALL TORSCL(TAU_ROSS,ROSS_MEAN,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
	 ELSE
	  TAU_ROSS(1:ND)=0.0D0
	 END IF
	 TA(1:ND_ATM)=6.65D-15*ED(1:ND_ATM)
	 CALL TORSCL(TAU_ES,TA,R,TB,TC,ND_ATM,METHOD,TYPE_ATM)
!
! 
!
! Loop to do many different species/ionization stages.
!
	SPECIES=' '
5000	CONTINUE
	WRITE(T_OUT,*)' '
	CALL GEN_IN(SPECIES,'Species (e.g., OIV, OVI) to guess d.c.''s for (or exit)')
	IF(UC(SPECIES(1:2)) .EQ. 'EX')STOP
	FILENAME=TRIM(SPECIES)//'_F_OSCDAT'
!
! Open Oscillator file.
!
	IOS=100
	DO WHILE(IOS .NE. 0)
	  CALL GEN_ASCI_OPEN(LU_HEAD,'HEAD_INFO','UNKNOWN',' ','WRITE',IZERO,IOS)
	  CALL GEN_IN(FILENAME,'Oscillator file')
	  CALL RD_ENERGY(NAME,G,ENERGY,FEDGE,NLEV,N_MAX,
	1       ION_EN,ZION,EN_DATE,FILENAME,LU_IN,LU_HEAD,IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error occurred reading Oscillator file: try again'
	  END IF
	CLOSE(LU_HEAD)
	END DO
	WRITE(T_OUT,*)'Successfully read oscillator file'
	WRITE(T_OUT,*)' '
!
! Check EDDFACTOR file extends to high enough frequencies.
!
	IF(FEDGE(1) .GT. NU(1))THEN
	  WRITE(T_OUT,*)' '
	  WRITE(T_OUT,*)'Error --- maximum frequency in EDDFACTOR too small '
	  WRITE(T_OUT,*)'EDDFACTOR NU_MAX=',NU(1)
	  WRITE(T_OUT,*)'Required  NU_MAX>',FEDGE(1)
	  WRITE(T_OUT,*)' '
	  STOP
	END IF
!
	ML=1
	DO WHILE(FEDGE(1) .LT. NU(ML))
	  ML=ML+1
	END DO
	EDGE_ML=ML-1
!
! Compute the recombination and photoionization rate to the ground state.
! For simplicity we assume the ground state population is set by a balance
! between photoionizations and recombinations. We ignore the frequency
! dependence of the photoionization cross-section, and note that its numerical
! value does not matter.
!
	IF(.NOT. ALLOCATED(PHOT_SUM))THEN
	  ALLOCATE (PHOT_SUM(1:ND))
	  ALLOCATE (RECOM_SUM(1:ND))
	END IF
	PHOT_SUM(1:ND)=0.0D0
	RECOM_SUM(1:ND)=0.0D0
	DO ML=1,EDGE_ML-1
	  DO I=1,ND
	    PHOT_SUM(I)=PHOT_SUM(I)+0.5D0*(NU(ML)-NU(ML+1))*(RJ(I,ML)+RJ(I,ML+1))
	    T1=(TWOHCSQ*(NU(ML)**3)+RJ(I,ML))*EXP(-HDKT*NU(ML)/T(I))
	    T2=(TWOHCSQ*(NU(ML+1)**3)+RJ(I,ML+1))*EXP(-HDKT*NU(ML+1)/T(I))
	    RECOM_SUM(I)=RECOM_SUM(I)+0.5D0*(NU(ML)-NU(ML+1))*(T1+T2)
	  END DO
	END DO
!
! Add contribution where NU(EDGE_ML) is very different from FEDGE(1)
!
	IF(NU(EDGE_ML) .GT. 1.00000001D0*FEDGE(1))THEN
	  DO ML=1,10
	    NU1=NU(EDGE_ML)-(ML-1)*(NU(EDGE_ML)-FEDGE(1))/10
	    NU2=NU(EDGE_ML)-ML*(NU(EDGE_ML)-FEDGE(1))/10
	    T1=(NU(EDGE_ML)-NU1)/(NU(EDGE_ML)-NU(EDGE_ML+1))
	    T2=(NU(EDGE_ML)-NU2)/(NU(EDGE_ML)-NU(EDGE_ML+1))
	    DO I=1,ND
	      RJ1=(1.0D0-T1)*RJ(I,EDGE_ML)+T1*RJ(I,EDGE_ML+1)
	      RJ2=(1.0D0-T2)*RJ(I,EDGE_ML)+T2*RJ(I,EDGE_ML+1)
	      PHOT_SUM(I)=PHOT_SUM(I)+0.5D0*(NU1-NU2)*(RJ1+RJ2)
	      T3=(TWOHCSQ*(NU1**3)+RJ1)*EXP(-HDKT*NU1/T(I))
	      T4=(TWOHCSQ*(NU2**3)+RJ2)*EXP(-HDKT*NU2/T(I))
	      RECOM_SUM(I)=RECOM_SUM(I)+0.5D0*(NU1-NU2)*(T3+T4)
	    END DO
	  END DO
	END IF
!
! Compute the ground-state departure coefficient, and the excitation temperature
! of the ground state.
!
	IF(.NOT. ALLOCATED(GS_DC))THEN
	  ALLOCATE (GS_DC(1:ND))
	  ALLOCATE (T_EXC(1:ND))
	END IF
	GS_DC(1:ND)=RECOM_SUM(1:ND)/PHOT_SUM(1:ND)
	T1=5.0D0
	DO I=1,ND
	  IF(ABS(GS_DC(I)-1.0D0) .LT. 0.1D0)THEN
	    T1=T(I)
	  ELSE
	    DO J=1,10
	      T1=LOG( GS_DC(I)*(T1/T(I))**1.5 )/HDKT/FEDGE(1)+1.0/T(I)
	      T1=1.0/T1
	    END DO
	  END IF
	  T_EXC(I)=T1
	  IF(I .EQ. 1 .OR. I .EQ. ND)THEN
	    WRITE(T_OUT,'(I4,2X,A,F8.3,3X,A,F8.3,3X,A,ES8.2)')I,'T=',T(I),'T_EXC=',T_EXC(I),'RR/PR=',GS_DC(I)
	  END IF
	END DO	
	WRITE(T_OUT,*)GS_DC(1),FEDGE(1)
!
! Compute departure coefficients of all levels. We assume that they have the same
! departure coefficient as the ground state.
!
	IF(ALLOCATED(DC))DEALLOCATE(DC)
	ALLOCATE (DC(NLEV,ND))
	DO I=1,ND
	   DO J=1,NLEV
	     DC(J,I)=((T(I)/T_EXC(I))**1.5 )*EXP(HDKT*FEDGE(J)*(1.0D0/T_EXC(I)-1.0D0/T(I)))
	   END DO
	END DO
!
! Read in ION file to get ground state population. For a lower ionization species,
! ion population does not matter, since the actal ion population gets used when
! species is read into CMFGEN.
!
	WRITE(6,*)' '
	WRITE(6,*)'If adding CV to model with C2, CIII, CIV enter CIVOUT '
	WRITE(6,*)'If adding CI to model with C2, CIII, CIV enter "" '
	WRITE(6,*)' '
	FILENAME=' '
	CALL GEN_IN(FILENAME,'File with ground state population ')
	  IF(FILENAME .NE. ' ')THEN
            CALL GEN_ASCI_OPEN(LU_IN,FILENAME,'OLD',' ','READ',IZERO,IOS)
            I=0
            STRING=' '
            DO WHILE(INDEX(STRING,'!Format date') .EQ. 0 .AND. I .LE. 10)
              I=I+1
              READ(LU_IN,'(A)')STRING
            END DO
            IF( INDEX(STRING,'!Format date') .EQ. 0)REWIND(LU_IN)
            READ(LU_IN,*)T1,T2,NLEV_RD,ND_RD
	    IF(ALLOCATED(DC_RUB))DEALLOCATE(DC_RUB)
	    ALLOCATE (DC_RUB(NLEV_RD))
            IF(ND_RD .NE. ND)THEN
	      WRITE(T_OUT,*)'ND in ion file must be same as current ND'
	      STOP
	    END IF
	    IF(.NOT. ALLOCATED(GS_ION_POP))ALLOCATE (GS_ION_POP(ND))
            DO J=1,ND_RD
              READ(LU_IN,*)T1,GS_ION_POP(J)
              READ(LU_IN,*)(DC_RUB(I),I=1,NLEV_RD)
	    END DO
	  CLOSE(LU_IN)
	  GION_UP=1; GION_LOW=1
	  CALL GEN_IN(GION_UP,'Statistical weight for upper ion level (i.e., 2 for CV [=G(CVI)])')
	  CALL GEN_IN(GION_LOW,'Statistical weight for lower ion level(i.e., 1 for CV)')
!
	  DO I=1,ND
	    T1=DLOG(2.07078D-22*ED(I)*DC(1,I))
	    T1=GION_LOW*EXP(T1+HDKT*FEDGE(1)/T(I))/(T(I)**1.5D0)/GION_UP
	    ION_POP(I)=GS_ION_POP(I)/T1
	  END DO
	ELSE
	  ION_POP(1:ND)=1.0D-10
	END IF
!
	FILENAME=TRIM(SPECIES)//'_IN'
	CALL GEN_IN(FILENAME,'Output file for DCs --- old file will be overwritten')
	CALL GEN_IN(NLEV,'Number of levels to be output to file')
	CALL GEN_ASCI_OPEN(LU_OUT,FILENAME,'UNKNOWN',' ','WRITE',IZERO,IOS)
	  WRITE(LU_OUT,'(/,1X,A,T40,A)')'07-Jul-1997','!Format date'
	  WRITE(LU_OUT,2120)R(ND),RLUM,NLEV,ND
	  DO I=1,ND
	    WRITE(LU_OUT,2122)R(I),ION_POP(I),ED(I),T(I),0.0,V(I),CLUMP_FAC(I)
	    WRITE(LU_OUT,'(1X,1P,5E17.7)')(DC(J,I),J=1,NLEV)
	  END DO
	CLOSE(LU_OUT)
	WRITE(T_OUT,*)
	WRITE(T_OUT,*)' Check output file to ensure d.c''s ---> 1 at depth'
	WRITE(T_OUT,*)
2120	FORMAT(/,1X,ES14.8,5X,1PE12.6,5X,0P,I4,5X,I4)
2122	FORMAT(/,1X,1P,E16.8,6E17.8)
!
! Loop back to do additional species and ionization stages.
!
	SPECIES='EXIT'
	GOTO 5000
	END
