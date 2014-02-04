	MODULE HYDRO_PARAM_MODULE
!
	  REAL*8 GAM_EDD
	  REAL*8 GAM_LIM
	  REAL*8 GAM_FULL
	  REAL*8 AMU			!Atomic mass unit
	  REAL*8 MU_ATOM		!Mean atomic mass in amu
	  REAL*8 BC			!Boltzmann constant
	  REAL*8 C_CMS
	  REAL*8 SIGMA_TH		!Thompson cross-section
	  REAL*8 STEFAN_BC		!Stephan-Boltzmann constant
	  REAL*8 LOGG			!Log (surface gravity) (cgs units)
	  REAL*8 GRAV_CON
	  REAL*8 REFERENCE_RADIUS
	  REAL*8 PREV_REF_RADIUS
	  REAL*8 RADIUS_AT_TAU_23
!
! The following quantities are used when integrating the hydrostatic equation.
! KAP is used to refer to mass absorption opacities.
!
	  REAL*8 R_EST
	  REAL*8 T_EST
	  REAL*8 P_EST
	  REAL*8 TAU_EST
	  REAL*8 ED_ON_NA_EST 
	  REAL*8 ATOM_EST
	  REAL*8 TEFF
	  REAL*8 ROSS_ON_ES
	  REAL*8 KAP_ROSS
	  REAL*8 KAP_ES
	  REAL*8 SOUND_SPEED
	  REAL*8 MDOT
!
	  REAL*8 dPdR
	  REAL*8 dTAUdR
!
! These are used for the old model, obtained from RVTJ.
!
          REAL*8, ALLOCATABLE :: OLD_R(:)
          REAL*8, ALLOCATABLE :: OLD_T(:)
          REAL*8, ALLOCATABLE :: OLD_V(:)
          REAL*8, ALLOCATABLE :: OLD_SIGMA(:)
          REAL*8, ALLOCATABLE :: OLD_ED(:)
          REAL*8, ALLOCATABLE :: OLD_POP_ATOM(:)
          REAL*8, ALLOCATABLE :: OLD_POPION(:)
          REAL*8, ALLOCATABLE :: OLD_MASS_DENSITY(:)
          REAL*8, ALLOCATABLE :: OLD_CLUMP_FAC(:)
	  REAL*8, ALLOCATABLE :: OLD_SF(:)
          REAL*8, ALLOCATABLE :: OLD_TAU(:)
          REAL*8, ALLOCATABLE :: OLD_ESEC(:)
          REAL*8, ALLOCATABLE :: OLD_ROSS_MEAN(:)
          REAL*8, ALLOCATABLE :: OLD_FLUX_MEAN(:)
          REAL*8, ALLOCATABLE :: OLD_KAP_ROSS(:)
          REAL*8, ALLOCATABLE :: OLD_KAP_FLUX(:)
          REAL*8, ALLOCATABLE :: OLD_KAP_ESEC(:)
!
	  REAL*8 OLD_TEFF
	  REAL*8 OLD_REF_RADIUS
	  INTEGER OLD_ND
	  LOGICAL PURE_LTE_EST
	  LOGICAL OLD_MODEL
	  LOGICAL WIND_PRESENT 
	  LOGICAL PLANE_PARALLEL_MOD
	  LOGICAL RESET_REF_RADIUS
!
	END MODULE HYDRO_PARAM_MODULE 	
!
!
! Simple program to estimate the hydrostatic structure for a plane-parallel
! atmosphere. It currently requires the RVTJ file of a previously converged
! model, but this could be easily modified.
!
	PROGRAM HYD
	USE GEN_IN_INTERFACE
	USE HYDRO_PARAM_MODULE
	IMPLICIT NONE
!
! Altered 29-Nov-2012 : Put in error check to make sure connection velocity < sound speed.
! Altered Mar2007: Minor bug fixes.
!
	INTEGER NP,NC
!
! The following vectors are used for the atmospheric structure resulting
! from the solution of the hydrostatic and tau equations.
! 
	INTEGER, PARAMETER :: ND_MAX=1000
	INTEGER ND
        REAL*8 R(ND_MAX)
        REAL*8 V(ND_MAX)
        REAL*8 SIGMA(ND_MAX)
        REAL*8 T(ND_MAX)
        REAL*8 ED(ND_MAX)
        REAL*8 TAU(ND_MAX)
        REAL*8 P(ND_MAX)
        REAL*8 ED_ON_NA(ND_MAX)
	REAL*8 dPdR_VEC(ND_MAX)
!
        REAL*8 ROSS_MEAN(ND_MAX)
        REAL*8 FLUX_MEAN(ND_MAX)
        REAL*8 POP_ATOM(ND_MAX)
        REAL*8 MASS_DENSITY(ND_MAX)
        REAL*8 CLUMP_FAC(ND_MAX)
        REAL*8 POPION(ND_MAX)
        REAL*8 CHI_ROSS(ND_MAX)
        REAL*8 GAMMA_FULL(ND_MAX)
!
	REAL*8 TA(ND_MAX)
	REAL*8 TB(ND_MAX)
	REAL*8 TC(ND_MAX)
!
! These vectors are output, and contain the atmospheric structure to be used
! by CMFGEN.
!
	INTEGER NEW_ND
	REAL*8, ALLOCATABLE :: REV_TAU(:)
        REAL*8, ALLOCATABLE :: REV_R(:)
        REAL*8, ALLOCATABLE :: REV_V(:)
        REAL*8, ALLOCATABLE :: REV_SIGMA(:)
        REAL*8, ALLOCATABLE :: REV_POP_ATOM(:)
        REAL*8, ALLOCATABLE :: REV_ED(:)
        REAL*8, ALLOCATABLE :: REV_T(:)
        REAL*8, ALLOCATABLE :: REV_CHI_ROSS(:)
        REAL*8, ALLOCATABLE :: REV_GAMMA_FULL(:)
        REAL*8, ALLOCATABLE :: COEF(:,:)
!
! Wind parameters:
!
	REAL*8 VINF
	REAL*8 BETA
	REAL*8 RMAX
	REAL*8 CONNECTION_RADIUS
	REAL*8 CONNECTION_VEL
	INTEGER CONNECTION_INDX
!
! Parameters for for cumputing the final R grid.
!	
	REAL*8 dLOG_TAU
	REAL*8 V_SCL_FAC
	REAL*8 OBND_PARS(20)
	INTEGER NUM_OBND_PARAMS
	CHARACTER(LEN=16) OUT_BND_OPT
!
	CHARACTER(LEN=20) TIME
	CHARACTER(LEN=80) FILENAME,DIR_NAME,STRING
	CHARACTER(LEN=20) FILE_EXTENT
	CHARACTER(LEN=11) FORMAT_DATE
	CHARACTER(LEN=10) NAME_CONVENTION
!
	INTEGER I,J,IOS,SCRATREC
	INTEGER FILE_OPT
	INTEGER LEN_DIR
	REAL*8 T1,T2,T3
	REAL*8 NI_ZERO
	REAL*8 OLD_TAU_MAX
	REAL*8 SCL_HT
	REAL*8 H
	REAL*8 TAU_MAX
	REAL*8 RLUM
	REAL*8 OLD_MDOT
	REAL*8 RBOUND
	REAL*8 MASS_LOSS_SCALE_FACTOR
!
! Runge-Kutta estimates
!
	REAL*8 dP1,dTAU1
	REAL*8 dP2,dTAU2
	REAL*8 dP3,dTAU3
	REAL*8 dP4,dTAU4
!
	LOGICAL ASK  		!Ask of filenames or use defaults.
	LOGICAL FILE_PRES
	LOGICAL SCRAT
	LOGICAL USE_OLD_VEL
	LOGICAL L_TEMP
	LOGICAL FILE_OPEN
!
! These have cgs units.
!
	REAL*8 ATOMIC_MASS_UNIT,BOLTZMANN_CONSTANT,STEFAN_BOLTZ
	REAL*8 GRAVITATIONAL_CONSTANT,MASS_SUN,SPEED_OF_LIGHT
	INTEGER GET_INDX_DP
	EXTERNAL ATOMIC_MASS_UNIT,BOLTZMANN_CONSTANT,STEFAN_BOLTZ,GET_INDX_DP
	EXTERNAL GRAVITATIONAL_CONSTANT,MASS_SUN,SPEED_OF_LIGHT
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: T_IN=5		!Terminal Input
	INTEGER, PARAMETER :: T_OUT=6		!Terminal Output
	INTEGER, PARAMETER :: LUIN=8
	INTEGER, PARAMETER :: LUSCR=9
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	CHARACTER*2, PARAMETER :: FORMFEED=' '//CHAR(12)
!
! Set constants.
!
	BC=1.0D+04*BOLTZMANN_CONSTANT()   	!erg/10^4 K
	AMU=ATOMIC_MASS_UNIT()          	!gm
	C_CMS=SPEED_OF_LIGHT()          	!cm/s
	GRAV_CON=1.0D-20*GRAVITATIONAL_CONSTANT()*MASS_SUN()
	SIGMA_TH=6.65D-15			!cm^{-2} x 10^10
	STEFAN_BC=STEFAN_BOLTZ()
	MASS_LOSS_SCALE_FACTOR=3.02286D+23
	PLANE_PARALLEL_MOD=.FALSE.
!
! These are the default settings if no old model.
!
	PURE_LTE_EST=.TRUE.
	OLD_TAU_MAX=100.0D0
	RBOUND=0.0D0
	dLOG_TAU=0.25D0
	V_SCL_FAC=0.67D00
	OBND_PARS(:)=0.0D0
	NUM_OBND_PARAMS=1
	OUT_BND_OPT='DEFAULT'
!
!
! *************************************************************************
!
! Read in parameters describing the new model.
!
	CALL GEN_ASCI_OPEN(LUIN,'HYDRO_PARAMS','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(T_OUT,*)'Error opening HYDRO_PARAMS in WIND_HYD, IOS=',IOS
	  STOP
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)

	CALL RD_STORE_DBLE(LOGG,'LOG_G',L_TRUE,'Log surface gravity (cgs units)')
	CALL RD_STORE_DBLE(MU_ATOM,'MU_ATOM',L_TRUE,'Mean atomic mass (amu)')
	CALL RD_STORE_DBLE(TEFF,'TEFF',L_TRUE,'New effective temperature (10^4K)')
	CALL RD_STORE_LOG(PLANE_PARALLEL_MOD,'PP_MOD',L_FALSE,'Plane parallel model?')
	CALL RD_STORE_LOG(OLD_MODEL,'OLD_MOD',L_TRUE,'Does a model exist')
	L_TEMP=.TRUE.; IF(OLD_MODEL)L_TEMP=.FALSE.
	CALL RD_STORE_LOG(WIND_PRESENT,'WIND_PRES',L_TRUE,'Is a wind present?')
	L_TEMP=.TRUE.; IF(WIND_PRESENT)L_TEMP=.FALSE.
	CALL RD_STORE_DBLE(NI_ZERO,'ATOM_DEN',L_TEMP,'Atom density at outer boundary (/cm^3)')
	CALL RD_STORE_LOG(USE_OLD_VEL,'OLD_V',WIND_PRESENT,
	1                'Use old wind velocity law above connection point?')
	L_TEMP=L_TRUE; IF(USE_OLD_VEL)L_TEMP=.FALSE.
	CALL RD_STORE_DBLE(MDOT,'MDOT',L_TEMP,'Mass loss in Msun/yr')
	MDOT=MDOT*MASS_LOSS_SCALE_FACTOR
	L_TEMP=WIND_PRESENT .AND. .NOT. USE_OLD_VEL
	CALL RD_STORE_DBLE(VINF,'VINF',L_TEMP,'Terminal velocity of wind in km/s')
	CALL RD_STORE_DBLE(BETA,'BETA',L_TEMP,'Beta exponent for velocity law')
	CALL RD_STORE_DBLE(CONNECTION_RADIUS,'CON_R',WIND_PRESENT,'Connection radius (10^10cm)')
	CALL RD_STORE_DBLE(REFERENCE_RADIUS,'REF_R',L_TRUE,'Reference radius (10^10cm)')
	CALL RD_STORE_DBLE(RMAX,'MAX_R',L_TEMP,'Maximum radius in terms of Connection radius')
	RMAX=RMAX*CONNECTION_RADIUS
	CALL RD_STORE_DBLE(CONNECTION_VEL,'CON_V',WIND_PRESENT,'Connection velocity (km/s)')
	RESET_REF_RADIUS=.TRUE.
	CALL RD_STORE_LOG(RESET_REF_RADIUS,'RES_REF',USE_OLD_VEL,'Reset reference radius if using old velocity law')
	GAM_LIM=0.98D0
	CALL RD_STORE_DBLE(GAM_LIM,'GAM_LIM',L_FALSE,'Limiting Eddington factor')
!
	CALL RD_STORE_DBLE(dLOG_TAU,'dLOG_TAU',L_FALSE,'Lograritihmic spacing in Tau for new R grid')
	CALL RD_STORE_DBLE(V_SCL_FAC,'VSCL_FAC',L_FALSE,'Maximum V(I-1)/V(I) for new R grid (<1)')
	I=10
	CALL RD_STORE_NCHAR(OUT_BND_OPT,'OB_OPT',I,L_FALSE,'Outer boundary option: POW, SPECIFY, DEFAULT, NONE')
	J=0; CALL RD_STORE_INT(J,'NOB_PARS',L_FALSE,'Number of outer boudary parameters')
	DO I=1,J
	  NUM_OBND_PARAMS=J
	  WRITE(STRING,'(I3)')I
	  STRING='OB_P'//ADJUSTL(STRING)
	  CALL RD_STORE_DBLE(OBND_PARS(I),TRIM(STRING),L_TRUE,'Paremeters for outer boundary condition')
	END DO  
	CALL CLEAN_RD_STORE()
C
        CLOSE(UNIT=LUIN)
        CLOSE(UNIT=LUSCR)
!
!
! Read in old model. This model will be used to help refine the
! hydrostatic structure. It used to provide corrections to the
! approximate formulae used in this program.
!
	IF(OLD_MODEL)THEN
	   PURE_LTE_EST=.FALSE.
	   ASK=.FALSE.
	   SCRAT=.FALSE.
!
10	  FILENAME=' '
	  DIR_NAME=' '
	  IOS=2			!Filename has to exits, blank not allowed.
	  WRITE(T_OUT,*)'Append (ask) to so that subsequent FILE names '//
	1         'are not defaulted'
	  WRITE(T_OUT,*)'Append (scrat) to get scratch output.'
	  WRITE(T_OUT,*)' '
	  FILENAME='RVTJ'
	  CALL GEN_IN(FILENAME,'Structure file')
	  STRING=FILENAME
	  CALL SET_CASE_UP(STRING,0,0)
	  IF( INDEX(STRING,'(SCRAT)') .NE. 0)SCRAT=.TRUE.
	  IF( INDEX(STRING,'(ASK)') .NE. 0)ASK=.TRUE.
	  I= INDEX(STRING,'(')
	  IF(I .NE. 0)FILENAME=FILENAME(1:I-1)//' '
	  INQUIRE(FILE=FILENAME,EXIST=FILE_PRES)
	  IF(.NOT. FILE_PRES)GOTO 10
C
	  WRITE(T_OUT,*)'Reading in dimensions and header from RVTJ file'
	  CALL RD_RVTJ_PARAMS_V2(OLD_MDOT,RLUM,T1,TIME,NAME_CONVENTION,
	1                     OLD_ND,NC,NP,FILENAME,LUIN)
	  ALLOCATE (OLD_R(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_V(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_SIGMA(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_T(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_ED(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_ROSS_MEAN(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_FLUX_MEAN(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_ESEC(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_KAP_ROSS(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_KAP_FLUX(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_KAP_ESEC(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_POP_ATOM(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_MASS_DENSITY(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_POPION(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_CLUMP_FAC(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_TAU(OLD_ND),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE (OLD_SF(OLD_ND),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(T_OUT,*)'Error in WIND_HYD -- error allocating atmospheric vectors'
	    WRITE(T_OUT,*)'STATUS=',IOS
	    STOP
	  END IF
	  WRITE(T_OUT,*)'Reading in vectors from RVTJ file'
	  CALL RD_RVTJ_VEC(OLD_R,OLD_V,OLD_SIGMA,OLD_ED,OLD_T,OLD_ROSS_MEAN,OLD_FLUX_MEAN,
	1         OLD_POP_ATOM,OLD_POPION,OLD_MASS_DENSITY,OLD_CLUMP_FAC,OLD_ND,LUIN)
!
! Compute mass-absorption coefficients (units cm^2/gm). These show much less variation than 
! do the conventional opacities used in CMFGEN.
!
	  OLD_ESEC(1:OLD_ND)=6.65D-15*OLD_ED(1:OLD_ND)
	  OLD_KAP_ROSS=1.0D-10*OLD_ROSS_MEAN/OLD_MASS_DENSITY
	  OLD_KAP_FLUX=1.0D-10*OLD_FLUX_MEAN/OLD_MASS_DENSITY
	  OLD_KAP_ESEC=1.0D-10*OLD_ESEC/OLD_MASS_DENSITY
!
! Compute the Rosseland optical depth scale. Clumping not yet taken into account.
!
	  WRITE(T_OUT,*)'Computing old Rosseland optical depth scale.'
	  TB(1:OLD_ND)=OLD_CLUMP_FAC(1:OLD_ND)*OLD_ROSS_MEAN(1:OLD_ND)
	  CALL TORSCL(OLD_TAU,TB,OLD_R,TA,TC,OLD_ND,'LOGMON','EXP')
	  WRITE(15,*)'OLD_TAU(1)=',OLD_TAU(1)
	  T1=2.0D0/3.0D0
	  I=GET_INDX_DP(T1,OLD_TAU,OLD_ND)
	  T2=(LOG(T1)-LOG(OLD_TAU(I)))/(LOG(OLD_TAU(I+1))-LOG(OLD_TAU(I)))
	  IF(PLANE_PARALLEL_MOD)THEN
	    OLD_REF_RADIUS=OLD_R(OLD_ND)
	    WRITE(T_OUT,*)'Reference radius of of old model is',OLD_REF_RADIUS
	  ELSE
	    OLD_REF_RADIUS=(1.0D0-T1)*OLD_R(I)+T1*OLD_R(I+1)
	    WRITE(T_OUT,*)'Reference radius (Tau=2/3) of of old model is',OLD_REF_RADIUS
	  END IF
!
	  OPEN(UNIT=23,FILE='OLD_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	    WRITE(23,'(/,A,ES14.4,/)')'Reference radius (Tau=2/3) of of old model is',OLD_REF_RADIUS
	    WRITE(23,'(A,17X,A,3(5X,A))')' Depth','R','        V','Tau(Ross)','    Gamma'
	    GAM_EDD=1.0D+06*SIGMA_TH*STEFAN_BC*(TEFF**4)/MU_ATOM/C_CMS/(10**LOGG)/AMU
	    DO I=1,OLD_ND
              T1=OLD_ED(I)/OLD_POP_ATOM(I)
	      WRITE(23,'(I6,2X,ES16.6,3ES14.4)')I,OLD_R(I),OLD_V(I),OLD_TAU(I),
	1                   T1*GAM_EDD*(OLD_FLUX_MEAN(I)/OLD_ESEC(I))
	    END DO
!
	    T1=OLD_REF_RADIUS/6.9599D0
	    OLD_TEFF=0.5770D0*(RLUM/T1/T1)**0.25           !in 10^4 K
	    WRITE(6,*)'Fudging OLD_TEFF'
	    OLD_TEFF=TEFF
!
	    OLD_SF(1:OLD_ND)=OLD_T(1:OLD_ND)/OLD_TEFF/(OLD_TAU(1:OLD_ND)+0.67D0)**0.25D0
	    OLD_TAU_MAX=OLD_TAU(OLD_ND)
!
	    WRITE(23,'(/,A,8(7X,A))')' Index','     R','   Tau','    Na',' Ne/Na','     T',
	1                  ' Kross','Kr/Kes','Kf/Kes'
	    DO I=1,OLD_ND
	      WRITE(23,'(I6,8ES13.4)')I,OLD_R(I),OLD_TAU(I),OLD_POP_ATOM(I),OLD_ED(I)/OLD_POP_ATOM(I),
	1        OLD_T(I),OLD_KAP_ROSS(I),OLD_ROSS_MEAN(I)/OLD_ESEC(I),OLD_FLUX_MEAN(I)/OLD_ESEC(I)
	    END DO
!
	    WRITE(23,'(/,A,/,A)')FORMFEED,' Old model mass absorption coeficients'
	    WRITE(23,'(A,3(8X,A))')' Index',' Kross',' Kflux','  Kes'
	    DO I=1,OLD_ND
	      WRITE(23,'(I6,3ES14.5)')I,OLD_KAP_ROSS(I),OLD_KAP_FLUX(I),OLD_KAP_ESEC(I)
	    END DO
	  CLOSE(UNIT=23)
!
	  WRITE(15,*)' '
	  WRITE(15,'(A,ES14.4)')' Old TEFF is',OLD_TEFF
	  WRITE(15,'(A,ES14.4)')' Optical depth at inner boundary (om):',OLD_TAU(OLD_ND)
	  WRITE(15,*)' '
!
	  IF(USE_OLD_VEL)THEN
	    WRITE(6,*)'            MU_ATOM from input file is',MU_ATOM
	    MU_ATOM=OLD_MASS_DENSITY(OLD_ND)/OLD_POP_ATOM(OLD_ND)/AMU
	    WRITE(6,*)'Adjusted MU_ATOM based on RVTJ file is',MU_ATOM
	    WRITE(6,*)'            Mass loss from input file is',MDOT/MASS_LOSS_SCALE_FACTOR
	    MDOT=MU_ATOM*OLD_V(OLD_ND)*OLD_POP_ATOM(OLD_ND)*OLD_R(OLD_ND)*OLD_R(OLD_ND)
	    WRITE(6,*)'Adjusted mass loss based on RVTJ file is',MDOT/MASS_LOSS_SCALE_FACTOR
	  END IF
!
	END IF			!Old model is present
!
!
!
!
	PREV_REF_RADIUS=-1.0
	DO WHILE(ABS(REFERENCE_RADIUS/PREV_REF_RADIUS-1.0) .GT. 1.0D-05)
	  WRITE(15,*)' Beginning new hydro loop'
	  WRITE(T_OUT,'(/,A)')' Beginning new hydro loop'
!
! Compute Eddington ratio. This formulae is set for one electron per ion.
! This formula holds at all radii, since g and Teff both scale as 1/r^2.
!
	  GAM_EDD=1.0D+06*SIGMA_TH*STEFAN_BC*(TEFF**4)/MU_ATOM/C_CMS/(10**LOGG)/AMU
!
	  WRITE(15,'(A,ES14.6)')'          Surface gravity is:',LOGG
	  WRITE(15,'(A,ES14.6)')'             Mass of star is:',10**(LOGG)*(REFERENCE_RADIUS**2)/GRAV_CON
	  WRITE(15,'(A,ES14.6)')'         Mean atomic mass is:',MU_ATOM
	  WRITE(15,'(A,ES14.6)')'             Atom density is:',NI_ZERO
	  WRITE(15,'(A,ES14.6)')'New effective temperature is:',TEFF
	  WRITE(15,'(A,ES14.6)')'      Eddington parameter is:',GAM_EDD
!
!
! Set parameters/initial conditions at the outer boundary of the
! hydrostatic structure.
!
! The first section assumes we have a wind present. Given this we have 3 choices:
!   (1) We have an old model and will use its wind
!   (2) We have an old model, but will input a new wind.
!   (3) We don't have an old model.
!
	  IF(WIND_PRESENT)THEN
!
! In this case will use exactly the same grid as for the old model beyond
! the connection velocity (i.e., at larger V in the wind).
!
	    IF(OLD_MODEL .AND. USE_OLD_VEL)THEN
	      DO I=1,OLD_ND
	        IF(OLD_V(I) .LT. CONNECTION_VEL)EXIT
	      END DO
	      IF( OLD_V(I+1)/CONNECTION_VEL .LT. CONNECTION_VEL/OLD_V(I))I=I+1
	      J=I
	      CONNECTION_INDX=I
	      WRITE(T_OUT,*)'Connection index is',CONNECTION_INDX
	      DO I=1,J
	        POP_ATOM(I)=OLD_POP_ATOM(I)
	        T(I)=OLD_T(I)                   !(TEFF/OLD_TEFF)*OLD_T(I)
	        P(I)=BC*T(I)*POP_ATOM(I)*(1.0D0+OLD_ED(I)/OLD_POP_ATOM(I))
	        R(I)=OLD_R(I)
	        TAU(I)=OLD_TAU(I)
	        ED_ON_NA(I)=OLD_ED(I)/OLD_POP_ATOM(I)
	        ROSS_ON_ES=OLD_ROSS_MEAN(I)/OLD_ESEC(I)
	        GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA(I)*(OLD_FLUX_MEAN(I)/OLD_ESEC(I)) )
	        ED(I)=ED_ON_NA(1)*POP_ATOM(I)
	        CHI_ROSS(I)=ED_ON_NA(I)*SIGMA_TH*POP_ATOM(I)*ROSS_ON_ES
	        GAMMA_FULL(I)=GAM_FULL
	       END DO
	       I=J
!
! In this case we simply adopt the parameters at the connection velocity to help
! define the initial conditions for integration of the hydrostatic equation.
!
! We also need to define the dp/dr so as to get the correct slope on the
! velocity law at the connection point.
!
	    ELSE IF(OLD_MODEL)THEN
	      DO I=1,OLD_ND
	        IF(OLD_V(I) .LT. CONNECTION_VEL)EXIT
	      END DO
	      IF( OLD_V(I+1)/CONNECTION_VEL .LT. CONNECTION_VEL/OLD_V(I))I=I+1
	      CONNECTION_INDX=I
	      ED_ON_NA_EST=OLD_ED(I)/OLD_POP_ATOM(I)
	      ROSS_ON_ES=OLD_ROSS_MEAN(I)/OLD_ESEC(I)
	      GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA_EST*(OLD_FLUX_MEAN(I)/OLD_ESEC(I)) )
	      T1=1.0D-10*BC*TEFF*(1+ED_ON_NA_EST)/( (10.0**LOGG)*(1.0D0-GAM_FULL)*MU_ATOM*AMU )
	      WRITE(T_OUT,*)'      Scale height is',T1
	      T2=CONNECTION_VEL*(1.0D0/T1-2.0D0/CONNECTION_RADIUS)
	      CALL WIND_VEL_LAW_V1(R,V,SIGMA,VINF,BETA,RMAX,
	1          CONNECTION_RADIUS,CONNECTION_VEL,T2,ITWO,J,ND_MAX)
	      DO I=1,J
	        POP_ATOM(I)=MDOT/MU_ATOM/R(I)/R(I)/V(I)
	      END DO
!
! To get other quantities we interpolate as a function of density.
! The atom density should be monotonic. At the outer boundary,
! we simply use the boundary value.
!
	      TA(1:OLD_ND)=LOG(OLD_POP_ATOM(1:OLD_ND))
	      TB(1:OLD_ND)=OLD_ED(1:OLD_ND)/OLD_POP_ATOM(1:OLD_ND)
	      TC(1:J)=LOG(POP_ATOM(1:J))
	      DO I=1,J
	        IF(TC(I) .LT. TA(1))TC(I)=TA(1)
	      END DO
	      CALL MON_INTERP(ED_ON_NA,J,IONE,TC,J,TB,OLD_ND,TA,OLD_ND)
	      WRITE(6,*)'MON1'
	      TB(1:OLD_ND)=OLD_T(1:OLD_ND)*TEFF/OLD_TEFF
	      CALL MON_INTERP(T,J,IONE,TC,J,TB,OLD_ND,TA,OLD_ND)
	      WRITE(6,*)'MON2'
	      TB(1:OLD_ND)=OLD_ROSS_MEAN(1:OLD_ND)/OLD_ESEC(1:OLD_ND)
	      CALL MON_INTERP(CHI_ROSS,J,IONE,TC,J,TB,OLD_ND,TA,OLD_ND)
	      WRITE(6,*)'MON3'
	      DO I=1,J
	        ED(I)=ED_ON_NA(I)*POP_ATOM(I)
	        GAMMA_FULL(I)=GAM_FULL
	        CHI_ROSS(I)=ED_ON_NA(I)*SIGMA_TH*POP_ATOM(I)*CHI_ROSS(I)
	        P(I)=BC*T(I)*POP_ATOM(I)*(1.0D0+ED(I)/POP_ATOM(I))
	      END DO
	      CALL TORSCL(TAU,CHI_ROSS,R,TB,TC,J,'LOGMON',' ')
	      I=J
	    ELSE
	      T(1)=0.75D0*TEFF
	      POP_ATOM(1)=MDOT/CONNECTION_RADIUS/CONNECTION_RADIUS/CONNECTION_VEL/MU_ATOM
	      CALL GET_LTE_ROSS_V2(KAP_ROSS,KAP_ES,T1,POP_ATOM(1),T(1))
	      ED_ON_NA(1)=T1/POP_ATOM(1)
	      P(1)=BC*T(1)*POP_ATOM(1)*(1.0D0+ED_ON_NA(1))
	      V(1)=CONNECTION_VEL
	      ROSS_ON_ES=KAP_ROSS/KAP_ES
	      GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA(1)*(KAP_ROSS/KAP_ES))
	      T1=1.0D-10*BC*TEFF*(1+ED_ON_NA_EST)/( (10.0**LOGG)*(1.0D0-GAM_FULL)*MU_ATOM*AMU )
	      WRITE(T_OUT,*)'Scale height is',T1
	      T2=CONNECTION_VEL*(1.0D0/T1-2.0D0/CONNECTION_RADIUS)
	      CALL WIND_VEL_LAW_V1(R,V,SIGMA,VINF,BETA,RMAX,
	1            CONNECTION_RADIUS,CONNECTION_VEL,T2,ITWO,J,ND_MAX)
	      DO I=1,J
	        WRITE(73,'(4ES14.4)')R(I),V(I)
	      END DO
	      DO I=1,J
	        POP_ATOM(I)=MDOT/MU_ATOM/R(I)/R(I)/V(I)
	        GAMMA_FULL(I)=GAM_FULL
	        CHI_ROSS(I)=ED_ON_NA(1)*SIGMA_TH*POP_ATOM(I)*ROSS_ON_ES
	        ED_ON_NA(I)=ED_ON_NA(1)
	        ED(I)=ED_ON_NA(I)*POP_ATOM(I)
	        T(I)=TEFF
	        P(I)=BC*T(I)*POP_ATOM(I)*(1.0D0+ED(I)/POP_ATOM(I))
	      END DO
	      DO I=1,J
	        WRITE(73,'(8ES14.4)')R(I),POP_ATOM(I),GAMMA_FULL(I),CHI_ROSS(I),ED_ON_NA(I),ED(I),T(I),P(I)
	      END DO
	      CALL TORSCL(TAU,CHI_ROSS,R,TB,TC,J,'LOGMON',' ')
	      WRITE(T_OUT,*)'Optical depth ar connection point is',TAU(J)
	      I=J
	    END IF
!
! The following section is for the case of no wind.
!
	  ELSE
!
! This option is for a pure spherical or plane-parallel model. The depth variation
! of gravity is taken into account.
! 
	    IF(OLD_MODEL)THEN
	      POP_ATOM(1)=NI_ZERO
	      T(1)=(TEFF/OLD_TEFF)*OLD_T(1)
	      P(1)=BC*T(1)*POP_ATOM(1)*(1.0D0+OLD_ED(1)/OLD_POP_ATOM(1))
	      IF(RBOUND .EQ. 0.0D0)RBOUND=OLD_R(1)
	      R(1)=RBOUND
	      TAU(1)=OLD_TAU(1)
	      ED_ON_NA(1)=OLD_ED(1)/OLD_POP_ATOM(1)
	      ROSS_ON_ES=OLD_ROSS_MEAN(1)/OLD_ESEC(1)
	      GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA(1)*(OLD_FLUX_MEAN(1)/OLD_ESEC(1)) )
	      ED(1)=ED_ON_NA(1)*POP_ATOM(1)
	      CHI_ROSS(1)=ED_ON_NA(1)*SIGMA_TH*POP_ATOM(1)*ROSS_ON_ES
	      GAMMA_FULL(1)=GAM_FULL
	    ELSE
!
! No old model. The guesses in the outer region may need improving.
!
	      WRITE(6,*)' '
	      WRITE(6,*)' Using LTE plane-parallel structure'
	      WRITE(6,*)' '
	      POP_ATOM(1)=NI_ZERO
	      T(1)=0.75D0*TEFF
	      CALL GET_LTE_ROSS_V2(KAP_ROSS,KAP_ES,T1,POP_ATOM(1),T(1))
	      ED_ON_NA(1)=T1/POP_ATOM(1)
	      P(1)=BC*T(1)*POP_ATOM(1)*(1.0D0+ED_ON_NA(1))
	      ROSS_ON_ES=KAP_ROSS/KAP_ES
	      GAM_FULL=MIN(GAM_LIM,GAM_EDD*ROSS_ON_ES*ED_ON_NA(1))
	      SCL_HT=(10**LOGG)*(1.0D0-GAM_FULL)*MU_ATOM*AMU/BC/(1.0D0+ED_ON_NA(1))/T(1)
	      SCL_HT=1.0D-10/SCL_HT
	      IF(RBOUND .EQ. 0.0D0)RBOUND=REFERENCE_RADIUS+12*SCL_HT
	      R(1)=RBOUND
	      TAU(1)=R(1)*SCL_HT*6.65D-15*POP_ATOM(1)*ED_ON_NA(1)*(KAP_ROSS/KAP_ES)
	      ED(1)=ED_ON_NA(1)*POP_ATOM(1)
	      CHI_ROSS(1)=ED_ON_NA(1)*SIGMA_TH*POP_ATOM(1)*ROSS_ON_ES
	      GAMMA_FULL(1)=GAM_FULL
	    END IF
	    I=1
	  END IF
!
	  SOUND_SPEED=1.0D+10
	  IF(WIND_PRESENT)THEN
	    SOUND_SPEED=1.0D+04*(1.0D0+ED_ON_NA(I))*BOLTZMANN_CONSTANT()*T(I)/MU_ATOM/ATOMIC_MASS_UNIT()
	    SOUND_SPEED=1.0D-05*SQRT(SOUND_SPEED)
	    WRITE(6,'(A)')' '
	    WRITE(6,'(A,3ES14.4)')' The sound speed at the iwind connection point in km/s is:',SOUND_SPEED
	    IF(SOUND_SPEED .LT. CONNECTION_VEL)THEN
	      WRITE(6,'(A)')' '
	      WRITE(6,*)'ERROR --- your connection velocity is larger than the sound speed'
	      WRITE(6,*)'The recommened connection velocity is 0.5 to 0.75 times the SOUND_SPEED'
	      WRITE(6,'(A)')' '
	      STOP
	    END IF
	  END IF
!
!
! Units 75 & 76 are used in NEW_ESTIMATES
!
	  OPEN(UNIT=75,STATUS='UNKNOWN',ACTION='WRITE',FILE='DIAGNOSTIC_EST_1')
	  WRITE(75,'(4X,9(4X,A))')'   R_EST','GAM_FULL',' GAM_EST','     TAU','      FM',
	1                         '      RM','KAP_ROSS','  KAP_ES','   ED/NA'
!
	  OPEN(UNIT=76,STATUS='UNKNOWN',ACTION='WRITE',FILE='DIAGNOSTIC_EST_2')
	  WRITE(76,'(4X,8(4X,A))')'   R_EST','  KR_OLD',' KES_OLD','   OLD_T',
	1                         ' OLD_TAU','     KES','       T','      KR'
!
	  CALL SET_LINE_BUFFERING(75)
	  CALL SET_LINE_BUFFERING(76)
!
! The boudary condition for the integration of the hydrostatic equation
! has been set, either at the outer boundary, or at the wind connection point.
! we can now perform the integration of the hydrostatic equation.
!
	  DO WHILE( TAU(I) .LT. MAX(100.0D0,OLD_TAU_MAX) )
	    I=I+1
!
! Compute the atmospheric pressure scale height. We use this to determine the
! step size. We use GAM_EDD as this provides a smaller and more consistent step
! size. It is only used to determine the step size.
!
!	    SCL_HT=(10**LOGG)*(1.0D0-GAM_FULL)*MU_ATOM*AMU/BC/(1.0D0+ED_ON_NA(I-1))/T(I-1)
	    SCL_HT=(10**LOGG)*(1.0D0-GAM_EDD)*MU_ATOM*AMU/BC/(1.0D0+ED_ON_NA(I-1))/T(I-1)
	    SCL_HT=1.0D-10/SCL_HT
!
! We set the step size to the pressure scale height on 5.
!
	    H=SCL_HT/10.0D0
!
!
	    WRITE(15,*)' '
	    WRITE(15,'(A,ES12.4,A)')'       Scale height is',SCL_HT,' 10^10 cm'
	    WRITE(15,*)'I=',I
	    WRITE(15,*)'P(I-1)=',P(I-1)
	    WRITE(15,*)'R(I-1)=',R(I-1)
	    WRITE(15,*)'T(I-1)=',T(I-1)
	    WRITE(15,*)'TAU(I-1)=',TAU(I-1)
	    WRITE(15,*)'ED_ON_NA(I-1)=',ED_ON_NA(I-1)
	    WRITE(15,*)'POP_ATOM(I-1)=',POP_ATOM(I-1)
	    WRITE(15,*)'GAMMA_FULL=',GAM_FULL
!
! Set estimates at current location. Then integrate hydrostatic
! equation using 4th order Runge-Kutta.
!
	    P_EST=P(I-1)
	    TAU_EST=TAU(I-1)
	    T_EST=T(I-1)
	    ED_ON_NA_EST=ED_ON_NA(I-1)
	    ATOM_EST=POP_ATOM(I-1)
	    WRITE(15,*)ATOM_EST,P_EST/BC/T_EST/(1+ED_ON_NA_EST)
	    R_EST=R(I-1)-H
	    CALL DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP1=H*dPdR
	    dTAU1=H*dTAUdR
	    WRITE(15,'(4ES14.5)')R_EST,P_EST,ATOM_EST,ED_ON_NA_EST
	    WRITE(15,'(5ES14.5)')dPdR,dTAUdR,T_EST,dP1,dTAU1
!
	    P_EST=P(I-1)+dP1/2
	    TAU_EST=TAU(I-1)+dTAU1/2
	    ATOM_EST=P_EST/BC/T_EST/(1+ED_ON_NA_EST)
	    R_EST=R(I-1)-0.5D0*H
	    CALL DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP2=H*dPdR
	    dTAU2=H*dTAUdR
	    WRITE(15,'(4ES14.5)')R_EST,P_EST,ATOM_EST,ED_ON_NA_EST
	    WRITE(15,'(5ES14.5)')dPdR,dTAUdR,T_EST,dP2,dTAU2
!
	    P_EST=P(I-1)+dP2/2
	    TAU_EST=TAU(I-1)+dTAU2/2
	    ATOM_EST=P_EST/BC/T_EST/(1+ED_ON_NA_EST)
	    R_EST=R(I-1)-0.5D0*H
	    CALL DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP3=H*dPdR
	    dTAU3=H*dTAUdR
	    WRITE(15,'(4ES14.5)')R_EST,P_EST,ATOM_EST,ED_ON_NA_EST
	    WRITE(15,'(5ES14.5)')dPdR,dTAUdR,T_EST,dP3,dTAU3
!
	    P_EST=P(I-1)+dP3
	    TAU_EST=TAU(I-1)+dTAU2
	    ATOM_EST=P_EST/BC/T_EST/(1+ED_ON_NA_EST)
	    R_EST=R(I-1)-H
	    CALL DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP4=H*dPdR
	    dTAU4=H*dTAUdR
	    WRITE(15,'(4ES14.5)')R_EST,P_EST,ATOM_EST,ED_ON_NA_EST
	    WRITE(15,'(5ES14.5)')dPdR,dTAUdR,T_EST,dP4,dTAU4
!
! Update values at next grid point.
!
	    TAU(I)=TAU(I-1)+(dTAU1+2*dTAU2+2*dTAU3+dTAU4)/6.0D0
	    P(I)=P(I-1)+(dP1+2*dP2+2*dP3+dP4)/6.0D0
!
	    CALL NEW_ESTIMATES(TAU(I),T_EST)
	    R(I)=R(I-1)-H
	    T(I)=T_EST
	    ED_ON_NA(I)=ED_ON_NA_EST
	    POP_ATOM(I)=P(I)/BC/T(I)/(1.0D0+ED_ON_NA(I))
	    ED(I)=ED_ON_NA(I)*POP_ATOM(I)
	    CHI_ROSS(I)=ED_ON_NA(I)*SIGMA_TH*POP_ATOM(I)*ROSS_ON_ES
	    GAMMA_FULL(I)=GAM_FULL
	    ND=I
!
	  END DO		!Loop over intire inner atmosphere
	  CLOSE(UNIT=75)
	  CLOSE(UNIT=76)
!
! Output estimates at the last depth.
!
	  I=ND
	  WRITE(15,*)'I=',I
	  WRITE(15,*)'P(I)=',P(I)
	  WRITE(15,*)'R(I)=',R(I)
	  WRITE(15,*)'T(I)=',T(I)
	  WRITE(15,*)'TAU(I)=',TAU(I)
	  WRITE(15,*)'ED_ON_NA(I)=',ED_ON_NA(I)
	  WRITE(15,*)'POP_ATOM(I)=',POP_ATOM(I)
!
! Output diagnostic files. These are on the calculate grid --- not the final
! grid.
!
	  ALLOCATE(COEF(ND,4))
	  WRITE(6,*)'H8'
	  DO I=1,ND
	  WRITE(170,*)I,P(I),R(I)
	  END DO
	  CLOSE(UNIT=170)
	  CALL MON_INT_FUNS_V2(COEF,P,R,ND)
	  WRITE(6,*)'H8'
	  DO I=1,ND
	    dPdR_VEC(I)=1.0D-10*COEF(I,3)/POP_ATOM(I)/AMU/MU_ATOM
	  END DO
	  DEALLOCATE(COEF)
!
	  INQUIRE(UNIT=18,OPENED=FILE_OPEN)
	  IF(FILE_OPEN)THEN
	    WRITE(18,'(/,A,/)')FORMFEED
	  ELSE
	    OPEN(UNIT=18,FILE='NEW_CALC_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	  END IF
	  WRITE(18,'(A,10(7X,A))')' Index','     R','  Vel','   Tau','    Na',' Ne/Na','     T',
	1                  ' Kross','Kr/Kes',' Gamma','dpdR/ROH'
	  DO I=1,ND
	    T1=1.0D+10*POP_ATOM(I)*AMU*MU_ATOM
	    T2=SIGMA_TH*POP_ATOM(I)*ED_ON_NA(I)
	    T3=MDOT/MU_ATOM/POP_ATOM(I)/R(I)/R(I)
	    WRITE(18,'(I6,10ES13.4)')I,R(I),T3,TAU(I),POP_ATOM(I),ED_ON_NA(I),
	1              T(I),CHI_ROSS(I)/T1,CHI_ROSS(I)/T2,GAMMA_FULL(I),dPdR_VEC(I)
	  END DO
!
! Adjust the grid so that we get the correct reference radius,
! defined as Tau(Ross)=2/3.
!
	  T1=2.0D0/3.0D0
	  I=GET_INDX_DP(T1,TAU,ND)
	  T2=(LOG(T1)-LOG(TAU(I)))/(LOG(TAU(I+1))-LOG(TAU(I)))
	  T2=(1.0D0-T1)*R(I)+T1*R(I+1)
	  T1=T2-REFERENCE_RADIUS
	  RADIUS_AT_TAU_23=T2
	  IF(WIND_PRESENT .AND. USE_OLD_VEL)THEN
	    PREV_REF_RADIUS=REFERENCE_RADIUS
	  ELSE IF(WIND_PRESENT)THEN
	    PREV_REF_RADIUS=T2
	    CONNECTION_RADIUS=CONNECTION_RADIUS-T1
	    WRITE(15,*)'    Old reference radius is',T2
	    WRITE(15,*)'Desired reference radius is',REFERENCE_RADIUS
	    WRITE(15,*)'Current sound speed (reference)',SOUND_SPEED
	    WRITE(T_OUT,*)'    Old reference radius is',T2
	    WRITE(T_OUT,*)'Desired reference radius is',REFERENCE_RADIUS
	  ELSE
	    PREV_REF_RADIUS=T2
	    RBOUND=RBOUND-T1
	    NI_ZERO=NI_ZERO*REFERENCE_RADIUS/PREV_REF_RADIUS
	    WRITE(15,*)'    Old reference radius is',T2
	    WRITE(15,*)'Desired reference radius is',REFERENCE_RADIUS
	    WRITE(T_OUT,*)'    Old reference radius is',T2
	    WRITE(T_OUT,*)'Desired reference radius is',REFERENCE_RADIUS
	  END IF

	END DO			!Loop to set R(Tau=2/3)=REFERENCE_RADIUS
!
	IF(WIND_PRESENT .AND. USE_OLD_VEL .AND. RESET_REF_RADIUS)THEN
	  REFERENCE_RADIUS=RADIUS_AT_TAU_23
	END IF
!
	CLOSE(UNIT=18)
	CLOSE(UNIT=19)
	WRITE(6,'(A)')' '
	CALL DP_CURVE(ND,R,POP_ATOM)
	IF(OLD_MODEL)CALL DP_CURVE(OLD_ND,OLD_R,OLD_POP_ATOM)
	CALL GRAMON_PGPLOT(' ',' ',' ',' ')
	WRITE(6,'(A)')' '
!
!
! We now create the revised grid, At present it is equally spaced in Log(tau) with
! 2 extra points at either end of the grid.
!
	NEW_ND=OLD_ND
	IF(.NOT. OLD_MODEL)NEW_ND=60
	CALL GEN_IN(NEW_ND,'Number of data points for output grid')
	ALLOCATE (REV_TAU(NEW_ND))
	ALLOCATE (REV_R(NEW_ND))
	ALLOCATE (REV_V(NEW_ND))
	ALLOCATE (REV_SIGMA(NEW_ND))
	ALLOCATE (REV_POP_ATOM(NEW_ND))
	ALLOCATE (REV_ED(NEW_ND))
	ALLOCATE (REV_T(NEW_ND))
	ALLOCATE (REV_CHI_ROSS(NEW_ND))
	ALLOCATE (REV_GAMMA_FULL(NEW_ND))
	ALLOCATE (COEF(NEW_ND,4))
!
	TAU_MAX=TAU(ND)
	IF(OLD_MODEL)THEN
	  WRITE(T_OUT,'(A,ES14.4)')'Maximum optical depth in old model was',OLD_TAU(OLD_ND)
	END IF
	CALL GEN_IN(TAU_MAX,'Maximum optical depth')
	IF(TAU_MAX .GT. TAU(ND))THEN
	  WRITE(T_OUT,*)'Error --- TAU_MAX cannot be greater than calculated grid Tau'
	  WRITE(T_OUT,*)'Setting TAU to maximum value'
	  TAU_MAX=TAU(ND)
	END IF 
!
	DO I=1,ND
	  V(I)=MDOT/MU_ATOM/POP_ATOM(I)/R(I)/R(I)
	END DO
	IF(WIND_PRESENT .AND. USE_OLD_VEL)THEN
	  J=CONNECTION_INDX
	  I=NEW_ND-J+1
	  CALL DET_R_GRID_V2(REV_TAU(J),I,ND_MAX,TAU_MAX,
	1           dLOG_TAU,V_SCL_FAC,'NONE',OBND_PARS,NUM_OBND_PARAMS,
	1           R(J),V(J),TAU(J),ND-J+1)
	  DO I=1,NEW_ND
	     WRITE(82,*)I,REV_TAU(I)
	  END DO
	ELSE
	  CALL DET_R_GRID_V2(REV_TAU,NEW_ND,ND_MAX,TAU_MAX,
	1           dLOG_TAU,V_SCL_FAC,OUT_BND_OPT,OBND_PARS,NUM_OBND_PARAMS,
	1           R,V,TAU,ND)
	END IF
!
! We now compute the revised R grid. We then interplate on Log (r^2.rho) which 
! is equivalent to interpolating on log V. This guarentees monotocity of V.
!
	WRITE(T_OUT,*)'Calling mon_interp'
	TAU(1:ND)=LOG(TAU(1:ND))
	REV_TAU(1:NEW_ND)=LOG(REV_TAU(1:NEW_ND))
	CALL MON_INTERP(REV_R,NEW_ND,IONE,REV_TAU,NEW_ND,R,ND,TAU,ND)
	POP_ATOM(1:ND)=LOG(POP_ATOM(1:ND)*R(1:ND)*R(1:ND))
	CALL MON_INTERP(REV_POP_ATOM,NEW_ND,IONE,REV_R,NEW_ND,POP_ATOM,ND,R,ND)
	REV_POP_ATOM(1:NEW_ND)=EXP(REV_POP_ATOM(1:NEW_ND))/REV_R(1:ND)/REV_R(1:ND)
	POP_ATOM(1:ND)=EXP(POP_ATOM(1:ND))/R(1:ND)/R(1:ND)
	WRITE(6,*)'Done final inteprolation of atom density'
!
! Compute revised velocity.
!
	IF(PLANE_PARALLEL_MOD)THEN
	  T1=REFERENCE_RADIUS-REV_R(NEW_ND)
	  R(1:ND)=R(1:ND)+T1
	  REV_R(1:NEW_ND)=REV_R(1:NEW_ND)+T1
	END IF
!
	IF(WIND_PRESENT .AND. USE_OLD_VEL)THEN
	  REV_V(1:CONNECTION_INDX)=OLD_V(1:CONNECTION_INDX)
	  DO I=CONNECTION_INDX+1,NEW_ND
	    REV_V(I)=MDOT/MU_ATOM/REV_POP_ATOM(I)/REV_R(I)/REV_R(I)
	  END DO
	ELSE
	  DO I=1,NEW_ND
	    REV_V(I)=MDOT/MU_ATOM/REV_POP_ATOM(I)/REV_R(I)/REV_R(I)
	  END DO
	END IF
!
! Compute SIGMA by performing a monotonic cubic fit to V as a function of R.
!
	CALL MON_INT_FUNS_V2(COEF,REV_V,REV_R,NEW_ND)
	DO I=1,NEW_ND
	  REV_SIGMA(I)=REV_R(I)*COEF(I,3)/REV_V(I)-1.0D0
	END DO
!
!
! Output revised hydrostatic structure.
!
	OPEN(UNIT=40,FILE='RVSIG_COL_NEW',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(40,'(A)')'!'
	  WRITE(40,'(A)')'! Note: The effective temperature and surface gravity are defined'
	  WRITE(40,'(A)')'! at the reference radius, which (except when using the old'
	  WRITE(40,'(A)')'! velocity) is the location where Tau=2/3.'
	  WRITE(40,'(A)')'!'
	  WRITE(40,'(A,ES16.6)')'! Effective temperature (10^4 K) is:',TEFF
	  WRITE(40,'(A,ES16.6)')'!      Log surface gravity (cgs) is:',LOGG
	  WRITE(40,'(A,ES16.6)')'!         Core radius (10^10 cm) is:',REV_R(NEW_ND)
	  WRITE(40,'(A,ES16.6)')'!    Reference radius (10^10 cm) is:',REFERENCE_RADIUS
	  WRITE(40,'(A,ES16.6)')'!              Luminosity (Lsun) is:',( (TEFF/0.5770D0)**4 )*( (REFERENCE_RADIUS/6.9599)**2 )
	  WRITE(40,'(A,ES16.6)')'!            Mass (Msun) of star is:',10**(LOGG)*REFERENCE_RADIUS*REFERENCE_RADIUS/GRAV_CON
	  WRITE(40,'(A,ES16.6)')'!       Mass loss rate (Msun/yr) is:',MDOT/MASS_LOSS_SCALE_FACTOR
	  WRITE(40,'(A,ES16.6)')'!         Mean atomic mass (amu) is:',MU_ATOM
	  WRITE(40,'(A,ES16.6)')'!            Eddington parameter is:',GAM_EDD
	  WRITE(40,'(A,ES16.6)')'!                   Atom density is:',NI_ZERO
	  WRITE(40,'(A,F14.8)') '! Ratio of inner to outer radius is:',REV_R(1)/REV_R(NEW_ND)
	  WRITE(40,'(A)')'!'
	  WRITE(40,'(3X,I5,10X,A)')NEW_ND,'!Number of depth points'
	  WRITE(40,'(A)')'!'
	  WRITE(40,'(A,4X,A,3(7X,A),3X,A)')'!','R(10^10cm)','V(km/s)','  Sigma','    Tau','  Index'
	  DO I=1,NEW_ND
	    WRITE(40,'(F15.8,3ES14.6,6X,I4)')REV_R(I),REV_V(I),REV_SIGMA(I),EXP(REV_TAU(I)),I
	  END DO
	CLOSE(UNIT=40)
!
!
! Output estimate data for comparison with new model data. R is used as the
! dependent variable.
!
	CALL MON_INTERP(REV_T,NEW_ND,IONE,REV_R,NEW_ND,T,ND,R,ND)
	CALL MON_INTERP(REV_ED,NEW_ND,IONE,REV_R,NEW_ND,ED,ND,R,ND)
	CALL MON_INTERP(REV_CHI_ROSS,NEW_ND,IONE,REV_R,NEW_ND,CHI_ROSS,ND,R,ND)
	CALL MON_INTERP(REV_GAMMA_FULL,NEW_ND,IONE,REV_R,NEW_ND,GAMMA_FULL,ND,R,ND)
!
	OPEN(UNIT=41,FILE='FIN_CAL_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(41,'(A,9(6X,A7))')'Index','     R','   TAU','    V','     T','    Na',
	1              ' Ne/Na',' Xross','Xr/Xes','   Gam'
	  DO I=1,NEW_ND
	    WRITE(41,'(I5,9ES13.4)')I,REV_R(I),EXP(REV_TAU(I)),REV_V(I),REV_T(I),
	1             REV_POP_ATOM(I),REV_ED(I)/REV_POP_ATOM(I),REV_CHI_ROSS(I),
	1             REV_CHI_ROSS(I)/SIGMA_TH/REV_ED(I),REV_GAMMA_FULL(I)
	  END DO
	CLOSE(UNIT=41)
!
	STOP
	END 
!	
! Subroutine to compute dPdR and dTAUdR for use with the
! Runge-Kutta integration.
!
	SUBROUTINE DERIVS(P,TAU,TEMP,ED_ON_NA,POP_ATOM)
	USE HYDRO_PARAM_MODULE
	IMPLICIT NONE
!
	REAL*8 P
	REAL*8 TAU
	REAL*8 TEMP
	REAL*8 ED_ON_NA
	REAL*8 POP_ATOM
	REAL*8 T1,T2,T3
!
	CALL NEW_ESTIMATES(TAU,TEMP)
!
	IF(PLANE_PARALLEL_MOD)THEN
	  T1=(10.0**LOGG)*(1.0D0-GAM_FULL)
	ELSE
	  T1=(10.0**LOGG)*(1.0D0-GAM_FULL)*(REFERENCE_RADIUS/R_EST)**2
	END IF
	T2=MDOT/MU_ATOM/POP_ATOM/R_EST/R_EST
	T3=(T2/SOUND_SPEED)**2
	dPdR=1.0D+10*T1*MU_ATOM*AMU*P/BC/TEMP/(1+ED_ON_NA)/(1.0D0-T3)
	dTAUdR=ED_ON_NA*POP_ATOM*SIGMA_TH*ROSS_ON_ES
!
	WRITE(17,'(A,ES14.4)')'      P_EST=',P
	WRITE(17,'(A,ES14.4)')'      T_EST=',TEMP
	WRITE(17,'(A,ES14.4)')'   ATOM_EST=',POP_ATOM
	WRITE(17,'(A,ES14.4)')'   ED_ON_NA=',ED_ON_NA
	WRITE(17,'(A,ES14.4)')'       dPdR=',dPdR
	WRITE(17,'(A,ES14.4)')'     dTAUdR=',dTAUdR
	WRITE(17,'(A,ES14.4)')'   GAM_FULL=',GAM_FULL
	WRITE(17,'(A,ES14.4)')' ROSS_ON_ES=',ROSS_ON_ES
!
	RETURN
	END
!
! Get new estimates of the model parameters. At present we
! can choose pure LTE estimates, or scaled LTE estimates.
!
	SUBROUTINE NEW_ESTIMATES(TAU,TEMP)
	USE HYDRO_PARAM_MODULE
	IMPLICIT NONE
!
	REAL*8 TEMP
	REAL*8 TAU
!
	REAL*8 T1
	REAL*8 SF
	REAL*8 FM
	REAL*8 ES
	REAL*8 RM
	REAL*8 LTE_ED
!
	REAL*8 ED_EST
	REAL*8 OLD_ATOM
	REAL*8 OLD_TEMP
	REAL*8 KAP_ROSS_OLD
	REAL*8 KAP_ES_OLD
!
	INTEGER, SAVE ::INDX=1
!
	IF(PURE_LTE_EST)THEN
	  TEMP=TEFF*(TAU+0.67D0)**0.25D0
	  CALL GET_LTE_ROSS_V2(KAP_ROSS,KAP_ES,LTE_ED,ATOM_EST,TEMP)
	  ED_ON_NA_EST=LTE_ED/ATOM_EST
	  GAM_FULL=GAM_EDD*(KAP_ROSS/KAP_ES)*ED_ON_NA_EST
	  ROSS_ON_ES=KAP_ROSS/KAP_ES
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	  RETURN
	END IF
!
! Determine parameters at current grid point. Special allowance
! has to be made because of the boundaries. As we are accesing
! the grid sequentially we can use INDX as the starting location
! of the search. We just need to check that we are not starting
! a new sequence.
!
	IF(TAU .LT. OLD_TAU(INDX))INDX=1
	DO WHILE(TAU .GT. OLD_TAU(INDX) .AND. INDX .LT. OLD_ND)
	  INDX=INDX+1
	END DO
	IF(INDX .EQ. 1)THEN
	  TEMP=TEMP				!Ue passed value !OLD_SF(1)*TEFF
	  ED_ON_NA_EST=OLD_ED(1)/OLD_POP_ATOM(1)
	  ROSS_ON_ES=OLD_KAP_ROSS(1)/OLD_KAP_ESEC(1)
	  GAM_FULL=GAM_EDD*(OLD_KAP_FLUX(1)/OLD_KAP_ESEC(1))*ED_ON_NA_EST
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	ELSE IF(TAU .GE. OLD_TAU(OLD_ND))THEN
	  TEMP=TEFF*OLD_SF(OLD_ND)*(TAU+0.67D0)**0.25D0
	  ED_ON_NA_EST=OLD_ED(OLD_ND)/OLD_POP_ATOM(OLD_ND)
	  ROSS_ON_ES=OLD_KAP_ROSS(OLD_ND)/OLD_KAP_ESEC(OLD_ND)
	  GAM_FULL=GAM_EDD*(OLD_KAP_FLUX(OLD_ND)/OLD_KAP_ESEC(OLD_ND))*ED_ON_NA_EST
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	ELSE
!
! Determine parameters at current TAU in old model.
!
	  T1=(TAU-OLD_TAU(INDX-1))/(OLD_TAU(INDX)-OLD_TAU(INDX-1))
	  SF=(1.0D0-T1)*OLD_SF(INDX-1)+T1*OLD_SF(INDX)
	  TEMP=SF*TEFF*(TAU+0.67D0)**0.25D0
	  FM=(1.0D0-T1)*OLD_KAP_FLUX(INDX-1)+T1*OLD_KAP_FLUX(INDX)
	  RM=(1.0D0-T1)*OLD_KAP_ROSS(INDX-1)+T1*OLD_KAP_ROSS(INDX)
	  ES=(1.0D0-T1)*OLD_KAP_ESEC(INDX-1)+T1*OLD_KAP_ESEC(INDX)
!
! Get parameters at interpolation point in TAU space. We use these to
! estimate Ne/N and the ratio of K(non-LTE)/K(LTE).
!
          ED_EST=(1.0D0-T1)*OLD_ED(INDX-1)+T1*OLD_ED(INDX)
          OLD_ATOM= (1.0D0-T1)*OLD_POP_ATOM(INDX-1)+T1*OLD_POP_ATOM(INDX)
          OLD_TEMP= (1.0D0-T1)*OLD_T(INDX-1)+T1*OLD_T(INDX)
          ED_ON_NA_EST=ED_EST/OLD_ATOM
!
! Some fiddling may be required here to choose the optimal density for
! switching.
!
	  IF(ATOM_EST .GT. 1.0D+09)THEN
	    CALL GET_LTE_ROSS_V2(KAP_ROSS_OLD,KAP_ES_OLD,LTE_ED,OLD_ATOM,OLD_TEMP)
	    CALL GET_LTE_ROSS_V2(T1,KAP_ES,LTE_ED,ATOM_EST,TEMP)
            KAP_ROSS=(T1-KAP_ES)*( (RM-ES)/(KAP_ROSS_OLD-KAP_ES_OLD))+KAP_ES
!            KAP_ROSS=(T1-KAP_ES)*(RM/KAP_ROSS_OLD)+KAP_ES
!            KAP_ROSS=T1*(RM/KAP_ROSS_OLD)
!	     GAM_FULL=GAM_EDD*(FM/RM)*(KAP_ROSS/KAP_ES)*ED_ON_NA_EST
	     KAP_ROSS=RM; KAP_ES=ES
	     GAM_FULL=GAM_EDD*(FM/RM)*(KAP_ROSS/KAP_ES)*ED_ON_NA_EST
	    ROSS_ON_ES=KAP_ROSS/KAP_ES
	   ELSE
	    GAM_FULL=GAM_EDD*(FM/ES)*ED_ON_NA_EST
	    ROSS_ON_ES=RM/ES
	   END IF
	   GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	   WRITE(76,'(ES16.8,8ES12.3)')R_EST,KAP_ROSS_OLD,KAP_ES_OLD,OLD_TEMP,T1,KAP_ES,TEMP,KAP_ROSS
	   WRITE(75,'(ES16.9,8ES12.3)')R_EST,GAM_FULL,ED_ON_NA_EST*GAM_EDD*FM/ES,TAU,FM,RM,KAP_ROSS,KAP_ES,ED_ON_NA_EST
	END IF
!  
	RETURN
	END 
