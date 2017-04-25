!
! Subroutine to calculate an updated hydrostatic structure. Routine is to be called from
! CMFGEN.
!
! Subroutine allows for:
!                      (a) A correction for the vdv/dr dynamical term below the sonic point.
!                      (b) A contributin by turbulent pressure.
!
	SUBROUTINE DO_CMF_HYDRO_V2(POPS,MOD_LUM,MOD_TEFF,MOD_LOGG,MOD_MASS,MOD_RSTAR,
	1		  MOD_RMAX,MOD_MDOT,MOD_VINF,MOD_BETA,
	1                 MOD_VTURB,PLANE_PARALLEL,PLANE_PARALLEL_NO_V,
	1                 MAIN_COUNTER,DONE_HYDRO,MOD_NC,MOD_ND,MOD_NP,NT)
	USE CMF_HYDRO_MODULE
	USE OLD_GRID_MODULE
	USE UPDATE_KEYWORD_INTERFACE
	IMPLICIT NONE
!
! Altered 05-Jun-2015 - Fixed bug; GAM_LIM_STORE was not being set to GAM_LIM when it was read
!                           in from the HDYRO_DEFAULTS file.
!                       Now save old RVSIG_COL file as RVSIG_COL_IT_#.
! Altered 30-Jul-2011 - Added VC_ON_SS as parameter
! Altered 17-Jun-2011 - Added check to make sure dVdR is not negative at the connection point.
!                         When -ve, we lower GAM_LIM.
! Altered 31-Mar-2011 - Added variable NO_ITS_DONE and corresponding KEYWORD to file.
! Altered 03-Mar-2011 - Revised check on GAM_EDD to properly account for ionization state of gas.
!                         Program value of GAM_EDD remains unchanged.
!                         Extra parameter BETA2 can now be read in to fiddle with velocity law.
!                         WIND_VEL_LAW_V2 is now called.
! Altered 31-Aug-2010 : TAU_REF can be a parameter (default is 2/3). For W-R stars.
! Altered 03-Aug-2010 : Match velocity at 0.75 x sound_speed (old vale was 0.5).
! Altered 18-May-2008 : Insert a limit as to the number of iterations (ITERATION_COUNT).
! Altered 11-May-2008 : Inserted MOD_RSTAR and changed to _V2.
! Created 15-Feb-2008 :
!
	INTEGER MOD_NC			!Model numbr of core rays
	INTEGER MOD_ND			!Model number of depth points
	INTEGER MOD_NP			!Model number of impact papramters
	INTEGER NT			!Number of levels
!
	REAL*8 POPS(NT,MOD_ND)
!
	REAL*8 MOD_MDOT			!Mass loss rate (N=MOD_MDOT/V(kms)/r(10^10cm)^^2)
	REAL*8 MOD_LUM			!In Lsun
	REAL*8 MOD_TEFF			!In unitos of 10^4 K
	REAL*8 MOD_LOGG			!In cgs units
	REAL*8 MOD_MASS			!Mass of star in Msun (returned and output to VADAT)
	REAL*8 MOD_RSTAR		!Returned
	REAL*8 MOD_RMAX			!R(ND) on input
	REAL*8 MOD_VINF			!km/s
	REAL*8 MOD_BETA			!Classic velocity law exponent
	REAL*8 MOD_VTURB
!
	INTEGER MAIN_COUNTER
	LOGICAL PLANE_PARALLEL
	LOGICAL PLANE_PARALLEL_NO_V
	LOGICAL DONE_HYDRO
!
! The following vectors are used for the atmospheric structure resulting
! from the solution of the hydrostatic and tau equations.
! 
	INTEGER, PARAMETER :: ND_MAX=4000
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
! Parameters for for cumputing the final R grid.
!       
	REAL*8 dLOG_TAU
	REAL*8 V_SCL_FAC
	REAL*8 OBND_PARS(20)
	INTEGER NUM_OBND_PARAMS
	CHARACTER(LEN=16) OUT_BND_OPT
	CHARACTER(LEN=80) STRING
	CHARACTER(LEN=20) HYDRO_OPT
!
! Wind parameters:
!
	REAL*8 VINF
	REAL*8 BETA
	REAL*8 BETA2
	REAL*8 RMAX
	REAL*8 CONNECTION_VEL
	REAL*8 CONNECTION_RADIUS
	INTEGER CONNECTION_INDX
!
	INTEGER ND			!Number of points in RUnge-Kutta integration
	INTEGER I,J,IOS
!
	REAL*8 PI			!
	REAL*8 NI_ZERO			!Density at outer boundary
	REAL*8 SCL_HT			!Atmopsheric scale height
	REAL*8 H			!Step size for Runge-Kutta integration
	REAL*8 PTURB_ON_NA		!Turbulent pressure due to turbulence/ (# of atoms)
	REAL*8 MASS_LOSS_SCALE_FACTOR	!Factor to convert MOD_MDOt to Msun/yr
	REAL*8 RBOUND
	REAL*8 TAU_MAX
	REAL*8 OLD_TAU_MAX
	REAL*8 BOLD,BNEW
	REAL*8 T1,T2,T3
	REAL*8 GAM_LIM_STORE
	REAL*8 VC_ON_SS
	REAL*8 TAU_REF
	REAL*8 MAX_ED_ON_NA
!
! Runge-Kutta estimates
!
	REAL*8 dP1,dTAU1
	REAL*8 dP2,dTAU2
	REAL*8 dP3,dTAU3
	REAL*8 dP4,dTAU4
!
	LOGICAL USE_OLD_VEL
	LOGICAL L_TEMP
	LOGICAL FILE_OPEN
	LOGICAL FILE_PRES
	LOGICAL VERBOSE_OUTPUT
	LOGICAL UPDATE_GREY_SCL
!
! These have cgs units.
!
	REAL*8 ATOMIC_MASS_UNIT,BOLTZMANN_CONSTANT,STEFAN_BOLTZ
	REAL*8 GRAVITATIONAL_CONSTANT,MASS_SUN,SPEED_OF_LIGHT,LUM_SUN
	INTEGER GET_INDX_DP,ERROR_LU
	EXTERNAL ATOMIC_MASS_UNIT,BOLTZMANN_CONSTANT,STEFAN_BOLTZ,GET_INDX_DP
	EXTERNAL GRAVITATIONAL_CONSTANT,MASS_SUN,SPEED_OF_LIGHT,ERROR_LU,LUM_SUN
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
	INTEGER, PARAMETER :: LUIN=7
	INTEGER, PARAMETER :: LUSCR=8
	INTEGER  STRT_HYDRO_ITS
	INTEGER  FREQ_HYDRO_ITS
	INTEGER  NO_HYDRO_ITS
	INTEGER  NO_ITS_DONE
	INTEGER  ITERATION_COUNTER
	INTEGER  VEL_LAW
!
	INTEGER  LU
	INTEGER  LUV
	INTEGER  LU_ERR
!
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	CHARACTER*2, PARAMETER :: FORMFEED=' '//CHAR(12)
!
! Set constants.
!
	LU_ERR=ERROR_LU()
	BC=1.0D+04*BOLTZMANN_CONSTANT()   			!erg/10^4 K
	AMU=ATOMIC_MASS_UNIT()          			!gm
	C_CMS=SPEED_OF_LIGHT()   	       			!cm/s
	GRAV_CON=1.0D-20*GRAVITATIONAL_CONSTANT()*MASS_SUN()
	SIGMA_TH=6.65D-15					!cm^{-2} x 10^10
	STEFAN_BC=STEFAN_BOLTZ()
	MASS_LOSS_SCALE_FACTOR=3.02286D+23
	PI=ACOS(-1.0D0)
	STRT_HYDRO_ITS=20
	FREQ_HYDRO_ITS=8
	HYDRO_OPT='DEFAULT'
!
	CALL STORE_OLD_GRID(MU_ATOM,MOD_ND)
!
	TEFF=MOD_TEFF
	LOGG=MOD_LOGG
	VINF=MOD_VINF
	BETA=MOD_BETA
	MDOT=MOD_MDOT
	RMAX=MOD_RMAX/OLD_R(OLD_ND)
	VTURB=MOD_VTURB
	VEL_LAW=ITWO
!
	VC_ON_SS=0.75D0
	TAU_REF=2.0D0/3.0D0
        dLOG_TAU=0.25D0
        V_SCL_FAC=0.75D00
        OBND_PARS(:)=0.0D0
        NUM_OBND_PARAMS=1
        OUT_BND_OPT='DEFAULT'
	UPDATE_GREY_SCL=.FALSE.
!
! Set defaults:
!
	IF(PLANE_PARALLEL)THEN
	  NI_ZERO=1.0D+06
	  PLANE_PARALLEL_MOD=.TRUE.
	  WIND_PRESENT=.TRUE.
	  RESET_REF_RADIUS=.FALSE.
	ELSE IF(PLANE_PARALLEL_NO_V)THEN
	  NI_ZERO=1.0D+06
	  PLANE_PARALLEL_MOD=.TRUE.
	  WIND_PRESENT=.FALSE.
	  RESET_REF_RADIUS=.FALSE.
	ELSE
	  PLANE_PARALLEL_MOD=.FALSE.
	  WIND_PRESENT=.TRUE.
	  RESET_REF_RADIUS=.TRUE.
	END IF
	USE_OLD_VEL=.FALSE.
	GAM_LIM=0.98D0
	GAM_LIM_STORE=GAM_LIM
!
! These are the default settings if no old model.
!
	PURE_LTE_EST=.FALSE.
	OLD_TAU_MAX=100.0D0
	RBOUND=0.0D0
	VERBOSE_OUTPUT=.TRUE.
!
! *************************************************************************
!
! Read in parameters describing the new model.
!
	CALL GEN_ASCI_OPEN(LUIN,'HYDRO_DEFAULTS','OLD',' ','READ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LU_ERR,*)'Error opening HYDRO_DEFAULTS in WIND_HYD, IOS=',IOS
	  RETURN 
	END IF
	CALL RD_OPTIONS_INTO_STORE(LUIN,LUSCR)

        CALL RD_STORE_INT(NO_HYDRO_ITS,'N_ITS',L_TRUE,'Number of hydro iterations remaining')
        CALL RD_STORE_INT(NO_ITS_DONE,'ITS_DONE',L_TRUE,'Number of hydro iterations completed')
        CALL RD_STORE_INT(STRT_HYDRO_ITS,'STRT_ITS',L_FALSE,'Iteration to start first hydro iteration')
        CALL RD_STORE_INT(FREQ_HYDRO_ITS,'FREQ_ITS',L_FALSE,'Frequency for hydro iterations')
	CALL RD_STORE_DBLE(NI_ZERO,'ATOM_DEN',L_FALSE,'Atom density at outer boundary (/cm^3)')
	CALL RD_STORE_LOG(USE_OLD_VEL,'OLD_V',L_FALSE,'Use old velocity law above connection point?')
	CALL RD_STORE_DBLE(RMAX,'MAX_R',L_FALSE,'Maximum radius in terms of Connection radius')
	CALL RD_STORE_LOG(RESET_REF_RADIUS,'RES_REF',L_FALSE,'Reset reference radius if using old velocity law')
	CALL RD_STORE_DBLE(GAM_LIM,'GAM_LIM',L_FALSE,'Limiting Eddington factor')
	GAM_LIM_STORE=GAM_LIM
	CALL RD_STORE_DBLE(VC_ON_SS,'VC_ON_SS',L_FALSE,'Connection velocity on sound speed')
	CALL RD_STORE_LOG(UPDATE_GREY_SCL,'UP_GREY_SCL',L_FALSE,'Update GREY_SCL_FAC_IN')
	CALL RD_STORE_DBLE(TAU_REF,'TAU_REF',L_FALSE,'Reference radius for g and Teff')
	CALL RD_STORE_INT(VEL_LAW,'VEL_LAW',L_FALSE,'Velocity law for wind region (2 or 3)')
	BETA2=BETA
	CALL RD_STORE_DBLE(BETA2,'BETA2',L_FALSE,'Second exponent for velocity law')
!
! Therse are the parameters used to define the new R grid to be output to RVSIG_COL.
! 
	CALL RD_STORE_DBLE(dLOG_TAU,'dLOG_TAU',L_FALSE,'Logarithmic spacing in Tau for new R grid')
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
	CALL RD_STORE_CHAR(HYDRO_OPT,'HYDRO_OPT',L_FALSE,'FIXED_R or FIXED_V_FLUX or DEFAULT')
	IF(HYDRO_OPT .EQ. 'FIXED_V_FLUX')THEN
	  CALL RD_STORE_DBLE(OLD_TEFF,'OLD_TEFF',L_TRUE,'Effective temperatre of input model')
	END IF
	CALL CLEAN_RD_STORE()
!
        CLOSE(UNIT=LUIN)
        CLOSE(UNIT=LUSCR)
!
! Decide here whether we will do an iteration or not.
!
	DONE_HYDRO=.FALSE.
	IF(NO_HYDRO_ITS .EQ. 0)RETURN
	IF(MAIN_COUNTER .LT. STRT_HYDRO_ITS)RETURN
	IF( MOD( (MAIN_COUNTER-STRT_HYDRO_ITS),FREQ_HYDRO_ITS ) .NE. 0)RETURN 
!
! Begin Hydro computation.
!
	IF(VERBOSE_OUTPUT)THEN
	  CALL GET_LU(LUV)
	  OPEN(UNIT=LUV,FILE='HYDRO_ITERATION_INFO',STATUS='UNKNOWN')
	  CALL SET_LINE_BUFFERING(LUV)
	END IF
	CALL GET_LU(LU)				!For files open/shut immediately
!
! In TORSCL_V3, TA is TAU, TB is dTAU, and TC is used fro dCHIdR. 
!
	WRITE(6,'(/,A)')' Updating hydrostatic structure of the model'
	IF(HYDRO_OPT .EQ. 'FIXED_R_REF')THEN
	  WRITE(6,'(A)')' Using FIXED_R_REF option in DO_CMF_HYDRO_V2'
	  CHI_ROSS(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	  I=7			!Use CHI(1) and CHI(I) to computed exponent.
	  CALL TORSCL_V3(TA,CHI_ROSS,OLD_R,TB,TC,MOD_ND,'LOGMON','PCOMP',I,L_FALSE)
	  DO I=1,MOD_ND
	    IF(TA(I) .GT. TAU_REF)THEN
	      T1=(TAU_REF-TA(I-1))/(TA(I)-TA(I-1))
	      REFERENCE_RADIUS=T1*OLD_R(I)+(1.0D0-T1)*OLD_R(I-1)
	      EXIT
	    END IF
	  END DO
!
	  T1=REFERENCE_RADIUS*TEFF*TEFF
	  MOD_LUM=4.0D+36*PI*STEFAN_BC*T1*T1/LUM_SUN()
	  CALL UPDATE_KEYWORD(MOD_LUM,'[LSTAR]','VADAT',L_TRUE,L_TRUE,LUIN)
	  CALL UPDATE_KEYWORD('DEFAULT','[HYDRO_OPT]','HYDRO_DEFAULTS',L_TRUE,L_TRUE,LUIN)
	  WRITE(6,'(A)')' DO_CMF_HYDRO_V2 has adjusted LSTAR in VADAT'
!
	ELSE IF(HYDRO_OPT .EQ. 'FIXED_V_FLUX')THEN
	  WRITE(6,'(A)')' Using FIXED_V_FLUX option in DO_CMF_HYDRO_V2'
	  I=7	!Use CHI(1) and CHI(I) to computed exponent.
	  CHI_ROSS(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	  CALL TORSCL_V3(TA,CHI_ROSS,OLD_R,TB,TC,MOD_ND,'LOGMON','PCOMP',I,L_FALSE)
	  DO I=1,MOD_ND
	    IF(TA(I) .GT. TAU_REF)THEN
	      T1=(TAU_REF-TA(I-1))/(TA(I)-TA(I-1))
	      REFERENCE_RADIUS=T1*OLD_R(I)+(1.0D0-T1)*OLD_R(I-1)
	      EXIT
	    END IF
	  END DO
          BOLD=1.0D0/(EXP(1.4388/0.55D0/OLD_TEFF)-1.0D0)
          BNEW=1.0D0/(EXP(1.4388/0.55D0/TEFF)-1.0D0)
	  REFERENCE_RADIUS=REFERENCE_RADIUS*SQRT(BOLD/BNEW)
	  T1=REFERENCE_RADIUS*TEFF*TEFF
	  MOD_LUM=4.0D+36*PI*STEFAN_BC*T1*T1/LUM_SUN()
	  CALL UPDATE_KEYWORD(MOD_LUM,'[LSTAR]','VADAT',L_TRUE,L_TRUE,LUIN)
	  CALL UPDATE_KEYWORD('DEFAULT','[HYDRO_OPT]','HYDRO_DEFAULTS',L_TRUE,L_FALSE,LUIN)
	  CALL UPDATE_KEYWORD(TEFF,'[OLD_TEFF]','HYDRO_DEFAULTS',L_FALSE,L_TRUE,LUIN)
	  WRITE(6,'(A)')' DO_CMF_HYDRO_V2 has adjusted LSTAR in VADAT'
!
	ELSE IF(HYDRO_OPT .EQ. 'DEFAULT')THEN
	  WRITE(6,'(A)')' Reference radius: based on effective temperature and luminosity of star'
	  REFERENCE_RADIUS=1.0D-18*SQRT(MOD_LUM*LUM_SUN()/TEFF**4/STEFAN_BC/4.0D0/PI)
	ELSE
	  WRITE(6,'(A)')' Error in do_cmf_hydro_v3.f: invlaid HYDRO_OPT option'
	  WRITE(6,'(2A)')' HYDRO_OPT=',TRIM(HYDRO_OPT)
	  STOP
	END IF
	WRITE(6,*)'Reference radius is',REFERENCE_RADIUS
!
!
!
! Compute the grey temperature distribution (returned in TC) and the Rosseland
! optical depth scale (returned in TA).
!
	TB(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	CALL COMP_GREY_V2(POPS,TC,TA,TB,LU_ERR,MOD_NC,MOD_ND,MOD_NP,NT)
	IF(MINVAL(OLD_ROSS_MEAN(1:MOD_ND)) .LE. 0.0D0)THEN
	   WRITE(LU_ERR,'(A)')' Bad Rosseland optical depth scale for T/TGREY output'
	   WRITE(LU_ERR,'(A)')' GREY_SCL_FAC_OUT not output'
	ELSE IF(UPDATE_GREY_SCL)THEN
	  OPEN(UNIT=LUIN,FILE='GREY_SCL_FAC_IN',STATUS='UNKNOWN')
	    WRITE(LUIN,'(A)')'!'
	    WRITE(LUIN,'(A,8X,A,7X,A,7X,A,6X,A)')'!','Log(Tau)','T/T(grey)','T(10^4 K)','L'
	    WRITE(LUIN,'(A)')'!'
	    WRITE(LUIN,*)MOD_ND
	    DO I=1,MOD_ND
	      WRITE(LUIN,'(2X,3ES16.6,4X,I3)')LOG10(TA(I)),OLD_T(I)/TC(I),OLD_T(I),I
	    END DO
	  CLOSE(LUIN)
	END IF
!
! Compute the Rosseland optical depth scale.
!
	TB(1:MOD_ND)=OLD_CLUMP_FAC(1:MOD_ND)*OLD_ROSS_MEAN(1:MOD_ND)
	CALL TORSCL(OLD_TAU,TB,OLD_R,TA,TC,MOD_ND,'LOGMON','EXP')
	IF(PLANE_PARALLEL_MOD)THEN
	  OLD_REF_RADIUS=OLD_R(OLD_ND)
	ELSE
	  T1=TAU_REF
	  I=GET_INDX_DP(T1,OLD_TAU,MOD_ND)
	  T2=(LOG(T1)-LOG(OLD_TAU(I)))/(LOG(OLD_TAU(I+1))-LOG(OLD_TAU(I)))
	  OLD_REF_RADIUS=(1.0D0-T1)*OLD_R(I)+T1*OLD_R(I+1)
	  IF(VERBOSE_OUTPUT)THEN
	    WRITE(LUV,*)'OLD_TAU(1)=',OLD_TAU(1)
	    WRITE(LUV,*)'Reference radius (Tau=',TAU_REF,') of of old model is',OLD_REF_RADIUS
	  END IF
	END IF
!
	T1=OLD_REF_RADIUS/6.9599D0
	OLD_TEFF=TEFF
	OLD_SF(1:MOD_ND)=OLD_T(1:MOD_ND)/OLD_TEFF/(OLD_TAU(1:MOD_ND)+0.67D0)**0.25D0
	OLD_TAU_MAX=OLD_TAU(MOD_ND)
!
	IF(VERBOSE_OUTPUT)THEN
	  OPEN(UNIT=LU,FILE='HYDRO_OLD_MODEL',STATUS='UNKNOWN')
	  WRITE(LU,'(/,A,8(7X,A))')' Index','     R','   Tau','    Na',' Ne/Na','     T',
	1                  ' Kross','Kr/Kes','Kf/Kes'
	  DO I=1,MOD_ND
	    WRITE(LU,'(I6,8ES13.4)')I,OLD_R(I),OLD_TAU(I),OLD_POP_ATOM(I),OLD_ED(I)/OLD_POP_ATOM(I),
	1        OLD_T(I),OLD_KAP_ROSS(I),OLD_ROSS_MEAN(I)/OLD_ESEC(I),OLD_FLUX_MEAN(I)/OLD_ESEC(I)
	  END DO
!
	  WRITE(LU,'(/,A,/,A)')FORMFEED,' Old model mass absorption coeficients'
	  WRITE(LU,'(A,3(8X,A))')' Index',' Kross',' Kflux','  Kes'
	  DO I=1,MOD_ND
	    WRITE(LU,'(I6,3ES14.5)')I,OLD_KAP_ROSS(I),OLD_KAP_FLUX(I),OLD_KAP_ESEC(I)
	  END DO
	  CLOSE(UNIT=LU)
!
	  WRITE(LUV,*)' '
	  WRITE(LUV,'(A,ES14.4)')' Old TEFF is',OLD_TEFF
	  WRITE(LUV,'(A,ES14.4)')' Optical depth at inner boundary (om):',OLD_TAU(MOD_ND)
	  WRITE(LUV,*)' '
	END IF
!
	IF(WIND_PRESENT)THEN
	  T1=1.0D-06*BOLTZMANN_CONSTANT()/MU_ATOM/ATOMIC_MASS_UNIT()
	  DO I=1,MOD_ND
	    ED_ON_NA_EST=OLD_ED(I)/OLD_POP_ATOM(I)
	    MOD_SOUND_SPEED=SQRT(T1*(1.0D0+ED_ON_NA_EST)*OLD_T(I)+0.5*VTURB*VTURB)
	    SOUND_SPEED=SQRT(T1*(1.0D0+ED_ON_NA_EST)*OLD_T(I))
	    IF(OLD_V(I) .LT. VC_ON_SS*MOD_SOUND_SPEED)THEN
	      CONNECTION_VEL=OLD_V(I)
	      CONNECTION_RADIUS=OLD_R(I)
	      CONNECTION_INDX=I
	      RMAX=RMAX*CONNECTION_RADIUS
	      WRITE(6,*)'  Connection velocity is',CONNECTION_VEL
	      WRITE(6,*)'    Connection radius is',CONNECTION_RADIUS
	      WRITE(6,*)'       Maximum radius is',RMAX
	      WRITE(6,*)'     Connection INDEX is',CONNECTION_INDX
	      WRITE(6,*)'          Sound speed is',SOUND_SPEED
	      WRITE(6,*)' Modified sound speed is',MOD_SOUND_SPEED
	      EXIT
	    END IF
	  END DO
	END IF
!
!
!
! Compute Eddington ratio, GAM_EDD. This formulae is set for one electron per ion.
! This formula holds at all radii, since g and Teff both scale as 1/r^2.
!
	MAX_ED_ON_NA=0.0D0
	DO I=1,MOD_ND
	  MAX_ED_ON_NA=MAX(MAX_ED_ON_NA,OLD_ED(I)/OLD_POP_ATOM(I))
	END DO
	GAM_EDD=1.0D+06*SIGMA_TH*STEFAN_BC*(TEFF**4)/MU_ATOM/C_CMS/(10**LOGG)/AMU
	IF(GAM_EDD*MAX_ED_ON_NA .GT. 1.0D0)THEN
	  WRITE(LU_ERR,*)'An invalid Eddington parameter has been computed in DO_CMF_HYDRO_V2'
	  WRITE(LU_ERR,*)'Check the validity of Teff and Log G'
	  WRITE(LU_ERR,*)'The computed (maximum) Eddington parameter is ',GAM_EDD*MAX_ED_ON_NA
	  WRITE(LU_ERR,*)'                            Teff(K)/1.0+04 is',TEFF
	  WRITE(LU_ERR,*)'                         LOG_G (cgs units) is',LOGG
	  WRITE(LU_ERR,*)'                                MAX(Ne/Na) is',MAX_ED_ON_NA
	  STOP
	END IF
!
	IF(VERBOSE_OUTPUT)THEN
	  WRITE(LUV,'(A,ES14.6)')'          Surface gravity is:',LOGG
	  WRITE(LUV,'(A,ES14.6)')'             Mass of star is:',10**(LOGG)*(REFERENCE_RADIUS**2)/GRAV_CON
	  WRITE(LUV,'(A,ES14.6)')'         Mean atomic mass is:',MU_ATOM
	  WRITE(LUV,'(A,ES14.6)')'             Atom density is:',NI_ZERO
	  WRITE(LUV,'(A,ES14.6)')'New effective temperature is:',TEFF
	  WRITE(LUV,'(A,ES14.6)')'      Eddington parameter is:',GAM_EDD
	END IF
!
	PREV_REF_RADIUS=-1.0D0
	ITERATION_COUNTER=0
	DO WHILE(ABS(REFERENCE_RADIUS/PREV_REF_RADIUS-1.0D0) .GT. 1.0D-05)
	  ITERATION_COUNTER=ITERATION_COUNTER+1
	  IF(VERBOSE_OUTPUT)WRITE(LUV,*)' Beginning new hydro loop'
!
! The turbulent pressure is taken to be 0.5. roh . VTURB^2
!
	  PTURB_ON_NA=0.5D+10*VTURB*VTURB*MU_ATOM*AMU
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
! NB: 1.0D-06 is  1.0E+04/(1.0E+05)**2. The first factor occurs since T is in
! units of 10^4K. The second is to convert V from cm/s to km/s, with allowance
! for the sqrt.
!
	      I=CONNECTION_INDX
	      ED_ON_NA_EST=OLD_ED(I)/OLD_POP_ATOM(I)
	      ROSS_ON_ES=OLD_ROSS_MEAN(I)/OLD_ESEC(I)
	      T2=-1.0D0
	      DO WHILE(T2 .LT. 0.0D0)
	        GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA_EST*(OLD_FLUX_MEAN(I)/OLD_ESEC(I)) )
	        T1=1.0D-10*(BC*TEFF*(1+ED_ON_NA_EST)+PTURB_ON_NA)/( (10.0D0**LOGG)*(1.0D0-GAM_FULL)*MU_ATOM*AMU )
	        T2=CONNECTION_VEL*(1.0D0/T1-2.0D0/CONNECTION_RADIUS)
	        WRITE(LU_ERR,*)'      Scale height is',T1
	        WRITE(LU_ERR,*)'      Connection dVdR',T2
                IF(T2 .LE. 0.0D0)THEN
	          WRITE(LU_ERR,*)'    Connection radius',CONNECTION_RADIUS
                  WRITE(LU_ERR,*)'Resetting GAM_LIM due to -ve velocity gradient'
                  WRITE(LU_ERR,*)'         GAM_FULL was',GAM_FULL
                  WRITE(LU_ERR,*)'          GAM_LIM was',GAM_LIM
	          GAM_LIM=GAM_LIM-0.01D0
                  WRITE(LU_ERR,*)'       New GAM_LIM is',GAM_LIM
	        END IF
	      END DO
	      CALL WIND_VEL_LAW_V2(R,V,SIGMA,VINF,BETA,BETA2,RMAX,
	1          CONNECTION_RADIUS,CONNECTION_VEL,T2,VEL_LAW,J,ND_MAX)
	      DO I=1,J
	        POP_ATOM(I)=MDOT/MU_ATOM/R(I)/R(I)/V(I)
	      END DO
!
! To get other quantities we interpolate as a function of density.
! The atom density should be monotonic. At the outer boundary,
! we simply use the boundary value.
!
	      TA(1:MOD_ND)=LOG(OLD_POP_ATOM(1:MOD_ND))
	      TB(1:MOD_ND)=OLD_ED(1:MOD_ND)/ OLD_POP_ATOM(1:MOD_ND)
	      TC(1:J)=LOG(POP_ATOM(1:J))
	      DO I=1,J
	        IF(TC(I) .LT. TA(1))TC(I)=TA(1)
	      END DO
	      CALL MON_INTERP(ED_ON_NA,J,IONE,TC,J,TB,MOD_ND,TA,MOD_ND)
	      TB(1:MOD_ND)=OLD_T(1:MOD_ND)*TEFF/OLD_TEFF
	      CALL MON_INTERP(T,J,IONE,TC,J,TB,MOD_ND,TA,MOD_ND)
	      TB(1:MOD_ND)=OLD_ROSS_MEAN(1:MOD_ND)/OLD_ESEC(1:MOD_ND)
	      CALL MON_INTERP(CHI_ROSS,J,IONE,TC,J,TB,MOD_ND,TA,MOD_ND)
	      DO I=1,J
	        ED(I)=ED_ON_NA(I)*POP_ATOM(I)
	        GAMMA_FULL(I)=GAM_FULL
	        CHI_ROSS(I)=ED_ON_NA(I)*SIGMA_TH*POP_ATOM(I)*CHI_ROSS(I)
	        P(I)=(BC*T(I)*(1.0D0+ED(I)/POP_ATOM(I))+PTURB_ON_NA)*POP_ATOM(I)
	      END DO
	      CALL TORSCL(TAU,CHI_ROSS,R,TB,TC,J,'LOGMON',' ')
	      I=J
!
! The following section is for the case of no wind.
!
	  ELSE
!
! This option is for a pure spherical or plane-parallel model. The depth variation
! of gravity is taken into account.
! 
	      POP_ATOM(1)=NI_ZERO
	      T(1)=(TEFF/OLD_TEFF)*OLD_T(1)
	      P(1)=(BC*T(1)*(1.0D0+OLD_ED(1)/OLD_POP_ATOM(1))+PTURB_ON_NA)*POP_ATOM(1)
	      IF(RBOUND .EQ. 0.0D0)RBOUND=OLD_R(1)
	      R(1)=RBOUND
	      TAU(1)=OLD_TAU(1)
	      ED_ON_NA(1)=OLD_ED(1)/OLD_POP_ATOM(1)
	      ROSS_ON_ES=OLD_ROSS_MEAN(1)/OLD_ESEC(1)
	      GAM_FULL=MIN( GAM_LIM,GAM_EDD*ED_ON_NA(1)*(OLD_FLUX_MEAN(1)/OLD_ESEC(1)) )
	      ED(1)=ED_ON_NA(1)*POP_ATOM(1)
	      CHI_ROSS(1)=ED_ON_NA(1)*SIGMA_TH*POP_ATOM(1)*ROSS_ON_ES
	      GAMMA_FULL(1)=GAM_FULL
	      V(1)=MDOT/MU_ATOM/POP_ATOM(1)/R(1)/R(1)
	      I=1
	  END IF
!
! Set sound speed at the point we match the wind to the hydrostatic structure.
! In no wind is present, we set the sound speed to a large number. This effectively
! set the correction to zero.
!
	  IF(WIND_PRESENT)THEN
	    SOUND_SPEED=1.0D+04*(1.0D0+ED_ON_NA(I))*BOLTZMANN_CONSTANT()*T(I)/MU_ATOM/ATOMIC_MASS_UNIT()
	    SOUND_SPEED=1.0D-05*SQRT(SOUND_SPEED)
	    WRITE(6,'(A,3ES14.4)')'SOUND_SPEED',SOUND_SPEED
	  ELSE
	    IF(VTURB .EQ. 0.0D0)THEN
	      SOUND_SPEED=1.0D+30
	    ELSE
	      SOUND_SPEED=0.0D+00
	    END IF
	  END IF
!
!
!
! The boudary condition for the integration of the hydrostatic equation
! has been set, either at the outer boundary, or at the wind connection point.
! we can now perform the integration of the hydrostatic equation.
!
	  DO WHILE( TAU(I) .LT. MAX(100.0D0,OLD_TAU_MAX) )
	    I=I+1
	    IF(I .GT. ND_MAX)THEN
	      WRITE(6,*)'Error id DO_CMF_HYDRO'
	      WRITE(6,*)'ND_MAX is too small'
	      STOP
	    END IF
!
! Compute the atmospheric pressure scale height. Close to the sonic point,
! Gamma can be very close to 1 which leads to a large scale height. But this
! will be come smaller as we move away from the sonic point. Thus to ensure
! that the step size is sufficiently small, we set the scale height using
! GAM_EDD (i.e., GAMMA computed the electron scattering opacity only) rather
! than GAM_FULL.
!
	    SCL_HT=(10.0D0**LOGG)*(1.0D0-GAM_EDD)*MU_ATOM*AMU/
	1            (BC*(1.0D0+ED_ON_NA(I-1))*T(I-1)+PTURB_ON_NA)
	    SCL_HT=1.0D-10/SCL_HT
!
! We set the step size to a fraction of the pressure scale height.
!
	    H=SCL_HT/10.0D0
!
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LUV,*)' '
	      WRITE(LUV,'(A,ES12.4,A)')'       Scale height is',SCL_HT,' 10^10 cm'
	      WRITE(LUV,*)'I=',I
	      WRITE(LUV,*)'P(I-1)=',P(I-1)
	      WRITE(LUV,*)'  Pgas=',BC*T(I-1)*POP_ATOM(I-1)*(1+ED_ON_NA(I-1))
	      WRITE(LUV,*)' Pturb=',0.5D+10*POP_ATOM(I-1)*MU_ATOM*AMU*VTURB*VTURB
	      WRITE(LUV,*)'R(I-1)=',R(I-1)
	      WRITE(LUV,*)'V(I-1)=',V(I-1)
	      WRITE(LUV,*)'T(I-1)=',T(I-1)
	      WRITE(LUV,*)'TAU(I-1)=',TAU(I-1)
	      WRITE(LUV,*)'ED_ON_NA(I-1)=',ED_ON_NA(I-1)
	      WRITE(LUV,*)'POP_ATOM(I-1)=',POP_ATOM(I-1)
	      WRITE(LUV,*)'GAMMA_FULL=',GAM_FULL
	    END IF
!
! Set estimates at current location. Then integrate hydrostatic
! equation using 4th order Runge-Kutta.
!
	    IF(VERBOSE_OUTPUT)WRITE(LUV,'(5(5X,A,5X))')'  dPdR','dTAUdR','  T_EST','   dPn',' dTAUn'
100	    P_EST=P(I-1)
	    TAU_EST=TAU(I-1)
	    T_EST=T(I-1)
	    ED_ON_NA_EST=ED_ON_NA(I-1)
	    ATOM_EST=POP_ATOM(I-1)
	    R_EST=R(I-1)-H
	    CALL CMF_HYDRO_DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP1=H*dPdR
	    dTAU1=H*dTAUdR
	    IF(VERBOSE_OUTPUT)WRITE(LUV,'(5ES16.8)')dPdR,dTAUdR,T_EST,dP1,dTAU1
!
	    P_EST=P(I-1)+dP1/2
	    TAU_EST=TAU(I-1)+dTAU1/2
	    ATOM_EST=P_EST/(BC*T_EST*(1+ED_ON_NA_EST)+PTURB_ON_NA)
	    R_EST=R(I-1)-0.5D0*H
	    CALL CMF_HYDRO_DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP2=H*dPdR
	    dTAU2=H*dTAUdR
	    IF(VERBOSE_OUTPUT)WRITE(LUV,'(5ES16.8)')dPdR,dTAUdR,T_EST,dP2,dTAU2
!
	    P_EST=P(I-1)+dP2/2
	    TAU_EST=TAU(I-1)+dTAU2/2
	    ATOM_EST=P_EST/(BC*T_EST*(1+ED_ON_NA_EST)+PTURB_ON_NA)
	    R_EST=R(I-1)-0.5D0*H
	    CALL CMF_HYDRO_DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP3=H*dPdR
	    dTAU3=H*dTAUdR
	    IF(VERBOSE_OUTPUT)WRITE(LUV,'(5ES16.8)')dPdR,dTAUdR,T_EST,dP3,dTAU3
!
	    P_EST=P(I-1)+dP3
	    TAU_EST=TAU(I-1)+dTAU2
	    ATOM_EST=P_EST/(BC*T_EST*(1+ED_ON_NA_EST)+PTURB_ON_NA)
	    R_EST=R(I-1)-H
	    CALL CMF_HYDRO_DERIVS(P_EST,TAU_EST,T_EST,ED_ON_NA_EST,ATOM_EST)
	    dP4=H*dPdR
	    dTAU4=H*dTAUdR
	    IF(VERBOSE_OUTPUT)WRITE(LUV,'(5ES16.8)')dPdR,dTAUdR,T_EST,dP4,dTAU4
!
! Update values at next grid point.
!
	    TAU(I)=TAU(I-1)+(dTAU1+2*dTAU2+2*dTAU3+dTAU4)/6.0D0
	    P(I)=P(I-1)+(dP1+2*dP2+2*dP3+dP4)/6.0D0
!
	    CALL CMF_HYDRO_NEW_EST(TAU(I),T_EST)
	    R(I)=R(I-1)-H
	    T(I)=T_EST
	    ED_ON_NA(I)=ED_ON_NA_EST
	    POP_ATOM(I)=P(I)/(BC*T(I)*(1.0D0+ED_ON_NA(I))+PTURB_ON_NA)
	    ED(I)=ED_ON_NA(I)*POP_ATOM(I)
	    CHI_ROSS(I)=ED_ON_NA(I)*SIGMA_TH*POP_ATOM(I)*ROSS_ON_ES
	    GAMMA_FULL(I)=GAM_FULL
	    V(I)=MDOT/MU_ATOM/POP_ATOM(I)/R(I)/R(I)
!	    IF(.NOT. PLANE_PARALLEL_NO_V)THEN
!	      IF(V(I) .GE. V(I-1))THEN
!	        GAM_LIM=GAM_LIM-0.01
!	        IF(VERBOSE_OUTPUT)WRITE(LUV,*)'Resetting GAM_LIM due to -ve velocity gradient'
!	        GOTO 100
!	      END IF 
!	    END IF
	    GAM_LIM=GAM_LIM_STORE
	    ND=I
!
	  END DO		!Loop over inner atmosphere
	  CLOSE(UNIT=75)
	  CLOSE(UNIT=76)
!
! Output estimates at the last depth.
!
	  I=ND
	  IF(VERBOSE_OUTPUT)THEN
	    WRITE(LUV,*)'I=',I
	    WRITE(LUV,*)'P(I)=',P(I)
	    WRITE(LUV,*)'R(I)=',R(I)
	    WRITE(LUV,*)'T(I)=',T(I)
	    WRITE(LUV,*)'TAU(I)=',TAU(I)
	    WRITE(LUV,*)'ED_ON_NA(I)=',ED_ON_NA(I)
	    WRITE(LUV,*)'POP_ATOM(I)=',POP_ATOM(I)
	  END IF
!
! Output diagnostic files. These are on the calculate grid --- not the final
! grid.
!
	  ALLOCATE(COEF(ND,4))
	  CALL MON_INT_FUNS_V2(COEF,P,R,ND)
	  DO I=1,ND
	    dPdR_VEC(I)=1.0D-10*COEF(I,3)/POP_ATOM(I)/AMU/MU_ATOM
	  END DO
	  DEALLOCATE(COEF)
!
	  IF(VERBOSE_OUTPUT)THEN
	    INQUIRE(UNIT=LU,OPENED=FILE_OPEN)
	    IF(FILE_OPEN)THEN
	      WRITE(LU,'(/,A,/)')FORMFEED
	    ELSE
	      OPEN(UNIT=LU,FILE='NEW_CALC_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	    END IF
	    WRITE(LU,'(A,10(7X,A))')' Index','     R','  Vel','   Tau','    Na',' Ne/Na','     T',
	1                  ' Kross','Kr/Kes',' Gamma','dpdR/ROH'
	    DO I=1,ND
	      T1=1.0D+10*POP_ATOM(I)*AMU*MU_ATOM
	      T2=SIGMA_TH*POP_ATOM(I)*ED_ON_NA(I)
	      T3=MDOT/MU_ATOM/POP_ATOM(I)/R(I)/R(I)
	      WRITE(LU,'(I6,10ES13.4)')I,R(I),T3,TAU(I),POP_ATOM(I),ED_ON_NA(I),
	1              T(I),CHI_ROSS(I)/T1,CHI_ROSS(I)/T2,GAMMA_FULL(I),dPdR_VEC(I)
	    END DO
	    CLOSE(LU)
	  END IF
!
! Adjust the grid so that we get the correct reference radius,
! defined as Tau(Ross)=TAU_REF.
!
	  T1=TAU_REF
	  I=GET_INDX_DP(T1,TAU,ND)
	  T2=(LOG(T1)-LOG(TAU(I)))/(LOG(TAU(I+1))-LOG(TAU(I)))
	  T2=(1.0D0-T2)*R(I)+T2*R(I+1)
	  T1=T2-REFERENCE_RADIUS
	  RADIUS_AT_TAU_23=T2
	  IF(WIND_PRESENT .AND. USE_OLD_VEL)THEN
	    PREV_REF_RADIUS=REFERENCE_RADIUS
	  ELSE IF(WIND_PRESENT)THEN
	    PREV_REF_RADIUS=T2
	    CONNECTION_RADIUS=CONNECTION_RADIUS-T1
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LUV,*)'    Old reference radius is',T2
	      WRITE(LUV,*)'Desired reference radius is',REFERENCE_RADIUS
	    END IF
	  ELSE IF(PLANE_PARALLEL_MOD)THEN
!
! For a plane-parallel model, the reference radius is the inner boundary,
! defined at TAU=100, not the radius at which tau=2/3.
!
	    T1=100.0D0
	    I=GET_INDX_DP(T1,TAU,ND)
	    T2=(LOG(T1)-LOG(TAU(I)))/(LOG(TAU(I+1))-LOG(TAU(I)))
	    T2=(1.0D0-T2)*R(I)+T2*R(I+1)
	    T1=T2-REFERENCE_RADIUS
	    PREV_REF_RADIUS=T2
	    RBOUND=RBOUND-T1
	    NI_ZERO=NI_ZERO*REFERENCE_RADIUS/PREV_REF_RADIUS
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LUV,*)'    Old reference radius is',T2
	      WRITE(LUV,*)'Desired reference radius is',REFERENCE_RADIUS
	    END IF
	  ELSE
	    PREV_REF_RADIUS=T2
	    RBOUND=RBOUND-T1
	    NI_ZERO=NI_ZERO*REFERENCE_RADIUS/PREV_REF_RADIUS
	    IF(VERBOSE_OUTPUT)THEN
	      WRITE(LUV,*)'    Old reference radius is',T2
	      WRITE(LUV,*)'Desired reference radius is',REFERENCE_RADIUS
	    END IF
	  END IF
!
	  IF(ITERATION_COUNTER .GE. 100)THEN
	    WRITE(LU_ERR,*)'Exceed iteration count in DO_CMF_HYDRO_V2.'
	    WRITE(LU_ERR,*)'Aborting update of the hydro structure.'
	    RETURN
	  END IF
	END DO			!Loop to set R(Tau=2/3)=REFERENCE_RADIUS
!
	IF(WIND_PRESENT .AND. USE_OLD_VEL .AND. RESET_REF_RADIUS)THEN
	  REFERENCE_RADIUS=RADIUS_AT_TAU_23
	END IF
!
	IF(VERBOSE_OUTPUT)THEN
	  CLOSE(UNIT=LUV)
	END IF
!
!
! We now create the revised grid, At present it is equally spaced in Log(tau) with
! 2 extra points at either end of the grid.
!
	NEW_ND=MOD_ND
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
	TAU_MAX=100.0D0
	IF(TAU_MAX .GT. TAU(ND))THEN
	  WRITE(LU_ERR,*)'Error --- TAU_MAX cannot be greater than calculated grid Tau'
	  WRITE(LU_ERR,*)'Setting TAU to maximum value in DO_CMF_HYDRO'
	  TAU_MAX=TAU(ND)
	END IF 
!
	DO I=1,ND
	  V(I)=MDOT/MU_ATOM/POP_ATOM(I)/R(I)/R(I)
	END DO
	IF(WIND_PRESENT .AND. USE_OLD_VEL)THEN
	  J=CONNECTION_INDX
	  I=NEW_ND-J+1
	  CALL DET_R_GRID_V1(REV_TAU(J),I,ND_MAX,TAU_MAX,L_FALSE,
	1                    R(J),V(J),TAU(J),ND-J+1)
	  REV_TAU(1:J)=TAU(1:J)
	ELSE
	  CALL DET_R_GRID_V2(REV_TAU,NEW_ND,ND_MAX,TAU_MAX,
	1        dLOG_TAU,V_SCL_FAC,OUT_BND_OPT,OBND_PARS,NUM_OBND_PARAMS,
	1        R,V,TAU,ND)
	END IF
!
! We now compute the revised R grid. We then interplate on Log (r^2.rho) which 
! is equivalent to interpolating on log V. This guarentees monotocity of V.
!
	TAU(1:ND)=LOG(TAU(1:ND))
	REV_TAU(1:NEW_ND)=LOG(REV_TAU(1:NEW_ND))
	CALL MON_INTERP(REV_R,NEW_ND,IONE,REV_TAU,NEW_ND,R,ND,TAU,ND)
	IF(PLANE_PARALLEL_MOD)THEN
	  POP_ATOM(1:ND)=LOG(POP_ATOM(1:ND))
	  CALL MON_INTERP(REV_POP_ATOM,NEW_ND,IONE,REV_R,NEW_ND,POP_ATOM,ND,R,ND)
	  REV_POP_ATOM(1:NEW_ND)=EXP(REV_POP_ATOM(1:NEW_ND))
	  POP_ATOM(1:ND)=EXP(POP_ATOM(1:ND))
	ELSE
	  POP_ATOM(1:ND)=LOG(POP_ATOM(1:ND)*R(1:ND)*R(1:ND))
	  CALL MON_INTERP(REV_POP_ATOM,NEW_ND,IONE,REV_R,NEW_ND,POP_ATOM,ND,R,ND)
	  REV_POP_ATOM(1:NEW_ND)=EXP(REV_POP_ATOM(1:NEW_ND))/REV_R(1:ND)/REV_R(1:ND)
	  POP_ATOM(1:ND)=EXP(POP_ATOM(1:ND))/R(1:ND)/R(1:ND)
	END IF
!
! This set the reference radius to R(ND), and keeps it the same as was read in.
!
	IF(PLANE_PARALLEL_MOD)THEN
	  DO I=1,ND
	    R(I)=R(I)+(REFERENCE_RADIUS-REV_R(NEW_ND))
	  END DO
	  DO I=1,NEW_ND
	    REV_R(I)=REV_R(I)+(REFERENCE_RADIUS-REV_R(NEW_ND))
	  END DO
	END IF
!
! Compute revised velocity.
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
!
! Saves the current RVSIG_FILE for recovery puposes. Except for the first
! iteration, this could be recovered from RVSIG_COL. For portability, we
! use only regular fortran commands.
!
	INQUIRE(FILE='RVSIG_COL',EXIST=FILE_PRES)
	IF(FILE_PRES)THEN
	  STRING=' '
	  WRITE(STRING,'(I4.4)')MAIN_COUNTER
	  STRING='RVSIG_COL_IT_'//TRIM(STRING)
   	  OPEN(UNIT=LUIN,FILE='RVSIG_COL',STATUS='OLD',ACTION='READ')
	  OPEN(UNIT=LU,FILE=TRIM(STRING),STATUS='UNKNOWN',ACTION='WRITE')
	  DO WHILE(1 .EQ. 1)
	    READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    WRITE(LU,'(A)')TRIM(STRING)
	    IF(IOS .NE. 0)EXIT
	  END DO
	  CLOSE(LUIN)
	  CLOSE(LU)
	END IF
!
! Output revised hydrostatic structure. This can be used to restart the current
! model from scratch.
!
	OPEN(UNIT=LU,FILE='RVSIG_COL',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(A)')'! Note: The effective temperature and surface gravity are defined'
	  WRITE(LU,'(A)')'! at the reference radius, which (except when using the old'
	  WRITE(LU,'(A)')'! velocity) is the location where Tau=2/3.'
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(A,ES16.6)')'! Effective temperature (10^4 K) is:',TEFF
	  WRITE(LU,'(A,ES16.6)')'!      Log surface gravity (cgs) is:',LOGG
	  WRITE(LU,'(A,ES16.6)')'!         Core radius (10^10 cm) is:',REV_R(NEW_ND)
	  WRITE(LU,'(A,ES16.6)')'!    Reference radius (10^10 cm) is:',REFERENCE_RADIUS
	  WRITE(LU,'(A,ES16.6)')'!              Luminosity (Lsun) is:',( (TEFF/0.5770D0)**4 )*( (REFERENCE_RADIUS/6.9599D0)**2 )
	  WRITE(LU,'(A,ES16.6)')'!            Mass (Msun) of star is:',10.0D0**(LOGG)*REFERENCE_RADIUS*REFERENCE_RADIUS/GRAV_CON
	  WRITE(LU,'(A,ES16.6)')'!       Mass loss rate (Msun/yr) is:',MDOT/MASS_LOSS_SCALE_FACTOR
	  WRITE(LU,'(A,ES16.6)')'!         Mean atomic mass (amu) is:',MU_ATOM
	  WRITE(LU,'(A,ES16.6)')'!            Eddington parameter is:',GAM_EDD
	  WRITE(LU,'(A,ES16.6)')'!                   Atom density is:',NI_ZERO
	  WRITE(LU,'(A,F14.8)') '! Ratio of inner to outer radius is:',REV_R(1)/REV_R(NEW_ND)
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(3X,I5,10X,A)')NEW_ND,'!Number of depth points'
	  WRITE(LU,'(A)')'!'
	  WRITE(LU,'(A,4X,A,3(7X,A),3X,A)')'!','R(10^10cm)','V(km/s)','  Sigma','    Tau','  Index'
	  IF(REV_R(1) .GT. 999999.0D0)THEN
	    DO I=1,NEW_ND
	      WRITE(LU,'(F18.8,3ES14.6,6X,I4)')REV_R(I),REV_V(I),REV_SIGMA(I),EXP(REV_TAU(I)),I
	    END DO
	  ELSE
	    DO I=1,NEW_ND
	      WRITE(LU,'(F15.8,3ES14.6,6X,I4)')REV_R(I),REV_V(I),REV_SIGMA(I),EXP(REV_TAU(I)),I
	    END DO
	  END IF
	CLOSE(UNIT=LU)
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
	OPEN(UNIT=LU,FILE='FIN_CAL_GRID',STATUS='UNKNOWN',ACTION='WRITE')
	  WRITE(41,'(A,9(6X,A7))')'Index','     R','   TAU','    V','     T','    Na',
	1              ' Ne/Na',' Xross','Xr/Xes','   Gam'
	  DO I=1,NEW_ND
	    WRITE(LU,'(I5,9ES13.4)')I,REV_R(I),EXP(REV_TAU(I)),REV_V(I),REV_T(I),
	1             REV_POP_ATOM(I),REV_ED(I)/REV_POP_ATOM(I),REV_CHI_ROSS(I),
	1             REV_CHI_ROSS(I)/SIGMA_TH/REV_ED(I),REV_GAMMA_FULL(I)
	  END DO
	CLOSE(UNIT=LU)
!
	CALL SET_NEW_GRID(REV_R,REV_V,REV_SIGMA,REV_ED,REV_CHI_ROSS,NEW_ND)
	CALL SET_ABUND_CLUMP(T1,T2,LU_ERR,NEW_ND)
	CALL ADJUST_POPS(POPS,LU_ERR,NEW_ND,NT)
!
	NO_HYDRO_ITS=NO_HYDRO_ITS-1
	NO_ITS_DONE=NO_ITS_DONE+1
	CALL UPDATE_KEYWORD(NO_HYDRO_ITS,'[N_ITS]','HYDRO_DEFAULTS',L_TRUE,L_FALSE,LUIN)
	CALL UPDATE_KEYWORD(NO_ITS_DONE,'[ITS_DONE]','HYDRO_DEFAULTS',L_FALSE,L_TRUE,LUIN)
	DONE_HYDRO=.TRUE.
!
! Make sure VADAT is consitent with revised RGRID & parameters in HYDRO_PARAMS.
! We also return the correct RSTAR, RMAX and stellar mass.
!
	MOD_RSTAR=REV_R(NEW_ND)
	MOD_RMAX=REV_R(1)
	MOD_MASS=10.0D0**(LOGG)*(REFERENCE_RADIUS**2)/GRAV_CON
	CALL UPDATE_KEYWORD(REV_R(NEW_ND),'[RSTAR]','VADAT',L_TRUE,L_FALSE,LUIN)
	T1=MOD_RMAX/REV_R(NEW_ND)
	CALL UPDATE_KEYWORD(T1,'[RMAX]','VADAT',L_FALSE,L_FALSE,LUIN)
	CALL UPDATE_KEYWORD(MOD_MASS,'[MASS]','VADAT',L_FALSE,L_TRUE,LUIN)
!	CALL UPDATE_KEYWORD(REV_V(1),'[VINF]','VADAT',L_FALSE,L_TRUE,LUIN)
	WRITE(LU_ERR,'(A)')' Revised hydrostatic structure and output new RVSIG_COL file'
	WRITE(LU_ERR,'(A)')' Updated RSTAR, RMAX and MASS in VADAT'
!
	RETURN
	END 
!	
! Subroutine to compute dPdR and dTAUdR for use with the
! Runge-Kutta integration.
!
	SUBROUTINE CMF_HYDRO_DERIVS(P,TAU,TEMP,ED_ON_NA,POP_ATOM)
	USE CMF_HYDRO_MODULE
	USE OLD_GRID_MODULE
	IMPLICIT NONE
!
	REAL*8 P
	REAL*8 TAU
	REAL*8 TEMP
	REAL*8 ED_ON_NA
	REAL*8 POP_ATOM
	REAL*8 T1,T2,T3
!
	CALL CMF_HYDRO_NEW_EST(TAU,TEMP)
!
	IF(PLANE_PARALLEL_MOD)THEN
	  T1=(10.0D0**LOGG)*(1.0D0-GAM_FULL)
	ELSE
	  T1=(10.0D0**LOGG)*(1.0D0-GAM_FULL)*(REFERENCE_RADIUS/R_EST)**2
	END IF
	T2=(MDOT/MU_ATOM/POP_ATOM/R_EST/R_EST)**2/(SOUND_SPEED**2+0.5D0*VTURB**2)
!	T2=(MDOT/MU_ATOM/POP_ATOM/R_EST/R_EST)**2/(SOUND_SPEED**2)                       !+0.5D0*VTURB**2)
	T3=BC*TEMP*(1.0D0+ED_ON_NA)/MU_ATOM/AMU+0.5D+10*VTURB*VTURB
	dPdR=1.0D+10*T1*P/T3/(1.0D0-T2)
!	dPdR=1.0D+10*T1*MU_ATOM*AMU*P/BC/TEMP/(1+ED_ON_NA)/(1.0D0-T2)
	dTAUdR=ED_ON_NA*POP_ATOM*SIGMA_TH*ROSS_ON_ES
!
	WRITE(17,'(A,ES14.4)')'      P_EST=',P
	IF(VTURB .NE. 0)THEN
	  WRITE(17,'(A,ES14.4)')'    P(TURB)=',0.5D+10*MU_ATOM*AMU*POP_ATOM*VTURB*VTURB
	  WRITE(17,'(A,ES14.4)')'     P(GAS)=',BC*POP_ATOM*TEMP*(1+ED_ON_NA)
	END IF
	WRITE(17,'(A,ES14.4)')'    (V/C)^2=',T2
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
	SUBROUTINE CMF_HYDRO_NEW_EST(TAU,TEMP)
	USE CMF_HYDRO_MODULE
	USE OLD_GRID_MODULE
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
	REAL*8 KAP_ROSS_LTE
	REAL*8 KAP_ROSS_EST
	REAL*8 KAP_ES_OLD
	REAL*8 KAP_ES_LTE
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
	  TEMP=T_EST     				!Use value passed
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
	    CALL GET_LTE_ROSS_V2(KAP_ROSS_LTE,KAP_ES_LTE,LTE_ED,ATOM_EST,TEMP)
!	    KAP_ROSS_EST=(T1-KAP_ES)*( (RM-ES)/(KAP_ROSS_OLD-KAP_ES_OLD))+KAP_ES
!           KAP_ROSS=(T1-KAP_ES)*(RM/KAP_ROSS_OLD)+KAP_ES
!           KAP_ROSS=T1*(RM/KAP_ROSS_OLD)
!	    GAM_FULL=GAM_EDD*(FM/RM)*(KAP_ROSS/KAP_ES)*ED_ON_NA_EST
            KAP_ROSS_EST=(KAP_ROSS_LTE-KAP_ES_LTE)*( (RM-ES)/(KAP_ROSS_OLD-KAP_ES_OLD))+ES
	    KAP_ROSS=RM; KAP_ES=ES
!	    GAM_FULL=GAM_EDD*(FM/RM)*(KAP_ROSS/KAP_ES)*ED_ON_NA_EST
	    GAM_FULL=GAM_EDD*(FM/RM)*(0.5D0*(KAP_ROSS_EST+RM)/KAP_ES)*ED_ON_NA_EST
	    ROSS_ON_ES=KAP_ROSS/KAP_ES
	    WRITE(76,'(ES16.8,3ES12.4,/,16X,3ES12.4,/,16X,ES12.4,/,16X,4ES12.4)')
	1                       R_EST,KAP_ROSS_OLD,KAP_ES_OLD,OLD_TEMP,
	1                       KAP_ROSS_LTE,KAP_ES_LTE,TEMP,
	1                       KAP_ROSS_EST,
	1                       RM,ES,FM,GAM_FULL
	  ELSE
	    GAM_FULL=GAM_EDD*(FM/ES)*ED_ON_NA_EST
	    ROSS_ON_ES=RM/ES
	  END IF
	  GAM_FULL=MIN(GAM_LIM,GAM_FULL)
	  WRITE(75,'(ES16.9,8ES12.3)')R_EST,GAM_FULL,ED_ON_NA_EST*GAM_EDD*FM/ES,TAU,FM,RM,KAP_ROSS,KAP_ES,ED_ON_NA_EST
	END IF
!  
	RETURN
	END 
