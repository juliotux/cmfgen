	MODULE CMF_HYDRO_MODULE
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
	  REAL*8 MDOT
	  REAL*8 GRAV_CON
	  REAL*8 REFERENCE_RADIUS
	  REAL*8 PREV_REF_RADIUS
	  REAL*8 RADIUS_AT_TAU_23
	  REAL*8 SOUND_SPEED
	  REAL*8 MOD_SOUND_SPEED
	  REAL*8 VTURB
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
	  REAL*8 V_EST
	  REAL*8 F_EST
!
	  REAL*8 dPdR
	  REAL*8 dTAUdR
!
	  REAL*8 OLD_TEFF
	  REAL*8 OLD_REF_RADIUS
	  LOGICAL PURE_LTE_EST
	  LOGICAL OLD_MODEL
	  LOGICAL WIND_PRESENT 
	  LOGICAL PLANE_PARALLEL_MOD
	  LOGICAL RESET_REF_RADIUS
	  INTEGER, SAVE :: LAST_HYDRO_ITERATION=0

	END MODULE CMF_HYDRO_MODULE 	
