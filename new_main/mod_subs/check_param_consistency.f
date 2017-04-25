!
! File to check the consistency of CMFGEN parameters
!
! Altered: 23-May-2010 - Fixed pontential compile problems with .NE. and logical operators.
! Created: 07-Jan-2008
!
	SUBROUTINE CHECK_PARAM_CONSISTENCY( )
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	CHARACTER(LEN=10) ER_LAB
!
	ER_LAB='Warning:'
	IF(STOP_IF_BAD_PARAM)ER_LAB='ERROR:'
	LUER=ERROR_LU()
!	
	IF(SN_MODEL)THEN
	  IF(.NOT. DO_FULL_REL_OBS)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' DO_FULL_RELL_OBS should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(.NOT. DO_FULL_REL_CMF)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' DO_FULL_REL_CMF should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(.NOT. USE_J_REL .AND. .NOT. USE_DJDT_RTE .AND. .NOT. USE_LAM_ES)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' USE_JREL or USE_DJDT_RTE should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	   END IF
	  IF(USE_J_REL .AND. .NOT. INCL_REL_TERMS)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' INCL_REL should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
	  IF(USE_J_REL .AND. .NOT. INCL_ADVEC_TERMS_IN_TRANS_EQ)THEN
	    WRITE(LUER,*)TRIM(ER_LAB),' Problem with control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),' INCL_ADV_TRANS should be true for SN MODEL'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
!
	  IF(USE_DJDT_RTE .AND. .NOT. INCL_ADIABATIC)THEN
	    WRITE(LUER,*)'Error: inconsistency in control parameters in VADAT'
	    WRITE(LUER,*)'SN model has USE_DJDT_RTE set to TRUE.'
	    WRITE(LUER,*)'INCL_ADIABATIC must also be set to TRUE.'
	    STOP
	  END IF
	  IF(USE_DJDT_RTE .AND. .NOT. DO_CO_MOV_DDT)THEN
	    WRITE(LUER,*)'Error: inconsistency in control parameters in VADAT'
	    WRITE(LUER,*)'SN model has USE_DJDT_RTE set to TRUE.'
	    WRITE(LUER,*)'DO_CO_MOV_DDT should also be set to TRUE.'
	    IF(STOP_IF_BAD_PARAM)STOP
	  END IF
!
! AS TIME_SEQ_NO is not actually used, this check is not necessary.
!
!	  IF(USE_DJDT_RTE .AND. (TIME_SEQ_NO .EQ. 1) )THEN
!	    WRITE(LUER,*)'Error: inconsistency in control parameters in VADAT'
!	    WRITE(LUER,*)'SN model has TS_NO=1 implying initial model'
!	    WRITE(LUER,*)TRIM(ER_LAB),'USE_JREL=TRUE and USE_DJDT_RTE=FALSE for an initial model'
!	    STOP
!	  END IF
!
	  IF(USE_J_REL .AND. USE_DJDT_RTE)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)TRIM(ER_LAB),'Error: USE_JREL or USE_DJDT_RTE cannot be true at the same time.'
	    WRITE(LUER,*)TRIM(ER_LAB),' For initial time seqnece model, only USE_JRL should be true.'
	    STOP
	  END IF
	  IF(PLANE_PARALLEL  .OR. PLANE_PARALLEL_NO_V)THEN
	    WRITE(LUER,*)'Error in control parameters in VADAT'
	    WRITE(LUER,*)'PP_NOV and PP_MOD cannot be TRUE for a SN model'
	    STOP
	  END IF
	END IF
!
	IF(PLANE_PARALLEL  .AND. PLANE_PARALLEL_NO_V)THEN
	  WRITE(LUER,*)'Error in control parameters in VADAT'
	  WRITE(LUER,*)'PP_NOV and PP_MOD cannot both be TRUE'
	  WRITE(LUER,*)'Use PP_NOV for a plane-parallel model with no velocity field'
	  WRITE(LUER,*)'Use PP_MOD for a plane-parallel model with a velocity field'
	  STOP
	END IF
!
	RETURN
	END	
