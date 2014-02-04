C
C Compute opacity and emissivity. This is a general include file
C provided program uses exactly the same variables. Can be achieved
C by copying declaration statements from CMFGEN. Always advisable
C to use ``IMPLICIT NONE''.
C
	IF(COMPUTE_NEW_CROSS)THEN
C
C Compute EXP(-hv/kT) and zero CHI, and ETA.
C
	  T1=-HDKT*CONT_FREQ
	  DO I=1,ND
	    EMHNUKT(I)=EXP(T1/T(I))
	    CHI(I)=0.0D0
	    ETA(I)=0.0D0
	  END DO
C
C Compute continuum intensity incident from the core assuming a TSTAR
C blackbody.
C
	  IC=TWOHCSQ*(CONT_FREQ**3)/(EXP(HDKT*CONT_FREQ/TSTAR)-1.0D0)
C
C Compute opacity and emissivity. ESOPAC must be call first since
C CHI is effectively zeroed in that routine.
C
	  IF(GLOBAL_LINE_SWITCH .EQ. 'BLANK')THEN
	    CALL ESOPAC(ESEC,ED,ND)		!Electron scattering emission factor.
	    CALL ESOPAC(CHI,ED,ND)		!Electron scattering opacity.
	  ELSE
!
! Add in electron opacity for this species.
!
	    TA(1:ND)=0.0D0 
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES .AND. ATM(ID)%XzV_TRANS_TYPE .EQ. 'BLANK')THEN
	        DO I=1,ND
	          T1=0.0D0
	          DO J=1,ATM(ID)%NxzV_F
	            T1=T1+ATM(ID)%XzV_F(J,I)
	          END DO
	          TA(I)=TA(I)+T1*(ATM(ID)%ZXzV-1.0D0)
	        END DO
	        IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
	          TA(1:ND)=TA(1:ND)+ATM(ID)%ZXzV*ATM(ID+1)%XzV(1,1:ND)
	        END IF
	      END IF
	    END DO
	    CALL ESOPAC(ESEC,TA,ND)
	    CALL ESOPAC(CHI,TA,ND)
	  END IF
C
C Free-free and bound-free opacities.
C
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES .AND. (GLOBAL_LINE_SWITCH .EQ. 'BLANK' .OR.
	1                 ATM(ID)%XzV_TRANS_TYPE .EQ. 'BLANK') )THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        CALL GENOPAETA_V8(ID,CHI,ETA,CONT_FREQ,
	1           ATM(ID)%XzV_F,     ATM(ID)%XzVLTE_F, ATM(ID)%EDGEXzV_F,
	1           ATM(ID)%GIONXzV_F, ATM(ID)%ZXzV,     ATM(ID)%NXzV_F,
	1           ATM(ID+1)%XzV,     ATM(ID+1)%XzVLTE, ATM(ID+1)%NXzV, J,
	1           ATM(ID)%XzV_ION_LEV_ID(J),
	1           ED,T,EMHNUKT,L_TRUE,ND,LST_DEPTH_ONLY)
	      END DO
	    END IF
	  END DO
C
C 
C
C Add in 2-photon emissivity and opacity.
C
	  IF(GLOBAL_LINE_SWITCH .EQ. 'BLANK')THEN
	    CALL TWO_PHOT_OPAC(ETA,CHI,POPS,T,CONT_FREQ,ND,NT)
	  END IF
C
C Compute X-ray opacities and emissivities (at local electron temp).
C
	  IF(XRAYS)THEN
	    INCLUDE 'XOPAC_V4.INC'
	  END IF
C
	  CHI_C_EVAL(:)=CHI(:)
	  ETA_C_EVAL(:)=ETA(:)
C
	END IF
C
C 
C
C Section to revise continuum opacities etc so that they are computed at
C the correct frequency. We have stored the orginal continuum opacity and
C emissivity in CHI_C_EVAL and ETA_C_EVAL, which were computed at CONT_FREQ.
C
	IF(FL .NE. CONT_FREQ)THEN
C
	  DO J=1,ND
	    EMHNUKT_CONT(J)=EXP(-HDKT*CONT_FREQ/T(J))
	    EMHNUKT(J)=EXP(-HDKT*FL/T(J))
	  END DO
C
C Compute continuum intensity incident from the core assuming a TSTAR
C blackbody.
C
	  IC=TWOHCSQ*(FL**3)/(EXP(HDKT*FL/TSTAR)-1.0D0)
C
C We assume that the photoionization cross-section has not changed since the
C last iteration. Using the result that the stimulated emission occurs in
C LTE and is given by
C                     ETA/(2hv^3/c^2)
C we can adjust CHI and ETA so that the condition of constant photoionization
C cross-section is met. This adjustment automatically ensures that ETA/CHI 
C gives the Planck function in LTE. 
C
	  T1=(FL/CONT_FREQ)**3
	  T2=TWOHCSQ*(CONT_FREQ**3)
	  T3=TWOHCSQ*(FL**3)
	  DO J=1,ND
	    T4=ETA_C_EVAL(J)*T1*(EMHNUKT(J)/EMHNUKT_CONT(J))
	    CHI(J)=CHI_C_EVAL(J)+(ETA_C_EVAL(J)/T2-T4/T3)
	    ETA(J)=T4
	  END DO
	ELSE
C
C We reset CHI and ETA in case shock X-ray emission has been added to ETA.
C
	  CHI(1:ND)=CHI_C_EVAL(1:ND)
	  ETA(1:ND)=ETA_C_EVAL(1:ND)
	END IF
C
C 
C
C The shock emission is added separately since it does not occur at the
C local electron temperature.
C
	IF(XRAYS)THEN
	  IF(FILL_VEC_SQ(1) .EQ. 0)THEN
	    DO I=1,ND
	      FILL_VEC_SQ(I)=4.0*FILL_FAC_XRAYS*FILL_FAC_XRAYS*
	1                   EXP(-V_SHOCK/V(I))
	    END DO
	  END IF
	  T1=CHIFF*TWOHCSQ*EXP(-HDKT*CONT_FREQ/T_SHOCK)/SQRT(T_SHOCK)
C
C Since T_SHOCK is depth indpendent, we can compute the free-free
C Gaunt factors are depth independent.
C
C ####
C
	  DO I=1,20
	    T2=I
	    TA(I)=GFF(CONT_FREQ,T_SHOCK,T2)
	  END DO
C
	  DO I=1,ND
	    T2=0.0D0
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)T2=T2+ATM(ID)%ZXzV*ATM(ID)%ZXzV*
	1                ATM(ID)%DXzV(I)*TA(NINT(ATM(ID)%ZXzV))
	    END DO
	    ZETA(I)=T1*T2*FILL_VEC_SQ(I)*ED(I)	!Zeta is temporary
            ETA(I)=ETA(I)+ZETA(I)
          END DO
	  IF(SECTION .EQ. 'CONTINUUM')THEN
	    IF(ML .EQ. 1)XRAY_LUM(1:ND)=0.0D0
	    XRAY_LUM(1:ND)=XRAY_LUM(1:ND)+ZETA(1:ND)*FQW(ML)
	  END IF
	END IF
C
C The continuum source function is defined by:
C                                              S= ZETA + THETA.J
	DO I=1,ND
	  ZETA(I)=ETA(I)/CHI(I)
	  THETA(I)=ESEC(I)/CHI(I)
	END DO
C
C Store TOTAL continuum line emissivity and opacity.
C
	ETA_CONT(:)=ETA(:)
	CHI_CONT(:)=CHI(:)
