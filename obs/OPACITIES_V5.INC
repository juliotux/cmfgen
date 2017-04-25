!
! Altered 17-Dec-2011: Changed opcity call to GENOPAETA_V10. This contains LOG_XzVLTE_F
!                        which is designed to allow lower wind temperatures.
!                        Altered computation of IC.
! Altered 15-Aug-2002: Minimum emissivity set to avoid floating point underflow.
!
!
! Compute opacity and emissivity. This is a general include file
! provided program uses exactly the same variables. Can be achieved
! by copying declaration statements from CMFGEN. Always advisable
! to use ``IMPLICIT NONE''.
!
	IF(COMPUTE_NEW_CROSS)THEN
!
! Compute EXP(-hv/kT) and zero CHI, and ETA.
!
	  T1=-HDKT*CONT_FREQ
	  DO I=1,ND
	    EMHNUKT(I)=EXP(T1/T(I))
	    CHI(I)=0.0D0
	    ETA(I)=0.0D0
	  END DO
!
! Compute continuum intensity incident from the core assuming a TSTAR
! blackbody.
!
	  T1=EXP(-HDKT*CONT_FREQ/TSTAR)
	  IC=TWOHCSQ*T1*(CONT_FREQ**3)/(1.0D0-T1)
!
! Compute opacity and emissivity. ESOPAC must be call first since
! CHI is effectively zeroed in that routine.
!
	  CALL ESOPAC(ESEC,ED,ND)		!Electron scattering emission factor.
!
! Add in Rayleigh scattering contribution. We set the Rayleigh vector to a negligble,
! but non-zero value, since it can then be interpolated in the LOG plane.
!
	  T1=1.0D+06/HUGE(T1)
          CHI_RAY(1:ND)=T1
          IF(ATM(1)%XzV_PRES .AND. INCL_RAY_SCAT)THEN
	    CALL RAYLEIGH_SCAT(CHI_RAY,ATM(1)%XzV_F,ATM(1)%AXzV_F,ATM(1)%EDGEXZV_F,
	1             ATM(1)%NXzV_F,CONT_FREQ,ND)
	  END IF
	  CHI_SCAT(1:ND)=ESEC(1:ND)+CHI_RAY(1:ND)
	  CHI(1:ND)=CHI_SCAT(1:ND)
!
! Free-free and bound-free opacities.
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) REDUCTION(+:CHI,ETA) PRIVATE(ID,J)
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        CALL GENOPAETA_V10(ID,CHI,ETA,CONT_FREQ,
	1           ATM(ID)%XzV_F,     ATM(ID)%XzVLTE_F,     ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%EDGEXzV_F,
	1           ATM(ID)%GIONXzV_F, ATM(ID)%ZXzV,         ATM(ID)%NXzV_F,
	1           ATM(ID+1)%XzV,     ATM(ID+1)%LOG_XzVLTE, ATM(ID+1)%NXzV, J,
	1           ATM(ID)%XzV_ION_LEV_ID(J),
	1           ED,T,EMHNUKT,L_TRUE,ND,LST_DEPTH_ONLY)
	      END DO
	    END IF
	  END DO
!$OMP END PARALLEL DO
!
! 
!
! Add in 2-photon emissivity and opacity.
!
	  IF(COMPUTE_EDDFAC .AND. TWO_PHOTON_METHOD .EQ. 'USE_RAD')THEN
	    CALL TWO_PHOT_OPAC_V3(ETA,CHI,POPS,T,CONT_FREQ,'OLD_DEFAULT',ND,NT)
	  ELSE
	    CALL TWO_PHOT_OPAC_V3(ETA,CHI,POPS,T,CONT_FREQ,TWO_PHOTON_METHOD,ND,NT)
	  END IF
!
! Compute X-ray opacities and emissivities (at local electron temp).
!
	  IF(XRAYS)THEN
	    INCLUDE 'XOPAC_V4.INC'
	  END IF
!
	  CHI_C_EVAL(:)=CHI(:)
	  ETA_C_EVAL(:)=ETA(:)
!
	END IF
!
! 
!
! Section to revise continuum opacities etc so that they are computed at
! the correct frequency. We have stored the orginal continuum opacity and
! emissivity in CHI_C_EVAL and ETA_C_EVAL, which were computed at CONT_FREQ.
!
	IF(FL .NE. CONT_FREQ)THEN
!
	  DO J=1,ND
	    EMHNUKT_CONT(J)=EXP(-HDKT*CONT_FREQ/T(J))
	    EMHNUKT(J)=EXP(-HDKT*FL/T(J))
	  END DO
!
! Compute continuum intensity incident from the core assuming a TSTAR
! blackbody.
!
	  T1=EXP(-HDKT*FL/TSTAR)
	  IC=TWOHCSQ*T1*(FL**3)/(1.0D0-T1)
!
! We assume that the photoionization cross-section has not changed since the
! last iteration. Using the result that the stimulated emission occurs in
! LTE and is given by
!                     ETA/(2hv^3/c^2)
! we can adjust CHI and ETA so that the condition of constant photoionization
! cross-section is met. This adjustment automatically ensures that ETA/CHI
! gives the Planck function in LTE.
!
	  T1=(FL/CONT_FREQ)**3
	  T2=TWOHCSQ*(CONT_FREQ**3)
	  T3=TWOHCSQ*(FL**3)
	  DO J=1,ND
	    T4=ETA_C_EVAL(J)*T1*EXP(-HDKT*(FL-CONT_FREQ)/T(J))
	    CHI(J)=CHI_C_EVAL(J)+(ETA_C_EVAL(J)/T2-T4/T3)
	    ETA(J)=T4
	  END DO
	ELSE
!
! We reset CHI and ETA in case shock X-ray emission has been added to ETA.
!
	  CHI(1:ND)=CHI_C_EVAL(1:ND)
	  ETA(1:ND)=ETA_C_EVAL(1:ND)
	END IF
!
! 
!
! The shock emission is added separately since it does not occur at the
! local electron temperature.
!
	IF(XRAYS)THEN

	  IF(FF_XRAYS)THEN

! Since T_SHOCK is depth indpendent, Z^2 * (the free-free Gaunt factors)
! are depth independent.
!
! We use T3 for the Electron density. We asume H, He, and C are fully ionized
! in the X-ray emitting plasma. All other species are assumed to have Z=6.0
!
	    DO I=1,ND
	      ZETA(I)=0.D0	 !Zeta is temporary
	    END DO
	    IF((T_SHOCK_1.GT.0.0D0).AND.(FILL_FAC_XRAYS_1.GT.0.D0))THEN
	      T1=CHIFF*TWOHCSQ*(FILL_FAC_XRAYS_1**2)*
	1                    EXP(-HDKT*CONT_FREQ/T_SHOCK_1)/SQRT(T_SHOCK_1)
	      T2=1.0D0 ; TA(1)=GFF(CONT_FREQ,T_SHOCK_1,T2)
	      T2=2.0D0 ; TA(2)=4.0D0*GFF(CONT_FREQ,T_SHOCK_1,T2)
	      T2=6.0D0 ; TA(3)=36.0D0*GFF(CONT_FREQ,T_SHOCK_1,T2)
	      DO I=1,ND
	        T2=TA(1)*POP_SPECIES(I,1)+TA(2)*POP_SPECIES(I,2) +
	1         TA(3)*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        T3=POP_SPECIES(I,1)+POP_SPECIES(I,2) +
	1          6.0D0*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        ZETA(I)=ZETA(I)+T1*T2*T3*EXP(-V_SHOCK_1/V(I))	 !Zeta is temporary
              END DO
	    END IF
	    IF((T_SHOCK_2 .GT. 0.0D0).AND.(FILL_FAC_XRAYS_2.GT.0.D0))THEN
	      T1=CHIFF*TWOHCSQ*(FILL_FAC_XRAYS_2**2)*
	1                    EXP(-HDKT*CONT_FREQ/T_SHOCK_2)/SQRT(T_SHOCK_2)
	      T2=1.0D0 ; TA(1)=GFF(CONT_FREQ,T_SHOCK_2,T2)
	      T2=2.0D0 ; TA(2)=4.0D0*GFF(CONT_FREQ,T_SHOCK_2,T2)
	      T2=6.0D0 ; TA(3)=36.0D0*GFF(CONT_FREQ,T_SHOCK_2,T2)
	      DO I=1,ND
	        T2=TA(1)*POP_SPECIES(I,1)+TA(2)*POP_SPECIES(I,2) +
	1       TA(3)*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        T3=POP_SPECIES(I,1)+POP_SPECIES(I,2) +
	1          6.0D0*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	        ZETA(I)=ZETA(I)+T1*T2*T3*EXP(-V_SHOCK_2/V(I))	 !Zeta is temporary
            END DO
	    END IF

	  ELSE

!
! Use X-ray emission as tubulated by a PLASMA code. Emission should  be
! tabulated per electron and per ion.
!
	    CALL GET_SCL_XRAY_FLUXES_V1(CONT_FREQ,
	1              XRAY_EMISS_1,XRAY_EMISS_2,
	1              NU_EVAL_CONT,NCF,FREQ_INDX,
	1              VSMOOTH_XRAYS,SECTION)
!
! We use T3 for the Electron density. We asume H, He, and C are fully ionized
! in the X-ray emitting plasma. All other species are assumed have Z=6.0
	    DO I=1,ND
	    T1=0.D0
	    T2=0.D0
	    ZETA(I)=0.D0
	    IF(FILL_FAC_XRAYS_1.GT.0.D0) T1=EXP(-V_SHOCK_1/V(I))*(FILL_FAC_XRAYS_1)**2
	    IF(FILL_FAC_XRAYS_2.GT.0.D0) T2=EXP(-V_SHOCK_2/V(I))*(FILL_FAC_XRAYS_2)**2
	    T3=POP_SPECIES(I,1)+2.0D0*POP_SPECIES(I,2)+
	1            6.0D0*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	    ZETA(I)=(T1*XRAY_EMISS_1+T2*XRAY_EMISS_2)*T3*POP_ATOM(I)
	    END DO

	  END IF

	  IF(XRAY_SMOOTH_WIND)ZETA(1:ND)=ZETA(1:ND)*CLUMP_FAC(1:ND)
	  ETA(1:ND)=ETA(1:ND)+ZETA(1:ND)

	  IF(SECTION .EQ. 'CONTINUUM')THEN
	    IF(ML .EQ. 1)XRAY_LUM(1:ND)=0.0D0
	    TA(1:ND)=ZETA(1:ND)*CLUMP_FAC(1:ND)*FQW(FREQ_INDX)
	    XRAY_LUM(1:ND)=XRAY_LUM(1:ND)+TA(1:ND)
	  END IF
	END IF
!
! Set a minimum emissivity to avoid floating point underflow.
! Mainly important when X-rays not present.
!
          DO I=1,ND
            IF(ETA(I) .LT. 1.0D-280)ETA(I)=1.0D-280
          END DO
!
! The continuum source function is defined by:
!     S= ZETA + THETA.J
	DO I=1,ND
	  ZETA(I)=ETA(I)/CHI(I)
	  THETA(I)=CHI_SCAT(I)/CHI(I)
	END DO
!
! Store TOTAL continuum line emissivity and opacity.
!
	ETA_CONT(:)=ETA(:)
	CHI_CONT(:)=CHI(:)
!
