!
! Subroutine to compute the opacities and emissivities for CMFGEN.
!
	SUBROUTINE COMP_OPAC(POPS,NU_EVAL_CONT,FQW,
	1                FL,CONT_FREQ,FREQ_INDX,NCF,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE OPAC_MOD
	IMPLICIT NONE
!
! Altered 05-Apr-2011 : Now call GENOPAETA_V10 (6-Feb-2011)
! Altered 11-Jun-2006: Installed CHI_NOSCAT and ETA_NOSCAT. Scattering
!                        opacity now computed after bound-free, free-free,
!                        two photon, and X-ray opacity have been computed.
! Altered 11-Jun-2005: Bug fix: X-ray opacity/emissivity was not being
!                        added correctly.
! Altered 13-Sep-2004: Installed ETA_MECH to keep track of mechanical
!                        energy input.
! Altered 03-Mar-2004: Bug fix. EMHNUKT was not being computed when
!                        COMPUTE_NEW_CROSS=.FALSE. and NU=NU_CONT.
!                        Value at earlier frequncy was being used.
!                        Computation ox X-ray cross-sections now
!                        directly included (not INCLUDE file).
! Created 16-Feb-2004: Based on OPACITIES_V4
!
	INTEGER ND
	INTEGER NT
	INTEGER NCF
!
	REAL*8 POPS(NT,ND)
	REAL*8 NU_EVAL_CONT(NCF)
	REAL*8 FQW(NCF)
	REAL*8 FL
	REAL*8 CONT_FREQ
!
	INTEGER FREQ_INDX
	CHARACTER*(*) SECTION
	LOGICAL LST_DEPTH_ONLY
!
! Constants for opacity etc.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Internally used variables
!
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL*8 XCROSS_V2
	REAL*8 GFF
	EXTERNAL XCROSS_V2,GFF
!
	REAL*8 TA(ND)
	REAL*8 T1,T2,T3,T4
	INTEGER I,J
	INTEGER ID
	INTEGER PHOT_ID
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
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
	    EMHNUKT_CONT(I)=EXP(T1/T(I))
	    CHI(I)=0.0D0
	    ESEC(I)=0.0D0
	    ETA(I)=0.0D0
	  END DO
!
! Compute continuum intensity incident from the core assuming a TSTAR
! blackbody.
!
          T1=EXP(-HDKT*CONT_FREQ/TSTAR)
	  IC=TWOHCSQ*T1*(CONT_FREQ**3)/(1.0D0-T1)
!
! Free-free and bound-free opacities.
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) REDUCTION(+:CHI,ETA) PRIVATE(ID,J,PHOT_ID)
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        PHOT_ID=J
	        CALL GENOPAETA_V10(ID,CHI,ETA,CONT_FREQ,
	1           ATM(ID)%XzV_F,      ATM(ID)%XzVLTE_F,     ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%EDGEXzV_F,
	1           ATM(ID)%GIONXzV_F,  ATM(ID)%ZXzV,         ATM(ID)%NXzV_F,
	1           ATM(ID+1)%XzV,      ATM(ID+1)%LOG_XzVLTE, ATM(ID+1)%NXzV, 
	1           PHOT_ID,            ATM(ID)%XzV_ION_LEV_ID(J),
	1           ED,T,EMHNUKT_CONT,L_TRUE,ND,LST_DEPTH_ONLY)
	      END DO
	    END IF
	  END DO
!$OMP END PARALLEL DO
!
	  IF(ADD_ADDITIONAL_OPACITY)THEN
	     DO I=1,ND
	       T1=ADD_OPAC_SCL_FAC*6.65D-15*POP_ATOM(I)
	       CHI(I)=CHI(I)+T1
	       ETA(I)=ETA(I)+T1*TWOHCSQ*(CONT_FREQ**3)*EMHNUKT_CONT(I)/(1.0D0-EMHNUKT_CONT(I))
	     END DO
	     IF(FREQ_INDX .EQ. 1)THEN
	       LUER=ERROR_LU()
	       WRITE(LUER,*)'Warning --- adding additional opacity to model'
	       WRITE(LUER,*)'Remember to remove additional opacity for final model'
	     END IF
	  END IF
!
! 
!
! Add in 2-photon emissivity and opacity.
!
	  IF(LST_DEPTH_ONLY)THEN
	    CALL TWO_PHOT_OPAC_V3(ETA,CHI,POPS,T,CONT_FREQ,'LTE',ND,NT)
	  ELSE IF(COMPUTE_EDDFAC .AND. TWO_PHOTON_METHOD .EQ. 'USE_RAD')THEN
	    CALL TWO_PHOT_OPAC_V3(ETA,CHI,POPS,T,CONT_FREQ,'OLD_DEFAULT',ND,NT)
	  ELSE
	    CALL TWO_PHOT_OPAC_V3(ETA,CHI,POPS,T,CONT_FREQ,TWO_PHOTON_METHOD,ND,NT)
	  END IF
!
! Compute X-ray opacities and emissivities due to K (& L) shell ionization. In all cases
! it is assumed that 2 electrons are ejected. The K (& L) shell cross-sections are
! assumed to be independent of the level of the valence electron. In practice,
! ionizations will generally be determined by the population of the ground
! configuration.
!
! Since the cross-sections are level independent, we can sum over the levels before
! we add the contribution to the opacity/emissivity.
!
! This section was originaly XOPAC_V4.INC. See that file for erlier corrections.
!
	  IF(XRAYS)THEN
!!$OMP PARALLEL DO SCHEDULE(DYNAMIC) REDUCTION(CHI,ETA) PRIVATE(ID,J,PHOT_ID)
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	        T2=AT_NO(SPECIES_LNK(ID))+1-ATM(ID)%ZXzV
	        T1=XCROSS_V2(CONT_FREQ,AT_NO(SPECIES_LNK(ID)),T2,IZERO,IZERO,L_FALSE,L_FALSE)
	        IF(T1 .NE. 0.0D0)THEN
	          J=1
	          IF(LST_DEPTH_ONLY)J=ND
	          DO I=J,ND
	            T2=0.0D0			!Temporary CHI
	            T3=0.0D0			!Temporary ETA
	            T4=(ATM(ID+1)%XzVLTE_F(1,I)*EMHNUKT_CONT(I))/ATM(ID+1)%XzV_F(1,I)
	            DO J=1,ATM(ID)%NXzV_F
		      T2=T2+ATM(ID)%XzV_F(J,I)
	              T3=T3+ATM(ID)%XzVLTE_F(J,I)
	            END DO
	            CHI(I)=CHI(I)+T1*(T2-T3*T4)
	            ETA(I)=ETA(I)+T1*T3*T4*TWOHCSQ*(CONT_FREQ**3)
	          END DO
	        END IF
	      END IF
	    END DO
!!$OMP END PARALLEL DO
	  END IF
!
	  CHI_NOSCAT(1:ND)=CHI(1:ND)
	  ETA_NOSCAT(1:ND)=ETA(1:ND)
!
! Compute scattering opacity. ESEC is zeroed in ESOPAC.
!
	  CALL ESOPAC(ESEC,ED,ND)		!Electron scattering emission factor.
!
! Add in Rayleigh scattering contribution.
!
	  CHI_RAY(1:ND)=0.0D0
	  IF(SPECIES_PRES(1) .AND. INCL_RAY_SCAT)THEN
	    CALL RAYLEIGH_SCAT(CHI_RAY,ATM(1)%XzV_F,ATM(1)%AXzV_F,ATM(1)%EDGEXZV_F,
	1             ATM(1)%NXzV_F,CONT_FREQ,ND)
	  END IF
	  CHI_SCAT(1:ND)=ESEC(1:ND)+CHI_RAY(1:ND)
!
! Now compute total opacity --- scattering + non scattering.
!
	  CHI(1:ND)=CHI(1:ND)+CHI_SCAT(1:ND)
!
	  CHI_C_EVAL(:)=CHI(:)
	  ETA_C_EVAL(:)=ETA(:)
	  CHI_NOSCAT_EVAL(:)=CHI_NOSCAT(:)
	  ETA_NOSCAT_EVAL(:)=ETA_NOSCAT(:)
!
	END IF
!
! 
!
! Evaluate EXP(-hv/kT) for current frequency. This is needed by routines 
! such as COMP_VAR_JREC etc.
!
	  DO J=1,ND
	    EMHNUKT(J)=EXP(-HDKT*FL/T(J))
	  END DO
!
! Section to revise continuum opacities etc so that they are computed at
! the correct frequency. We have stored the orginal continuum opacity and
! emissivity in CHI_C_EVAL and ETA_C_EVAL, which were computed at CONT_FREQ.
!
	IF(FL .NE. CONT_FREQ)THEN
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
	    CHI_NOSCAT(J)=CHI_NOSCAT_EVAL(J)+(ETA_C_EVAL(J)/T2-T4/T3)
	    ETA_NOSCAT(J)=T4
	  END DO
!
	ELSE
!
! We reset CHI and ETA in case shock X-ray emission has been added to ETA,
! or CONT_FREQ was not the first frequency.
!
	  CHI(1:ND)=CHI_C_EVAL(1:ND)
	  ETA(1:ND)=ETA_C_EVAL(1:ND)
	  CHI_NOSCAT(1:ND)=CHI_NOSCAT_EVAL(1:ND)
	  ETA_NOSCAT(1:ND)=ETA_NOSCAT_EVAL(1:ND)
	END IF
!
! 
!
! The shock emission is added separately since it does not occur at the
! local electron temperature.
!
	IF(XRAYS)THEN
!
	  IF(FF_XRAYS)THEN
!
! Since T_SHOCK is depth indpendent, Z^2 * (the free-free Gaunt factors) 
! are depth independent.
!
! We use T3 for the Electron density. We asume H, He, and C are fully ionized
! in the X-ray emitting plasma. All other species are assumed to have Z=6.0
!
	    IF(T_SHOCK_1 .NE. 0.0D0)THEN
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
	        ZETA(I)=T1*T2*T3*EXP(-V_SHOCK_1/V(I))	 !Zeta is temporary
              END DO
	    END IF
	    IF(T_SHOCK_2 .NE. 0.0D0)THEN
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
	        ZETA(I)=T1*T2*T3*EXP(-V_SHOCK_2/V(I))	 !Zeta is temporary
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
!
	    DO I=1,ND
	      T1=EXP(-V_SHOCK_1/V(I))*(FILL_FAC_XRAYS_1)**2
	      T2=EXP(-V_SHOCK_2/V(I))*(FILL_FAC_XRAYS_2)**2
	      T3=POP_SPECIES(I,1)+2.0D0*POP_SPECIES(I,2)+
	1              6.0D0*(POP_ATOM(I)-POP_SPECIES(I,1)-POP_SPECIES(I,2))
	      ZETA(I)=(T1*XRAY_EMISS_1+T2*XRAY_EMISS_2)*T3*POP_ATOM(I)
	    END DO
	  END IF
!
	  IF(XRAY_SMOOTH_WIND)ZETA(1:ND)=ZETA(1:ND)*CLUMP_FAC(1:ND)
          ETA(1:ND)=ETA(1:ND)+ZETA(1:ND)
	  ETA_MECH(1:ND)=ZETA(1:ND)
!
! Changed 06-Aug-2003: Clumping was not beeing allowed for when computing
! the shock luminosity.
!
	  T1=0.241838D0		!eV to 10^15Hz
	  IF(SECTION .EQ. 'CONTINUUM')THEN
	    IF(FREQ_INDX .EQ. 1)THEN
	       XRAY_LUM_TOT(1:ND)=0.0D0
	       XRAY_LUM_0P1(1:ND)=0.0D0
	       XRAY_LUM_1KEV(1:ND)=0.0D0
	    END IF
	    TA(1:ND)=ZETA(1:ND)*CLUMP_FAC(1:ND)*FQW(FREQ_INDX)
	    XRAY_LUM_TOT(1:ND)=XRAY_LUM_TOT(1:ND)+TA(1:ND)
	    IF(FL .GE. 100.0D0*T1)XRAY_LUM_0P1(1:ND)=XRAY_LUM_0P1(1:ND)+TA(1:ND)
	    IF(FL .GE. 1000.0D0*T1)XRAY_LUM_1KEV(1:ND)=XRAY_LUM_1KEV(1:ND)+TA(1:ND)
	  END IF
	ELSE
	  ETA_MECH(1:ND)=0.0D0
	END IF
!
! Set a minimum emissivity. Mainly important when X-rays are not present.
!
	DO I=1,ND
	  IF(ETA(I) .LT. 1.0D-280)THEN
	    ETA(I)=1.0D-280
	    ETA_NOSCAT(I)=1.0D-280
	  END IF
	END DO
!
! The continuum source function is defined by:
!                                              S= ZETA + THETA.J
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
	RETURN
	END
