!
! Subroutine to compute how the continuum opacities and emissivities vary with
! the populations.
!
	SUBROUTINE COMP_VAR_OPAC(POPS,RJ,FL,CONT_FREQ,FREQ_INDX,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE OPAC_MOD
	IMPLICIT NONE
!
! Altered 05-Apr-2011 : Now call VAR_OP_V10 (6-Feb-2011)
! Created 26-Feb-2004 : Based on VAROPAC_V4
!
	INTEGER ND
	INTEGER NT
!
	REAL*8 POPS(NT,ND)
	REAL*8 RJ(ND)
	REAL*8 FL
	REAL*8 CONT_FREQ
	INTEGER FREQ_INDX
	CHARACTER*(*) SECTION
	LOGICAL LST_DEPTH_ONLY
!
! Constants for opacity etc.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables.
!
	REAL*8 T1,T2,T3,T4
	INTEGER J
	INTEGER K
	INTEGER L
	INTEGER ID
	INTEGER PHOT_ID
!
! Compute the opacity AND emissivity variation as a function of the changes
! in population levels. As the variation is linear, they can be added
! independently. Note that the VCHI and VETA are not zeroed when the
! variation routines are called.
!
	IF(COMPUTE_NEW_CROSS)THEN
!
! Zero VCHI and VETA, then allow for the variation of the electron
! scattering opacity and emissivity with ED.
!
! TA is used as a work vector.
!
	  CALL DP_ZERO(VCHI,NT*ND)
	  CALL DP_ZERO(VETA,NT*ND)
	  CALL DP_ZERO(VCHI_ALL,NT*ND)
	  CALL DP_ZERO(VETA_ALL,NT*ND)
!
! We use EMHNUKT_CONT since we are evaluating the variation in the opacities
! and emissivities at CONT_FREQ, not FL
!
	  DO L=1,ND
	    EMHNUKT_CONT(L)=EXP(-HDKT*CONT_FREQ/T(L))
	  END DO
!
! Parallelizing over ID appears to be inefficient (using reduction) since
! the whole arrays (VCHI, etc) must be added at the end of each call to VAR_op_V8.
!
	  DO ID=1,NUM_IONS
	    IF(ATM(ID)%XzV_PRES)THEN
	      DO J=1,ATM(ID)%N_XzV_PHOT
	        PHOT_ID=J
	        CALL VAR_OP_V10(VCHI,VETA,VCHI_ALL,VETA_ALL,
	1         ATM(ID)%XzV, ATM(ID)%XzVLTE, ATM(ID)%LOG_XzVLTE,  ATM(ID)%dlnXzVLTE_dlnT, ATM(ID)%NXzV,
	1         ATM(ID)%XzVLTE_F_ON_S, ATM(ID)%EDGEXzV_F, ATM(ID)%NXzV_F, ATM(ID)%F_TO_S_XzV,
	1         ATM(ID+1)%XzV,    ATM(ID+1)%LOG_XzVLTE,  ATM(ID+1)%dlnXzVLTE_dlnT,
	1         ATM(ID+1)%NXzV,   PHOT_ID,           ATM(ID)%XzV_ION_LEV_ID(J),
	1         ED,T,EMHNUKT_CONT,IMP_VAR,
	1         CONT_FREQ,ATM(ID)%ZXzV,ID,L_TRUE,
	1         ATM(ID)%EQXzV,ATM(ID+1)%EQXzV,NT,ND,LST_DEPTH_ONLY)
	      END DO
	    END IF
	 END DO
!
	IF(ADD_ADDITIONAL_OPACITY)THEN
	  DO L=1,ND
	    T1=ADD_OPAC_SCL_FAC*6.65D-15*POP_ATOM(L)
	    T2=T1*HDKT*CONT_FREQ/T(L)/T(L)
	    VETA(NT,L)=VETA(NT,L)+T2*TWOHCSQ*(CONT_FREQ**3)*EMHNUKT_CONT(L)/
	1                   (1.0D0-EMHNUKT_CONT(L))**2
	  END DO
	END IF
!
!
! 
!
! Add in 2-photon opacities and emissivities.
!
	  CALL TWO_PHOT_VAR_OPAC(VETA,VCHI,POPS,T,CONT_FREQ,ND,NT)
!
! Altered 04-Mar-2004: Call inserted directly into subroutine. No longer done
!                        as include file.
! Altered 21-Jun-2002: Changed to V4 and inserted LST_DEPTH_ONLY in call.
!
	  IF(XRAYS)THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	        CALL VAR_X_OPA_ETA_V4(VCHI,VETA,
	1          ATM(ID)%XzV,ATM(ID)%XzVLTE,ATM(ID)%dlnXzVLTE_dlnT,
	1          ATM(ID)%NXzV,
	1          ATM(ID+1)%XzV, ATM(ID+1)%XzVLTE,
	1          ATM(ID+1)%dlnXzVLTE_dlnT,
	1          ATM(ID+1)%NXzV,ED, ATM(ID+1)%DXzV,
	1          T,IMP_VAR,ATM(ID)%EQXzV, ATM(ID+2)%EQXzv,
	1          AT_NO(SPECIES_LNK(ID)), ATM(ID)%ZXzV,
	1          CONT_FREQ,EMHNUKT_CONT,NT,ND,LST_DEPTH_ONLY)
	      END IF
	    END DO
	  END IF
!
! _ALL will contain the variation of ETA/CHI for ALL species.
!
	  IF(LST_DEPTH_ONLY)THEN
!
!$OMP PARALLEL WORKSHARE
!
	    VCHI_ALL(:,ND)=VCHI_ALL(:,ND)+VCHI(:,ND)
	    VETA_ALL(:,ND)=VETA_ALL(:,ND)+VETA(:,ND)
	    VCHI_SAV(:,ND)=VCHI(:,ND)
	    VETA_SAV(:,ND)=VETA(:,ND)
	    VCHI_ALL_SAV(:,ND)=VCHI_ALL(:,ND)
	    VETA_ALL_SAV(:,ND)=VETA_ALL(:,ND)
!
!$OMP END PARALLEL WORKSHARE
!
	  ELSE
!
!$OMP PARALLEL WORKSHARE
!
	    VCHI_ALL=VCHI_ALL+VCHI
	    VETA_ALL=VETA_ALL+VETA
	    VCHI_SAV=VCHI
	    VETA_SAV=VETA
	    VCHI_ALL_SAV=VCHI_ALL
	    VETA_ALL_SAV=VETA_ALL
!
!$OMP END PARALLEL WORKSHARE
!
	  END IF
	ELSE IF(CONT_FREQ .EQ. FL)THEN
!
!$OMP PARALLEL WORKSHARE
!
	    VCHI=VCHI_SAV
	    VETA=VETA_SAV
	    VCHI_ALL=VCHI_ALL_SAV
	    VETA_ALL=VETA_ALL_SAV
!
!$OMP END PARALLEL WORKSHARE
!
	END IF
!
! 
!
! This section of code can be utilized by both the DTDR and CONTINUUM code
! sections.
!
  	IF(CONT_FREQ .NE. FL)THEN
!
! Scale dCHI and dETA to allow for the slight variation in frequency.
! Correction should generally be small.
!
! NB: Previously T4 was defined by T4= EMHNUKT(L)/EMHNUKT_CONT(L). In the
! presence of X-rays, and with NU=1000 or larger, EMHNUKT_CONT could be
! zero, thus causing a divide by zero.
!
	  T1=(FL/CONT_FREQ)**3
	  T2=TWOHCSQ*(CONT_FREQ**3)
	  T3=TWOHCSQ*(FL**3)
!
	  IF(LST_DEPTH_ONLY)THEN
	    L=ND
	    T4=T1*EXP(-HDKT*(FL-CONT_FREQ)/T(L))
            DO K=1,NT
               VETA(K,L)=VETA_SAV(K,L)*T4
               VETA_ALL(K,L)=VETA_ALL_SAV(K,L)*T4
            END DO
!
! Need to correct VETA for the variation in T in the factor T4.
! Over correction at present because impurity species included.
!
	    VETA(NT,L)=VETA(NT,L) + 
	1                ETA_C_EVAL(L)*T4*HDKT*(FL-CONT_FREQ)/T(L)/T(L)
	    VETA_ALL(NT,L)=VETA_ALL(NT,L) + 
	1                ETA_C_EVAL(L)*T4*HDKT*(FL-CONT_FREQ)/T(L)/T(L)
!
! Can now correct the opacity.
!
	     DO K=1,NT
               VCHI(K,L)=VCHI_SAV(K,L)+(VETA_SAV(K,L)/T2-VETA(K,L)/T3)
               VCHI_ALL(K,L)=VCHI_ALL_SAV(K,L)+(VETA_ALL_SAV(K,L)/T2-VETA_ALL(K,L)/T3)
	     END DO
	  ELSE
!
! NB: We also need to correct VETA for the variation in T in the factor T4.
! Over correction at present because impurity species included.
!
!$OMP PARALLEL DO PRIVATE(T4)
	    DO L=1,ND
	      T4=T1*EXP(-HDKT*(FL-CONT_FREQ)/T(L))
              DO K=1,NT
                VETA(K,L)=VETA_SAV(K,L)*T4
              END DO
	      VETA(NT,L)=VETA(NT,L) + 
	1                  ETA_C_EVAL(L)*T4*HDKT*(FL-CONT_FREQ)/T(L)/T(L)
	    END DO
!
!$OMP PARALLEL DO PRIVATE(T4)
	    DO L=1,ND
	      T4=T1*EXP(-HDKT*(FL-CONT_FREQ)/T(L))
              DO K=1,NT
                VETA_ALL(K,L)=VETA_ALL_SAV(K,L)*T4
              END DO
	      VETA_ALL(NT,L)=VETA_ALL(NT,L) + 
	1                  ETA_C_EVAL(L)*T4*HDKT*(FL-CONT_FREQ)/T(L)/T(L)
	    END DO
!
! Can now correct the opacity.
!
!$OMP  PARALLEL WORKSHARE
!
            VCHI=VCHI_SAV+(VETA_SAV/T2-VETA/T3)
            VCHI_ALL=VCHI_ALL_SAV+(VETA_ALL_SAV/T2-VETA_ALL/T3)
!
!$OMP END PARALLEL WORKSHARE
!
	  END IF  
!
	END IF			!Correct cross-section?
!
! 
!
!
! In the CONTINUUM code section (with blanketing) the variation of the opacity 
! and emissivity due to the electron scattering term is handled separately.
! It is required in the DTDR section, DIELECTRONIC and LINE sections.
! In DTDR the cross-sections can also be held fixed. In order that our
! correction process implemented above works, we must keep the electron
! scattering cross-section separate (primarily in the emissivity).
!
	IF(SECTION .NE. 'CONTINUUM')THEN
	  IF(ATM(1)%XzV_PRES)THEN
	    DO K=1,ND
	      T1=CHI_RAY(K)/ATM(1)%XzV_F(1,K)
	      VCHI(1,K)=VCHI(1,K)+T1
	      VETA(1,K)=VETA(1,K)+T1*RJ(K)
	  END DO
	  END IF
	  DO K=1,ND
	    VCHI(NT-1,K)=VCHI(NT-1,K)+ESEC(K)/ED(K)
	    VETA(NT-1,K)=VETA(NT-1,K)+ESEC(K)*RJ(K)/ED(K)
	  END DO
	END IF
!
	RETURN
	END
