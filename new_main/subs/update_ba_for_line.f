	SUBROUTINE UPDATE_BA_FOR_LINE(FL,FQW,FREQ_INDX,
	1              POPS,JREC,dJRECdT,JPHOT,
	1              ND,NT,NUM_BNDS,NION,DIAG_INDX,
	1              TX_OFFSET,MAX_SIM,NM,NCF,NLF,
	1              LUER,FINAL_CONSTANT_CROSS,LST_ITERATION)
	USE MOD_CMFGEN
	USE CONTROL_VARIABLE_MOD
	USE VAR_RAD_MOD
	USE OPAC_MOD
	USE STEQ_DATA_MOD
	USE LINE_VEC_MOD
	USE LINE_MOD
	USE RADIATION_MOD
	IMPLICIT NONE
!
! Altered: 04-Oxt-2016 : Now call VSEBYJ_MULTI_V8 and VSEBYJ_X_V7.
! Altered: 06-Feb-2010 : Now call VSEBYJ_MULTI_V7.
!
        INTEGER ND
        INTEGER NT
        INTEGER NUM_BNDS
	INTEGER NION
	INTEGER DIAG_INDX
        INTEGER NM
	INTEGER MAX_SIM
	INTEGER NCF
	INTEGER NLF
        INTEGER TX_OFFSET
	INTEGER FREQ_INDX
	INTEGER LUER
	LOGICAL FINAL_CONSTANT_CROSS
	LOGICAL LST_ITERATION
	LOGICAL UPDATE_dZ
!
        REAL*8 FL
	REAL*8 NU
	REAL*8 FQW
	REAL*8 POPS(NT,ND)
	REAL*8 JREC(ND)
	REAL*8 dJRECdT(ND)
	REAL*8 JPHOT(ND)
!
        INTEGER NL,NUP
        INTEGER MNL,MNUP
        INTEGER MNL_F,MNUP_F
        INTEGER DPTH_INDX
        INTEGER VAR_INDX
        INTEGER INDX_BA_METH
!
	REAL*8 VB(ND),VC(ND)
        REAL*8 T1,T2,T3,T4
	REAL*8 SCL_FAC
	REAL*8, SAVE :: SUM_BA
	REAL*8, SAVE :: FL_OLD
!
	INTEGER I,J,K,L
	INTEGER JJ
	INTEGER ID
	INTEGER ID_SAV
	INTEGER DST,DEND
!
! Constants for opacity etc. These are set in CMFGEN.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 OPLIN,EMLIN
!
!***********************************************************************
!
!*******************FUNCTION DEFINITIONS********************************
!
        INTEGER GET_DIAG
        INTEGER BNDST
        INTEGER BNDEND
        INTEGER BND_TO_FULL
!
! This function takes a band-index and converts it the equivalent index
! in the full matrix. L=BND_TO_FULL(J,K) is equivalent to the statements:
!     IF(NUM_BNDS .EQ. ND)THEN L=J ELSE L=K+J-DIAG_INDX END IF
! The second indice is the equation depth.
!
	BND_TO_FULL(J,K)=(NUM_BNDS/ND)*(DIAG_INDX-K)+K+J-DIAG_INDX
!
! This function computes the index L on BA( , ,?,K) corresponding
! to the local depth variable (i.e that at K). It is equivalent
! to IF (NUM_BNDS .EQ. ND)THEN L=K ELSE L=DIAG END IF
!
	GET_DIAG(K)=(NUM_BNDS/ND)*(K-DIAG_INDX)+DIAG_INDX
!
! These two functions compute the start and end indices when updating
! VJ. eg. we do not wish to update VJ( ,1, ) if we are using the banded
! matrix since this refers to a variable beyond the outer atmosphere.
!
	BNDST(K)=MAX( (NUM_BNDS/ND)*(K-DIAG_INDX)+1+DIAG_INDX-K, 1 )
	BNDEND(K)=MIN( (NUM_BNDS/ND)*(K-DIAG_INDX)+ND+DIAG_INDX-K, NUM_BNDS )
!
! 
!
	UPDATE_dZ=.FALSE.
	DPTH_INDX=6
	DPTH_INDX=MIN(DPTH_INDX,ND)             !Thus no problem if 84 > ND
	VAR_INDX=NT
	VAR_INDX=MIN(VAR_INDX,NT)
	IF(FREQ_INDX .EQ. 1)THEN
	  FL_OLD=FL
	  SUM_BA=0.0D0
	END IF
!
!
! Modify the BA matrix for terms in the statistical equilibrium
! equations which are multiplied by RJ. NB. This is not for the
! variation of RJ - rather the multiplying factors. This section
! must be done for a LAMBDA iteration.
!
! We use DST and DEND to avoid reading in the entire diagonal of the BA
! array for each call to BA. Calling the routines ND times may take
! too long --- therefore try alternative method.
!
! WE have replaced BA and BAION by BA_PAR and BAION_PAR in calls to
! VSEBYJ. Since these are diagonal components, we have had to replace 
! NUM_BANDS by IONE in the calls.
!
	IF(COMPUTE_BA)THEN
	  CALL TUNE(IONE,'COMPUTE_BA')
!
!	  DO K=1,ND
!	    DST=K
!	    DEND=K
	  DO K=1,1
	    DST=1       
	    DEND=ND
!
! NB: In this equation the matrices are NOT passed for WS, dWS and NU.
!
	    IF(FINAL_CONSTANT_CROSS)THEN
	      CALL TUNE(1,'BA_CONT_UP')
	      DO ID=1,NUM_IONS-1
	        ID_SAV=ID
	        IF(ATM(ID)%XzV_PRES)THEN
	          DO J=1,ATM(ID)%N_XzV_PHOT
	            CALL VSEBYJ_MULTI_V8(ID_SAV,
	1             ATM(ID)%WSXzV(1,1,J), ATM(ID)%dWSXzVdT(1,1,J),
	1             ATM(ID)%XzV, ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT, 
	1             ATM(ID)%NXzV,
	1             ATM(ID+1)%XzV, ATM(ID+1)%LOG_XzVLTE,
	1             ATM(ID+1)%dlnXzVLTE_dlnT, ATM(ID+1)%NXzV,
	1             ATM(ID)%XzV_ION_LEV_ID(J),ED,T,
	1             JREC,dJRECdt,JPHOT,FIXED_T,NUM_BNDS,ND,DST,DEND)
	          END DO
	        END IF
	      END DO
	      CALL TUNE(2,'BA_CONT_UP')
	    END IF
!
! 
! Note that ATM(ID+2)%EQXzV is the ion equation. Since Auger ionization,
! 2 electrons are ejected.
!
	    IF(XRAYS .AND. FINAL_CONSTANT_CROSS)THEN
	      CALL TUNE(1,'BA_XRAY_UP')
	      DO ID=1,NUM_IONS-1
	        IF(ATM(ID)%XzV_PRES .AND. ATM(ID+1)%XzV_PRES)THEN
	          CALL VSEBYJ_X_V7(ID,ATM(ID)%WSE_X_XzV,
	1             ATM(ID)%XzV, ATM(ID)%XzVLTE, ATM(ID)%dlnXzVLTE_dlnT, ATM(ID)%NXzV,
	1             ATM(ID+1)%XzV_F, ATM(ID+1)%XzVLTE_F, ATM(ID+1)%EDGEXzV_F,
	1             ATM(ID+1)%NXzV_F,ATM(ID+1)%DXzV,    ATM(ID+2)%EQXzV,
	1             ED,T,JREC,dJRECdT,JPHOT,FIXED_T,
	1             ND,NION,DST,DEND)
	        END IF
	      END DO
	      CALL TUNE(2,'BA_XRAY_UP')
	    END IF		!End X-rays
!
! 
!
! Increment the large simultaneous perturbation matrix due to the
! variation in J. We generally update the _PAR matrices. This is done to
! ensure better cancelation of terms similar in size as occurs when
! the optical depth is large.
!
! Every N_PAR frequencies we update the the full matrices, and zero the
! ?_PAR arrays. Note that the non-diagonal part of BA and BAION
! (i.e BA(I,J,K,L) with L .NE. DIAG_INDX) are presently updated for
! every frequency.
!
! We now pass the continuum emissivity and opacity without any scattering 
! contribution. We use TA as a zeroed vec for ESEC, which was subtracted
! from CHI_CONT inside BA_UPDATE_V7.
!
	    INDX_BA_METH=ND+1
	    IF(NEW_LINE_BA)INDX_BA_METH=MIN(MAX(INDX_BA_METH_RD,1),ND)
	    IF(NEW_LINE_BA .AND. .NOT. LAMBDA_ITERATION)THEN
	      TA(1:INDX_BA_METH-1)=ETA_NOSCAT(1:INDX_BA_METH-1); TA(INDX_BA_METH:ND)=ETA(INDX_BA_METH:ND)
	      TB(1:INDX_BA_METH-1)=CHI_NOSCAT(1:INDX_BA_METH-1); TB(INDX_BA_METH:ND)=CHI(INDX_BA_METH:ND)
	      TC(1:INDX_BA_METH-1)=0.0D0;                        TC(INDX_BA_METH:ND)=CHI_SCAT(INDX_BA_METH:ND)
              CALL BA_UPDATE_V7(VJ,VCHI_ALL,VETA_ALL,
	1             TA,TB,TC,T,POPS,RJ,FL,FQW,
	1             COMPUTE_NEW_CROSS,FINAL_CONSTANT_CROSS,DO_SRCE_VAR_ONLY,
	1             BA_CHK_FAC,NION,NT,NUM_BNDS,ND,DST,DEND)
	    ELSE IF(.NOT. LAMBDA_ITERATION)THEN
	      TA(1:ND)=0.0D0
              CALL BA_UPDATE_V7(VJ,VCHI_ALL,VETA_ALL,
	1             ETA_NOSCAT,CHI_NOSCAT,TA,T,POPS,RJ,FL,FQW,
	1             COMPUTE_NEW_CROSS,FINAL_CONSTANT_CROSS,DO_SRCE_VAR_ONLY,
	1             BA_CHK_FAC,NION,NT,NUM_BNDS,ND,DST,DEND)
	    END IF
!
	    IF(LST_ITERATION .AND. VERBOSE_OUTPUT .AND. .NOT. LAMBDA_ITERATION)THEN
	        T1=VCHI_ALL(VAR_INDX,DPTH_INDX)*RJ(DPTH_INDX)
                T2=(CHI(DPTH_INDX)-CHI_SCAT(DPTH_INDX))*VJ(VAR_INDX,DIAG_INDX,DPTH_INDX)
                T3=T1-VETA_ALL(VAR_INDX,DPTH_INDX)+T2
!                SUM_BA=SUM_BA+FQW*(T1+CHI_NOSCAT(DPTH_INDX)*VJ(VAR_INDX,DIAG_INDX,DPTH_INDX)-VETA_ALL(VAR_INDX,DPTH_INDX))
                SUM_BA=SUM_BA+FQW*T3
	        WRITE(200,'(I10,12ES16.8)')FREQ_INDX,FL,RJ(DPTH_INDX),STEQ_T(DPTH_INDX),T1,VETA_ALL(VAR_INDX,DPTH_INDX),T2,T3,FQW*T3,SUM_BA
	        WRITE(223,'(I10,12ES16.8)')FREQ_INDX,FL,RJ(DPTH_INDX),ETA(DPTH_INDX)/(CHI(DPTH_INDX)-CHI_SCAT(DPTH_INDX)),
	1                                     ETA_CONT(DPTH_INDX),CHI_CONT(DPTH_INDX),
	1                                     VETA_ALL(VAR_INDX,DPTH_INDX)*T(DPTH_INDX)/ETA_CONT(DPTH_INDX),
	1                                     VCHI_ALL(VAR_INDX,DPTH_INDX)*T(DPTH_INDX)/CHI_CONT(DPTH_INDX)
	    END IF
!
	  END DO		!K=1,ND : DST,DEND
	  CALL TUNE(ITWO,'COMPUTE_BA')
!
	  CALL TUNE(IONE,'ADD_PAR')
	  IF( MOD(FREQ_INDX,N_PAR) .EQ. 0 .OR. FREQ_INDX .EQ. NCF )THEN
            CALL ADD_PAR_TO_FULL_V2(NION,DIAG_INDX)
  	  END IF
	  CALL TUNE(ITWO,'ADD_PAR')
!
	END IF			!End compute_ba
!
! 
!
! We now allow for the variation of the LINE terms in the S.E. equations.
! In this first section we perform a LAMBDA iteration for the LINE terms.
! We do this if we are doing a FULL Lambda iteration, or if the line is
! being treated as a WEAK_LINE.
!
	CALL TUNE(1,'BA_LAM')
	DO SIM_INDX=1,MAX_SIM
	  IF( COMPUTE_BA .AND. (LAMBDA_ITERATION .OR. WEAK_LINE(SIM_INDX)) )THEN
	    IF (END_RES_ZONE(SIM_INDX) )THEN
!
! If we are doing a lambda iteration, we are assuming that JBAR is fixed.
! Thus, dZ/dCHIL and dZ/dETAL is given by the following.
!
! NB: We do not set VB(I) to if NEG_OPACITY(I)=.TRUE. as we have not altered
!     CHIL : We have only changed CHI which effects JBAR only.
!
	      OPAC_FAC=LINE_OPAC_CON(SIM_INDX)
	      STIM_FAC=GLDGU(SIM_INDX)*OPAC_FAC
	      EMIS_FAC=LINE_EMIS_CON(SIM_INDX)
	      NUP=SIM_NUP(SIM_INDX)
 	      NL=SIM_NL(SIM_INDX)
!
	      I=SIM_LINE_POINTER(SIM_INDX)
	      ID=VEC_ID(I)
	      MNL_F=VEC_MNL_F(I)
	      MNUP_F=VEC_MNUP_F(I)
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
! Altered 15-Feb-2005: No longer need VC. We now xplicitly cancel NU/ETAL depdence with Nu.
! Altered 9-Dec-2009: Changed sign of VB.
!
	      DO I=1,ND
	        VC(I)=JBAR_SIM(I,SIM_INDX)*(CHIL_MAT(I,SIM_INDX)/ETAL_MAT(I,SIM_INDX))
	        VB(I)=JBAR_SIM(I,SIM_INDX)*(ATM(ID)%XzV_F(MNUP_F,I)/ETAL_MAT(I,SIM_INDX))
	      END DO
!
! Because we write the heating/cooling term interms of EINA, we need
! to do the scaling when either SCL_LINE_COOL_RATES or SCL_SL_LINE_OPAC is true.
!
	      SCL_FAC=1.0D0
	      IF(SCL_LINE_COOL_RATES)THEN
	        SCL_FAC=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
	        IF(ABS(SCL_FAC-1.0D0) .GT. SCL_LINE_HT_FAC)SCL_FAC=1.0D0
	      END IF
	      DO K=1,ND
	        L=GET_DIAG(K)
	        dRATE_dUP=EINA(SIM_INDX)*U_STAR_RATIO(K,SIM_INDX)*
	1                   (ZNET_SIM(K,SIM_INDX)+VC(K)+STIM_FAC*VB(K))
	        dRATE_dLOW=-EINA(SIM_INDX)*OPAC_FAC*L_STAR_RATIO(K,SIM_INDX)*VB(K)
	        SE(ID)%BA(MNUP,MNUP,L,K)=SE(ID)%BA(MNUP,MNUP,L,K)-dRATE_dUP
	        SE(ID)%BA(MNUP,MNL,L,K) =SE(ID)%BA(MNUP,MNL,L,K) -dRATE_dLOW
	        SE(ID)%BA(MNL,MNUP,L,K) =SE(ID)%BA(MNL,MNUP,L,K) +dRATE_dUP
	        SE(ID)%BA(MNL,MNL,L,K)  =SE(ID)%BA(MNL,MNL,L,K)  +dRATE_dLOW
!
	        T4=SCL_FAC
	        IF(POP_ATOM(K) .GE. SCL_LINE_DENSITY_LIMIT)T4=1.0D0
	        dRATE_dUP=T4*U_STAR_RATIO(K,SIM_INDX)*(LINE_EMIS_CON(SIM_INDX)*LINE_QW_SUM(K,SIM_INDX)+STIM_FAC*JBAR_SIM(K,SIM_INDX))
	        dRATE_dLOW=T4*L_STAR_RATIO(K,SIM_INDX)*LINE_OPAC_CON(SIM_INDX)*JBAR_SIM(K,SIM_INDX)
	        BA_T(NL,L,K) =BA_T(NL,L,K) + dRATE_dLOW
	        BA_T(NUP,L,K)=BA_T(NUP,L,K) - dRATE_dUP
	      END DO
	    END IF	!Resonance zone
	  END IF	!Lambda/Weak line
	END DO		!Loop over SIM_INDX
	CALL TUNE(2,'BA_LAM')
!
! 
	CALL TUNE(IONE,'LINE_BA')
	IF(COMPUTE_BA .AND. .NOT. LAMBDA_ITERATION)THEN
	  IF(NUM_BNDS .EQ. ND)THEN
	     WRITE(LUER,*)'Error --- NUM_BNDS .EQ. ND not implemented'
	     STOP
	  END IF
!
! Update the dZ matrix which gives the variation of the Net Rate
! (=1-JBAR/SL) with respect to eta, chi, ED, ... etc.
!
! dZ( X, X depth [1:NUM_BNDS], Z depth, N_SIM)
!
! Scaling MAX_SIM*ND*NUM_BNDS*NM
!
! NB: J goes from 1 to NUM_BNDS (and not BNDST(K),BNDEND(K)) since we have no
! reference to elements in a full matrix (such as FCHI).
!
! NB: OPAC_FAC and STIM_FAC are not multiplied by NEG_OPAC_FAC since we use
!     the uncorrected line opacity in the source function and hence the
!     expression for ZNET.
!
	  IF(UPDATE_dZ)THEN
	    DO L=1,ND
	      TA(L)=CHI_NOSCAT_PREV(L)/CHI_NOSCAT(L)
	      IF(ETA_CONT(L) .EQ. 0)THEN
	        TB(L)=1.0D0
	      ELSE
	        TB(L)=ETA_PREV(L)/ETA_CONT(L)
	      END IF
	      IF(TA(L) .GT. 5.0)THEN
	        TA(L)=0.0D0; TB(L)=0.0D0
	      END IF
	    END DO
	  END IF
!
	  CALL TUNE(IONE,'dZ_LINE')
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(T1,MUL_FAC,OPAC_FAC,EMIS_FAC,STIM_FAC,J,K,L,NL,NUP)
	  DO SIM_INDX=1,MAX_SIM
	    IF(RESONANCE_ZONE(SIM_INDX) .AND. .NOT. WEAK_LINE(SIM_INDX))THEN
!
! Adjust variation of line with continuum to current continuum frequency.
!
	      IF(UPDATE_dZ)THEN
	        DO K=1,ND
	          DO J=BNDST(K),BNDEND(K)
 	            L=BND_TO_FULL(J,K)
	            T1=ETA_CONT(L)*HDKT*(FL_OLD-FL)/T(L)/T(L)
	            dZ(3,J,K,SIM_INDX)=dZ(3,J,K,SIM_INDX)*TA(L)
	            dZ(4,J,K,SIM_INDX)=dZ(4,J,K,SIM_INDX)*TB(L)
	            dZ(6,J,K,SIM_INDX)=dZ(6,J,K,SIM_INDX)+dZ(4,J,K,SIM_INDX)*T1
	          END DO
	        END DO
	      END IF
!
	      IF(INCLUDE_dSLdT)THEN
	        NL=SIM_NL(SIM_INDX)
	        NUP=SIM_NUP(SIM_INDX)
	        J=DIAG_INDX
	        DO K=1,ND
	          MUL_FAC=LINE_QW_SIM(K,SIM_INDX)*CHIL_MAT(K,SIM_INDX)/ETAL_MAT(K,SIM_INDX)
	          OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*dL_RAT_dT(K,SIM_INDX)*POPS(NL,K)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*dU_RAT_dT(K,SIM_INDX)*POPS(NUP,K)
	          STIM_FAC=LINE_OPAC_CON(SIM_INDX)*dU_RAT_dT(K,SIM_INDX)*GLDGU(SIM_INDX)*POPS(NUP,K)
	          dZ(6,J,K,SIM_INDX)=dZ(6,J,K,SIM_INDX)- 
	1               MUL_FAC*( dJ_LOC(6,J,K) + RJ(K)*(OPAC_FAC-STIM_FAC)/CHIL_MAT(K,SIM_INDX) 
	1              - RJ(K)*EMIS_FAC/ETAL_MAT(K,SIM_INDX) )
	        END DO
	      END IF
!
	      DO K=1,ND
	        MUL_FAC=LINE_QW_SIM(K,SIM_INDX)*CHIL_MAT(K,SIM_INDX)/ETAL_MAT(K,SIM_INDX)
	        OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*L_STAR_RATIO(K,SIM_INDX)
	        EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*U_STAR_RATIO(K,SIM_INDX)
	        STIM_FAC=LINE_OPAC_CON(SIM_INDX)*U_STAR_RATIO(K,SIM_INDX)*GLDGU(SIM_INDX)
	        DO J=BNDST(K),BNDEND(K)
	          IF(J .NE. DIAG_INDX)THEN
	            DO L=3,NM
	              IF(DO_THIS_TX_MATRIX(L))THEN
	                dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) - MUL_FAC*dJ_LOC(L,J,K)
	              END IF
	            END DO
	          ELSE 
	            DO L=3,NM
	              IF(DO_THIS_TX_MATRIX(L))THEN
	                IF( L .EQ. LOW_POINTER(SIM_INDX) )THEN
	                  dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) -
	1                    MUL_FAC*( dJ_LOC(L,J,K) + RJ(K)*OPAC_FAC/CHIL_MAT(K,SIM_INDX) )
	                ELSE IF( L .EQ. UP_POINTER(SIM_INDX) )THEN
	                  dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) -
	1                    MUL_FAC*( dJ_LOC(L,J,K)-
	1                    RJ(K)*STIM_FAC/CHIL_MAT(K,SIM_INDX)-
	1                    RJ(K)*EMIS_FAC/ETAL_MAT(K,SIM_INDX) )
	                ELSE IF(L .EQ. 6 .AND. INCLUDE_dSLdT)THEN
	                ELSE
	                  dZ(L,J,K,SIM_INDX)=dZ(L,J,K,SIM_INDX) - MUL_FAC*dJ_LOC(L,J,K)
	                END IF
	              END IF
	            END DO
	          END IF
	        END DO
	      END DO
	    END IF 
	  END DO
!$OMP END PARALLEL DO
	  CALL TUNE(ITWO,'dZ_LINE')
!
! 
!
! Check whether any lines have been finalized. If so we can update the 
! variation matrices.
!
! NB: Weak lines are handled in the LAMBDA ITERATION section.
!
! We use ZENT_VAR_LIMIT to reduce the number of operations. When ZNET
! is approximatley unity, there is little point in linearizing ZNET.
! NB: An ideal vale for ZNET_VAR_LIMIT is probably 0.01 or 0.001. If 
! ZNET_VAR_LIMIT is zero, all depths will be included in the linearization, 
! independent of ZNET. A very large value of ZNET (i.e. 10^4), will imply
! an interation on the NET_RATES, with no linearization.
!
	  DO SIM_INDX=1,MAX_SIM
!
	    IF( END_RES_ZONE(SIM_INDX) .AND. .NOT. WEAK_LINE(SIM_INDX))THEN
!
	      I=SIM_LINE_POINTER(SIM_INDX)
	      ID=VEC_ID(I)
!
! Can now update rate equation for variation in Z.
!
	      I=NT*NUM_BNDS*ND
	      CALL DP_ZERO(dZ_POPS,I)
!
! Check which matrices get used to compute dZ_POPS. We include:
!   (i) all transitions designated to be important.
!  (ii) all transitions for levels in the same ionization stage.
!
	      DO I=TX_OFFSET+1,NM
	        USE_THIS_VAR_MAT(I)=.FALSE.
	        IF(VAR_IN_USE_CNT(I) .GT. 0 .AND. IMP_TRANS_VEC(I))THEN
	          USE_THIS_VAR_MAT(I)=.TRUE.
	        ELSE IF(VAR_IN_USE_CNT(I) .GT. 0)THEN
	          NL=VAR_LEV_ID(I)
	          IF(NL .GE. ATM(ID)%EQXzV .AND.
	1                      NL .LE. ATM(ID)%EQXzV+ATM(ID)%NXzV)
	1         USE_THIS_VAR_MAT(I)=.TRUE.
	        END IF
	      END DO
!
! NOPS= ND*NUM_BNDS*( 4NT + 2 + 7NUM_SIM )
!
!
! Set up a temporary vector to handle Rayleigh scattering.
!
	      TA(1:ND)=0.0D0
	      IF(SPECIES_PRES(1))THEN			!If H present!
	        DO L=1,ND
	          TA(L)=CHI_RAY(L)/ATM(1)%XzV_F(1,L)	          
	        END DO
	      END IF
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(T4,I,J,JJ,L,NL)
	      DO K=1,ND
		T4=ABS(ZNET_SIM(K,SIM_INDX)-1.0D0)
	          IF(T4 .GT. ZNET_VAR_LIMIT)THEN
	          DO J=BNDST(K),BNDEND(K)	  	  !Since refer to VCHI etc.
 	            L=BND_TO_FULL(J,K)
	            DO JJ=1,SE(ID)%N_IV			!Bad notation
	              I=SE(ID)%LNK_TO_F(JJ)
	              dZ_POPS(I,J,K)=dZ_POPS(I,J,K) +
	1                ( VCHI(I,L)*dZ(3,J,K,SIM_INDX) + 
	1                    VETA(I,L)*dZ(4,J,K,SIM_INDX) )
	            END DO
	            dZ_POPS(NT-1,J,K)=dZ_POPS(NT-1,J,K) +
	1                            ESEC(L)*dZ(5,J,K,SIM_INDX)/ED(L)
	            dZ_POPS(1,J,K)=dZ_POPS(1,J,K) +
	1                            TA(L)*dZ(5,J,K,SIM_INDX)
	            dZ_POPS(NT,J,K)=dZ_POPS(NT,J,K) + dZ(6,J,K,SIM_INDX)
!
! Now must do line terms.
!
	            DO I=TX_OFFSET+1,NM
	              IF(USE_THIS_VAR_MAT(I))THEN
	                NL=VAR_LEV_ID(I)
	                dZ_POPS(NL,J,K)=dZ_POPS(NL,J,K) + dZ(I,J,K,SIM_INDX)
	              END IF
	            END DO
	          END DO	        !Over Variable depth (1:NUM_BNDS)
	        END IF		!ABS|ZNET-1| > 0.01
	      END DO			!Over J depth.
!$OMP END PARALLEL DO
!
! Update variation equations. NB. We do not need to update the Radiative
! Equilibrium equation, since this is automatically updated with the
! continuum.
!
! Note that there is no update of the NT equation since this is 
! automatically included in the direct CHI*J-ETA term.
!
! NOPS = 6*ND*NT*NUM_BNDS
!
	      CALL TUNE(IONE,'dBA_LINE')
	      NL=SIM_NL(SIM_INDX)
	      NUP=SIM_NUP(SIM_INDX)
	      I=SIM_LINE_POINTER(SIM_INDX)
              MNL_F=VEC_MNL_F(I)
	      MNUP_F=VEC_MNUP_F(I)
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
! NB: We compute T2 as U_STAR_RATIO may have been scaled, and this should
! not effect the rate equations.
! 
	      OPAC_FAC=LINE_OPAC_CON(SIM_INDX)
	      STIM_FAC=LINE_OPAC_CON(SIM_INDX)*GLDGU(SIM_INDX)
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(T1,T2,T4,J,K,JJ)
	      DO L=1,ND					!Equation depth
!
	        K=DIAG_INDX
	        T2=EINA(SIM_INDX)*U_STAR_RATIO(L,SIM_INDX)*ZNET_SIM(L,SIM_INDX)
	        SE(ID)%BA(MNL,MNUP,K,L)=SE(ID)%BA(MNL,MNUP,K,L) + T2
	        SE(ID)%BA(MNUP,MNUP,K,L)=SE(ID)%BA(MNUP,MNUP,K,L) - T2
!
		J=SE(ID)%LNK_TO_IV(NT)
	        T2=EINA(SIM_INDX)*dU_RAT_dT(L,SIM_INDX)*POPS(NUP,L)*ZNET_SIM(L,SIM_INDX)
	        SE(ID)%BA(MNL,J,K,L)=SE(ID)%BA(MNL,J,K,L) + T2
	        SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L) - T2
!
	        T1=EINA(SIM_INDX)*ATM(ID)%XzV_F(MNUP_F,L)
		T4=ABS(ZNET_SIM(L,SIM_INDX)-1.0D0)
	        IF(T4 .GT. ZNET_VAR_LIMIT)THEN
	          DO K=1,NUM_BNDS			!Variable depth
	            DO J=1,SE(ID)%N_IV
	              JJ=SE(ID)%LNK_TO_F(J)
	              SE(ID)%BA(MNL,J,K,L) =SE(ID)%BA(MNL,J,K,L) +dZ_POPS(JJ,K,L)*T1
	              SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L)-dZ_POPS(JJ,K,L)*T1
	            END DO
	          END DO
	        END IF
	      END DO
!$OMP END PARALLEL DO
!
	      IF(SCL_LINE_COOL_RATES)THEN
	        SCL_FAC=(AVE_ENERGY(NL)-AVE_ENERGY(NUP))/FL_SIM(SIM_INDX)
	        IF(ABS(SCL_FAC-1.0D0) .GT. SCL_LINE_HT_FAC)SCL_FAC=1.0D0
	      ELSE
	        SCL_FAC=1.0D0
	      END IF
!
	      IF(.NOT. FIXED_T)THEN
	        DO L=INDX_BA_METH,ND
	          K=GET_DIAG(L)
	          T3=SCL_FAC
	          IF(POP_ATOM(L) .GE. SCL_LINE_DENSITY_LIMIT)T3=1.0D0
	          dRATE_dUP=T3*U_STAR_RATIO(L,SIM_INDX)*(LINE_EMIS_CON(SIM_INDX)*LINE_QW_SUM(L,SIM_INDX)+STIM_FAC*JBAR_SIM(L,SIM_INDX))
	          dRATE_dLOW=T3*L_STAR_RATIO(L,SIM_INDX)*LINE_OPAC_CON(SIM_INDX)*JBAR_SIM(L,SIM_INDX)
	          BA_T(NUP,K,L)=BA_T(NUP,K,L) - dRATE_dUP
	          BA_T(NL,K,L)=BA_T(NL,K,L) + dRATE_dLOW
	          T1=LINE_OPAC_CON(SIM_INDX)*(dL_RAT_dT(L,SIM_INDX)*POPS(NL,L)-GLDGU(SIM_INDX)*dU_RAT_dT(L,SIM_INDX)*POPS(NUP,L))
	          T2=LINE_EMIS_CON(SIM_INDX)*dU_RAT_dT(L,SIM_INDX)*POPS(NUP,L)
	          BA_T(NT,K,L)=BA_T(NT,K,L) + T3*(T1*JBAR_SIM(L,SIM_INDX) - T2*LINE_QW_SUM(L,SIM_INDX)) 
	        END DO
	      END IF
!
	      IF(.NOT. FIXED_T)THEN
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(T1,T3,T4,J,K,JJ)
	        DO L=1,INDX_BA_METH-1
	          K=GET_DIAG(L)
	          T3=SCL_FAC
	          IF(POP_ATOM(L) .GE. SCL_LINE_DENSITY_LIMIT)T3=1.0D0
	          BA_T(NUP,K,L)=BA_T(NUP,K,L)-T3*ETAL_MAT(L,SIM_INDX)*ZNET_SIM(L,SIM_INDX)/POPS(NUP,L)
	          T1=T3*LINE_EMIS_CON(SIM_INDX)*dU_RAT_dT(L,SIM_INDX)*POPS(NUP,L)*ZNET_SIM(L,SIM_INDX)
	          BA_T(NT,K,L)=BA_T(NT,K,L) - T1
		  T4=ABS(ZNET_SIM(L,SIM_INDX)-1.0D0)
	          IF(T4 .GT. ZNET_VAR_LIMIT)THEN
	            DO K=1,NUM_BNDS			!Variable depth
	              DO J=1,SE(ID)%N_IV
	                JJ=SE(ID)%LNK_TO_F(J)
                        BA_T(JJ,K,L)=BA_T(JJ,K,L)-T3*dZ_POPS(JJ,K,L)*ETAL_MAT(L,SIM_INDX)
	              END DO
	            END DO
	          END IF
	        END DO
!$OMP END PARALLEL DO
	      END IF
!
	      CALL TUNE(ITWO,'dBA_LINE')
!
! Must now zero dZ since next time it is used it will be for a new line.
!
	      I=NM*NUM_BNDS*ND
	      CALL DP_ZERO(dZ(1,1,1,SIM_INDX),I)
!
	    END IF	!End_res_zone .and. .not. weak_line
	  END DO	!SIM_INDX
	END IF		!Compute_ba .and. .not. lambda_iteration
	CALL TUNE(ITWO,'LINE_BA')
	FL_OLD=FL
!
	RETURN
	END
