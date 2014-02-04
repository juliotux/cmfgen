!
! Subroutine to compute JBAR and update the variation matrix (BA) for lines treated using 
! the solution Sobolev approximation. Routine is based on a section of code in CMFGEN_SUB.F
!
! Routine was tested and appears to work. This option is becoming obsolete in CMFGEN.
!
	SUBROUTINE SUB_SOB_LINE(SECTION,POPS,CHIL,ETAL,NEG_OPACITY,
	1                    FL,CONT_FREQ,AMASS,
	1                    EDDINGTON,IMPURITY_CODE,VAR_SOB_JC,LST_ITERATION,
	1                    EW,CONT_INT,
	1                    NL,NUP,NT,ND,NC,NP,
	1                    NDEXT,NCEXT,NPEXT,
	1                    NLF,DIAG_INDX,NUM_BNDS)
!
	USE ANG_QW_MOD
	USE CONTROL_VARIABLE_MOD
	USE CMF_SOB_MOD
	USE LINE_MOD
	USE LINE_VEC_MOD
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE RADIATION_MOD
	USE STEQ_DATA_MOD
 	USE VAR_RAD_MOD
	IMPLICIT NONE
!
! Cleaned 17-Dec-2004
!
	REAL*8 FL
	REAL*8 CONT_FREQ
	REAL*8 AMASS
	REAL*8 EW
	REAL*8 CONT_INT
!
	INTEGER NL
	INTEGER NUP
	INTEGER NT
	INTEGER ND,NC,NP
	INTEGER NDEXT,NCEXT,NPEXT
	INTEGER NLF
	INTEGER DIAG_INDX
	INTEGER NUM_BNDS
	LOGICAL EDDINGTON
	LOGICAL LST_ITERATION
	LOGICAL VAR_SOB_JC
	LOGICAL IMPURITY_CODE
	CHARACTER(LEN=*) SECTION
!
	REAL*8 POPS(NT,ND)
	REAL*8 CHIL(ND)
	REAL*8 ETAL(ND)
	LOGICAL NEG_OPACITY(ND)
!
! Constants for opacity etc.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 OPLIN,EMLIN
!
! Local variables.
!
	REAL*8 VB(ND)
	REAL*8 VC(ND)
	REAL*8 T1,T2,T3,T4
	REAL*8 SRAT
	INTEGER I,J,K,L
	INTEGER L1,L2,U1,U2,IT
	INTEGER MNL,MNL_F
	INTEGER MNUP,MNUP_F
	INTEGER MNT
	INTEGER LS
	INTEGER JJ
	INTEGER ID
	INTEGER PHOT_ID
	INTEGER LUER,ERROR_LU
	INTEGER TX_OFFSET
	INTEGER FREQ_INDX
	INTEGER NM
	INTEGER NM_KI
	INTEGER MAX_SIM
	EXTERNAL ERROR_LU
!
	INTEGER GET_DIAG
	INTEGER BNDST
	INTEGER BNDEND
	INTEGER BND_TO_FULL
	LOGICAL LST_DEPTH_ONLY 
	LOGICAL COMPUTE_LAM
	LOGICAL FULL_ES
	LOGICAL FIRST_FREQ
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
	BNDEND(K)=MIN( (NUM_BNDS/ND)*(K-DIAG_INDX)+ND+DIAG_INDX-K,NUM_BNDS )
!
!
	COMPUTE_LAM=.FALSE.
	FULL_ES=.TRUE.
	LUER=ERROR_LU()
	FIRST_FREQ=.FALSE.		!Only used if SECTION='CONTINUUM'
	FREQ_INDX=0                     !Not used although passed
	TX_OFFSET=0                     !Not used
	NM=-1
	NM_KI=-1
	MAX_SIM=-1
C
C Use the escape probability approximation for lines originating
C in all levels.
C
	  CALL TUNE(IONE,'SOBOLEV')
	    CALL SOBJBAR_SIM(SOURCE,CHI,CHIL,ETAL,V,SIGMA,R,P,AQW,
	1                  JBAR,BETA,
	1                  ZNET_SIM,VB_SIM,VC_SIM,VB_2,VC_2,BETAC_SIM,
	1                  CHIL_MAT,ETAL_MAT,BB_COR,
	1                  FL,DIF,DBB,IC,THK_CONT,
	1                  NLF,NC,NP,ND,NUM_SIM_LINES,METHOD)
	  CALL TUNE(ITWO,'SOBOLEV')
C
C Estimate the line EW using a Modified SObolev approximation.
C
	  IF(LST_ITERATION)THEN
C
C We use TA as a temporary vector which indicates the origin
C of the line emission. Not required in this code as used only
C for display purposes.
C
	    CALL SOBEW(SOURCE,CHI,CHI_SCAT,CHIL,ETAL,
	1              V,SIGMA,R,P,AQW,HQW,TA,EW,CONT_INT,
	1              FL,DIF,DBB,IC,THK_CONT,L_FALSE,NC,NP,ND,METHOD)
C
	  END IF
C
C If the line opacity was negative, we set the variation of JBAR
C with CHIL to zero.
C
	  DO SIM_INDX=1,NUM_SIM_LINES
	    DO I=1,ND
	      IF(NEG_OPACITY(I))THEN
	        VB_SIM(I,SIM_INDX)=0.0D0
	        VB_2(I)=0.0D0
	      END IF
	    END DO
	  END DO
C
C Allow for the variation of the continuous radiation field.
C
	  IF(COMPUTE_BA .AND. VAR_SOB_JC .AND. .NOT. 
	1                LAMBDA_ITERATION .AND. .NOT. IMPURITY_CODE)THEN
C
C	    INCLUDE 'VARCONT.INC'
            CALL DO_VAR_CONT(POPS,SECTION,EDDINGTON,
	1                  FL,CONT_FREQ,FREQ_INDX,FIRST_FREQ,TX_OFFSET,
	1                  ND,NC,NP,NUM_BNDS,DIAG_INDX,NT,NM,
	1                  NDEXT,NCEXT,NPEXT,MAX_SIM,NM_KI)
C
	    DO SIM_INDX=1,NUM_SIM_LINES
	      DO I=1,ND
	        BETAC_SIM(I,SIM_INDX)=BETAC_SIM(I,SIM_INDX)/RJ(I)
	      END DO
	    END DO
C
C Increment the large simultaneous perturbation matrix. Note that terms are
C
C EINA*POPS(NUP,L)*BETAC*VJ(J,K,L) 
C
C        and
C
C ETAL(L)*BETAC(L)*VJ(J,K,L) but 
C
C ETAL=FL*EINA(1)*EMLIN*POPS. Therefore use FL*EMLIN as constant.
C
	    CALL TUNE(IONE,'SOBCONTBA')
	    DO SIM_INDX=1,NUM_SIM_LINES
	      NL=SIM_NL(SIM_INDX)
	      NUP=SIM_NUP(SIM_INDX)
	      I=SIM_LINE_POINTER(SIM_INDX)
	      ID=VEC_ID(I)
	      MNL_F=VEC_MNL_F(I)
	      MNUP_F=VEC_MNUP_F(I)
	      MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	      MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
!
	      T1=FL_SIM(SIM_INDX)*EMLIN
	      DO L=1,ND	  		  	!S.E. equation depth
	        T2=ETAL_MAT(L,SIM_INDX)*BETAC_SIM(L,SIM_INDX)
	        T3=1.0D-06*RJ(L)
	        DO K=BNDST(L),BNDEND(L)	 	!Variable depth.
	          LS=BND_TO_FULL(K,L)
!**********************************
!*********************************
   	          DO  J=1,SE(ID)%N_IV	    		!Variable
	            JJ=SE(ID)%LNK_TO_F(J)
!   	          DO  J=1,NT	    		!Variable
!	            IF( DABS(VJ(J,K,L)*POPS(J,LS)) .GE. T3 )THEN
	              T4=T2*VJ(JJ,K,L)
	              SE(ID)%BA(MNL,J,K,L) =SE(ID)%BA(MNL,J,K,L) - T4/T1
	              SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L) + T4/T1
	              BA_T(JJ,K,L)=BA_T(JJ,K,L) + T4
!	            END IF
	          END DO
	        END DO
	      END DO
	    END DO
	    CALL TUNE(ITWO,'SOBCONTBA')
	  END IF			!BA Matrix computed (compute_ba).
C
C If we are doing a lambda iteration, we are assuming that JBAR is fixed.
C Thus, dZ/dCHIL and dZ/dETAL is given by the following.
C
	IF(COMPUTE_BA .AND. LAMBDA_ITERATION)THEN
	  DO SIM_INDX=1,NUM_SIM_LINES
	    DO I=1,ND
	      VB_SIM(I,SIM_INDX)=-JBAR(I)/ETAL_MAT(I,SIM_INDX)
	      IF(NEG_OPACITY(I))VB_SIM(I,SIM_INDX)=0.0D0
	      VC_SIM(I,SIM_INDX)=JBAR(I)*CHIL_MAT(I,SIM_INDX)/
	1                 ETAL_MAT(I,SIM_INDX)/ETAL_MAT(I,SIM_INDX)
	    END DO
	  END DO
	  DO I=1,ND
	    VB_2(I)=0.0D0
	    VC_2(I)=0.0D0
	  END DO
	END IF
C
C Evaluate contribution to statistical equilibrium equation, and
C and increment variation matrices.
C
	DO SIM_INDX=1,NUM_SIM_LINES
	  T1=FL_SIM(SIM_INDX)*EMLIN
	  NUP=SIM_NUP(SIM_INDX)
	  NL=SIM_NL(SIM_INDX)
	  I=SIM_LINE_POINTER(SIM_INDX)
	  ID=VEC_ID(I)
	  MNL_F=VEC_MNL_F(I);     MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	  MNUP_F=VEC_MNUP_F(I);   MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	  DO K=1,ND
	    T2=ETAL_MAT(K,SIM_INDX)*ZNET_SIM(K,SIM_INDX)
	    SE(ID)%STEQ(MNUP,K)=SE(ID)%STEQ(MNUP,K) - T2/T1
	    SE(ID)%STEQ(MNL,K) =SE(ID)%STEQ(MNL,K) + T2/T1
	    STEQ_T(K)=STEQ_T(K) - T2
	  END DO
	END DO
C
C NB: We now include the temperature variation of L_STAR_RATIO and U_STAR_RATIO
C     in the variation of the net rate. 
C
C Definitions:
C            RATE_FAC : Factor which multiplies the net rate to form the total
C                       rate from the upper level.
C            SRAT:      Ratio of total source function to individual source
C                       function corrected for any difference in frequency.
C            SIM_INDX   is used to refer to the transition of interest.
C            K          Depth index
C            J          Line opacity/emissivity identifier
C                                    
	IF(COMPUTE_BA)THEN
	  DO SIM_INDX=1,NUM_SIM_LINES
	    NL=SIM_NL(SIM_INDX)
	    NUP=SIM_NUP(SIM_INDX)
	    I=SIM_LINE_POINTER(SIM_INDX)
	    ID=VEC_ID(I)
	    MNL_F=VEC_MNL_F(I)
	    MNUP_F=VEC_MNUP_F(I)
	    MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	    MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	    IT=SE(ID)%N_IV
!
	    T4=FL_SIM(SIM_INDX)*EMLIN
	    DO K=1,ND
	      L=GET_DIAG(K)
	      SRAT=(ETAL(K)/ETAL_MAT(K,SIM_INDX)/BB_COR(K,SIM_INDX))*
	1          (CHIL_MAT(K,SIM_INDX)/CHIL(K))
	      RATE_FAC=EINA(SIM_INDX)*POPS(NUP,K)*U_STAR_RATIO(K,SIM_INDX)
	      DO J=1,NUM_SIM_LINES
	        OPAC_FAC=OSCIL(J)*OPLIN
	        STIM_FAC=-GLDGU(J)*OPAC_FAC
	        EMIS_FAC=EINA(J)*FL_SIM(J)*EMLIN
	        L1=SIM_NL(J);    L2=SE(ID)%LNK_TO_IV(L1)
	        U1=SIM_NUP(J);   U2=SE(ID)%LNK_TO_IV(U1)
!
	        IF(L2 .NE. 0 .AND. U2 .NE. 0)THEN
	          IF(J .EQ. SIM_INDX)THEN
	            dRATE_dUP=RATE_FAC*( ZNET_SIM(K,SIM_INDX)/POPS(NUP,K) +
	1                 U_STAR_RATIO(K,J)*
	1                (STIM_FAC*VB_SIM(K,J)+EMIS_FAC*VC_SIM(K,J)) )
	            dRATE_dLOW=RATE_FAC*OPAC_FAC*VB_SIM(K,J)*
	1                 L_STAR_RATIO(K,J)
	            dRATE_dT=RATE_FAC*
	1              (  OPAC_FAC*POPS(NL,K)*VB_SIM(K,J)*dL_RAT_dT(K,J)+
	1                  POPS(NUP,K)*dU_RAT_dT(K,J)*(
	1                    EMIS_FAC*VC_SIM(K,J)+
	1                        STIM_FAC*VB_SIM(K,J) ) +
	1                 ZNET_SIM(K,SIM_INDX)*dU_RAT_dT(K,J)
	1                            /U_STAR_RATIO(K,J)
	1              )
	          ELSE
	            dRATE_dUP=RATE_FAC*SRAT*U_STAR_RATIO(K,J)*
	1                     ( STIM_FAC*VB_2(K)+EMIS_FAC*VC_2(K)*BB_COR(K,J) )
	            dRATE_dLOW=RATE_FAC*OPAC_FAC*VB_2(K)*SRAT*L_STAR_RATIO(K,J)
	            dRATE_dT=RATE_FAC*SRAT*
	1               ( OPAC_FAC*POPS(SIM_NL(J),K)*VB_2(K)*dL_RAT_dT(K,J)+
	1                 POPS(SIM_NUP(J),K)*dU_RAT_dT(K,J)*(
	1                 EMIS_FAC*VC_2(K)*BB_COR(K,J)+STIM_FAC*VB_2(K)) )
	          END IF
!
	          SE(ID)%BA(MNUP,U2,L,K) = SE(ID)%BA(MNUP,U2,L,K)  - dRATE_dUP
	          SE(ID)%BA(MNUP,L2,L,K) = SE(ID)%BA(MNUP,L2,L,K) - dRATE_dLOW
	          SE(ID)%BA(MNUP,IT,L,K) = SE(ID)%BA(MNUP,IT,L,K) - dRATE_dT
	          SE(ID)%BA(MNL,U2,L,K ) = SE(ID)%BA(MNL,U2,L,K)  + dRATE_dUP
	          SE(ID)%BA(MNL,L2,L,K)  = SE(ID)%BA(MNL,L2,L,K)  + dRATE_dLOW
	          SE(ID)%BA(MNL,IT,L,K)  = SE(ID)%BA(MNL,IT,L,K)  + dRATE_dT
	          BA_T(SIM_NL(J),L,K) =BA_T(SIM_NL(J),L,K) - T4*dRATE_dLOW
	          BA_T(SIM_NUP(J),L,K)=BA_T(SIM_NUP(J),L,K) - T4*dRATE_dUP
	          BA_T(NT,L,K)=BA_T(NT,L,K) - T4*dRATE_dT
	       END IF
	      END DO
	    END DO
	  END DO
	END IF
!
	RETURN
	END
