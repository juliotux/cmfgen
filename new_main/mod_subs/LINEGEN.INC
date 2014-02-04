C
C This section of the program computes the variation of Jmn as
C a function of the population levels and updates both the statistical
C equilibrium equations and the variation matrix.
C
C Altered 07-Dec-2004 - VAROPAC replaced by call to COMP_VAR_OPAC.
C                       Bug fixed: NL and NUP were being used to access
C                         SE(ID)%BA matrix instead of MNL & MNUP (on RHS).
C Altered 06-Jan-1997 - Super levels allowed for in linearization.
C
C Altered 12-Jun-91 - PERTJFEAUNEW call changed to PERTJFEAU_IBC. This
C                     routine has improved handling of the outer boundary
C                     condition when optically thick condition is required.
C
C Altered 13-Feb-88 - Changed to allow Banded-Newton Raphson Operators.
C                     The number of bands is given by NUM_BNDS and must
C                     be odd. If NUM_BNDS=ND, the matrix is assumed
C                     to represent the full Newton-Raphson Operator.
C                     The index, DIAG_INDX must also be specified, but is
C                     only used for the BANDED case. Note 
C                     DIAG_INDX=( NUM_BNDS/2+1 )
C Altered 25-Jan-87 - Update of STEQ is now located in JBAR.  Consequently
C                     this section can be advoided if we only wish to 
C                     increment STEQ and not change the large variation matrix.
C
C Altered 2-MAY-85 (Bug fixed in variation of radiative equilibrium equation)
C
C BA and VJ have the following dimensions.
C
C                     BA(NT,NT,NUM_BNDS,ND)
C                     VJ(NT,NUM_BANDS,ND)
C
C                     In addittion, the variables HYD_PRES, HEI_PRES
C                     and HE2_PRES have been installed. These have
C                     been included so the routine is completely general.
C                     These LOGICAL variables should be specified by
C                     paramter statements - they can then be optimized
C                     out of the code at compile time.
C
C
C Altere              - Carbon routines instaled etc.
C
C Altered 11-Apr-1989 - Pertubation of JCONT made identical to that in
C                       VARCONT. Routine now calls NEWPERTJD rather
C                       than PERTJD. In addition, dZdCHIL is now set to
C                       zero if the NEG_OPACITY flag is set.
C
C Altered 31-May-1989 - Line varaition section made into a subroutine
C                       (VAR_FORMSOL). VAR_MOMHAM included.
C
C Altered 19-Oct-1989 - IMPURITY_CODE variable installed. If this variable is
C                       TRUE, the continuous radiation field is held fixed.
C 
C
C Routine to compute dJdETA, dJdCHI etc.
C
	IF(EDDINGTON)THEN
C
C We use WM for WORKMAT and F2DA for RHS_dHdCHI. These (ND,ND) matrices
C are only used as work matrices in VAR_MOMHAM
C
	  CALL TUNE(1,'VARMOM')
	  CALL VAR_MOMJBAR(ETA,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,
	1                  FQAF,FQAFD,TX,TVX,KI,WM,F2DA,
	1                  F_LINE,G_LINE,HBC_LINE,IN_HBC_LINE,NBC_LINE,
	1                  PF,PROF,LFQW,FL,DIF,DBB,IC,METHOD,
	1                  NLF,ND,NNM)
	  CALL TUNE(2,'VARMOM')
    	ELSE
	  CALL TUNE(1,'VARFOR')
	  CALL VAR_FORMSOL(ETA,CHI,ESEC,CHIL,ETAL,V,SIGMA,R,P,
	1                      FQAF,FQAFD,TX,TVX,KI,AQW,
	1                      PF,PROF,LFQW,ERF,FL,DIF,DBB,IC,AMASS,
	1                      THK_LINE,THK_CONT,
	1                      NLF,NC,NP,ND,NNM,METHOD)
	  CALL TUNE(2,'VARFOR')
	END IF
C 
C
C Zero VZNET matrix
C
	CALL ZERO(VZNET,NT*NUM_BNDS*ND)
C
C Convert FQAF from dJBAR to dZNET. We note that dZnet = -dJ/S
C -(J/etal)dchil + (J/S/etal)detal
C
	DO K=1,ND				!Intensity depth
	  TA(K)=-CHIL(K)/ETAL(K)
	END DO
C
	DO J=1,NNM
	  DO L=1,ND				!Variable depth
	    DO K=1,ND
	      FQAF(K,L,J)=FQAF(K,L,J)*TA(K)
	    END DO
	  END DO
	END DO
C
	DO K=1,ND
	  FQAF(K,K,1)=FQAF(K,K,1)-JBAR(K)/ETAL(K)
	  FQAF(K,K,2)=FQAF(K,K,2)-JBAR(K)*TA(K)/ETAL(K)
          FQAFD(K)=FQAFD(K)*TA(K)
 	END DO
C
C If the negative opacity flag is set, SET dZ/dCHIL to zero for
C for all depths. NB First subscript refers to ZNET, second to
C depth of CHIL.
C
	DO I=1,ND
	  IF(NEG_OPACITY(I))THEN
	    DO K=1,ND
	      F2DA(K,I)=0.0D0
	    END DO
	  END IF
	END DO
C
C Increment variation matrix due to variation in line opacity,
C Line emissivity and the line profile .
C
	T1=OSCIL(1)*OPLIN
	T2=-GLDGU(1)*T1
	T3=EINA(1)*EMLIN*FL
	DO K=1,ND
	  DO J=BNDST(K),BNDEND(K)
	    L=BND_TO_FULL(J,K)
	    VZNET(NL,J,K)=VZNET(NL,J,K)+FQAF(K,L,1)*
	1                        T1*L_STAR_RATIO(L,1)
	    VZNET(NUP,J,K)=VZNET(NUP,J,K)+
	1      (FQAF(K,L,1)*T2+FQAF(K,L,2)*T3)*U_STAR_RATIO(L,1)
	    VZNET(NT,J,K)=VZNET(NT,J,K)+
	1      FQAF(K,L,1)*T1*dL_RAT_dT(L,1)*POPS(NL,L)+
	1      (FQAF(K,L,1)*T2+FQAF(K,L,2)*T3)*dU_RAT_dT(L,1)*POPS(NUP,L)
	  END DO
	END DO
C 
C
C Increment VZNET matrix due to continuum opacity and emissivity
C variation if the variation matrices have been computed. This is
C only done if we are not running the impurity code -- if running
C impurity code the continuous radiation field is fixed.
C
	IF(NNM .EQ. 4 .AND. .NOT. IMPURITY_CODE)THEN
C**************************************************************************
C
C Solve for the perturbations to J in terms of the perturbations
C to CHI and ETA. F2DA is the dCHI matrix ; FC the dETA matrix
C and FA is the d(diffusion approx) vector for the boundary
C condition at the core.
C
	  IF(THIS_FREQ_EXT)THEN
	    CALL TUNE(1,'DJFEAUEXT')
	      CALL PERTJFEAU_IBC(F2DAEXT,FCEXT,FAEXT,
	1            DTAU,CHIEXT,REXT,ZETAEXT,
	1            THETAEXT,RJEXT,QEXT,FEXT,dCHIdR,
	1            TA,TB,TC,HBC_J,HBC_S,INBC,DBB,DIF,
	1            THK_CONT,NDEXT,METHOD)
C
C Put variation matrices on old grid. Note that FA is used for the
C diffusion approximation.
C
	      CALL REGRID_dCHI(F2DA,CHI,ND,POS_IN_NEW_GRID,
	1                          F2DAEXT,CHIEXT,NDEXT,COEF,INDX)
	      CALL REGRID_dCHI(FC,ETA,ND,POS_IN_NEW_GRID,
	1                          FCEXT,ETAEXT,NDEXT,COEF,INDX)
	      DO I=1,ND
	        FA(I)=FAEXT(POS_IN_NEW_GRID(I))
	      END DO
	    CALL TUNE(2,'DJFEAUEXT')
	  ELSE IF(EDDINGTON)THEN
	    CALL TUNE(1,'PERTJFEAU')
	      CALL PERTJFEAU_IBC(F2DA,FC,FA,
	1            DTAU,CHI,R,ZETA,
	1            THETA,RJ,QEDD,FEDD,dCHIdR,
	1            TA,TB,TC,HBC_J,HBC_S,INBC,DBB,DIF,
	1            THK_CONT,ND,METHOD)
	    CALL TUNE(2,'PERTJFEAU')
	  ELSE
	    CALL TUNE(1,'PERTJD')
	      CALL NEWPERTJD(F2DA,FC,FA,FB,VK,WM,AQW
	1       ,DTAU,CHI,dCHIdR,R,Z,P,THETA,SOURCE,TA,TB,TC,XM
	1       ,DIF,DBB,IC,ESEC,THK_CONT,NC,ND,NP,METHOD)
	    CALL TUNE(2,'PERTJD')
	  END IF
C
C Update %CHI and %ETA matrices to allow for electron scatering.
C
	  CALL FQNEW(FQAF,F2DA,FC,WM,ESEC,ND,NNM)
C
	  IF(DIF)THEN
	    DO I=1,ND
	      DO J=1,ND
	        FQAFD(I)=FQAFD(I)+FQAF(I,J,4)*FA(J)*ESEC(J)
	      END DO
	    END DO
	  END IF
C 
C
C Compute the opacity AND emisisvity variation as a function of the changes 
C in the population levels. 
C
C	  INCLUDE 'VAROPAC_V4.INC'
	  CALL COMP_VAR_OPAC(POPS,RJ,FL,CONT_FREQ,FREQ_INDX,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
C 
C
C	  Increment the J variation matrices due to opacity variation
C
	  DO K=1,ND	  	  			  !Depth of intensity
	    DO J=BNDST(K),BNDEND(K)	  	  	  !Depth of variable
	      L=BND_TO_FULL(J,K)
	      DO I=1,NT	  	  		  !Variable
	        VZNET(I,J,K)=VZNET(I,J,K)+
	1          ( FQAF(K,L,3)*VCHI(I,L)+FQAF(K,L,4)*VETA(I,L) )
	      END DO
	    END DO
	  END DO
C
	END IF		!NNM IF
C 
C
C Update VZNET for perturbations in diffusion approximation.
C 
	IF(DIF .AND. NUM_BNDS .EQ. ND)THEN
	  T1=DBB/DTDR
	  DO J=1,ND
	    DO I=1,NT-1
	      VZNET(I,NUM_BNDS,J)=VZNET(I,NUM_BNDS,J)+FQAFD(J)*T1*DIFFW(I)
	    END DO
	    VZNET(NT,NUM_BNDS,J)=VZNET(NT,NUM_BNDS,J) +
	1                           FQAFD(J)*(DDBBDT+T1*DIFFW(NT))
	  END DO
	ELSE
	  T1=DBB/DTDR
	  DO J=DIAG_INDX,NUM_BNDS
	    K=ND+DIAG_INDX-J
	    DO I=1,NT-1
	      VZNET(I,J,K)=VZNET(I,J,K)+FQAFD(K)*T1*DIFFW(I)
	    END DO
	    VZNET(NT,J,K)=VZNET(NT,J,K)+FQAFD(K)*(DDBBDT+T1*DIFFW(NT))
	  END DO
	END IF
C
C Increment the BA matrix due to Jmn line terms. The variation of
C U_STAR with T is now allowed for.
C
	I=SIM_LINE_POINTER(1)
	ID=VEC_ID(I)
	MNL_F=VEC_MNL_F(I);     MNL=ATM(ID)%F_TO_S_XzV(MNL_F)
	MNUP_F=VEC_MNUP_F(I);   MNUP=ATM(ID)%F_TO_S_XzV(MNUP_F)
	DO L=1,ND					!Equation depth
	  K=GET_DIAG(L)
	  T2=EINA(1)*ZNET(L)*U_STAR_RATIO(L,1)
	  SE(ID)%BA(MNL,MNUP,K,L)=SE(ID)%BA(MNL,MNUP,K,L) + T2
	  SE(ID)%BA(MNUP,MNUP,K,L)=SE(ID)%BA(MNUP,MNUP,K,L) - T2
	  BA_T(NUP,K,L)=BA_T(NUP,K,L) - EMLIN*FL*T2
	  T2=EINA(1)*ZNET(L)*dU_RAT_dT(L,1)*POPS(NUP,L)
          J=SE(ID)%N_IV
	  SE(ID)%BA(MNL,J,K,L)=SE(ID)%BA(MNL,J,K,L) + T2
	  SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L) - T2
	  BA_T(NT,K,L)=BA_T(NT,K,L) - EMLIN*FL*T2
	  T1=EINA(1)*U_STAR_RATIO(L,1)*POPS(NUP,L)
	  DO K=1,NUM_BNDS				!Variable depth
	    DO J=1,SE(ID)%N_IV
	      JJ=SE(ID)%LNK_TO_F(J)
	      SE(ID)%BA(MNL,J,K,L) =SE(ID)%BA(MNL,J,K,L) +VZNET(JJ,K,L)*T1
	      SE(ID)%BA(MNUP,J,K,L)=SE(ID)%BA(MNUP,J,K,L)-VZNET(JJ,K,L)*T1
	      BA_T(JJ,K,L)=BA_T(JJ,K,L)-VZNET(JJ,K,L)*ETAL(L)
	    END DO
	  END DO
	END DO
