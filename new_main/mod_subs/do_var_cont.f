!
! Subroutine to compute the variation of J (continuum) with respect to the populations.
! Line blanketing and velocity fields can be treated. This routine replaces VARCONT.INC
! Routine computes
!                   dJ(I,B,K) = dJ(K)/dN(I)
!         where K is the J depth, I the variable, and B the variable depth.
!
! Available options for computing dJ. Thes should be the same as for COMP_J_BLANK:
!
!       CONT_VEL=.TRUE.			!Velocity field taken into account.
!            ACURATE=.TRUE.		!Enhanced spatial grid
!            ACCURATE=.FALSE.		!Regular spatial grid
!
!If the following options hold, we use Eddington factors and the enhanced grid.
!
!       CONT_VEL=.FALSE.  & THIS_FREQ_EXT=.TRUE.
!
!If the following options hold, we use Eddington factors and the regular grid.
!
!	CONT_VEL=.FALSE., EDDINGTON=.TRUE. & THIS_FREQ_EXT=.FALSE.
!
!If the following options hold we use ray-by ray solution for the full computaion.
!Spherical model, and no velocity terms.
!
!	CONT_VEL=.FALSE., EDDINGTON=.FALSE. & THIS_FREQ_EXT=.FALSE.
!
	SUBROUTINE DO_VAR_CONT(POPS,SECTION,EDDINGTON,
	1                    FL,CONT_FREQ,FREQ_INDX,FIRST_FREQ,TX_OFFSET,
	1                    ND,NC,NP,NUM_BNDS,DIAG_INDX,NT,NM,
	1                    NDEXT,NCEXT,NPEXT,MAX_SIM,NM_KI)
	USE ANG_QW_MOD
	USE MOD_CMFGEN
	USE OPAC_MOD
	USE LINE_MOD 
	USE RADIATION_MOD
	USE VAR_RAD_MOD
	USE CONTROL_VARIABLE_MOD
	IMPLICIT NONE
!
! Altered 17-Oct-2016 : Variation routines changed to: VAR_MOM_J_CMF_V12.F, VAR_MOM_J_DDT_V4.F, and VAR_JREL_V4.F .
! Incorporated 02-Jan-2014: Changed to allow depth dependent profiles.
! Altered 01-Dec-2012 : Added loop variables to OMP calls.
! Altered 31-Jan-2010 : Changed to VAR_MOM_J_DDT_V2 from VAR_MOM_J_DDT_V1.
! Altered 09-Jan-2010 : Included variation of line opacity & emissivity with T. This variation arises
!                           from the use of super-levels.
! Finalized 17-Dec-2004
!
	INTEGER NUM_BNDS
	INTEGER DIAG_INDX
	INTEGER ND,NC,NP
	INTEGER NT
	INTEGER NM
!
! NM_KI is the 3rd dimension of KI matrix (dCHI,dETA). Only needs to be 2 for this routine.
!
	INTEGER NM_KI
	INTEGER MAX_SIM
	INTEGER NDEXT,NCEXT,NPEXT
	INTEGER FREQ_INDX
	INTEGER TX_OFFSET
!
	REAL*8 POPS(NT,ND)
	REAL*8 FL				!Current frequency in units of 10^15 Hz
	REAL*8 CONT_FREQ			!frequency at which opacity was evaluated.
!
	CHARACTER(LEN=*) SECTION
	LOGICAL FIRST_FREQ
!
! Use Eddington factors to compute J. This option is needed since continuum and lines, could,
! in principal, use different options.
!
	LOGICAL EDDINGTON
!
! Local variables
!
	REAL*8 T1
	INTEGER X_INDX
	INTEGER I,J,L,K
	INTEGER NL,NUP
	INTEGER LOW,UP
	LOGICAL RAT_TOO_BIG
	LOGICAL LST_DEPTH_ONLY
	REAL*8, SAVE :: FL_OLD
!
! Constants for opacity etc. These are set in CMFGEN.
!
        COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
        COMMON/LINE/ OPLIN,EMLIN
        REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
        REAL*8 OPLIN,EMLIN
!
! These two functions compute the start and end indices when updating
! VJ. eg. we do not wish to update VJ( ,1, ) if we are using the banded
! matrix since this refers to a variable beyond the outer atmosphere.
!
	INTEGER BND_TO_FULL
	INTEGER BNDST
	INTEGER BNDEND
!
        BNDST(K)=MAX( (NUM_BNDS/ND)*(K-DIAG_INDX)+1+DIAG_INDX-K, 1 )
        BNDEND(K)=MIN( (NUM_BNDS/ND)*(K-DIAG_INDX)+ND+DIAG_INDX-K,NUM_BNDS )
!
! This function takes a band-index and converts it the equivalent index
! in the full matrix. L=BND_TO_FULL(J,K) is equivalent to the statements:
!     IF(NUM_BNDS .EQ. ND)THEN L=J ELSE L=K+J-DIAG_INDX END IF
! The second indice is the equation depth.
!
        BND_TO_FULL(J,K)=(NUM_BNDS/ND)*(DIAG_INDX-K)+K+J-DIAG_INDX
!
!**************************************************************************
!
! Include file to calculate the variation of J as a function of
! the population values. For use with the Sobolev approximation,
! and the Dielectronic recombination section.
!
! Altered 17-Jun-2003 - RJEXT_ES (not RJ_ES) used for non-coherent electron 
!                         scattering when updating the emissivity in the ACCURATE
!                         and CONT_VEL section.
! Altered 11-Dec-1997 - Section for computing J in CONT_VEL section  with
!                         additional data points inserted. Designed to give
!                         higher accuracy, especially near ionization fronts.
!                         Additional points are inserted for ALL frequencies.
!                         For this option, ACCURATE must be set to TRUE.
!
! Altered 12-Jun-91 - Call to PERTJFEAUNEW replaced by call to
!                     PERTJFEAU_IBC. This has improved handling of
!                     the outer boundary condition.
! Altered 13-Feb-88 - Changed to allow Banded-Newton Raphson Operators.
!                     The number of bands is given by NUM_BNDS and must
!                     be odd. If NUM_BNDS=ND, the matrix is assumed
!                     to represent the full Newton-Raphson Operator.
!                     The index, DIAG_INDX must also be specified, but is
!                     only used for the BANDED case. Note 
!                     DIAG_INDX=( NUM_BNDS/2+1 )
! Altered 6-May-86 -  PERTJFEAU installed. (now varcontfeau)
! Created 26-Jan-88 - (Final version).
!
! BA and VJ have the following dimensions.
!
!                     BA(NT,NT,NUM_BNDS,ND)
!                     VJ(NT,NUM_BNDS,ND)
!
!                     In addition, the variables HYD_PRES, HEI_PRES
!                     and HE2_PRES have been installed. These have
!                     been included so the routine is completely general.
!                     These LOGICAL variables should be specified by
!                     parameter statements - they can then be optimized
!                     out of the code at compile time.
!
! Altered 21-Mar-1989 - VETA routine installed. VETA and VCHI are now
!                       computed simultaneously.
!
! Altered 04-Apr-89 - EDDINGTON logical variable installed.
!                     THIS_FREQ_EXT restored (was THIS_FREQ).
!                     ACCURATE variable removed.
!
! Altered 31-May-1989 - THICK variable replaced by THK_CONT.
!
! 
!**************************************************************************
!
! Indicates to  COMP_VAR_OPAC that we are computing the opacity at ALL depths.
!
	LST_DEPTH_ONLY=.FALSE.
!
! Solve for the perturbations to J in terms of the perturbations
! to CHI and ETA. F2DA is the dCHI matrix ; FC the dETA matrix
! and FA is the d(diffusion approx) vector for the boundary
! condition at the core.
!
	  IF(THIS_FREQ_EXT .AND. .NOT. CONT_VEL)THEN
	    CALL TUNE(1,'DJFEAUEXT')
	      CALL PERTJFEAU_IBC(F2DAEXT,FCEXT,FAEXT,
	1            DTAU,CHIEXT,REXT,ZETAEXT,
	1            THETAEXT,RJEXT,QEXT,FEXT,dCHIdR,
	1            TA,TB,TC,HBC_J,HBC_S,INBC,DBB,DIF,
	1            THK_CONT,NDEXT,METHOD)
!
! Put variation matrices on old grid. Note that FA is used for the
! diffusion approximation.
!
	      CALL REGRID_dCHI(F2DA,CHI,ND,POS_IN_NEW_GRID,
	1                          F2DAEXT,CHIEXT,NDEXT,COEF,INDX)
	      CALL REGRID_dCHI(FC,ETA,ND,POS_IN_NEW_GRID,
	1                          FCEXT,ETAEXT,NDEXT,COEF,INDX)
	      DO I=1,ND
	        FA(I)=FAEXT(POS_IN_NEW_GRID(I))
	      END DO
	    CALL TUNE(2,'DJFEAUEXT')
!
! 
!
	  ELSE IF(CONT_VEL .AND. .NOT. ACCURATE)THEN
	    IF(FIRST_FREQ)THEN
	      TX(:,:,:)=0.0D0
	      TVX(:,:,:)=0.0D0
	      dJ_DIF_d_T(:)=0.0D0
	      dJ_DIF_d_dTdR(:)=0.0D0
	      dRSQH_DIF_d_T=0.0D0
	      dRSQH_DIF_d_dTdR=0.0D0
	      FL_OLD=FL
	    ELSE
	      RAT_TOO_BIG=.FALSE.
	      DO L=1,ND
	        TA(L)=CHI_NOSCAT_PREV(L)/CHI_NOSCAT(L)
	        IF(ETA_CONT(L) .EQ. 0.0D0)THEN
                  TB(L)=1.0D0
	        ELSE
	          TB(L)=ETA_PREV(L)/ETA_CONT(L)
	        END IF
	        IF(TA(L) .GT. 5.0D0)THEN
	          TA(L)=0.0D0; TB(L)=0.0D0
	        END IF
!	        IF(TA(L) .GT. 1.5)RAT_TOO_BIG=.TRUE.
	      END DO
	      IF(RAT_TOO_BIG)THEN
	        DO L=1,ND
	          TA(L)=0.0D0
	          TB(L)=0.0D0
	        END DO
	      END IF
!$OMP PARALLEL DO PRIVATE(J,K,T1)
	      DO J=1,ND
	        T1=ETA_CONT(J)*HDKT*(FL_OLD-FL)/T(J)/T(J)
	        DO K=1,ND
	          TX(K,J,3)=TX(K,J,3) * TA(J)
	          TX(K,J,4)=TX(K,J,4) * TB(J)
	          TX(K,J,6)=TX(K,J,6) + TX(K,J,4)*T1
	        END DO
	      END DO
!$OMP PARALLEL DO PRIVATE(J,K,T1)
	      DO J=1,ND
	        T1=ETA_CONT(J)*HDKT*(FL_OLD-FL)/T(J)/T(J)
	        DO K=1,ND-1
	          TVX(K,J,3)=TVX(K,J,3) * TA(J)
	          TVX(K,J,4)=TVX(K,J,4) * TB(J)
	          TVX(K,J,6)=TVX(K,J,6) + TVX(K,J,4)*T1
	        END DO
	      END DO
	    END IF
	    DO I=1,NM
	      DO_THIS_TX_MATRIX(I)=.TRUE.
	    END DO
	    DO I=TX_OFFSET+1,NM
	      IF(VAR_IN_USE_CNT(I) .EQ. 0)THEN
	        DO_THIS_TX_MATRIX(I)=.FALSE.
	      END IF
	    END DO
!
! Use TA as temporary storage for the emissivity.
!
	    IF(COHERENT_ES)THEN
	      TA(1:ND)=ETA_CLUMP(1:ND)
	      ES_COH_VEC(1:ND)=CHI_SCAT_CLUMP(1:ND)/CHI_CLUMP(1:ND)
	    ELSE
!
! Two scenarios: 
!    (i) We use a lambda iteration to allow for the variation of J in
!           the electron scattering term. ES_COH_VEC (==THETA) is zero,
!           and the electron scattering emissivity (related to RJ_ES)
!           is included with the emissivity directly.
!   (ii) We assume that the variation of J in the electron scattering
!           term can be treated coherently. ES_COH_VEC is non-zero.
!           A correction to the emissivity is made to allow for the
!           difference between RJ and RJ_ES.
!
!	      ES_COH_VEC(1:ND)=0.0D0
!	      IF(MIXED_ES_VAR)THEN
!	        DO I=ND,1,-1
!	          IF( ABS(RJ(I)-RJ_ES(I))/(RJ(I)+RJ_ES(I)) .GT.
!	1                          0.5D0*ES_VAR_FAC)EXIT
!	          ES_COH_VEC(I)=ESEC_CLUMP(I)/CHI_CLUMP(I)
!	        END DO
!	      END IF
!	      DO I=1,ND
!	        IF(ES_COH_VEC(I) .EQ. 0.0D0)THEN
!	          TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*RJ_ES(I)
!	        ELSE
!	          TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*(RJ_ES(I)-RJ(I))
!	        END IF
!	      END DO
!	    END IF
!
	      IF(MIXED_ES_VAR)THEN
	        T1=2.0D0
	        DO I=1,ND
	          IF(RJ_ES(I) .GT. RJ(I))THEN
	            TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*
	1                   (RJ_ES(I)-RJ(I)/T1)
	            ES_COH_VEC(I)=ESEC_CLUMP(I)/CHI_CLUMP(I)/T1
	          ELSE
	            TA(I)=ETA_CLUMP(I)+ESEC_CLUMP(I)*RJ_ES(I)/T1
	            ES_COH_VEC(I)=ESEC_CLUMP(I)/CHI_CLUMP(I)*
	1                    (RJ_ES(I)/RJ(I))/T1
	          END IF
	        END DO
	      ELSE
	        ES_COH_VEC(1:ND)=0.0D0
	        TA(1:ND)=ETA_CLUMP(1:ND)+ESEC_CLUMP(1:ND)*RJ_ES(1:ND)
	      END IF
	   END IF
!
	   CALL TUNE(1,'VAR_MOM_J')
	   IF(PLANE_PARALLEL_NO_V)THEN
	     CALL VAR_MOM_PP_V1(R,TA,CHI_CLUMP,CHI_SCAT_CLUMP,FEDD,
	1           TX,dJ_DIF_d_T,dJ_DIF_d_dTdR,DO_THIS_TX_MATRIX,
	1           HBC_CMF,NBC_CMF,INBC,
	1           DIF,DBB,dDBBdT,dTdR,IC,METHOD,COHERENT_ES,ND,NM)
	   ELSE IF(PLANE_PARALLEL)THEN
	     CALL PP_VAR_MOM_CMF_V1(TA,CHI_CLUMP,CHI_SCAT_CLUMP,V,SIGMA,R,
	1           TX,TVX,dJ_DIF_d_T,dJ_DIF_d_dTdR, 
	1           dRSQH_DIF_d_T,dRSQH_DIF_d_dTdR,FEDD,GEDD,N_ON_J,
	1           INBC,HBC_CMF(1),HBC_CMF(2),NBC_CMF(1),NBC_CMF(2),
	1           FIRST_FREQ,L_FALSE,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1           DO_THIS_TX_MATRIX,METHOD,COHERENT_ES,ND,NM)
           ELSE IF(USE_J_REL)THEN
             CALL VAR_JREL_V4(TA,CHI_CLUMP,CHI_SCAT_CLUMP,ES_COH_VEC,V,SIGMA,R,
	1                  TX,TVX,dJ_DIF_d_T,dJ_DIF_d_dTdR,
	1                  dRSQH_DIF_d_T,dRSQH_DIF_d_dTdR,KI,WM,RHS_dHdCHI,
	1                  FIRST_FREQ,FL,dLOG_NU,H_CHK_OPTION,
	1                  INNER_BND_METH,OUTER_BND_METH,IB_STAB_FACTOR,
	1                  dTdR,DBB,dDBBdT,IC,
	1                  INCL_ADVEC_TERMS_IN_TRANS_EQ,INCL_REL_TERMS,
	1                  DO_THIS_TX_MATRIX,METHOD,ND,NM,NM_KI)
	   ELSE IF (USE_LAM_ES)THEN
	   ELSE IF(USE_DJDT_RTE)THEN
	     CALL VAR_MOM_J_DDT_V4(TA,CHI_CLUMP,CHI_SCAT_CLUMP,ES_COH_VEC,V,R,
	1           TX,TVX,dJ_DIF_d_T,dJ_DIF_d_dTdR,
	1           dRSQH_DIF_d_T,dRSQH_DIF_d_dTdR,
	1           KI,WM,RHS_dHdCHI,FEDD,
	1           INCL_DJDT_TERMS,USE_DR4JDT,DJDT_RELAX_PARAM,FIRST_FREQ,FL,dLOG_NU,
	1           dTdR,DBB,dDBBdT,DO_THIS_TX_MATRIX,METHOD,
	1           H_CHK_OPTION,INNER_BND_METH,OUTER_BND_METH,
	1           ND,NM,NM_KI)
	   ELSE
	    CALL VAR_MOM_J_CMF_V12(TA,CHI_CLUMP,CHI_SCAT_CLUMP,
	1           ES_COH_VEC,V,SIGMA,R,
	1           TX,TVX,dJ_DIF_d_T,dJ_DIF_d_dTdR, 
	1           dRSQH_DIF_d_T,dRSQH_DIF_d_dTdR, 
	1           KI,WM,RHS_dHdCHI,
	1           FIRST_FREQ,dLOG_NU,
	1           INNER_BND_METH,dTdR,DBB,dDBBdT,IC,IB_STAB_FACTOR,
	1           FL,H_CHK_OPTION,OUT_BC_TYPE,DO_THIS_TX_MATRIX,
	1           METHOD,ND,NM,NM_KI)
	   END IF
	   CALL TUNE(2,'VAR_MOM_J')
!
! Correcting for clumping this way does it for both the continuum and lines.
!
	    IF(DO_CLUMP_MODEL)THEN
!$OMP PARALLEL DO
	      DO J=1,ND
	        DO K=1,ND
	          TX(K,J,1)=TX(K,J,1)*CLUMP_FAC(J)
	          TX(K,J,2)=TX(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
!$OMP PARALLEL DO
	      DO J=1,ND
	        DO K=1,ND-1
	          TVX(K,J,1)=TVX(K,J,1)*CLUMP_FAC(J)
	          TVX(K,J,2)=TVX(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
	    END IF
!
! For many ALO like calculations, only the variation of J with the source function
! is allowed for. The following allows us to do the same thing.
!
	   IF(DO_SRCE_VAR_ONLY)THEN
	      DO J=1,ND
	       T1=-ETA_CLUMP(J)/CHI_CLUMP(J)
	        DO K=1,ND
	          TX(K,J,1)=TX(K,J,2)*T1
	        END DO
	      END DO
	      DO J=1,ND
	       T1=-ETA_CLUMP(J)/CHI_CLUMP(J)
	        DO K=1,ND-1
	          TVX(K,J,1)=TVX(K,J,2)*T1
	        END DO
	      END DO
	   END IF
!
!	   K=70
!	   T1=CHI(K)-CHI_SCAT(K)
!	   WRITE(221,'(12ES16.8)')FL,ETA(K),CHI(K),T1,RJ(K),ETA(K)/T1,
!	1             TX(K,K,1),TX(K,K,3),TX(K,K,2),TX(K,K,4),
!	1             1.0D0+(TX(K,K,1)+TX(K,K,3))*T1/RJ(K), 
!	1             (TX(K,K,2)+TX(K,K,4))*T1-1.0D0
!	   WRITE(224,'(12ES16.8)')FL,ETA_CONT(K),CHI_CONT(K),T1,RJ(K),TX(K,K,1:6)
!
! We use TB as a temporary vector for J in the electron scattering emissivity.
! Its value depends on whether we have coherent or incoherent e.s.
!
	    IF(COHERENT_ES)THEN
	      TB(1:ND)=RJ(1:ND)
	    ELSE
	      TB(1:ND)=RJ_ES(1:ND)
	    END IF
!$OMP PARALLEL DO
	    DO J=1,ND
	      DO K=1,ND
	        TX(K,J,3)=TX(K,J,3) + TX(K,J,1)
	        TX(K,J,4)=TX(K,J,4) + TX(K,J,2)
	        TX(K,J,5)=TX(K,J,5) + TX(K,J,1) + TX(K,J,2)*TB(J)
	      END DO
	    END DO
!$OMP PARALLEL DO
	    DO J=1,ND
	      DO K=1,ND-1
	        TVX(K,J,3)=TVX(K,J,3) + TVX(K,J,1)
	        TVX(K,J,4)=TVX(K,J,4) + TVX(K,J,2)
	        TVX(K,J,5)=TVX(K,J,5) + TVX(K,J,1) + TVX(K,J,2)*TB(J)
	      END DO
	    END DO
!
!	   K=70
!	   T1=CHI(K)-CHI_SCAT(K)
!	   WRITE(225,'(12ES16.8)')FL,ETA(K),CHI(K),T1,RJ(K),ETA(K)/T1,
!	1             1.0D0+TX(K,K,3)*T1/RJ(K),TX(K,K,4)*T1-1.0D0,TX(K,K,5)*ESEC(K)/RJ(K)
!
! Update line variation matrices. Note that the matrices now refer to the
! variation with respect to levels (e.g. the lower and upper level) and
! not CHIL and ETAL. This should save arrays when we are dealing with
! overlapping transitions between the same level(s).
!
! For simplicity we have ignored the T dependance of L_STAR_RATIO and
! U_STAR_RATIO.
!
! NB: We paralleize over the second loop, rather than SIM_INDX, as the
!     variables LOW and UP may be the same for different SIM_INDX values.
!
	    CALL TUNE(1,'TX_TVX_VC')
	    DO SIM_INDX=1,MAX_SIM
	      LOW=LOW_POINTER(SIM_INDX);    UP=UP_POINTER(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
!$OMP PARALLEL DO PRIVATE(J,K,OPAC_FAC, EMIS_FAC, STIM_FAC)
	        DO J=1,ND
	          OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*NEG_OPAC_FAC(J)
	          STIM_FAC=OPAC_FAC*U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          OPAC_FAC=OPAC_FAC*L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*U_STAR_RATIO(J,SIM_INDX)
	          DO K=1,ND
	            TX(K,J,LOW)=TX(K,J,LOW) + OPAC_FAC*TX(K,J,1)
	          END DO
	          DO K=1,ND
	            TX(K,J,UP)=TX(K,J,UP) + ( EMIS_FAC*TX(K,J,2) - STIM_FAC*TX(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
!
! We now do the update for dH (i.e. TVX)
!
	    DO SIM_INDX=1,MAX_SIM
	      LOW=LOW_POINTER(SIM_INDX);    UP=UP_POINTER(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
!$OMP PARALLEL DO PRIVATE(J,K,OPAC_FAC, EMIS_FAC, STIM_FAC)
	        DO J=1,ND
	          OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*NEG_OPAC_FAC(J)
	          STIM_FAC=OPAC_FAC*U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          OPAC_FAC=OPAC_FAC*L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*U_STAR_RATIO(J,SIM_INDX)
	          DO K=1,ND-1
	            TVX(K,J,LOW)=TVX(K,J,LOW) + OPAC_FAC*TVX(K,J,1)
	          END DO
	          DO K=1,ND-1
	            TVX(K,J,UP)=TVX(K,J,UP) + ( EMIS_FAC*TVX(K,J,2) - STIM_FAC*TVX(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
	    CALL TUNE(2,'TX_TVX_VC')
!
	    IF(INCLUDE_dSLdT)THEN
	      DO SIM_INDX=1,MAX_SIM
	        NL=SIM_NL(SIM_INDX)
	        NUP=SIM_NUP(SIM_INDX)
	        IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
!$OMP PARALLEL DO PRIVATE(J,K,OPAC_FAC, EMIS_FAC)
	          DO J=1,ND
	            OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*NEG_OPAC_FAC(J)*
	1              (dL_RAT_dT(J,SIM_INDX)*POPS(NL,J)-GLDGU(SIM_INDX)*dU_RAT_dT(J,SIM_INDX)*POPS(NUP,J))
	            EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*dU_RAT_dT(J,SIM_INDX)*POPS(NUP,J)
	            DO K=1,ND
	              TX(K,J,6)=TX(K,J,6) + (OPAC_FAC*TX(K,J,1)+EMIS_FAC*TX(K,J,2))
	            END DO
	            DO K=1,ND-1
	              TVX(K,J,6)=TVX(K,J,6) + (OPAC_FAC*TVX(K,J,1)+EMIS_FAC*TVX(K,J,2))
	            END DO
	          END DO
	        END IF
	      END DO
	    END IF
!
! Now zero dCHI and dETA storage locations. These refer to the TOTAL 
! opacity and emissivity.
!
	    TX(:,:,1:2)=0.0D0
	    TVX(:,:,1:2)=0.0D0
!
! 
!
	  ELSE IF(CONT_VEL .AND. ACCURATE)THEN
	    IF(FIRST_FREQ)THEN
	      TX_EXT(:,:,:)=0.0D0
	      TVX_EXT(:,:,:)=0.0D0
	    ELSE
	      RAT_TOO_BIG=.FALSE.
	      DO L=1,ND
	        TA(L)=CHI_NOSCAT_PREV(L)/CHI_NOSCAT(L)
	        TB(L)=ETA_PREV(L)/ETA_CONT(L)
	        IF(TA(L) .GT. 1.5)RAT_TOO_BIG=.TRUE.
	      END DO
	      IF(RAT_TOO_BIG)THEN
	        DO L=1,ND
	          TA(L)=0.0D0
	          TB(L)=0.0D0
	        END DO
	      END IF
!$OMP PARALLEL DO
	      DO J=1,ND
	        DO K=1,NDEXT
	          TX_EXT(K,J,3)=TX_EXT(K,J,3) * TA(J)
	          TX_EXT(K,J,4)=TX_EXT(K,J,4) * TB(J)
	        END DO
	      END DO
!$OMP PARALLEL DO
	      DO J=1,ND
	        DO K=1,NDEXT-1
	          TVX_EXT(K,J,3)=TVX_EXT(K,J,3) * TA(J)
	          TVX_EXT(K,J,4)=TVX_EXT(K,J,4) * TB(J)
	        END DO
	      END DO
	    END IF
	    DO I=1,NM
	      DO_THIS_TX_MATRIX(I)=.TRUE.
	    END DO
	    DO I=TX_OFFSET+1,NM
	      IF(VAR_IN_USE_CNT(I) .EQ. 0)THEN
	        DO_THIS_TX_MATRIX(I)=.FALSE.
	      END IF
	    END DO
!
! Use TA as temporary storage for the emissivity.
!
	    IF(COHERENT_ES)THEN
	      TA(1:NDEXT)=ETAEXT(1:NDEXT)
	    ELSE
	      TA(1:NDEXT)=ETAEXT(1:NDEXT)+
	1                 ESECEXT(1:NDEXT)*RJEXT_ES(1:NDEXT)
	    END IF
	    CALL VAR_MOM_JEXT_CMF_V3(TA,CHIEXT,ESECEXT,
	1                  VEXT,SIGMAEXT,REXT,
	1                  ETA_CLUMP,CHI_CLUMP,ESEC_CLUMP,
	1                  INDX,COEF,INTERP_TYPE,ND,
	1                  TX_EXT,TVX_EXT,
	1                  dJ_DIF_d_T_EXT,dJ_DIF_d_dTdR_EXT,
	1                  dRSQH_DIF_d_T,dRSQH_DIF_d_dTdR, 
	1                  KI,RHS_dHdCHI,
	1                  FIRST_FREQ,dLOG_NU,DIF,dTdR,DBB,dDBBdT,IC,
	1                  DO_THIS_TX_MATRIX,METHOD,COHERENT_ES,NDEXT,NM,NM_KI)
!
! Correcting for clumping this way does it for both the continuum and lines.
!
	    IF(DO_CLUMP_MODEL)THEN
!$OMP PARALLEL DO
	      DO J=1,ND
	        DO K=1,NDEXT
	          TX_EXT(K,J,1)=TX_EXT(K,J,1)*CLUMP_FAC(J)
	          TX_EXT(K,J,2)=TX_EXT(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
!$OMP PARALLEL DO
	      DO J=1,ND
	        DO K=1,NDEXT-1
	          TVX_EXT(K,J,1)=TVX_EXT(K,J,1)*CLUMP_FAC(J)
	          TVX_EXT(K,J,2)=TVX_EXT(K,J,2)*CLUMP_FAC(J)
	        END DO
	      END DO
	    END IF
!
! We use TB as a temporary vector for J in the electron scattering emissivity.
! Its value depends on whether we have coherent or incoherent e.s.
!
	    IF(COHERENT_ES)THEN
	      TB(1:ND)=RJ(1:ND)
	    ELSE
	      TB(1:ND)=RJ_ES(1:ND)
	    END IF
!$OMP PARALLEL DO
	    DO J=1,ND
	      DO K=1,NDEXT
	        TX_EXT(K,J,3)=TX_EXT(K,J,3) + TX_EXT(K,J,1)
	        TX_EXT(K,J,4)=TX_EXT(K,J,4) + TX_EXT(K,J,2)
	        TX_EXT(K,J,5)=TX_EXT(K,J,5) + TX_EXT(K,J,1) + 
	1                                         TX_EXT(K,J,2)*TB(J)
	      END DO
	    END DO
!$OMP PARALLEL DO
	    DO J=1,ND
	      DO K=1,NDEXT-1
	        TVX_EXT(K,J,3)=TVX_EXT(K,J,3) + TVX_EXT(K,J,1)
	        TVX_EXT(K,J,4)=TVX_EXT(K,J,4) + TVX_EXT(K,J,2)
	        TVX_EXT(K,J,5)=TVX_EXT(K,J,5) + TVX_EXT(K,J,1) + 
	1                                        TVX_EXT(K,J,2)*TB(J)
	      END DO
	    END DO
!
! Update line variation matrices. Note that the matrices now refer to the
! variation with respect to levels (e.g. the lower and upper level) and
! not CHIL and ETAL. This should save arrays when we are dealing with
! overlapping transitions between the same level(s).
!
! For simplicity we have ignored the T dependance of L_STAR_RATIO and
! U_STAR_RATIO.
!
	    DO SIM_INDX=1,MAX_SIM
	      LOW=LOW_POINTER(SIM_INDX);    UP=UP_POINTER(SIM_INDX)
	      NL=SIM_NL(SIM_INDX);          NUP=SIM_NUP(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
!$OMP PARALLEL DO PRIVATE(J,K,OPAC_FAC, EMIS_FAC, STIM_FAC)
	        DO J=1,ND
	          OPAC_FAC=lINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*
	1                      U_STAR_RATIO(J,SIM_INDX)
	          STIM_FAC=LINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          DO K=1,NDEXT
	            TX_EXT(K,J,LOW)=TX_EXT(K,J,LOW) + OPAC_FAC*TX_EXT(K,J,1)
	            TX_EXT(K,J,UP)=TX_EXT(K,J,UP) + ( EMIS_FAC*TX_EXT(K,J,2) - STIM_FAC*TX_EXT(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
!
! We now do the update for dH (i.e. TVX_EXT)
!
	    DO SIM_INDX=1,MAX_SIM
	      LOW=LOW_POINTER(SIM_INDX);    UP=UP_POINTER(SIM_INDX)
	      NL=SIM_NL(SIM_INDX);          NUP=SIM_NUP(SIM_INDX)
	      IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
!$OMP PARALLEL DO PRIVATE(J,K,OPAC_FAC, EMIS_FAC, STIM_FAC)
	        DO J=1,ND
	          OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*L_STAR_RATIO(J,SIM_INDX)
	          EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*
	1                      U_STAR_RATIO(J,SIM_INDX)
	          STIM_FAC=LINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*
	1                      NEG_OPAC_FAC(J)*U_STAR_RATIO(J,SIM_INDX)*GLDGU(SIM_INDX)
	          DO K=1,NDEXT-1
	            TVX_EXT(K,J,LOW)=TVX_EXT(K,J,LOW) + OPAC_FAC*TVX_EXT(K,J,1)
	            TVX_EXT(K,J,UP)=TVX_EXT(K,J,UP) + ( EMIS_FAC*TVX_EXT(K,J,2) - STIM_FAC*TVX_EXT(K,J,1) )
	          END DO
	        END DO
	      END IF
	    END DO
!
	    IF(INCLUDE_dSLdT)THEN
	      DO SIM_INDX=1,MAX_SIM
	        NL=SIM_NL(SIM_INDX)
	        NUP=SIM_NUP(SIM_INDX)
	        IF(.NOT. WEAK_LINE(SIM_INDX) .AND. RESONANCE_ZONE(SIM_INDX))THEN
!$OMP PARALLEL DO PRIVATE(J,K,OPAC_FAC, EMIS_FAC)
	          DO J=1,ND
	            OPAC_FAC=LINE_OPAC_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*NEG_OPAC_FAC(J)*
	1               (dL_RAT_dT(J,SIM_INDX)*POPS(NL,J)-GLDGU(SIM_INDX)*dU_RAT_dT(J,SIM_INDX)*POPS(NUP,J))
	            EMIS_FAC=LINE_EMIS_CON(SIM_INDX)*LINE_PROF_SIM(J,SIM_INDX)*dU_RAT_dT(J,SIM_INDX)*POPS(NUP,J)
	            DO K=1,NDEXT
	              TX_EXT(K,J,6)=TX_EXT(K,J,6) + (OPAC_FAC*TX_EXT(K,J,1)+EMIS_FAC*TX_EXT(K,J,2))
	            END DO
	            DO K=1,NDEXT-1
	              TVX_EXT(K,J,6)=TVX_EXT(K,J,6) + (OPAC_FAC*TVX_EXT(K,J,1)+EMIS_FAC*TVX_EXT(K,J,2))
	            END DO
	          END DO
	        END IF
	      END DO
	    END IF
!
! Now zero dCHI and dETA storage locations. These refer to the TOTAL
! opacity and emissivity.
!
	    TX_EXT(:,:,1:2)=0.0D0		!dCHI,dETA
	    TVX_EXT(:,:,1:2)=0.0D0		!dCHI,dETA
!
! Compute dJ/d? just on the small grid. We don't need dH/d? since the
! statistical and radiative equilibrium equations depend only on J.
!
	    DO I=3,NM
	      DO J=1,ND
	        DO K=1,ND
	          TX(K,J,I)=TX_EXT(POS_IN_NEW_GRID(K),J,I)
	        END DO
	      END DO
	    END DO
!
	    DO K=1,ND
	      dJ_DIF_d_T(K)=dJ_DIF_d_T_EXT(POS_IN_NEW_GRID(K))
	      dJ_DIF_d_dTdR(K)=dJ_DIF_d_dTdR_EXT(POS_IN_NEW_GRID(K))
	    END DO
!
! 
!
	  ELSE IF(EDDINGTON)THEN
	    CALL TUNE(1,'PERTJFEAU')
	      CALL PERTJFEAU_IBC(F2DA,FC,FA,
	1            DTAU,CHI_CLUMP,R,ZETA,
	1            THETA,RJ,QEDD,FEDD,dCHIdR,
	1            TA,TB,TC,HBC_J,HBC_S,INBC,DBB,DIF,
	1            THK_CONT,ND,METHOD)
	    CALL TUNE(2,'PERTJFEAU')
	  ELSE
	    CALL TUNE(1,'PERTJD')
	      CALL MULTVEC(SOURCE,ZETA,THETA,RJ,ND)
	      CALL NEWPERTJD(F2DA,FC,FA,FB,VK,WM,AQW
	1       ,DTAU,CHI_CLUMP,dCHIdR,R,Z,P,THETA,SOURCE,TA,TB,TC,XM
	1       ,DIF,DBB,IC,CHI_SCAT,THK_CONT,NC,ND,NP,METHOD)
	    CALL TUNE(2,'PERTJD')
	  END IF
! 
!
! Compute the opacity AND emissivity variation as a function of the changes
! in population levels.
!
!
	CALL TUNE(1,'VAROPAC')
!	INCLUDE 'VAROPAC_V4.INC'
	CALL COMP_VAR_OPAC(POPS,RJ,FL,CONT_FREQ,FREQ_INDX,
	1                SECTION,ND,NT,LST_DEPTH_ONLY)
	CALL TUNE(2,'VAROPAC')
! 
!
! Zero VJ array.
!
	CALL DP_ZERO(VJ,NT*NUM_BNDS*ND) 
!
	IF(CONT_VEL)THEN
!
! NB: We no longer include the variation ESEC variation with CHI, but treat it
! separately. this correction has now been specifically included in VAROPAC.
!
!
! The matrix TX gives dJ ( J depth, X depth, X) where X represent some 
! fundamental parameter such as CHI_C, ETA_C, ESEC, ETAL_ etc
!
! We now convert TX  to a smaller matrix, taking into account only the 
! variation of X at Js depth, and neighbouring depths.
!
! dJ_LOC( X , X depth [1=3:NUM_BNDS], J depth )
!
! We only treat Z_INDX from 3 onwards as the first refer to CHI and ETA which
! have been zeroed.
!
!	  IF(FIRST_FREQ)WRITE(6,*)'ERROR REMOVED dVCHI in DO_VAR_CONT'
!	  DO X_INDX=4,NM
!
	  dJ_LOC=0.0D0
	  DO X_INDX=3,NM
	    IF(DO_THIS_TX_MATRIX(X_INDX))THEN
	      DO K=1,ND
	        DO J=BNDST(K),BNDEND(K)
	          L=BND_TO_FULL(J,K)
	          dJ_LOC(X_INDX,J,K)=TX(K,L,X_INDX)
	        END DO
	      END DO
	    END IF
	  END DO
!
! Compute VJ which gives the variation of J with respect to the atomic
! populations. NB: Electron scattering cross-section (6.65D-15) was replaced 
! by ESEC(L)/ED(L) 24-Sep-1997.
!
! NOPS= ND*NUM_BNDS*( 6NT + 2 + 7NUM_SIM )
!
	  CALL TUNE(1,'VC_VCHI')
!
!$OMP PARALLEL DO PRIVATE(I,L,NL)
!
	  DO K=1,ND
	    DO J=BNDST(K),BNDEND(K)
	      L=BND_TO_FULL(J,K)
	      DO I=1,NT
	        VJ(I,J,K)=VJ(I,J,K) +
	1              ( VCHI(I,L)*dJ_LOC(3,J,K) + VETA(I,L)*dJ_LOC(4,J,K) )
	      END DO
	      VJ(NT-1,J,K)=VJ(NT-1,J,K) + ESEC(L)*dj_LOC(5,J,K)/ED(L)
	      IF(SPECIES_PRES(1))THEN
	        VJ(1,J,K)=VJ(1,J,K) + CHI_RAY(L)*dj_LOC(5,J,K)/ATM(1)%XzV_F(1,L)
	      END IF
!                    
! Now must do line terms.
!
	      DO I=TX_OFFSET+1,NM
	        IF(VAR_IN_USE_CNT(I) .GT. 0 .AND. IMP_TRANS_VEC(I))THEN
	          NL=VAR_LEV_ID(I)
	          VJ(NL,J,K)=VJ(NL,J,K) + dj_LOC(I,J,K)
	        END IF
	      END DO
	    END DO	!Over Variable depth (1:NUM_BNDS)
	  END DO		!Over J depth.
!
!$OMP END PARALLEL DO
!
	  DO K=1,ND
	    DO J=BNDST(K),BNDEND(K)
	      VJ(NT,J,K)=VJ(NT,J,K)+dj_LOC(6,J,K)
	    END DO
	  END DO
!
	  CALL TUNE(2,'VC_VCHI')
!
!	  DO L=1,ND
!	    VCHI(EQNE,L)=VCHI(EQNE,L)+ESEC(L)/ED(L)
!	    VETA(EQNE,L)=VETA(EQNE,L)+ESEC(L)/ED(L)*RJ(L)
!	  END DO
!
! Update VJ for perturbations in diffusion approximation. The case ND=NUM_BNDS
! is no longer treated.
!
	  IF(DIF)THEN
	    T1=DBB/DTDR
	    DO J=DIAG_INDX,NUM_BNDS
	      K=ND+DIAG_INDX-J
	      DO I=1,NT-1
	        VJ(I,J,K)=VJ(I,J,K) + dJ_DIF_d_dTdR(K)*DIFFW(I)
	      END DO  
	      VJ(NT,J,K)=VJ(NT,J,K) +
	1                       dJ_DIF_d_T(K)+dJ_DIF_d_dTdR(K)*DIFFW(NT)
	    END DO
	  END IF
!
	ELSE
!
! Section to increment the variation matrices.
!
	  DO K=1,ND			!Depth of intensity
	    DO J=BNDST(K),BNDEND(K)		!Depth of variable index
	      L=BND_TO_FULL(J,K)
	      DO I=1,NT
	        VJ(I,J,K)=VJ(I,J,K)+
	1          (F2DA(K,L)*VCHI(I,L)+FC(K,L)*VETA(I,L))
	      END DO
	    END DO
	  END DO
!
! Update VJ for perturbations in diffusion approximation.
! 18-Dec-1991 replaced ND in VJ( ,ND,K) by VJ( ,NUM_BNDS,K) to avoid
!             compilations errors when NUM_BNDS .NE. ND
!             (only in first clause)
!
	  IF(DIF .AND. ND .EQ. NUM_BNDS)THEN
	    T1=DBB/DTDR
	    DO K=1,ND
	      DO I=1,NT-1
	        VJ(I,NUM_BNDS,K)=VJ(I,NUM_BNDS,K)+FA(K)*T1*DIFFW(I)
	      END DO
	      VJ(NT,NUM_BNDS,K)=VJ(NT,NUM_BNDS,K)+ 
	1                       FA(K)*(DDBBDT+T1*DIFFW(NT))
	    END DO
	  ELSE IF(DIF)THEN
	    T1=DBB/DTDR
	    DO J=DIAG_INDX,NUM_BNDS
	      K=ND+DIAG_INDX-J
	      DO I=1,NT-1
	        VJ(I,J,K)=VJ(I,J,K)+FA(K)*T1*DIFFW(I)
	      END DO
	      VJ(NT,J,K)=VJ(NT,J,K)+FA(K)*(DDBBDT+T1*DIFFW(NT))
	    END DO
	  END IF
	END IF
!
	FL_OLD=FL
	RETURN
	END
