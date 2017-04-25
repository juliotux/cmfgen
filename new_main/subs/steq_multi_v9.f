!
! Subroutine to compute the value of the statistical equilibrium
! equations and the variation of the statistical equilibrium matrix for
! terms which are radiation field independent.
!
! This routine is specifically designed for the handling of super levels.
! That is, we treat the process in a large atom but assume that the populations
! can be described by a smaller set of levels.
!
! Routine also handles level dissolution.
!
! Notation:
!
!         We use _F to denote populations and variables for the FULL atom,
!            with all terms and levels treated separately.
!	  We use _S to denote populations and variables for the SMALL model
!            atom, with many terms and levels treated as one (i.e using
!            SUPER levels).
!
!
! The COLLISION routine that is called has a special FORM, which is distinct
! from that in STEQGEN_V2.
!
! NB - ZION is the charge on the ion - thus ZHYD=1.0D0
!
! Routine also increments the ionization equilibrium equations. Routine no
! longer works for NUM_BNDS=ND.
!
! At present only collisional ionizations to ground state are considered.
!
! NION is the the first dimension of STEQ[ION]. In general we
! NION would be the total number of ionic species.
!
	SUBROUTINE STEQ_MULTI_V9(CNM,DCNM,ED,T,
	1       HN_S,HNST_S,dlnHNST_S_dlnT,N_S,DI_S,
	1       HN_F,HNST_F_ON_S,W_F,A_F,FEDGE_F,G_F,LEVNAME_F,N_F,
	1       F_TO_S_MAPPING,POP,NEXT_PRES,ZION,
	1       ID,COL_FILE,OMEGA_GEN,
	1       EQGS,NUM_BNDS,ND,NION,
	1       COMPUTE_BA,FIXED_T,LAST_ITERATION,
	1       DST,DEND)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 04-Oct-2016 - Changed to V9 (FIXED_T,LAST_ITERATION added to call)
!                         Now call SUBCOL_MULTI_v6.
! Altered 05-Apr-2011 - Changed to V8.
!                       HNST_F_ON_S (rather than HNST_F) is passed in call.
!                       HNST_F/HNST_S replaced by HNST_F_ON_S - done to faciliate
!                         modifications allowing lower temperatures.
!                       Most of editing done early 2011 
! Altered 18-Feb-2010 : Changed order of times in (HNST_S(J,I)*CNM(J,J))*ED(I)/DI_S(I).
!                         Using HNST*CNM instead of HNST/DI_S*CNM should prevent NANs
!                           over a larger dynamic range of T
! Altered 30-Jan-2004 : Replaced by 0.0 by 0.0D0 everywhere.
!                         Important in cool objects where collisions are more 
!                         important than the radiation field.
! Altered 01-Arp-2001 : Changed to V7
!                       Use STEQ_DATA_MOD. Call changed.
! Altered 14-Mar-2001 : Changed to V6
!                       Based on STEQ_MULT_V3 (&V4)
!                       NIV,LNK_F_TO_IV,COMPUTE_BA passed in call
!                       BA no longer altered by this routine, if BA is not
!                       being computed.
!                       Allow for treatment of important variables only.
!                       
!
! Altered 20-Sep-1999 : TMP_VEC_ED and TMP_VEC_COOL used in call to 
!                                                         SUBCOL_MULTI_V3
! Altered 14-Dec-1996 : SUB_PHOT replaces PHOT_FUN (superficial).
! Altered 15-Jun-1996 : T1 initialized before being passed to SUMBCOL_MULTI_V3.
! Altered 26-May-1996 : N_F_MAX removed. Now use dynamic memory allocation
!                         for OMEGA_F etc.
! Altered 03-Nov-1995 : Version changed to _V3
!                       HN_F inserted in call to SUBCOL_MULTI_V3 (prev. _V2)
!
! Altered 10-Nov-1995 : Call to CUBCOL_MULTI_V2 updated.
! Altered 27-Oct-1995 : Call altered to handle new SUBCOL routine.
!                       Now _V2.
! Altered 07-Jun-1995 : Bug fix. Wrong values of HNLTE_S etc being
!                        passed to SUBCOL (effectively those at d=1).
! Created 16-May-1995 : Based on STEQGEN_V2
!
	EXTERNAL OMEGA_GEN
!
	INTEGER NT
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER NION
	INTEGER EQGS
	INTEGER DST,DEND
!
! CNM, and DCNM are used as work arrays. DCNM refers to dCNMdT
!
	INTEGER ID
	INTEGER N_S,N_F
	REAL*8 CNM(N_S,N_S),DCNM(N_S,N_S)
!
	REAL*8 T(ND)
	REAL*8 ED(ND)
	REAL*8 DI_S(ND)
!
	REAL*8 HN_S(N_S,ND)
	REAL*8 HNST_S(N_S,ND)
	REAL*8 dlnHNST_S_dlnT(N_S,ND)
!
	REAL*8 HN_F(N_F,ND)
	REAL*8 HNST_F_ON_S(N_F,ND)
	REAL*8 W_F(N_F,ND)
	REAL*8 A_F(N_F,N_F)
	REAL*8 FEDGE_F(N_F)
	REAL*8 G_F(N_F)
	CHARACTER*(*) LEVNAME_F(N_F),COL_FILE
	INTEGER F_TO_S_MAPPING(N_F)
	REAL*8 ZION
!
	REAL*8 POP(ND)		!Population of species.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	LOGICAL NEXT_PRES
	LOGICAL COMPUTE_BA,FIXED_T,LAST_ITERATION
!
! Local variables.
!
	INTEGER EQION,IONE
	INTEGER I,J,K,L,M
	INTEGER JJ,KK
	INTEGER VION,VED,VT
	INTEGER EQ_ION_BAL
	INTEGER EQ_NUM_CONV
	REAL*8 T1,T2
	REAL*8 TMP_VEC_ED(1)
	REAL*8 TMP_VEC_COOL(1)
	PARAMETER (IONE=1)
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	REAL*8 OMEGA_F(N_F,N_F)
	REAL*8 dln_OMEGA_F_dlnT(N_F,N_F)
!
        EQION=N_S+1			!Ion equation
        VION=N_S+1
        VED=SE(ID)%N_IV-1
        VT=SE(ID)%N_IV
	EQ_NUM_CONV=SE(ID)%NUMBER_BAL_EQ
	M=(NUM_BNDS/2)+1
!
	DO I=DST,DEND			!Which depth
!
! Compute collisional cross-sections (and their T derivatives)
! We call this routine ND times so the CNM and DCM arrays can be
! smaller (i.e. no ND dimension).
!
! OMEGA_F,dln_OMEGA_dlnT are work arrays only.
! T1 is returned with the total cooling rate. Not used in this routine.
! We use arrays (even though of length 1) so that some F90 compilers 
! don't give an error message because a scaler is passed a vector.
!
	  TMP_VEC_ED(1)=1.0D0		!Electron density
	  TMP_VEC_COOL(1)=0.0D0		!Initialize cooling rate even 
!                                                        though not used here.
          CALL TUNE(1,'SUBCOL')
	  CALL SUBCOL_MULTI_V6(
	1         OMEGA_F,dln_OMEGA_F_dlnT,
	1         CNM,DCNM,
	1         HN_S(1,I),HNST_S(1,I),dlnHNST_S_dlnT(1,I),N_S,
	1         HN_F(1,I),HNST_F_ON_S(1,I),W_F(1,I),FEDGE_F,
	1         A_F,G_F,LEVNAME_F,N_F,
	1         ZION,ID,COL_FILE,OMEGA_GEN,
	1         F_TO_S_MAPPING,TMP_VEC_COOL,T(I),TMP_VEC_ED,IONE,
	1         COMPUTE_BA,FIXED_T,LAST_ITERATION)
          CALL TUNE(2,'SUBCOL')
!
! 
	  DO J=1,N_S			!Which S.E. equation
	    T1=0.0D0
	    T2=0.0D0
	    DO L=1,N_S
	      T1=T1+( HN_S(L,I)*CNM(L,J)-HN_S(J,I)*CNM(J,L) )
	      T2=T2+( HN_S(L,I)*DCNM(L,J)-HN_S(J,I)*DCNM(J,L) )
	    END DO
	    SE(ID)%STEQ(J,I)=SE(ID)%STEQ(J,I)+(T1+(HNST_S(J,I)-HN_S(J,I))*CNM(J,J))*ED(I)
	  END DO
!
	  IF(COMPUTE_BA)THEN
	    DO J=1,N_S			!Which S.E. equation
	      DO K=1,N_S			!Which variable
	        IF(K.EQ.J)THEN
	          T1=0.0D0
	          DO L=1,N_S
		    T1=T1+CNM(J,L)
	          END DO
		  SE(ID)%BA(J,K,M,I)=SE(ID)%BA(J,K,M,I)-T1*ED(I)
	        ELSE
	          SE(ID)%BA(J,K,M,I)=SE(ID)%BA(J,K,M,I)+ED(I)*CNM(K,J)
	        END IF
	      END DO
!
	      T1=0.0D0
	      T2=0.0D0
	      DO L=1,N_S
	        T1=T1+( HN_S(L,I)*CNM(L,J)-HN_S(J,I)*CNM(J,L) )
	        T2=T2+( HN_S(L,I)*DCNM(L,J)-HN_S(J,I)*DCNM(J,L) )
	      END DO
!
	      SE(ID)%BA(J,VION,M,I)=SE(ID)%BA(J,VION,M,I) + (HNST_S(J,I)*CNM(J,J))*ED(I)/DI_S(I)
	      SE(ID)%BA(J,VED,M,I) =SE(ID)%BA(J,VED,M,I) + T1+CNM(J,J)*(2*HNST_S(J,I)-HN_S(J,I))
	      SE(ID)%BA(J,VT,M,I)  =SE(ID)%BA(J,VT,M,I) +
	1                              ED(I)*( T2+(HNST_S(J,I)-HN_S(J,I))*DCNM(J,J)+
	1                              CNM(J,J)*HNST_S(J,I)*dlnHNST_S_dlnT(J,I)/T(I) )
	    END DO
	  END IF
!
! EQION is the ion equation 
!
	  T1=0.0D0
	  T2=0.0D0
	  DO J=1,N_S
	    T1=T1+(HNST_S(J,I)-HN_S(J,I))*CNM(J,J)
	    T2=T2+HNST_S(J,I)*CNM(J,J)
	  END DO
	  SE(ID)%STEQ(EQION,I)=SE(ID)%STEQ(EQION,I)-T1*ED(I)
!
	  IF(COMPUTE_BA)THEN
	    DO J=1,N_S
	      SE(ID)%BA(EQION,J,M,I)=SE(ID)%BA(EQION,J,M,I)+CNM(J,J)*ED(I)
	      SE(ID)%BA(EQION,VT,M,I)=SE(ID)%BA(EQION,VT,M,I) -
	1          ED(I)*(  (HNST_S(J,I)-HN_S(J,I))*DCNM(J,J) +
	1          CNM(J,J)*HNST_S(J,I)*dlnHNST_S_dlnT(J,I)/T(I)  )
	    END DO
	    SE(ID)%BA(EQION,VION,M,I)=SE(ID)%BA(EQION,VION,M,I) - T2*ED(I)/DI_S(I)
	    SE(ID)%BA(EQION,VED,M,I)=SE(ID)%BA(EQION,VED,M,I) - T1 - T2
	  END IF
!
!
          IF(COMPUTE_BA)THEN
	    DO J=1,N_S
              JJ=EQGS+J-1;
	      BA_ED(JJ,M,I)=BA_ED(JJ,M,I)+(ZION-1.0D0)
	      SE(ID)%BA(EQ_NUM_CONV,J,M,I)=SE(ID)%BA(EQ_NUM_CONV,J,M,I)+1.0D0
	    END DO
	  END IF
!
	  T1=0.0D0
	  DO L=1,N_S
	    T1=T1+HN_S(L,I)
	  END DO
	  SE(ID)%STEQ(EQ_NUM_CONV,I)=SE(ID)%STEQ(EQ_NUM_CONV,I)+T1
	  STEQ_ED(I)=STEQ_ED(I)+(ZION-1.0D0)*T1
!
! We only include DI in the population and charge conservation equations
! if the higher ionization species is not present. Necessary to do this as
! DI is the ground state of the next species. We also correct the
! conservation equation for POP if the higher ionization species is not
! present.
!
! Note the charge on DI is ZION. Since we will ALWAYS declare the
! ground state to be an important variable, JJ should never be zero.
!
	  IF(.NOT. NEXT_PRES)THEN
	    SE(ID)%STEQ(EQ_NUM_CONV,I)=SE(ID)%STEQ(EQ_NUM_CONV,I)+DI_S(I)-POP(I)
	    STEQ_ED(I)=STEQ_ED(I)+DI_S(I)*ZION
	  END IF
	  IF(COMPUTE_BA .AND. .NOT. NEXT_PRES)THEN
	    JJ=EQ_NUM_CONV
	    SE(ID)%BA(JJ,VION,M,I)=SE(ID)%BA(JJ,VION,M,I)+1.0D0
	    JJ=EQGS+N_S
	    BA_ED(JJ,M,I)=BA_ED(JJ,M,I)+ZION
	  END IF
!
	END DO		!Loop over depth (I)
!
	RETURN
	END
