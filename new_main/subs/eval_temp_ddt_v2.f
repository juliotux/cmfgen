!
! Routine to adjust the radiative equilibrium equation for the exchange of
! energy with the gas, and the work done on the gas. This code is designed for
! SN models, and assumes a Lagrangian formulation. V is the comoving variable. 
!
!  (1) Increments STEQ if adiabatic cooling is to be allowed for
!  (2) Increments the variation matrix [BA] if it is being computed, and
!       if adiabatic cooling is being included.
!  (3) Evaluates the adiabatic cooling terms, splitting them into 2 terms ---
!       the specific energy term, and the density term. This is done for diagnostic
!       purposes.
!
	SUBROUTINE EVAL_TEMP_DDT_V2(WORK,AD_CR_V,AD_CR_DT,
	1               POPS,AVE_ENERGY,HDKT,COMPUTE_BA,INCL_ADIABATIC,
	1               TIME_SEQ_NO,DIAG_INDX,NUM_BNDS,NT,ND)
	USE MOD_CMFGEN
 	USE STEQ_DATA_MOD
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
! Altered 01-Sep-2016: TIME_SEQ_NO changed from integer to real.
! Altered 22-Nov-2011 : Changed call to GET_POPS_AT_PREV_TIME_STEP from V4 to V5
!                          (POPS was added to call).
! Altered 11-Nov-2009 : Fixed bug in BA calculation --- was skipping depth ND.
!                          No longer calculate ION_EN. 
! Altered 12-Feb-2008 : VOL_EXP_FAC used. Now valid for all expansion laws.
!                          Returned populations are corrected for decays and adiabatic expansion.
! Altered 22-Mar-2007 : Call GET_POPS_AT_PREV_TIME_STEP_V4
! Altered 04-Feb-2007 : Bug fixed in d(STEQ_T)/dNe
! Altered 23-Jun-2006 : R, V removed from call to GET_POPS_AT_PREV_TIME_STEP_V2 since
!                         passed by module MOD_CMFGEN.
! Altered 05-Jun-2006 : Bug fixed with evaluation of terms for electron cooling
!                           equation.
! Created 13-Dec-2005 : Based on EVAL_ADIABATIC_V3
!
	INTEGER NT
	INTEGER ND
	INTEGER DIAG_INDX
	INTEGER NUM_BNDS
!
! Output:
!
	REAL*8 AD_CR_V(ND)
	REAL*8 AD_CR_DT(ND)
!
! Input:
!
	REAL*8 POPS(NT,ND)
	REAL*8 AVE_ENERGY(NT)
	REAL*8 WORK(ND)
	REAL*8 HDKT
	REAL*8 TIME_SEQ_NO
!
! Local vectors.
!
	REAL*8 EK_VEC(ND)
	REAL*8 EI_VEC(ND)
	REAL*8 P_VEC(ND)
	REAL*8 GAMMA(ND)
	REAL*8 INT_EN(ND)
!
	REAL*8 OLD_POPS(NT,ND)
	REAL*8 OLD_R(ND)
	REAL*8 OLD_T(ND)
	REAL*8 OLD_ED(ND)
	REAL*8 OLD_GAMMA(ND)
	REAL*8 OLD_POP_ATOM(ND)
	REAL*8 OLD_INT_EN(ND)
!
	REAL*8 TOT_ENERGY(NT)
!
! Local variables.
!
	LOGICAL, PARAMETER :: L_TRUE=.TRUE.
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL COMPUTE_BA,INCL_ADIABATIC
!
	INTEGER ERROR_LU
	REAL*8 BOLTZMANN_CONSTANT,FUN_PI,STEFAN_BOLTZ,SPEED_OF_LIGHT
	EXTERNAL BOLTZMANN_CONSTANT,FUN_PI,ERROR_LU,STEFAN_BOLTZ,SPEED_OF_LIGHT
!
	REAL*8 SCALE
	REAL*8 T1,T2,T3,T4,PI
	REAL*8 DELTA_T_SECS
	INTEGER I,J,K,L
	INTEGER LUER
	INTEGER ISPEC
	INTEGER ID
	INTEGER LU
	LOGICAL WRITE_CHK
!
! A full linearization is now obsolete, but check to make sure.
!
	LUER=ERROR_LU()
	IF(NUM_BNDS .EQ. ND)THEN
	  WRITE(LUER,*)'Error --- EVAL_TEMOP_DDT_V1 can''t handle a full linearization'
	  STOP
	END IF
!
! Compute the total excitation energy of each level.
!
	TOT_ENERGY(1:NT)=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  T1=0.0D0
	  T2=0.0D0
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    T2=T2+AVE_ENERGY(ATM(ID)%EQXzV)
	    DO I=1,ATM(ID)%NXzV
	      J=ATM(ID)%EQXzV+I-1
	      TOT_ENERGY(J)=(AVE_ENERGY(ATM(ID)%EQXzV)-AVE_ENERGY(J))+T1
	    END DO
	    J=ATM(ID)%EQXzV
	    T1=T1+AVE_ENERGY(J)			!Adding on ionization energy
	  END DO
	  ID=SPECIES_END_ID(ISPEC)-1
	  IF(ID .GT. 0)THEN
	    J=ATM(ID)%EQXzV
	    TOT_ENERGY(J+ATM(ID)%NXzV)=T1
	  END IF
	END DO
!
! Get the populations at the previous time step. These are put onto the same
! V grid as the current model. The returned populations are CORRECTED for 
! advection, radioactive decays, and are normalized.
!
	LU=7
	CALL GET_POPS_AT_PREV_TIME_STEP_V5(POPS,OLD_POPS,OLD_R,
	1      L_TRUE,L_TRUE,L_TRUE,TIME_SEQ_NO,ND,NT,LU)
	OLD_ED(:)=OLD_POPS(NT-1,:)
	OLD_T(:)=OLD_POPS(NT,:)
	OLD_POP_ATOM=0.0D0
	DO I=1,ND
	  DO J=1,NT-2
	    OLD_POP_ATOM(I)=OLD_POP_ATOM(I)+OLD_POPS(J,I)
	  END DO
	END DO
!
	DO ISPEC=1,NUM_SPECIES
	  DO J=1,ND
	    T1=0.0D0
	    T2=0.0D0
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV
	        K=ATM(ID)%EQXzV+I-1
	        T1=T1+POPS(K,J)
	        T2=T2+OLD_POPS(K,J)
	      END DO
	    END DO
	    IF(SPECIES_BEG_ID(ISPEC) .GT. 0)THEN
	      ID=SPECIES_END_ID(ISPEC)-1
	      I=ATM(ID)%EQXzV+ATM(ID)%NXzV
	      T1=T1+POPS(I,J)
	      T2=T2+OLD_POPS(I,J)
	      WRITE(155,'(I5,I5,2ES16.8)')ISPEC,J,T1,T2
	    END IF
	  END DO
	END DO
!
! Compute time step. The factor of 10^5 arises because R is in units of 10^10 cm, and
! V is in units of km/s.
!
	DELTA_T_SECS=1.0D+05*(R(ND)-OLD_R(ND))/V(ND)
!
! Compute the mean energy per atom. At first it is units of 10^15Hz.
!
	INT_EN(:)=0.0D0
	OLD_INT_EN(:)=0.0D0
	DO I=1,ND
	  DO J=1,NT-2
	     INT_EN(I)=INT_EN(I)+POPS(J,I)*TOT_ENERGY(J)
	     OLD_INT_EN(I)=OLD_INT_EN(I)+OLD_POPS(J,I)*TOT_ENERGY(J)
	     IF(I .EQ. 1)THEN
	       T1=POPS(J,I)*TOT_ENERGY(J)/POP_ATOM(I)
	       T2=OLD_POPS(J,I)*TOT_ENERGY(J)/OLD_POP_ATOM(I)
	       WRITE(220,'(I5,5ES14.6)')J,POPS(J,I),OLD_POPS(J,I),TOT_ENERGY(J),T1,T2
	     END IF
	  END DO
	END DO
	INT_EN=HDKT*INT_EN/POP_ATOM
	OLD_INT_EN=HDKT*OLD_INT_EN/OLD_POP_ATOM
!
! We now compute constants for each of the 4 terms. These make
! it simpler and cleaner for the evaluation of the linearization.
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the adiabatic cooling rate by 4pi.
! Also, since chi, and eta are a factor of 10^10 too large, we need to
! scale by 10^10. We need to introduce another factor of 10^4 since T is 
! in units of 10^4K.
!
	PI=FUN_PI()
	SCALE=1.0D+14*BOLTZMANN_CONSTANT()/4.0D0/PI
	DO I=1,ND
	  EK_VEC(I)=1.5D0*SCALE*POP_ATOM(I)/DELTA_T_SECS
	  EI_VEC(I)=SCALE*POP_ATOM(I)/DELTA_T_SECS
          P_VEC(I)=SCALE*(POP_ATOM(I)+ED(I))*T(I)/DELTA_T_SECS
	END DO
!
	DO I=1,ND
	  GAMMA(I)=ED(I)/POP_ATOM(I)
	  OLD_GAMMA(I)=OLD_ED(I)/OLD_POP_ATOM(I)
	END DO
!
	IF(INCL_ADIABATIC)THEN
	  DO I=1,ND
 	    WORK(I)=EK_VEC(I)*( (1.0D0+GAMMA(I))*T(I)- (1.0D0+OLD_GAMMA(I))*OLD_T(I) ) +
	1           EI_VEC(I)*(INT_EN(I)-OLD_INT_EN(I))      +
	1           P_VEC(I)*LOG(VOL_EXP_FAC(I))
	  END DO
	  DO I=1,ND
	    STEQ_T(I)=STEQ_T(I)-WORK(I)
	  END DO
	END IF
!
! Diagonal terms.
!
	IF(INCL_ADIABATIC .AND. COMPUTE_BA)THEN
	  DO I=1,ND
	    L=DIAG_INDX
	    BA_T(NT,L,I)=BA_T(NT,L,I)-EK_VEC(I)*(1+GAMMA(I)) -
	1                    P_VEC(I)*LOG(VOL_EXP_FAC(I))/T(I)
	    BA_T(NT-1,L,I)=BA_T(NT-1,L,I)-T(I)*EK_VEC(I)/POP_ATOM(I) -
	1                    P_VEC(I)*LOG(VOL_EXP_FAC(I))/(POP_ATOM(I)+ED(I))
	    T1=HDKT*EI_VEC(I)/POP_ATOM(I)
	    DO J=1,NT-2
	      BA_T(J,L,I)=BA_T(J,L,I)-T1*TOT_ENERGY(J)
	    END DO
	  END DO	!Loop of depth.
	END IF            !End COMPUTE_BA
!
! Now compute the adiabatic cooling rate (in ergs/cm^3/sec) for diagnostic
! purposes. The rate is output to the GENCOOL file.
!
! The factor of 4.0D-10*PI arises from the fact that A, B etc were computed
! for the radiative equilibrium equation. That equation has units a factor of 10^10/4Pi
! larger (10^10 because of opacity definition, 4Pi as we don't scale J).
!
! We split the adiabatic terms into 2 parts: The energy term, and the
! density (velocity) term. This split was useful for diagnostic purposes,
! but has now been kept for simplicity so that GENCOOL does not need to be changed.
!
	T1=4.0D-10*PI
	DO I=1,ND
	  AD_CR_V(I) =P_VEC(I)*LOG(VOL_EXP_FAC(I))
	  AD_CR_DT(I)=EK_VEC(I)*( (1.0D0+GAMMA(I))*T(I)- (1.0D0+OLD_GAMMA(I))*OLD_T(I) )
	END DO
	AD_CR_V=AD_CR_V*T1
	AD_CR_DT=AD_CR_DT*T1
!
! Output diagnostic information. As OLD_POP_ATOM has already been corrected for
! advection, we don't need to scale it by the change in volume.
!
	WRITE_CHK=.TRUE.
	IF(WRITE_CHK)THEN
	  OPEN(UNIT=7,FILE='ADIABAT_CHK',STATUS='UNKNOWN')
	   WRITE(7,'(A)')'!'
	   WRITE(7,'(A)')'! The terms (DEk/Dt, DEI/Dt, and DP/dt) listed below are included in the CMFGEN'
	   WRITE(7,'(A)')'! radiative equilibrium equation. No additional scaling is needed. The d terms'
	   WRITE(7,'(A)')'! represent the contribution to D... by current values'
	   WRITE(7,'(A)')'!'
	   WRITE(7,'(A,ES14.4)')'DELTA_T_SECS=',DELTA_T_SECS
	   WRITE(7,'(A,ES14.4)')'SCALE_FAC=',SCALE
	   WRITE(7,'(A,7(7X,A7),2X,A)')'Index','      V','      R','  OLD_R','      T','  OLD_T',
	1               '    GAM','OLD_GAM','(Nn/No)(Rn/Ro)**3-1'
	   WRITE(7,'(A,7(7X,A7))')'     ',' Ek(ev)',' Ei(ev)','    dEk','    dEI',
	1                    ' DEk/Dt',' DEI/Dt','?DP/Dt'
	   DO I=1,ND
	    T1=EK_VEC(I)*(1.0D0+GAMMA(I))*T(I)
	    T2=EI_VEC(I)*INT_EN(I)
	    T3=EK_VEC(I)*( (1.0D0+GAMMA(I))*T(I)-(1.0D0+OLD_GAMMA(I))*OLD_T(I) )
	    T4=EI_VEC(I)*(INT_EN(I)-OLD_INT_EN(I))
	    WRITE(7,'(I5,8ES14.5)')I,V(I),R(I),OLD_R(I),T(I),OLD_T(I),GAMMA(I),OLD_GAMMA(I),
	1                    (POP_ATOM(I)/OLD_POP_ATOM(I))-1.0D0
	    WRITE(7,'(5X,7ES14.5)')1.5D0*T(I)*8.6174D-01,INT_EN(I)*8.6174D-01,T1,T2,T3,T4,
	1                    P_VEC(I)*LOG(VOL_EXP_FAC(I))
	   END DO
	   CALL WR2D_V2(POPS,NT,ND,'POPS','*',.TRUE.,7)
	   CALL WR2D_V2(OLD_POPS,NT,ND,'OLD_POPS','*',.TRUE.,7)
	  CLOSE(UNIT=7)
!
	  OPEN(UNIT=7,FILE='ENERGY_COMP',STATUS='UNKNOWN')
	    WRITE(7,'(A)')'!'
	    WRITE(7,'(A)')'! Energy summary (ergs/cm^3)'
	    WRITE(7,'(A,E10.4)')'! Delta t=',DELTA_T_SECS
	    WRITE(7,'(A)')'!'
	    WRITE(7,'(A)')'! NB: Radiative energy density assumes a BB at T(elec).'
	    WRITE(7,'(A)')'!'
	    WRITE(7,'(A,3(14X,A),6(3X,A))')'Depth','R','V','T',
	1                 ' Ek(ejecta)','Ek(thermal)','E(internal)',
	1                 '       Erad','     Enuc/s','       Enuc'
	    T1=4*ACOS(-1.0D0)*1.0D-10
	    DO I=1,ND
	      T2=T1*EK_VEC(I)*(1.0D0+GAMMA(I))*T(I)*DELTA_T_SECS
	      T3=T1*EI_VEC(I)*INT_EN(I)*DELTA_T_SECS
	      T4=4.0D+16*STEFAN_BOLTZ()*(T(I)**4)/SPEED_OF_LIGHT()
	      WRITE(7,'(I5,3ES15.5,6ES14.5)')I,R(I),V(I),T(I),0.5D+10*DENSITY(I)*(V(I)**2),T2,T3,T4,
	1		RADIOACTIVE_DECAY_ENERGY(I),RADIOACTIVE_DECAY_ENERGY(I)*DELTA_T_SECS      
	    END DO
	  CLOSE(UNIT=7)
	END IF
!
	RETURN
	END
