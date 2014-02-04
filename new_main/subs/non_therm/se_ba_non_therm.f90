	SUBROUTINE SE_BA_NON_THERM(dE_RAD_DECAY,COMPUTE_BA,NT,ND)
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
        USE MOD_NON_THERM
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
! Altered 24-Sep-2011 : No longer multiply  RADIOACTIVE_DECAY_ENERGY by FRAC_ELEC_HEATING.
!                          This bug fixed much earlier in Luc's and CHendong's modeling.
! Created 16-Sep-2010
!
        INTEGER NT
        INTEGER ND
!
! Output: STEQ_T in MOD_CMFGEN is modified.
!
        REAL*8 dE_RAD_DECAY(ND)
	REAL*8 RADIOACTIVE_DECAY_ENERGY_eV(ND)
	REAL*8 LOCAL_ION_HEATING(ND)
	REAL*8 LOCAL_EXC_HEATING(ND)
	LOGICAL COMPUTE_BA
!
	REAL*8 PI
	REAL*8 ELECTRON_VOLT
	INTEGER GET_INDX_DP
	EXTERNAL GET_INDX_DP, ELECTRON_VOLT
!
	REAL*8 RATE
	REAL*8 T1,T2,T3
	REAL*8 ION_EXC_EN
	REAL*8 SCALE
	REAL*8 GUPPER
!
	INTEGER ID
	INTEGER ISPEC
	INTEGER IT
	INTEGER IKT
	INTEGER IKT_ST
	INTEGER DPTH_INDX
	INTEGER SE_ION_LEV
	INTEGER I,J
!
	INTEGER NL
	INTEGER NUP
	INTEGER NL_F
	INTEGER NUP_F
!
	INTEGER, PARAMETER :: LU_TH=8
	INTEGER, PARAMETER :: LU_ER=6
	REAL*8, SAVE :: SCALE_FACTOR=1.0D0
	REAL*8, PARAMETER :: Hz_to_eV=13.60D0/3.2897D0
!
        OPEN(UNIT=LU_TH,FILE='NON_THERM_SPEC_INFO',STATUS='UNKNOWN',POSITION='APPEND')
        CALL SET_LINE_BUFFERING(LU_TH)
!
	LOCAL_ION_HEATING=0.0D0
	LOCAL_EXC_HEATING=0.0D0
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the radiactive heating by 4pi.
! Also, since chi, and eta are a factor of 10^10 too large, we need to
! scale by 10^10. The BA matrix does not need to be altered.
!
        PI=ACOS(-1.0D0)
        SCALE=1.0D+10/4.0D0/PI
        STEQ_T=STEQ_T+SCALE*RADIOACTIVE_DECAY_ENERGY                                      !*FRAC_ELEC_HEATING
        dE_RAD_DECAY=RADIOACTIVE_DECAY_ENERGY
	RADIOACTIVE_DECAY_ENERGY_eV=SCALE_FACTOR*RADIOACTIVE_DECAY_ENERGY/ELECTRON_VOLT()
!
	DO DPTH_INDX=1,ND
!
! Estimate number of non-thermal electrons.
!
	  T1=0.0D0
	  DO IKT=1,NKT
	    T1=T1+YE(IKT,DPTH_INDX)*dXKT(IKT)/SQRT(XKT(IKT))
	  END DO
	  T1=T1*SQRT(0.5D0*9.109389D-28/1.602177D-12)*RADIOACTIVE_DECAY_ENERGY_eV(DPTH_INDX)
	  WRITE(LU_TH,'(A,I3,A,2ES12.4)')' Non-thermal and thermal electron densities at depth',  &
	          DPTH_INDX,' are:',T1,ED(DPTH_INDX)
!
! Compute the ionization contribution
!
	  IF (INCLUDE_NON_THERM_IONIZATION) THEN
!
! Loop over different ionization routes.
!
	    DO IT=1,NUM_THD
!
! Get the total excitation route. We ignore slight differences in IP.
!
	      ID=THD(IT)%LNK_TO_ION
	      ISPEC=THD(IT)%LNK_TO_SPECIES
	      IKT_ST=GET_INDX_DP(THD(IT)%ION_POT,XKT,NKT)
	      RATE = 0.0D0
	      DO IKT=IKT_ST,NKT
	        RATE = RATE + YE(IKT,DPTH_INDX)*dXKT(IKT)*THD(IT)%CROSS_SEC(IKT)
	      END DO
	      RATE=RATE*RADIOACTIVE_DECAY_ENERGY_eV(DPTH_INDX)
!
	      DO J=1,THD(IT)%N_ION_ROUTES
	         IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
                   SE_ION_LEV=ATM(ID)%NXzV+1
	           GUPPER=THD(IT)%SUM_GION
	           ION_EXC_EN=0.0D0
	         ELSE
	           NUP_F=THD(IT)%ION_LEV(J); NUP=ATM(ID+1)%F_TO_S_XzV(NUP_F)
                   SE_ION_LEV=SE(ID)%ION_LEV_TO_EQ_PNT(NUP)
	           GUPPER=ATM(ID+1)%GXzV_F(NUP_F)
                   SE_ION_LEV=ATM(ID)%NXzV+1               !NUP -- assume all to ground state at  present.
	           ION_EXC_EN=0.0D0                        !=ATM(ID+1)%EDGEXzV_F(1)-ATM(ID+1)%EDGEXzV_F(NUP_F)
	         END IF
	         DO I=1,THD(IT)%N_STATES
	           NL_F=THD(IT)%ATOM_STATES(I); NL=ATM(ID)%F_TO_S_XzV(NL_F)
	           T1=RATE*ATM(ID)%XzV_F(NL_F,DPTH_INDX)*GUPPER/THD(IT)%SUM_GION
	           SE(ID)%STEQ(NL,DPTH_INDX)=SE(ID)%STEQ(NL,DPTH_INDX)-T1
	           SE(ID)%STEQ(SE_ION_LEV,DPTH_INDX)=SE(ID)%STEQ(SE_ION_LEV,DPTH_INDX)+T1
	           LOCAL_ION_HEATING(DPTH_INDX)=LOCAL_ION_HEATING(DPTH_INDX)+ T1*(ATM(ID)%EDGEXzV_F(NL_F)+ION_EXC_EN)
	         END DO
	      END DO
!
! NB:  XzV_F= XzV  (XzVLTE_F/XzVLTE)
!
	      IF(COMPUTE_BA)THEN
	        DO J=1,THD(IT)%N_ION_ROUTES
	          IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
                    SE_ION_LEV=ATM(ID)%NXzV+1
	            GUPPER=THD(IT)%SUM_GION
	          ELSE
	            NUP_F=THD(IT)%ION_LEV(J); NUP=ATM(ID+1)%F_TO_S_XzV(NUP)
                    SE_ION_LEV=SE(ID)%ION_LEV_TO_EQ_PNT(NUP)
	            GUPPER=ATM(ID+1)%GXzV_F(NUP_F)
	          END IF
                  SE_ION_LEV=ATM(ID)%NXzV+1               !NUP -- assume all to ground state at  present.
	          DO I=1,THD(IT)%N_STATES
	            NL_F=THD(IT)%ATOM_STATES(I); NL=ATM(ID)%F_TO_S_XzV(NL_F)
	            T1=RATE*GUPPER/THD(IT)%SUM_GION * &
	                         (ATM(ID)%XzVLTE_F(NL_F,DPTH_INDX)/ATM(ID)%XzVLTE(NL,DPTH_INDX))
	            SE(ID)%BA_PAR(NL,NL,DPTH_INDX)=SE(ID)%BA_PAR(NL,NL,DPTH_INDX)-T1
	            SE(ID)%BA_PAR(SE_ION_LEV,NL,DPTH_INDX)=SE(ID)%BA_PAR(SE_ION_LEV,NL,DPTH_INDX)+T1
	          END DO
	        END DO
	      END IF
!
	    END DO			!Ionization species/route
	  END IF			!Include ioizations?
!
	  IF(INCLUDE_NON_THERM_EXCITATION)THEN
	    DO ID=1,NUM_IONS
	      IF(ATM(ID)%XzV_PRES)THEN
	        DO J=2,ATM(ID)%NXzV_F
	          DO I=1,MIN(10,J-1)
	            NL_F=I; NUP_F=J
	            CALL TOTAL_BETHE_RATE(RATE,NL_F,NUP_F,YE,XKT,dXKT,NKT,FAST_BETHE_METHOD,ID,DPTH_INDX,ND)
	            NUP=ATM(ID)%F_TO_S_XzV(NUP_F)
!
! NB:  XzV_F= XzV  (XzVLTE_F/XzVLTE)
!
	            NL=ATM(ID)%F_TO_S_XzV(NL_F)
	            T1=RATE*RADIOACTIVE_DECAY_ENERGY_eV(DPTH_INDX)
	            IF(COMPUTE_BA)THEN
	              T2=T1*(ATM(ID)%XzVLTE_F(NL_F,DPTH_INDX)/ATM(ID)%XzVLTE(NL,DPTH_INDX))
	              SE(ID)%BA_PAR(NL,NL,DPTH_INDX)=SE(ID)%BA_PAR(NL,NL,DPTH_INDX)-T2
	              SE(ID)%BA_PAR(NUP,NL,DPTH_INDX)=SE(ID)%BA_PAR(NUP,NL,DPTH_INDX)+T2
	            END IF
	            T1=T1*ATM(ID)%XzV_F(NL_F,DPTH_INDX)
	            SE(ID)%STEQ(NL,DPTH_INDX)=SE(ID)%STEQ(NL,DPTH_INDX)-T1
	            SE(ID)%STEQ(NUP,DPTH_INDX)=SE(ID)%STEQ(NUP,DPTH_INDX)+T1
	            LOCAL_EXC_HEATING(DPTH_INDX)=LOCAL_EXC_HEATING(DPTH_INDX)+T1* &
	                                (ATM(ID)%EDGEXzV_F(NL_F)-ATM(ID)%EDGEXzV_F(NUP_F))
	          END DO
	        END DO
	      END IF
	    END DO			!Loop over ionization stage
	  END IF			!Include excitations?
!
	END DO 		!loop over depth
!
	WRITE(LU_TH,'(//,A)')'Comparison of heating fractions (SE as evaluated in SE_BA_NON_THERM)'
	WRITE(LU_TH,'(//,A,6(3X,A))')       &
	            '  Fe_nuc(eV)','    Felec','Felec(SE)','     Fion',' Fion(SE)','     Fexc',' Fexc(SE)'
	DO I=1,ND
	  T2=Hz_TO_eV*LOCAL_ION_HEATING(I)/RADIOACTIVE_DECAY_ENERGY_eV(I)
	  T3=Hz_TO_eV*LOCAL_EXC_HEATING(I)/RADIOACTIVE_DECAY_ENERGY_eV(I)
	  T1=(1.0D0-T2-T3)
	  WRITE(LU_TH,'(ES12.4,6(ES12.4))')RADIOACTIVE_DECAY_ENERGY_eV(I),        &
	               FRAC_ELEC_HEATING(I),T1,         &
	               FRAC_ION_HEATING(I),T2,          &
	               FRAC_EXCITE_HEATING(I),T3
	END DO
	CLOSE(LU_TH)
!
	SCALE_FACTOR=MIN(SCALE_FACTOR*10.0D0, 1.0D0)
	WRITE(6,*)'New scale factor is',SCALE_FACTOR
!
	RETURN
	END
