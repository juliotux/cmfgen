!
! Quick routine designed to set the CMGFEN populations to their LTE values.
! This includeds the electron density.
!
! At presnet code probably only works when hydrogen id dominant, and ionized.
! Needs to be updated to use a modified form of DET_LTE_ED where an accurate
! ED and accuracy is passed.
!
	SUBROUTINE SET_POPS_TO_PURE_LTE(POPS,NT,ND)
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	USE CONTROL_VARIABLE_MOD, ONLY : DO_LEV_DISSOLUTION
!
! Created 17-Feb-2016
!
	INTEGER NT
	INTEGER ND
	REAL*8 POPS(NT,ND)
!
	REAL*8 REVISED_ED(ND)
	REAL*8 VEC_SUM(ND)
	REAL*8 Z_POP(NT)
	REAL*8 T1
!
	INTEGER ISPEC
	INTEGER ID
	INTEGER I,K,L
!
	Z_POP=0.0D0
	DO ID=1,NUM_IONS-1
          CALL SET_Z_POP(Z_POP, ATM(ID)%ZXzV, ATM(ID)%EQXzV,ATM(ID)%NXzV, NT, ATM(ID)%XzV_PRES)
        END DO
!
! This get ths ion populations.
!
	DO ID=1,NUM_IONS-1
	  CALL POPTOION(POPS, ATM(ID)%XzV, ATM(ID)%DXzV,ED,T,
	1    ATM(ID)%EQXzV, ATM(ID)%NXzV, NT, ND, ATM(ID)%XzV_PRES)
	END DO
	CALL SUP_TO_FULL_V4(POPS,Z_POP,DO_LEV_DISSOLUTION,ND,NT)
!
	DO K=1,5
	  REVISED_ED=0.0D0
	  DO ISPEC=1,NUM_SPECIES
	    IF(SPECIES_PRES(ISPEC))THEN
	      ID=SPECIES_END_ID(ISPEC)-1
	      VEC_SUM=ATM(ID)%DXzV_F
	      DO ID=SPECIES_END_ID(ISPEC)-1,SPECIES_BEG_ID(ISPEC),-1
	        CALL LTEPOP_WLD_V2(
	1            ATM(ID)%XzVLTE_F,  ATM(ID)%LOG_XzVLTE_F, ATM(ID)%W_XzV_F,
	1            ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,  ATM(ID)%ZXzV,
	1            ATM(ID)%GIONXzV_F, ATM(ID)%NXzV_F,
	1            ATM(ID)%DXzV_F,    ED,T,ND)
	        DO L=1,ND
	          VEC_SUM(L)=VEC_SUM(L)+SUM(ATM(ID)%XzVLTE_F(:,L))
	        END DO
	        IF(ID .NE. SPECIES_BEG_ID(ISPEC))ATM(ID-1)%DXzV_F(:)=ATM(ID)%XzVLTE_F(1,:)
	        ATM(ID)%XzV_F=ATM(ID)%XzVLTE_F
	      END DO
!
! Scale LTE pops to match population constraint.
!
	      VEC_SUM(:)=POP_SPECIES(:,ISPEC)/VEC_SUM(:)
	      DO ID=SPECIES_END_ID(ISPEC)-1,SPECIES_BEG_ID(ISPEC),-1
	        DO L=1,ND
	          ATM(ID)%XzV_F(:,L)=ATM(ID)%XzV_F(:,L)*VEC_SUM(L)
	          REVISED_ED(L)=REVISED_ED(L)+(ATM(ID)%ZXzV-1.0D0)*SUM(ATM(ID)%XzV_F(:,L))
	        END DO
	        IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)ATM(ID)%DXzV_F=ATM(ID)%DXzV_F*VEC_SUM
	        IF(ID .NE. SPECIES_BEG_ID(ISPEC))ATM(ID-1)%DXzV_F(:)=ATM(ID)%XzV_F(1,:)
	      END DO
	      ID=SPECIES_END_ID(ISPEC)-1
	      IF(ID .EQ. 55555)THEN
	        WRITE(6,*)ID,ISPEC
	        WRITE(6,*)VEC_SUM(1),VEC_SUM(ND)
	        WRITE(6,*)REVISED_ED(1),REVISED_ED(ND)
	        WRITE(6,*)ATM(ID)%DXzV_F(1),ATM(ID)%DXzV_F(ND)
	        WRITE(6,*)ATM(ID)%ZXzV
	        WRITE(6,*)ATM(ID)%DXzV_F
	        WRITE(6,*)REVISED_ED
	      END IF
	      DO L=1,ND
	        REVISED_ED(L)=REVISED_ED(L)+ATM(ID)%ZXzV*ATM(ID)%DXzV_F(L)
	      END DO
!	      REVISED_ED=REVISED_ED+ATM(ID)%ZXzV*ATM(ID)%DXzV_F
	    END IF
	  END DO
	  WRITE(6,*)'Doing iteration in SET_POPS_TO_PURE_LTE',K
	  T1=MAXVAL( ABS( (ED-REVISED_ED)/ED ) )
	  WRITE(6,*)'Convergence accuracy for the electron number (in %) is',100.0D0*T1
	  ED=REVISED_ED
	END DO
!
        DO ID=NUM_IONS-1,1,-1
	   CALL FULL_TO_SUP(
	1      ATM(ID)%XzV,   ATM(ID)%NXzV,       ATM(ID)%DXzV,      ATM(ID)%XzV_PRES,
	1      ATM(ID)%XzV_F, ATM(ID)%F_TO_S_XzV, ATM(ID)%NXzV_F,    ATM(ID)%DXzV_F,
	1      ATM(ID+1)%XzV, ATM(ID+1)%NXzV,     ATM(ID+1)%XzV_PRES, ND)
	END DO
!
! Store all quantities in POPS array. This is done here (rather than
! after final iteration) as it enable POPION to be readily computed.
!
	DO ID=1,NUM_IONS-1
	  CALL IONTOPOP(POPS, ATM(ID)%XzV, ATM(ID)%DXzV, ED,T,
	1         ATM(ID)%EQXzV, ATM(ID)%NXzV, NT,ND,
	1         ATM(ID)%XzV_PRES)
	END DO
!
	RETURN
	END
