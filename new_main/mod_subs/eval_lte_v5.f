!
! This subroutine (which replaces EVAL_LTE_V4.INC) evaluates the LTE populations with respect to
! the ground state for 
!
!  (a) The full atoms.
!  (B) The Super-Level model atom.
!
	SUBROUTINE EVAL_LTE_V5(DO_LEV_DISSOLUTION,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 18-Oct-2011 : Changed to V5 as call LTEPOP_WLD_V2 etc.
!                       V4 reverts to old routine.
! Altered 05-Apr-2011 : Now call LTEPOP_WLD_V2 (instead of V1) and
!                                LTE_POP_SL_V2 (instead of V1) (28-Nov-2010).
!                         Changes done to give a wider dynamic rangs in LTE populations. 
! Created 19-Dec-2004
!
	INTEGER ND
	LOGICAL DO_LEV_DISSOLUTION
!
! Local variables.
!
	INTEGER I
	INTEGER ID
!
! Revise vector constants for evaluating the level dissolution. These
! constants are the same for all species. These are stored in a common block,
! and are required by SUP_TO_FULL and LTE_POP_WLD.
!
	CALL COMP_LEV_DIS_BLK(ED,POPION,T,DO_LEV_DISSOLUTION,ND)
!
! The final statements set the population of the ground state of the next
! ionizations stages. These must be set since they are used in determining
! opacities and photoinization terms.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL LTEPOP_WLD_V2(
	1          ATM(ID)%XzVLTE_F,  ATM(ID)%LOG_XzVLTE_F, ATM(ID)%W_XzV_F,
	1          ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,  ATM(ID)%ZXzV,
	1          ATM(ID)%GIONXzV_F, ATM(ID)%NXzV_F,
	1          ATM(ID)%DXzV_F,    ED,T,ND)
	    CALL LTE_POP_SL_V2(
	1          ATM(ID)%XzVLTE,         ATM(ID)%LOG_XzVLTE,    ATM(ID)%dlnXzVLTE_dlnT,
	1          ATM(ID)%NXzV,           ATM(ID)%XzVLTE_F,      ATM(ID)%LOG_XzVLTE_F,
	1          ATM(ID)%XzVLTE_F_ON_S,  ATM(ID)%EDGEXzV_F,     ATM(ID)%F_TO_S_XzV,
	1          ATM(ID)%NXzV_F,         ATM(ID)%XzV_PRES, T,ND)
	    IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
	      DO I=1,ND
	        ATM(ID+1)%XzV(1,I)=ATM(ID)%DXzV_F(I)			!True if not present.
	        ATM(ID+1)%XzVLTE(1,I)=ATM(ID)%DXzV_F(I)			!True if not present.
	        ATM(ID+1)%LOG_XzVLTE(1,I)=LOG(ATM(ID)%DXzV_F(I))	!True if not present.
	        ATM(ID+1)%dlnXzVLTE_dlnT(1,I)=0.0D0
	      END DO
	    END IF
	  END IF
	END DO
!
	RETURN
	END
