!
! This include file evaluates the LTE populations with respect to
! the ground state for 
!
!  (a) The full atoms.
!  (B) The Super-Level model atom.
!
! Altered 05-Mar-2015 : Bug fix. Added ATM(ID+1)%LOG_XzVLTE(1,I)=LOG(ATM(ID)%DXzV_F(I)) - RSS
! Altered 17-Dec-2011 : Changed to LTEPOP_WLD_V2 and LTE_POP_SL_V2. These contain LOG_XzVLTE etc
!                           and are designed to faciltate lower wind temperatures.
! Altered 01-Mar-1999 - Variables now accessed by pointer reference.
! Altered 26-Oct-1995
! Created 1-Oct-1990 (Section taken from CMFGEN)>
!
! The final statements set the population of the ground state of the next
! ionizations stages. These must be set since they are used in determining
! opacities and photoinization terms.
!
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL LTEPOP_WLD_V2(
	1          ATM(ID)%XzVLTE_F,  ATM(ID)%LOG_XzVLTE_F,  ATM(ID)%W_XzV_F,
	1          ATM(ID)%EDGEXzV_F, ATM(ID)%GXzV_F,        ATM(ID)%ZXzV,
	1          ATM(ID)%GIONXzV_F, ATM(ID)%NXzV_F,
	1          ATM(ID)%DXzV_F,    ED,T,ND)
	    CALL LTE_POP_SL_V2(
	1          ATM(ID)%XzVLTE,          ATM(ID)%LOG_XzVLTE,   ATM(ID)%dlnXzVLTE_dlnT,
	1          ATM(ID)%NXzV,            ATM(ID)%XzVLTE_F,     ATM(ID)%LOG_XzVLTE_F,
	1          ATM(ID)%XzVLTE_F_ON_S,   ATM(ID)%EDGEXzV_F,    ATM(ID)%F_TO_S_XzV,
	1          ATM(ID)%NXzV_F,          ATM(ID)%XzV_PRES, T,ND)
	    IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
	      DO I=1,ND
	        ATM(ID+1)%XzV(1,I)=ATM(ID)%DXzV_F(I)	 	        !True if not present.
	        ATM(ID+1)%XzVLTE(1,I)=ATM(ID)%DXzV_F(I)  		!True if not present.
	        ATM(ID+1)%LOG_XzVLTE(1,I)=LOG(ATM(ID)%DXzV_F(I))        !True if not present (RSS fix)
	        ATM(ID+1)%dlnXzVLTE_dlnT(1,I)=0.0D0
	      END DO
	    END IF
	  END IF
	END DO
