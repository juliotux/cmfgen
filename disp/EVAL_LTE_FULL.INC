C
C This include file evaluates the LTE populations with respect to
C the ground state for 
C
C  (a) The full atoms.
C  (B) The Super-Level model atom.
C
C Altered 26-Oct-1995
C Created 1-Oct-1990 (Section taken from CMFGEN)>
C
C The final stements set the population of the ground state of the next
C ionizations stages. These must be set since they are used in determining
C opacities and photoinization terms.
C
	DO ID=1,NUM_IONS
	  IF(ATM(ID)%XzV_PRES)THEN
	    CALL LTEPOP_WLD_V2(
	1          ATM(ID)%XzVLTE_F,    ATM(ID)%LOG_XZVLTE_F,     ATM(ID)%W_XzV_F,
	1          ATM(ID)%EDGEXzV_F,   ATM(ID)%GXzV_F,           ATM(ID)%ZXzV,
	1          ATM(ID)%GIONXzV_F,   ATM(ID)%NXzV_F,           ATM(ID)%DXzV_F,
	1          ED,T,ND)
	    IF(.NOT. ATM(ID+1)%XzV_PRES)THEN
	      DO I=1,ND
	        ATM(ID+1)%XzV_F(1,I)=ATM(ID)%DXzV_F(I)
	        ATM(ID+1)%XzVLTE_F(1,I)=ATM(ID)%DXzV_F(I)
	      END DO
	    END IF
	  END IF
	END DO
