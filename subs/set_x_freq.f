C
C Routine to place bound-free edge frequencies into a vector.
C Option CHECK_CROSS allows bound-free edges with zero cross
C section to be deleted (must be set to FALSE for H, and HeII.
C
	SUBROUTINE SET_X_FREQ(FREQ,NCF,NCF_MAX,ZCORE,ZION,
	1                      NI_PRES,N2_PRES)
	IMPLICIT NONE
C
C Altered 16-Oct-2000 : Bug fix when cross-section is zero.
!                       An EDGE frequency of zero was being returned.
C Altered 26-May-1996 : ERROR_LU installed.
C                       Warning for non CNO species installed.
C Created 18-Jul-1994
C
	INTEGER NCF,NCF_MAX,J
	REAL*8 FREQ(NCF_MAX)
	REAL*8 ZCORE,ZION,CON_FAC
	LOGICAL NI_PRES,N2_PRES
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
C Data variables. C, N, AND O
C
	REAL*8 EDGE_C(6),EDGE_N(7),EDGE_O(8)
C
	DATA EDGE_C/490,392,347,317,296,280/
	DATA EDGE_N/666,552,496,459,432,412,395/
	DATA EDGE_O/870,739,672,627,595,570,550,533/
C
	IF( .NOT. NI_PRES)RETURN
	IF( .NOT. N2_PRES)RETURN
        IF(NCF+1 .GT. NCF_MAX)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in SET_X_FREQ --- NCF_MAX too small'
	  STOP
	END IF
	NCF=NCF+1
	CON_FAC=0.24191				!ev to 10^15 Hz
	J=ZCORE-ZION+1
	IF(ZCORE .EQ. 6)THEN
	  FREQ(NCF)=EDGE_C(J)*CON_FAC
        ELSE IF(ZCORE .EQ. 7)THEN
	  FREQ(NCF)=EDGE_N(J)*CON_FAC
        ELSE IF(ZCORE .EQ. 8)THEN
	  FREQ(NCF)=EDGE_O(J)*CON_FAC
	ELSE
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning: No X-ray edge set in SET_X_FREQ for:'
	  WRITE(LUER,*)'ZCORE=',ZCORE
	  NCF=NCF-1			!Since added one earlier.
        END IF
C
	END
