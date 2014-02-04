!
! Routine to write out cooling/heating for charge exchange reactions. 
!
	SUBROUTINE WR_CHG_COOL_V3(NETCR,TOTCR,COUNTER,ND,LU)
	USE CHG_EXCH_MOD_V3
	IMPLICIT NONE
!
! Altered 03-Oct-2011 : Checks on non-zero N_CHG before writing header.
! Created 26-Aug-1998 : based on WR_AD_COOL
!
	INTEGER COUNTER
	INTEGER ND,LU
	REAL*8 NETCR(ND)	!Accumlated net cooling rate
	REAL*8 TOTCR(ND)	!Accumulated sum of absoulte cooling rates.
!                 
	INTEGER I,J
	INTEGER MS,MF
!
	MS=(COUNTER-1)*10+1	
	MF=COUNTER*10
	IF(MF .GT. ND)MF=ND
!
	IF(.NOT. DO_CHG_EXCH)RETURN
	IF(N_CHG .EQ. 0)RETURN
!
	WRITE(LU,'(/3X,A)')'Charge exchange cooling rate [ergs/cm**3/sec]'
	DO J=1,N_CHG
	  IF(CHG_REACTION_AVAILABLE(J))THEN
	    DO I=MS,MF
	      NETCR(I)=NETCR(I)+COOL_CHG(J,I)
	      TOTCR(I)=TOTCR(I)+ABS(COOL_CHG(J,I))
	    END DO
	    WRITE(LU,999)(COOL_CHG(J,I),I=MS,MF)
	  END IF
	END DO
!
999	FORMAT(1X,1P,10E12.4)
	RETURN
	END
