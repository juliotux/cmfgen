C
C Routine to write out the adiabatic cooling. A flag indicates whether
C adiabatic cooling has been included in the T calculation.
C Rate is still output, but in different format.
C
	SUBROUTINE WR_AD_COOL(AD_COOL_V,AD_COOL_DT,NETCR,TOTCR,AD_INCL,
	1                       COUNTER,ND,LU)
	IMPLICIT NONE
C
C Created 26-Jul-1994
C
	INTEGER COUNTER,ND,LU
	REAL*8 AD_COOL_V(ND)
	REAL*8 AD_COOL_DT(ND)
	REAL*8 NETCR(ND)
	REAL*8 TOTCR(ND)
	LOGICAL AD_INCL
C
	INTEGER MS,MF,I
C
	MS=(COUNTER-1)*10+1	
	MF=COUNTER*10
	IF(MF .GT. ND)MF=ND
C
	IF(AD_INCL)THEN
	  DO I=MS,MF
	    NETCR(I)=NETCR(I)+(AD_COOL_V(I)+AD_COOL_DT(I))
	    TOTCR(I)=TOTCR(I)+ABS(AD_COOL_V(I))+ABS(AD_COOL_DT(I))
	  END DO
	  WRITE(LU,'(/3X,A)')'Adiabatic cooling rate (V term)'//
	1                        ' [ergs/cm**3/sec]'
	  WRITE(LU,999)(AD_COOL_V(I),I=MS,MF)
	  WRITE(LU,'(/3X,A)')'Adiabatic cooling rate (dTdR term)'//
	1                         ' [ergs/cm**3/sec]'
	  WRITE(LU,999)(AD_COOL_DT(I),I=MS,MF)
	ELSE
	  WRITE(LU,'(/3X,A)')'Ratio of adiabatic to total cooling'//
	1          ' rate (V term)[Not Incl.]'
	  WRITE(LU,999)(2.0D0*AD_COOL_V(I)/TOTCR(I),I=MS,MF)
	  WRITE(LU,'(/3X,A)')'Ratio of adiabatic to total cooling'//
	1          ' rate (dTdR term)[Not Incl.]'
	  WRITE(LU,999)(2.0D0*AD_COOL_DT(I)/TOTCR(I),I=MS,MF)
	END IF
C
999	FORMAT(1X,1P,10E12.4)
	RETURN
	END
