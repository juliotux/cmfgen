C
C General routine to write out the cooling rates. Routine must be
C called for each species separately.
C
C The variable FIRST indicates that this is the first time the routine
C has been called for the current value of counter.
C
C The variable LAST indicates that this is the last time the routine
C will be called for the current value of counter.
C
	SUBROUTINE WRCOOLGEN_V2(BFCRC2,FFC2,COLC2,DIECOOL,X_COOL,NT_COOL,
	1                      C2_PRES,NC2,DESC,NETRR,TOTRR,
	1                      COUNTER,ND,LU)
	IMPLICIT NONE
C
C Altered 25-Sep-2011 - Changed to _V2 from _NT (Chendong routine).
C                       Changed from originalroutine by addition of NT_COOL to call.
C Altered 26-Jun-1996 - Call to GEN_ASCI_OPEN installed.
C Altered 04-Jun-1996 - FORMFEED defined as a parameter.
C Altered 26-AUg-1991 - FILE used in OPEN statement. Variable format
C                       expressions (i.e. <L> deleted.) FORMFEED used instead
C                       of 1H1.
C                       Comma inseted between A and '' in WRITE statement for
C                         FFC2. Crashed in F90 but not F77.
C Created 24-Mar-1989 - Based on WRC2TOCIVCOOL
C
C
	CHARACTER*(*)DESC
	LOGICAL C2_PRES
	INTEGER NC2,ND,COUNTER,LU,IOS
	INTEGER, PARAMETER :: IZERO=0
	REAL*8 BFCRC2(NC2,ND),FFC2(ND),COLC2(ND),DIECOOL(ND),X_COOL(ND),NT_COOL(ND)
	REAL*8 R(ND),T(ND),ED(ND),TOTRR(ND),NETRR(ND)
C
	INTEGER ICHRLEN,ERROR_LU,LUER
	EXTERNAL ICHRLEN,ERROR_LU
C
C Local variables.
C
	INTEGER MS,MF,I,J,L
	CHARACTER*2, PARAMETER :: FORMFEED=' '//CHAR(12)
C
	IF(.NOT. C2_PRES)RETURN
C
C Set limits for this write.
C
	MS=(COUNTER-1)*10+1	
	MF=COUNTER*10
	IF(MF .GT. ND)MF=ND
C
C Write bound-free cooling (and dielectronic)
C
	L=ICHRLEN(DESC)
	WRITE( LU,
	1  '(/,3X,A,'' Bound-Free Cooling [ergs/cm**3/s]'')')DESC(1:L)
	DO I=1,NC2
	  WRITE(LU,999)(BFCRC2(I,J),J=MS,MF)
	END DO
C
	IF(DIECOOL(1) .NE. 0 .OR. DIECOOL(ND) .NE. 0)THEN
	  WRITE( LU,'(/,3X,A,'' Dielectronic and Implicit '//
	1         'Recombination Cooling [ergs/cm**3/s]'')' )DESC(1:L)
	  WRITE(LU,999)(DIECOOL(J),J=MS,MF)
	END IF
C
C Output Collisional cooling rate.
C
	WRITE( LU,'(/,3X,A,'' Collisional Cooling'')' )DESC(1:L)
	WRITE(LU,999)(COLC2(J),J=MS,MF)
C
C Output free- free cooling rate.
C
	WRITE( LU,'(/,3X,A,'' (ion) Free-Free Cooling'')' )DESC(1:L)
	WRITE(LU,999)(FFC2(J),J=MS,MF)
C
C Output Net K-shell (Auger ionzation) cooling rate (Normally negative
C and hence a heating term).
C
	IF(X_COOL(1) .NE. 0 .OR. X_COOL(ND) .NE. 0)THEN
	  WRITE( LU,'(/,3X,A,'' K-shell cooling'')' )DESC(1:L)
	  WRITE(LU,999)(X_COOL(J),J=MS,MF)
	END IF
C
C Output Non-thermal cooling
C
	IF(SUM(NT_COOL) .NE. 0.0D0)THEN
	  WRITE( LU,'(/,3X,A,'' Non-thermal cooling'')' )DESC(1:L)
	  WRITE(LU,999)(NT_COOL(J),J=MS,MF)
	END IF
C
	DO J=MS,MF
	  DO I=1,NC2
	    NETRR(J)=NETRR(J)+BFCRC2(I,J)
	    TOTRR(J)=TOTRR(J)+ABS(BFCRC2(I,J))
	  END DO
	  NETRR(J)=NETRR(J)+FFC2(J)+COLC2(J)+DIECOOL(J)+X_COOL(J)+NT_COOL(J)
	  TOTRR(J)=TOTRR(J)+ABS(FFC2(J))+ABS(COLC2(J))+
	1            ABS(DIECOOL(J))+ABS(X_COOL(J))+ABS(NT_COOL(J))
	END DO
C
	RETURN
999	FORMAT(1X,1P,10E12.4)
C
C This initializes arrays, and writes the R, T and ED arrays out.
C
	 ENTRY FSTCOOL(R,T,ED,NETRR,TOTRR,COUNTER,ND,LU)
C
	  MS=(COUNTER-1)*10+1	
	  MF=COUNTER*10
	  IF(MF .GT. ND)MF=ND
C
	  IF(COUNTER .EQ. 1)THEN
	    CALL GEN_ASCI_OPEN(LU,'GENCOOL','UNKNOWN',' ',' ',IZERO,IOS)
	    IF(IOS .NE. 0)THEN
	      LUER=ERROR_LU()
	      WRITE(LUER,*)'Error opening GENCOOL file'
	      WRITE(LUER,*)'IOSTAT=',IOS
	      RETURN
	    END IF
	  END IF
C
	  DO I=MS,MF
	    NETRR(I)=0.0D0
	    TOTRR(I)=0.0D0
	  END DO
	  IF(MS .NE. 1)WRITE(LU,'(A)')FORMFEED
	  WRITE(LU,'(/,3X,''Radius [1.0E+10cm] '')')
	  WRITE(LU,999)(R(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Temperature [1.0E+4K] '')')
	  WRITE(LU,999)(T(J),J=MS,MF)
	  WRITE(LU,'(/,3X,''Electron Density'')')
	  WRITE(LU,999)(ED(J),J=MS,MF)
C
	RETURN
C
C Write out the summed arrays.
C
	ENTRY ENDCOOL(NETRR,TOTRR,COUNTER,ND,LU)
C
C
	  MS=(COUNTER-1)*10+1	
	  MF=COUNTER*10
	  IF(MF .GT. ND)MF=ND
C
C Times *200 since we have added cooling and heating rates (absolute values).
C
	  DO J=MS,MF
	    TOTRR(J)=NETRR(J)/TOTRR(J)*200.0
	  END DO
C
	  WRITE(LU,'(//3X,''Net Cooling Rate [ergs/cm**3/sec]'') ')
	  WRITE(LU,999)(NETRR(J),J=MS,MF)
	  WRITE(LU,'(/3X,''Net Cooling as percentage of total.'') ')
	  WRITE(LU,999)(TOTRR(J),J=MS,MF)
	  IF( COUNTER .EQ. (ND+9)/10 )THEN
	    CLOSE(UNIT=LU)
	  END IF
C
	RETURN
C
	END
