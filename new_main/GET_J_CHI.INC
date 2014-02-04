C
C Routine to read in OPACITY and MEAN INTENSITY data for calculation
C with impurity code. This routine should replace
C                                             OPACITIES and
C                                             COMP_JCONT.
C
C Created 20-Oct-1989.
C
C Compute EXP(-hv/kT) and zero CHI, and ETA.
C
	T1=-HDKT*FL
	DO I=1,ND
	  EMHNUKT(I)=EXP(T1/T(I))
	END DO
C
C Compute continuum intensity incident from the core assuming a TSTAR
C blackbody.
C
	IC=TWOHCSQ*(FL**3)/(EXP(HDKT*FL/TSTAR)-1.0D0)
C
C Compute DBB and DDBBDT for diffusion approximation. DBB=dB/dR
C and DDBBDT= dB/dTR .
C
	T1=HDKT*FL/T(ND)
	T2=1.0D0-EMHNUKT(ND)
	DBB=TWOHCSQ*( FL**3 )*T1*DTDR/T(ND)*EMHNUKT(ND)/(T2**2)
	DDBBDT=DBB*(T1*(1.0D0+EMHNUKT(ND))/T2-2.0D0)/T(ND)
C
	READ(35,REC=ACCESS_F)(CHI(I),I=1,ND),(ESEC(I),I=1,ND),
	1                (ETA(I),I=1,ND),(RJ(I),I=1,ND),T1
	IF( .NOT. EQUAL(T1,FL,1.0D-14) )THEN
	  WRITE(LUER,*)'Error - incorrect reading of J and CHI '
	  WRITE(LUER,*)'Frequency is ',FL,'Old Frequency is ',T1
	  WRITE(LUER,*)'Error occurred in '//SECTION
	  STOP
	END IF
C
C Update record for next frequency
C
	ACCESS_F=ACCESS_F+1
C
	DO I=1,ND
	  ZETA(I)=ETA(I)/CHI(I)
	  THETA(I)=ESEC(I)/CHI(I)
	END DO
C
