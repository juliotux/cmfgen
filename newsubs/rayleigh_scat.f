!
! Subroutine to compute the opacity due to Rayleigh scattering at
! ND depth points. This routine has good accuracy (5%) redward of the
! Lyman lines. Within 1000 km/s of Lyman alpha, cros-section is continuous
! and constant.
!
! References:
!         Gavrila (Phys Rev., 1967, 163, p147)i
!         Hee-Won Lee and Hee Il Kim, 2004, MNRAS, 347, 802
!         Nussbaumer, H., Schmid, H. M., Vogel, M., 1989, A&A, 211, L27 
!
	SUBROUTINE RAYLEIGH_SCAT(RKI,HYD,AHYD,EDGE_HYD,NHYD,FREQ,ND)
	IMPLICIT NONE
!
! Altered 11-Apr-2005 : Major bug fix.
! Created 17-Mar-2004
!
	INTEGER NHYD
	INTEGER ND
!
	REAL*8 AHYD(NHYD,NHYD)		!No longer used
	REAL*8 HYD(NHYD,ND)
	REAL*8 EDGE_HYD(NHYD)		!Ground state value used only.
	REAL*8 RKI(ND)
	REAL*8 FREQ
!
! We include the oscillator strength in the file. That way we don't
! have to worry about whether the l states are split.
!
	REAL*8, SAVE ::  FOSC(20)
	DATA FOSC(2:20)/4.162D-01,7.910D-02,2.899D-02,1.394D-02,7.799D-03,
	1         4.814D-03,3.183D-03,2.216D-03,1.605D-03,1.201D-03,
	1         9.214D-04,7.227D-04,5.774D-04,4.686D-04,3.856D-04,
	1         3.211D-04,2.702D-04,2.296D-04,1.967D-04/
!
	REAL*8 EDGE_FREQ
	REAL*8 T1,T2
	REAL*8 DELF
	INTEGER I,IRES
!
! Use the results of Gavrila (Phys Rev., 1967, 163, p147) below
! the Lyman limit. Below 0.647 we use the fitting formula due to
! Ferland (Cloudy). Blueward of Lyman edge did a simple fit to
! results of Gavrila.
!
	IF(FREQ .GT. EDGE_HYD(1))THEN
	  T1=FREQ/EDGE_HYD(1)
	  RKI(1:ND)=RKI(1:ND)+6.65D-15*HYD(1,1:ND)*(1.0D0+1.66D0/SQRT(T1))
	  RETURN
	ELSE IF(FREQ .LT. 0.647D0*EDGE_HYD(1))THEN
	  T1=(FREQ/EDGE_HYD(1))**2
	  RKI(1:ND)=RKI(1:ND)+6.65D-15*HYD(1,1:ND)*
	1               T1*T1*(1.26D0+5.068D0*T1+708.3D0*(T1**5))
	  RETURN
	END IF
!
! Now sum up over all the resonances.
! See: Nussbaumer, H., Schmid, H. M., Vogel, M., 1989, A&A, 211, L27 
!
	T1=0.0D0
	IRES=1
	DELF=100
	DO I=2,20
	  EDGE_FREQ=EDGE_HYD(1)*(1.0D0-1.0D0/I/I)
	  T2=(EDGE_FREQ/FREQ)**2-1.0D0
	  IF(DELF .GT. ABS(T2))THEN
	    DELF=ABS(T2)
	    IRES=I
	  END IF
	  T2=MAX(ABS(T2),0.0001D0)		!Prevent overflow
	  IF(EDGE_FREQ .LT. FREQ)T2=-T2
	  T1=T1+FOSC(I)/T2
	END DO
!
! The 0.35 gives a match to fitting formula above, and prevents resonances
! from going to zero.
!
	T1=T1*T1+0.35D0
	IF(IRES .GT. 4)IRES=4
	T1=MIN(1.0D4*FOSC(IRES)*FOSC(IRES),T1)
	RKI(1:ND)=RKI(1:ND)+6.65D-15*T1*HYD(1,1:ND)
!
	RETURN
	END
