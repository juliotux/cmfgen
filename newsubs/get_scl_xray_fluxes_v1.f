!
! Routine to return X-ray EMISSIVITIES for a set of NFREQ frequencies,
! and for 2 different shock temperatures (in units of 10^4 K).
! The returned emissivities have units of 10^{-10} ergs/cm^3/s/Hz/steradian.
! At present, we assume that the X-ray emissivity is independent of density.
!
	SUBROUTINE GET_SCL_XRAY_FLUXES_V1(CUR_CONT_FREQ,EMISS1,EMISS2,
	1              CONT_FREQ,NFREQ,ML,
	1              VSMOOTH,SECTION)
	USE MOD_XRAY_FLUXES
	IMPLICIT NONE
!
	INTEGER NFREQ		!Number of frequencies
	INTEGER ML
	REAL*8 VSMOOTH
!
! NB: CONT_FREQ, NFREQ, and ML are ignored if SECTION =/ 'CONTINUUM'
!
	REAL*8 CUR_CONT_FREQ
	REAL*8 CONT_FREQ(NFREQ)
!
!The X-ray emissivity has units of ergs/cm^3/s/steradian/Hz, multiplied
! by a factor of 10^10 so as in program units.
!
	REAL*8 EMISS1
	REAL*8 EMISS2
	CHARACTER*(*) SECTION
!
! Local variables:
!
	REAL*8 PREV_CONT_FREQ
	REAL*8 NEXT_CONT_FREQ
!
	REAL*8 NU_LF_BIN,NU_HF_BIN
	REAL*8 HIGH_F
	REAL*8 LOW_F
	REAL*8 T1
	INTEGER J
	INTEGER LOC_LF,LOC_HF
!
! ******************************************************************
! ******************************************************************
!
!       Rebin the data onto the new frequency grid.
!
! For the continuum calculation, we bin the data such that the integral of 
! the fluxes on the new frequency grid agrees with that of the original data 
! on its grid. We assume the intgeral is done by the trapazoidal rule.
!
	IF(SECTION .EQ. 'CONTINUUM')THEN
!
! NB: CONT_FREQ is a vector containing the frequency at which the continuum
!     opacities and emissivities are evaluated for each frequency. Designed
!     for CMFGEN, where we don't evaluate the continuum opacities/emissivities
!     at every frequency. Evaluation at every frequecny, is a trivial
!     case of the above, and is thus automatically handled.
!
	  J=ML
	  DO WHILE(CONT_FREQ(J) .EQ. CUR_CONT_FREQ .AND. J .GT. 1)
	    J=J-1
	  END DO
	  PREV_CONT_FREQ=CONT_FREQ(J)
!
	  J=ML
	  DO WHILE(CONT_FREQ(J) .EQ. CUR_CONT_FREQ .AND. J .LT. NFREQ)
	    J=J+1
	  END DO
	  NEXT_CONT_FREQ=CONT_FREQ(J)
!	  
	  HIGH_F=0.5D0*(PREV_CONT_FREQ+CUR_CONT_FREQ)
	  LOW_F=0.5D0*(NEXT_CONT_FREQ+CUR_CONT_FREQ)
	  IF(HIGH_F .LT. LOW_F)THEN
	    T1=HIGH_F; HIGH_F=LOW_F; LOW_F=T1
	  END IF
	ELSE
!
! For the fluxes at a single frequency, as required by the Sobolev 
! approximation (for example) we simply average the fluxes over a
! bin VSMOOTH km/s broad.
!
	  HIGH_F=CUR_CONT_FREQ*(1.0D0+0.5D0*VSMOOTH/2.99D+05)
	  LOW_F=CUR_CONT_FREQ*(1.0D0-0.5D0*VSMOOTH/2.99D+05)
	END IF
!
! Now compute the fluxes.
!
	EMISS1=0.0D0
	EMISS2=0.0D0
	LOC_LF=NINT( (LOW_F-BIN_MIN)/BIN_SIZE ) +1
	LOC_HF=NINT( (HIGH_F-BIN_MIN)/BIN_SIZE ) +1
	IF(LOC_LF .GE. 1 .AND. LOC_HF .LE. N_BINS)THEN
	  IF(LOC_HF .EQ. LOC_LF)THEN
	    EMISS1=X_EMISS1(LOC_LF)
	    EMISS2=X_EMISS2(LOC_LF)
	  ELSE
!
! We integrate the fluxes over the half interval centered on the
! current frequency. We then normalize this integral by the half 
! interval.
!
	    NU_LF_BIN=BIN_MIN+(LOC_LF-1)*BIN_SIZE
	    NU_HF_BIN=BIN_MIN+(LOC_HF-1)*BIN_SIZE
	    EMISS1=X_EMISS1(LOC_LF)*(0.5D0*BIN_SIZE+NU_LF_BIN-LOW_F)
	    EMISS1=EMISS1+X_EMISS1(LOC_HF)*
	1                   (0.5D0*BIN_SIZE+HIGH_F-NU_HF_BIN)
	    EMISS2=X_EMISS2(LOC_LF)*(0.5D0*BIN_SIZE+NU_Lf_BIN-LOW_F)
	    EMISS2=EMISS2+X_EMISS2(LOC_HF)*
	1                   (0.5D0*BIN_SIZE+HIGH_F-NU_HF_BIN)
	    DO J=LOC_LF+1,LOC_HF-1
	      EMISS1=EMISS1+X_EMISS1(J)*BIN_SIZE
	      EMISS2=EMISS2+X_EMISS2(J)*BIN_SIZE
	    END DO
	    EMISS1=EMISS1/(HIGH_F-LOW_F)
	    EMISS2=EMISS2/(HIGH_F-LOW_F)
	  END IF
	END IF
!
	RETURN
	END
