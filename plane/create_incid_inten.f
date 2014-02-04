!
! Simple program to create an INCIDENT INTENSITY file for use with cmfgen and
! the plane-parallel subroutine:
!                              pp_form_cmf_v2.f
! The incident intensity is outout in table format, and thus arbitrary
! distributions are possible. At present, the blackbody flux is described
! a blackbody distribution.
!
	PROGRAM CREATE_INCID_INTEN
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
	REAL*8, ALLOCATABLE :: MU(:)
	REAL*8, ALLOCATABLE :: DIST(:)
	REAL*8 NU
	REAL*8 BNU
!
	REAL*8 TEFF_STAR
	REAL*8 TEFF_INCID
	REAL*8 REL_LUM
	REAL*8 DILUTION_FACTOR
!
	REAL*8, PARAMETER :: HDKT=4.7994145D0
	REAL*8, PARAMETER :: TWOHCSQ=0.0147452575D0
	REAL*8, PARAMETER :: NU_MAX=100.0D0
	REAL*8, PARAMETER :: NU_MIN=0.0001D0
!
	REAL*8 T1
	REAL*8 FLUX
	REAL*8 X
	REAL*8 PI
	REAL*8 NORMALIZED_ANGLE_INTEGRAL
!
	INTEGER I
	INTEGER NANG
	INTEGER NCF
	INTEGER NPTS_PER_DECADE
!
	CHARACTER*132 FILENAME
!
	NANG=11; NPTS_PER_DECADE=20; REL_LUM=0.1D0
	CALL GEN_IN(TEFF_STAR,'Teff for star (in 10^4K)')
	CALL GEN_IN(TEFF_INCID,'Teff describing incident radiation (in 10^4 K)')
	CALL GEN_IN(REL_LUM,'Relative fluxes (incident/stellar)')
	CALL GEN_IN(NPTS_PER_DECADE,'Number of points per decade in frequency')
	CALL GEN_IN(NANG,'Number of points in angle')
!
	ALLOCATE (MU(NANG))
	ALLOCATE (DIST(NANG))
!
! Set angles used to describe the incident radiation field.
!
	MU(1)=0.0D0
	MU(NANG)=1.0D0
	DO I=2,NANG-1
	  MU(I)=(I-1.0D0)/(NANG-1.0D0)
	END DO
!
! Set the distribution of flux as a function of angle.
!
	DO I=1,NANG
	  DIST(I)=1.0D0-4.0D0*(MU(I)-0.5D0)**2
	END DO
!
! Estimate the integral I(mu).mu dmu. For isotropic radiation, the integral
! is 0.5. We thus normalize the integral by 0.5.
!
	T1=0.0D0
	DO I=1,NANG-1
	  T1=T1+0.5D0*(MU(I+1)-MU(I))*(MU(I)*DIST(I)+MU(I+1)*DIST(I+1))
	END DO
	NORMALIZED_ANGLE_INTEGRAL=2.0D0*T1
        WRITE(6,*)'Normalized angle weighting integral is',NORMALIZED_ANGLE_INTEGRAL
!
! Determine the dilution factor to give the requested relative luminosity fir
! the requested incident Teff (fluxes) and the stellar model Teff.
!
	DILUTION_FACTOR=(REL_LUM/NORMALIZED_ANGLE_INTEGRAL)*(TEFF_STAR/TEFF_INCID)**4
!
! Determine number of frequency points.
!
	NCF=LOG10(NU_MAX/NU_MIN)*NPTS_PER_DECADE+1
!
! INCID_INTEN is the name of the output file to contain the incident fluxes as a
! function of angle and frequency.
!
	FILENAME='INCID_INTEN'
	CALL GEN_IN(FILENAME,'Output file')
	OPEN(UNIT=15,STATUS='UNKNOWN',FILE=TRIM(FILENAME))
!
	WRITE(15,'(A)')'! '
	WRITE(15,'(A)')'! The effective temperature of the incident radiation is defined as the'
	WRITE(15,'(A)')'! temperature of the blackbody that gives the equivalent radiation flux.'
	WRITE(15,'(A)')'! Teff(star) is used to set the relative luminosity and dilution factor.'
	WRITE(15,'(A)')'! Thus the dilution factor (W) is defined by:'
	WRITE(15,'(A)')'! '
	WRITE(15,'(A)')'!         Rel. Lum = W . NAI . (Teff/Teff[*])*0.25'
	WRITE(15,'(A)')'! '
	WRITE(15,'(A)')'! where NIA is the normalized angle integral. The actual dilution factor'
	WRITE(15,'(A)')'! used to scale the blackbody is output just prior to the angle grid.'
	WRITE(15,'(A)')'! This can be easily changed without recreating a new incident flux file.'
	WRITE(15,'(A)')'! Frequency is in units of 10^15 Hz, and I is in CGS units'
	WRITE(15,'(A)')'! '
	WRITE(15,'(A,ES14.6)')'! Teff(K) of incident radiation      :',TEFF_INCID*1.0D+04
	WRITE(15,'(A,ES14.6)')'! Normalized angle weighting integral:',NORMALIZED_ANGLE_INTEGRAL
	WRITE(15,'(A,ES14.6)')'! '
	WRITE(15,'(A,ES14.6)')'! Assumed Teff(K) of star is:         ',TEFF_STAR*1.0D+04
	WRITE(15,'(A,ES14.6)')'! Relative luminosity:                ',REL_LUM
	WRITE(15,'(A,ES14.6)')'! Dilution factor is:                 ',DILUTION_FACTOR
	WRITE(15,*)' '
!
	WRITE(15,'(A,T30,A)')'   10-Apr-2006','!Format date'
	WRITE(15,'(I14,T30,A)')NCF,'!Number of continuum frequencies'
	WRITE(15,'(I14,T30,A)')NANG,'!Number of angles'
	WRITE(15,'(ES14.6,T30,A)')DILUTION_FACTOR,'!Dilution factor used for scaling'
	WRITE(15,*)' '
	WRITE(15,*)'Angle grid'
	WRITE(15,*)(MU(I),I=1,NANG)
	WRITE(15,*)' '
	WRITE(15,*)'Intensity variation with angle'
	WRITE(15,*)(DIST(I),I=1,NANG)
!
	WRITE(15,*)' '
	T1=EXP(LOG(10.0D0)/NPTS_PER_DECADE)
	NU=NU_MAX*T1
	FLUX=0.0D0
	DO I=1,NCF
	  NU=NU/T1
	  X=EXP(-HDKT*NU/TEFF_INCID)
	  BNU=TWOHCSQ*NU*NU*NU*X/(1.0D0-X)
	  FLUX=FLUX+BNU*(NU*T1-NU/T1)/2.0D0
	  WRITE(15,'(X,2ES14.6)')NU,BNU
	END DO
	PI=4.0D0*ATAN(1.0D0)
	T1=(FLUX*PI*1.0D+15/5.6705D-05)**(0.25D0)
	WRITE(6,'(A,ES14.6,A)')'Equivalent effective temperature is',T1,' K'
!
	STOP
	END
