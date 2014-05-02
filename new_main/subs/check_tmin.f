!
! Subroutine designed to facilitate (or check for inconsistencies) in
! the use of CMFGEN.
!
! Routine indicataes to the user whether ions can be deleted from the
! model, or whether it may be necesary to add additional high
! ionization stages.
!
	SUBROUTINE CHECK_TMIN()
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered: 23-May-2009 : Bug fix
! Created: 14-Apr-2009
!
	REAL*8, PARAMETER :: HDKT=4.7994145D0
	REAL*8 T1,T2,TMIN
	INTEGER ISPEC
	INTEGER ID
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
! The following table gives IP in eV. Grabbed off the web.
! As its only used to guide the user, the IP's need not be very accurate.
!
! Determine ionzation fractions.
!
	LUER=ERROR_LU()
	TMIN=0.0D0
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    T1=HDKT*ATM(ID)%EDGEXzV_F(1)
	    T2=T1/LOG(1.0D0+300.0D0)
	    TMIN=MAX(TMIN,T2)
	  END DO
	END DO
	WRITE(6,*)' '
	WRITE(6,'(A)')' If T is too small may get expoential overflow when computing LTE populations'
	WRITE(6,'(A,F6.3)')' The recommended minimum value of T (T_MIN) in units of 10^4 K is:',TMIN
	WRITE(6,*)' '
!
	RETURN
	END
