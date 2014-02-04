	SUBROUTINE GET_LELEC(LELEC,XKT,NKT,N_ELEC)
	IMPLICIT NONE
!
	INTEGER NKT
	REAL*8 LELEC(NKT)
	REAL*8 XKT(NKT)
	REAL*8 N_ELEC
!
	real*8, parameter :: ELECTRON_VOLT=1.60217733D-12  ! erg
	real*8, parameter :: PLANCKS_CONSTANT=6.626075D-27 ! erg sec
	real*8, parameter :: SPEED_OF_LIGHT=2.99792458D+10 ! cm / sec
	real*8, parameter :: PI = 3.141592653589793238462643D0
	real*8, parameter :: ELECTRON_MASS=9.109389D-28    !gm
	real*8, parameter :: ELECTRON_CHARGE = 4.803206814D-10 ! esu
	real*8, parameter :: a0 = 0.529189379D-8    ! Bohr radius in cm
!
	INTEGER I
	INTEGER IKT
	REAL*8 T1,T2,XV,XE_CGS
	REAL*8 XI_E
	REAL*8 PLASMA_FREQ
	REAL*8 GAMMA_EULER
!
! plasma_freq is in 1/s, i.e. cgs
! xi_e is put in erg !!! [e^2] has dimensions of energy x length
! we use xe_cgs to work with cgs units
!
	  plasma_freq = sqrt( 4.0D0 * PI * n_elec * ELECTRON_CHARGE**2 / ELECTRON_MASS ) ! 1/s
	  xi_e = PLANCKS_CONSTANT/2.0D0/PI * plasma_freq                                 ! erg
!
! write(6,'(A25,ES15.5)') 'Plasma Frequency [Hz]: ',plasma_freq
! write(6,'(A25,ES15.5)') 'xi electron [eV]',xi_e/ELECTRON_VOLT
	
	  gamma_euler = 0.577215665D0
	
	  do ikt=1,nkt
	     xe_cgs = xkt(ikt) * ELECTRON_VOLT
!
! Not sure if this velcoity shoud be that of thermal electrons.
! Here I take that of the non-thermal electrons we are treating at xe_cgs
!
	     xv = sqrt(2.0D0*xe_cgs/ELECTRON_MASS) ! in cm/s!!!
	
	     if (xkt(ikt).lt.14.0D0) then
	        t1 = n_elec * 2.D0 * PI * ELECTRON_CHARGE**4 / xe_cgs * &
	             dlog( ELECTRON_MASS * xv**3/ gamma_euler / ELECTRON_CHARGE**2 / plasma_freq  )
	        lelec(ikt) = t1
	     else
	        t2 = n_elec * 2.D0 * PI * ELECTRON_CHARGE**4 / xe_cgs * &
	             dlog(2.0D0 * xe_cgs / xi_e)
	        lelec(ikt) =  t2
	     endif
	  enddo
	  lelec(1:nkt) = lelec(1:nkt) / ELECTRON_VOLT ! now in eV/cm
	
	return
	end subroutine get_lelec
