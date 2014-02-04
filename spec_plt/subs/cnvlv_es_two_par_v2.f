!
! Routine to convolve the mean intensity J at a given depth, with the
! electron scattering redistribution function. E.S.R.F is approximated
! with a TWO parameter fit.
!
	SUBROUTINE CNVLV_ES_TWO_PAR_V2(NU,RJ,J_ES,T_ELEC,T_OUT,NCF)
	IMPLICIT NONE
!
	INTEGER NCF
	REAL*8 NU(NCF)
	REAL*8 RJ(NCF)
	REAL*8 J_ES(NCF)
	REAL*8 T_ELEC
	INTEGER T_OUT
!
	REAL*8 A(NCF)
	REAL*8 H(NCF)
	REAL*8 C(NCF)
	REAL*8 D(NCF)
	REAL*8 A_STORE(NCF)
	REAL*8 C_STORE(NCF)
	REAL*8 BETA
	REAL*8 T1
	REAL*8 D1,D2,DH
	INTEGER ML,L
	INTEGER, PARAMETER :: IONE=1
!
	REAL*8 OLD_FLUX
	REAL*8 NEW_FLUX
	REAL*8 OLD_N_PHOT
	REAL*8 NEW_N_PHOT
!
! Parameters for fit to Electrons Scattering redistribution function
! (dipole form). From Rybicki and Hummer (A&A, 290,553)
!
	INTEGER, PARAMETER :: NCOEF=2
	REAL*8 ACOEF(2),BCOEF(2)
	DATA ACOEF/1.690703717290D0,-0.690703717290D0/
	DATA BCOEF/1.614249968779D0,2.154326524957D0/
!
! Compute those parts of the TRIDIAGONAL vectors which are independent of
! depth, and the fitting parameters.
!
	A_STORE(1)=0
	C_STORE(1)=-2.0D0/( LOG(NU(1)/NU(2)) )**2
	DO ML=2,NCF-1
	  D1=LOG(NU(ML-1)/NU(ML))
	  D2=LOG(NU(ML)/NU(ML+1))
	  DH=0.5*(D1+D2)
	  A_STORE(ML)=-1.0D0/D1/DH
	  C_STORE(ML)=-1.0D0/D2/DH
	END DO
	A_STORE(NCF)=-2.0D0/( LOG(NU(NCF-1)/NU(NCF)) )**2
	C_STORE(NCF)=0.0D0
!
! Compute the triadiagonal quantities for performing the convolution, and 
! perform the convolution. Due to the depth dependence of BETA, the vectors 
! are depth dependent. We convolve both J and the Planck Function.
!
	BETA=1.84E-03*SQRT(T_ELEC)
	J_ES(:)=0.0D0
	DO L=1,NCOEF
	  T1=BETA*BETA/BCOEF(L)/BCOEF(L)
	  A(:)=T1*A_STORE(:)			!Over frequency
	  H(:)=-1.0D0
	  C(:)=T1*C_STORE(:)
	  D(:)=RJ(:)
	  CALL THOMAS_RH(A,H,C,D,NCF,IONE)
	  J_ES(:)=J_ES(:)+ACOEF(L)*D(:)
	END DO
!
! Compute integrals as a function of depth to check flux conservation.
!
	OLD_FLUX=0.0D0
	NEW_FLUX=0.0D0
	DO ML=2,NCF-1
	  OLD_FLUX=OLD_FLUX+(NU(ML)-NU(ML+1))*(RJ(ML)+RJ(ML+1))
	  NEW_FLUX=NEW_FLUX+(NU(ML)-NU(ML+1))*(J_ES(ML)+J_ES(ML+1))
	END DO
	T1=200.0D0*(OLD_FLUX-NEW_FLUX)/(OLD_FLUX+NEW_FLUX)
	WRITE(T_OUT,'(A,1X,1P,E13.5)')'  %Flux error:',T1
!
! Now check photon number conservation
!
	OLD_N_PHOT=0.0D0
	NEW_N_PHOT=0.0D0
	DO ML=2,NCF-1
	  T1=LOG(NU(ML)/NU(ML+1))
	  OLD_N_PHOT=OLD_N_PHOT+T1*(RJ(ML)+RJ(ML+1))
	  NEW_N_PHOT=NEW_N_PHOT+T1*(J_ES(ML)+J_ES(ML+1))
	END DO
	T1=200.0D0*(OLD_N_PHOT-NEW_N_PHOT)/(OLD_N_PHOT+NEW_N_PHOT)
	WRITE(T_OUT,'(A,1X,1P,E13.5)')'  %Photon error:',T1
!
	RETURN
	END
