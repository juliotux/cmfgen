!
! Routine to convolve the mean intensity J at a given depth, with the
! electron scattering redistribution function. E.S.R.F is approximated
! with a ONE parameter fit.
!
	SUBROUTINE CNVLV_ES_ONE_PAR_V2(NU,RJ,J_ES,T_ELEC,T_OUT,NCF)
	IMPLICIT NONE
!
! Modified 28-Apr-1999: Changged to V2
!                       Photon checks inserted.
!
	INTEGER NCF
	REAL*8 NU(NCF)
	REAL*8 RJ(NCF)
	REAL*8 J_ES(NCF)
	REAL*8 T_ELEC
	INTEGER T_OUT
!
	REAL*8 A(NCF)
	REAL*8 B(NCF)
	REAL*8 C(NCF)
	REAL*8 D(NCF)
	REAL*8 BETA
	REAL*8 T1
	REAL*8 D1,D2,DH
	INTEGER ML
	INTEGER, PARAMETER :: IONE=1
!
	REAL*8 OLD_FLUX
	REAL*8 NEW_FLUX
	REAL*8 OLD_N_PHOT
	REAL*8 NEW_N_PHOT
!
	BETA=1.84E-03*SQRT(T_ELEC)
	T1=0.5D0*BETA*BETA
	C(1)=0.0D0
	A(1)=0.0D0
	B(1)=-1.0D0
	D(1)=RJ(1)
	DO ML=2,NCF-1
	  D1=LOG(NU(ML-1)/NU(ML))
	  D2=LOG(NU(ML)/NU(ML+1))
	  DH=0.5*(D1+D2)
	  A(ML)=-T1/D1/DH
	  B(ML)=-1.0D0
	  C(ML)=-T1/D2/DH
	  D(ML)=RJ(ML)
	END DO
	A(NCF)=0.0
	C(NCF)=0.0
	B(NCF)=-1.0D0
	D(NCF)=RJ(NCF)
!
	DO ML=1,NCF
	 T1=4.7994*NU(ML)/T_ELEC 
	 D(ML)=D(ML)/T1/T1
!	 IF(T1 .LT. 1)D(ML)=D(ML)/T1/T1
	END DO
	CALL THOMAS_RH(A,B,C,D,NCF,IONE)
	DO ML=1,NCF
	 T1=4.7994*NU(ML)/T_ELEC 
	 D(ML)=D(ML)*T1*T1
	END DO
!
	J_ES(:)=D(:)
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

