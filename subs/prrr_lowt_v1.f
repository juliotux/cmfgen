!
! Routine to increment the photoionization and recombination rates
! for an arbitrary ion. The bound-free cooling rate (in ergs/cm**3/s)
! is also computed. 
!
! This rountine is specifically designed for the case when the
! LTE populatons have been set to zero (because they are so
! large). In practice, the term appearing in the opacities/
! and rate equations is HNST*EXP(-hv/kT) and hence remains finite.
!
	SUBROUTINE PRRR_LOWT_V1(RJ,NU,FQW,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Altered 05-APr-2011 - Bug fixed and cleaned.
! Created 01-Feb-2011 - Based on prrr_sl_v6.f 
!
	INTEGER ND			!Number of depth points.
	REAL*8 RJ(ND)
	REAL*8 NU
	REAL*8 FQW
!
	INTEGER I,J,ID,IPR
	INTEGER ION_LEV			!Final (destination) level in ion.
	REAL*8 T1,T2
	REAL*8 LOG_JB_RAT
	REAL*8 H
!
! ConstantS for opacity etC.
!
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF
	REAL*8 HDKT,TWOHCSQ
!
        H=6.6261965D-12                                 !ergs/s (*1.0E+15 due to *nu)
!
! Note that JREC     = Int [ (2hv^3/c^2 +J) exp(-hv/kT)/v dv ]
!           JREC_CR  = Int [ (2hv^3/c^2 +J) exp(-hv/kT)   dv ]
!
        DO ID=1,NUM_IONS
          IF(ATM(ID)%XzV_PRES)THEN
            DO IPR=1,ATM(ID)%N_XzV_PHOT
	      ION_LEV=ATM(ID)%XzV_ION_LEV_ID(IPR)
!
	      LOG_JB_RAT=0.0D0
	      DO J=1,ND
	        T1=( TWOHCSQ*(NU**3)+RJ(J) )*FQW/NU
	        IF(ION_LEV .NE. 1)THEN
	          LOG_JB_RAT=LOG(ATM(ID+1)%XzV(ION_LEV,J)/ATM(ID+1)%XZV(1,J)) +
	1          (ATM(ID+1)%LOG_XzVLTE(1,J)-ATM(ID+1)%LOG_XzVLTE(ION_LEV,J))
	        END IF
	        DO I=1,ATM(ID)%NXzV
	          IF(ATM(ID)%WSXzV(I,J,IPR) .NE. 0 .AND. ATM(ID)%XzVLTE(I,J) .EQ. 0.0D0)THEN
	            T2=T1*EXP(LOG_JB_RAT+ATM(ID)%LOG_XzVLTE(I,J)-HDKT*NU/T(J))
	            ATM(ID)%ARRXzV(I,J)=ATM(ID)%ARRXzV(I,J)+ATM(ID)%WSXzV(I,J,IPR)*T2
	            ATM(ID)%BFCRXzV(I,J)=ATM(ID)%BFCRXzV(I,J) + 
	1                H*T2*( ATM(ID)%WCRXzV(I,J,IPR) + ATM(ID)%WSXzV(I,J,IPR)*NU)
	          END IF
	        END DO		!Level
	      END DO		!depth
!
	    END DO		!Phot route
	  END IF		!species present
	END DO			!ion
!
	RETURN
	END
