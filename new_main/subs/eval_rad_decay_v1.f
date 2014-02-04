!
! Routine to adjust the radiative equilibrium equation for energy added to the
! gas by radiactive decay. At present routine simply assumes all the
! energy is deposited locally. 
!
! This routine will need to be modified as we use more sophisticated prescriptions
! for the energy deposition.
!
	SUBROUTINE EVAL_RAD_DECAY_V1(dE_RAD_DECAY,NT,ND)
	USE MOD_CMFGEN
 	USE STEQ_DATA_MOD
	USE NUC_ISO_MOD
 	IMPLICIT NONE
!
! Created 13-Dec-2005 
!
	INTEGER NT
	INTEGER ND
!
! Output: STEQ_T in MOD_CMFGEN is modified.
!
	REAL*8 de_RAD_DECAY(ND)
!
	REAL*8 SCALE
	REAL*8 PI
!
! For historical reasons STEQ contains Int[chi.J - eta]dv. Rather than multiply
! this term everywhere by 4pi, we divide the radiactive heating by 4pi.
! Also, since chi, and eta are a factor of 10^10 too large, we need to
! scale by 10^10. The BA matrix does not need to be altered.
!
	PI=ACOS(-1.0D0)
	SCALE=1.0D+10/4.0D0/PI
	STEQ_T=STEQ_T+SCALE*RADIOACTIVE_DECAY_ENERGY
	de_RAD_DECAY=RADIOACTIVE_DECAY_ENERGY
!
	RETURN
	END
