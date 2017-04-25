!
! Self contained routine to compute the luminosity from radioactive decays.
! Thes routine was designed specifically to be inserted into do_species_decays_v2.f so
! that energy from each chain could be checked.
!
	SUBROUTINE GET_DECAY_LUM_V1(RAD_DECAY_LUM,dE_RAD_DECAY,ND)
	USE MOD_CMFGEN
	IMPLICIT NONE
!
! Created: 30-May-2016
!
	INTEGER ND
	REAL*8 dE_RAD_DECAY(ND)
	REAL*8 RAD_DECAY_LUM
!
	REAL*8 WORK(ND)
	INTEGER I
!
! Note:  3.2845D-03=(4pi*1.0D+30)/Lsun
! 10^30 arises since R is in units of 10^10 cm.
!
	DO I=1,ND
	  WORK(I)=dE_RAD_DECAY(I)*R(I)*R(I)*CLUMP_FAC(I)
	END DO
	CALL LUM_FROM_ETA_V2(WORK,R,'LINMON',ND)
	RAD_DECAY_LUM=3.2845D-03*SUM(WORK)
!
	RETURN
	END
