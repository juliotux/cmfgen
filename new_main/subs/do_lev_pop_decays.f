	SUBROUTINE DO_LEV_POP_DECAYS(OLD_POPS,ND,NT)
	USE MOD_CMFGEN
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NT
        REAL*8 OLD_POPS(NT,ND)
!
! Local variables.
!
	REAL*8 RATIO(ND)
!
	INTEGER ID		!Ion index
	INTEGER IP		!Parent index
	INTEGER ISPEC		!Species index
	INTEGER I		!Level index
	INTEGER K		!Depth index
!
	INTEGER IVAR		!Level in OLD_POPS
	INTEGER ION_IVAR	!Ion level in OLD_POPS
!
! Below we assume that all levels change in the same way,
! which is proportional to the change in the parent species.
!
	IF(.NOT. DO_RAD_DECAYS)RETURN
	DO IP=1,NUM_PARENTS
	  ISPEC=PAR(IP)%ISPEC
	  IF(SPECIES_PRES(ISPEC) .AND. PAR(IP)%DECAY_CHAIN_AVAILABLE)THEN
	    RATIO=PAR(IP)%OLD_POP_DECAY/PAR(IP)%OLD_POP
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)
	      IVAR=ATM(ID)%EQXzV
	      DO K=1,ND
	        DO I=1,ATM(ID)%NXzV
	          OLD_POPS(IVAR+I-1,K)=OLD_POPS(IVAR+I-1,K)*RATIO(K)
	        END DO
	      END DO
	    END DO
	  END IF
	END DO
!
	RETURN
	END
