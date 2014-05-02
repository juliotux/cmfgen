!
! Subroutine designed to correct species and isotope populations for 
! radioactive decays. Upon entry:
!     (1) radioactive data must be set
!     (2) ISO%OLD_POP must be set.
!
! Upon exit ISO%OLD_POP_DECAY is set.
!           PAR_ISO           is set.
!
! This routine currently only works for 1 and 2-step decay processes.
!
! e.g., Ni56 --> Co56 -- Fe56
!
	SUBROUTINE DO_SPECIES_DECAYS(DELTA_T,ND) 
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
! Altered: 02-Sep-2013 : Minor bug fix -- added excess energy (for those Ni
!                            that decay all the way to 56Fe in a single time step).
! Altered: 23-Nov-2012 : Added 1 step decay process
!
	INTEGER ND
	REAL*8 VEC1(ND)
	REAL*8 VEC2(ND)
	REAL*8 VEC3(ND)
	REAL*8 DELTA_T
!
	INTEGER IN
	INTEGER IP
	INTEGER JN
!
	INTEGER IS
	INTEGER JS
	INTEGER LS
!
	INTEGER K,L,M
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
	DO IS=1,NUM_ISOTOPES
	  ISO(IS)%OLD_POP_DECAY=ISO(IS)%OLD_POP
	END DO
	RADIOACTIVE_DECAY_ENERGY=0.0D0
!
	WRITE(6,*)'   In DO_SPECIES_DECAYS, Delta T=',DELTA_T
	DO IN=1,NUM_DECAY_PATHS
!
! Do decay chains with a single decay route.
!
	  IS=NUC(IN)%LNK_TO_ISO
	  JS=NUC(IN)%DAUGHTER_LNK_TO_ISO
	  IF(NUC(IN)%SEQUENCE ==  'E' .AND. IS*JS .NE. 0)THEN
	    VEC1=EXP(-NUC(IN)%DECAY_CONST*DELTA_T)
!	    WRITE(6,*)'Found one step nuclear reaction chain:',IN
	    ISO(IS)%OLD_POP_DECAY=ISO(IS)%OLD_POP*VEC1
	    ISO(JS)%OLD_POP_DECAY=ISO(JS)%OLD_POP + ISO(IS)%OLD_POP*(1.0D0-VEC1)
	    RADIOACTIVE_DECAY_ENERGY=RADIOACTIVE_DECAY_ENERGY +  
	1            ISO(IS)%OLD_POP*(1.0D0-VEC1)*NUC(IN)%ENERGY_PER_DECAY
!	    WRITE(6,*)RADIOACTIVE_DECAY_ENERGY
!
! Do two step sequences.
!
	  ELSE IF(NUC(IN)%SEQUENCE ==  'F' .AND. IS .NE. 0)THEN
	    VEC1=EXP(-NUC(IN)%DECAY_CONST*DELTA_T)
	    DO JN=1,NUM_DECAY_PATHS
	      IF(NUC(JN)%SEQUENCE ==  'S' .AND. NUC(JN)%BARYON_NUMBER .EQ.  NUC(IN)%BARYON_NUMBER .AND.
	1        NUC(JN)%SPECIES .EQ. NUC(IN)%DAUGHTER)THEN
!	        WRITE(6,*)'Found nuclear reaction chain:',IN,JN
	        JS=NUC(JN)%LNK_TO_ISO
	        IF(JS .EQ. 0)EXIT
	        VEC2=EXP(-NUC(JN)%DECAY_CONST*DELTA_T)
	        VEC3=NUC(IN)%DECAY_CONST*( EXP(-NUC(JN)%DECAY_CONST*DELTA_T) -
	1                   EXP(-NUC(IN)%DECAY_CONST*DELTA_T) )/(NUC(IN)%DECAY_CONST-NUC(JN)%DECAY_CONST)
	        LS=NUC(JN)%DAUGHTER_LNK_TO_ISO
		ISO(IS)%OLD_POP_DECAY=ISO(IS)%OLD_POP*VEC1
		ISO(JS)%OLD_POP_DECAY=ISO(JS)%OLD_POP*VEC2 + ISO(IS)%OLD_POP*VEC3
	        ISO(LS)%OLD_POP_DECAY=ISO(LS)%OLD_POP +
	1                             ISO(JS)%OLD_POP*(1.0D0-VEC2) + 
	1                             ISO(IS)%OLD_POP*(1.0D0-VEC1-VEC3)
	        RADIOACTIVE_DECAY_ENERGY=RADIOACTIVE_DECAY_ENERGY +  
	1            ISO(IS)%OLD_POP*(1.0D0-VEC1)*NUC(IN)%ENERGY_PER_DECAY +
	1            ISO(JS)%OLD_POP*(1.0D0-VEC2)*NUC(JN)%ENERGY_PER_DECAY +
	1            ISO(IS)%OLD_POP*(1.0D0-VEC1-VEC3)*NUC(JN)%ENERGY_PER_DECAY
	        EXIT
	      END IF
	    END DO
	  END IF
	END DO
	IF(DELTA_T .NE. 0.0D0)RADIOACTIVE_DECAY_ENERGY=RADIOACTIVE_DECAY_ENERGY/DELTA_T		!ergs/sec
!
	DO IP=1,NUM_PARENTS
	  PAR(IP)%OLD_POP=0.0D0 
	  PAR(IP)%OLD_POP_DECAY=0.0D0 
	END DO
	DO IS=1,NUM_ISOTOPES
	  IP=ISO(IS)%LNK_TO_PAR
	  PAR(IP)%OLD_POP=PAR(IP)%OLD_POP+ISO(IS)%OLD_POP
	  PAR(IP)%OLD_POP_DECAY=PAR(IP)%OLD_POP_DECAY+ISO(IS)%OLD_POP_DECAY
	END DO
!
	RETURN
	END
