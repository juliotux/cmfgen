!
! Subroutine to determine correspondence between SPECIES in the CHARGE exchange
! reactions, and the corresponding program variables. It is advised 
! (and program checks) that levels involved in charge exchange reactions should
! be distinct super-levels. 
!
! Program assumes (and checks) that each charge reaction involves full terms.
!
! Charge exchange reactions must have the form (and ordering)
!
!     Y(n+) + X([m-1]+)  <--> Y([n-1]+) + X(M+)
!
	SUBROUTINE SET_CHG_LEV_ID_V3(
	1            ID,SPECIES,LEVEL_NAMES,F_TO_S,N_F,N_S,ND,
	1            ZION,EQSPEC,EQION,EQHYD)
	USE CHG_EXCH_MOD_V3
	IMPLICIT NONE
!
! ALetered 09-Oct-1999 : Error reporting inproved.
! Created  20-Aug-1998 : BASED on V1. Very different calls and a change in
!                          philosiphy of how super levels are managed.
!
	INTEGER ID
	INTEGER N_S
	INTEGER N_F
	INTEGER ND
	INTEGER EQSPEC
	INTEGER EQION
	INTEGER EQHYD
!
	INTEGER F_TO_S(N_F)
!
	CHARACTER*(*) SPECIES
	CHARACTER*(*) LEVEL_NAMES(N_F)
!
	REAL*8 ZION
!
	INTEGER I,J,K,L
	INTEGER I_S,I_F
	CHARACTER*(30) LOC_NAME
	LOGICAL LEVEL_SET
!
	IF(.NOT. DO_CHG_EXCH)RETURN
	IF( .NOT. ALLOCATED(LEV_IN_POPS_CHG))THEN
	  ALLOCATE (CHG_REACTION_AVAILABLE(N_CHG))
	  ALLOCATE (ID_ION_CHG(N_CHG,4))
	  ALLOCATE (LEV_IN_POPS_CHG(N_CHG,4))
	  ALLOCATE (LEV_IN_ION_CHG(N_CHG,4))
	  ALLOCATE (Z_CHG(N_CHG,4))
          ALLOCATE (AI_AR_CHG(N_CHG,ND))
          ALLOCATE (dlnAI_AR_CHG_dlnT(N_CHG,ND))
          ALLOCATE (COOL_CHG(N_CHG,ND))
!
! Integer variables
!
	  ID_ION_CHG(:,:)=0
	  LEV_IN_POPS_CHG(:,:)=0
	  LEV_IN_ION_CHG(:,:)=0
!
! Real variables
!
	  Z_CHG(:,:)=0.0D0
          AI_AR_CHG(:,:)=0.0D0
          dlnAI_AR_CHG_dlnT(:,:)=0.0D0
          COOL_CHG(:,:)=0.0D0
	END IF
	INITIALIZE_ARRAYS=.TRUE.
!
! At present, CHARGE reactions cannot refer to individual j states.
!
	DO J=1,N_CHG
	  DO K=1,4
	    IF(INDEX(LEV_NAME_CHG(J,K),'[') .NE. 0 .OR.
	1          INDEX(ALT_LEV_NAME_CHG(J,K),'[') .NE. 0)THEN
	      WRITE(LUER,*)'Error in SET_CHG_EXCH'
	      WRITE(LUER,*)'Level names in CHG reaction file cannot',
	1            'refer to split levels'
	      STOP
	    END IF
	  END DO
	END DO
!
! Determine ionization stage and levels for species involved in charge 
! exchange reactions. Now determine whether the present species is in 
! the CHARGE exchange reaction list.
!
	DO J=1,N_CHG
	  LEVEL_SET=.FALSE.
	  DO K=1,4
	    IF(SPEC_ID_CHG(J,K) .EQ. SPECIES)THEN
	      I_S=0
	      DO I_F=1,N_F
	        LOC_NAME=LEVEL_NAMES(I_F)
	        L=INDEX(LOC_NAME,'[')
	        IF(L .NE. 0)LOC_NAME=LOC_NAME(1:L-1)
	        IF(LOC_NAME .EQ. LEV_NAME_CHG(J,K) .OR.
	1               LOC_NAME .EQ. ALT_LEV_NAME_CHG(J,K))THEN
	          IF(I_S .EQ. 0)THEN
	            ID_ION_CHG(J,K)=ID     
	            I_S=F_TO_S(I_F)
	            LEV_IN_ION_CHG(J,K)=I_S
	            LEV_IN_POPS_CHG(J,K)=EQSPEC+I_S-1		!In POPS
	            Z_CHG(J,K)=ZION-1.0D0
	            LEVEL_SET=.TRUE.
	          END IF              
	          IF(I_S .NE. F_TO_S(I_F))THEN
	            WRITE(LUER,*)'Inconsistent level IDs in SET_CHG_EXCH'
		    WRITE(LUER,*)' Charge exchange reaction:',J
		    WRITE(LUER,*)' Species:',K
	            STOP
	          END IF
	        END IF
	      END DO
!
	      IF( LEVEL_SET)THEN
!
! If the charge reaction involves the final ionization stage of an atomic
! species, we need to set the DATA when we have passed data for the lower
! ionization stage. Because of our conventions, we need only check when
! K=2 or 3.
!
		IF(K .EQ. 2 .AND. EQSPEC+N_S .EQ. EQHYD)THEN
	          Z_CHG(J,4)=ZION
	          ID_ION_CHG(J,4)=ID+1
	          LEV_IN_POPS_CHG(J,4)=EQSPEC+N_S
	          LEV_IN_ION_CHG(J,4)=1
	        END IF
	        IF( K .EQ. 3 .AND. EQSPEC+N_S .EQ. EQHYD )THEN
	          Z_CHG(J,1)=ZION
	          ID_ION_CHG(J,1)=ID+1
	          LEV_IN_POPS_CHG(J,1)=EQSPEC+N_S
	          LEV_IN_ION_CHG(J,1)=1
	        END IF
	      ELSE
	        WRITE(LUER,*)'***********************************'
	        WRITE(LUER,*)' *** WARNNING in SET_CHG_EXCH ***'
	        WRITE(LUER,*)' Species match, but no name match'
		WRITE(LUER,*)' Charge exchange reaction:',J
		WRITE(LUER,*)' Species:',K
		WRITE(LUER,*)' Check for level naming consistency'
	        WRITE(LUER,*)'***********************************'
	      END IF			!Level set
!
	    END IF			!Species verification.
	  END DO			!K
!
500	  CONTINUE
	END DO		!J: Which charge reaction
!
	RETURN
	END
