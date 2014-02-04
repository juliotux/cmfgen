!
! Subroutine to determine correspondence between SPECIES in the CHARGE exchange
! reactions, and the corresponding program variables. It is advised 
! (and program checks) that levels involved in charge exchange reactions should
! be distinct super-levels. Routine also defines the vector to compute the 
! inverse reaction rate.
!
! Program assumes (and checks) that each charge reaction involves full terms.
!
! Charge exchange reactions are assumed to have the form (and ordering)
!
!     Y(n+) + X([m-1]+)  <--> Y([n-1]+) + X(M+)
!
	SUBROUTINE SET_CHG_EXCH_V2(
	1            SPECIES,LEVEL_NAMES,EDGE_F,G_F,F_TO_S,N_F,N_S,ND,
	1            T,ZION,GION,EQSPEC,EQION,EQHYD)
	USE CHG_EXCH_MOD
	IMPLICIT NONE
!
! ALetered 09-Oct-1999 : Error reporting inproved.
! Created  20-Aug-1998 : BASED on V1. Very different calls and a change in
!                          philosiphy of how super levels are managed.
!
	INTEGER*4 N_S
	INTEGER*4 N_F
	INTEGER*4 ND
	INTEGER*4 EQSPEC
	INTEGER*4 EQION
	INTEGER*4 EQHYD
!
	REAL*8 EDGE_F(N_F)
	REAL*8 G_F(N_F)
	INTEGER*4 F_TO_S(N_F)
!
	CHARACTER*(*) SPECIES
	CHARACTER*(*) LEVEL_NAMES(N_F)
!
	REAL*8 T(ND)
	REAL*8 GION
	REAL*8 ZION
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local vectors and variables
!
	REAL*8 G_CHG_VEC(ND)
	REAL*8 dG_CHG_VEC(ND)
	REAL*8 T_VEC(ND)
	REAL*8 T1
!
	INTEGER*4 I,J,K,L
	INTEGER*4 I_S,I_F
	CHARACTER*(30) LOC_NAME
	LOGICAL LEVEL_SET
!
	IF(.NOT. DO_CHG_EXCH)RETURN
	IF( .NOT. ALLOCATED(LEV_CHG))THEN
	  ALLOCATE (Z_CHG(N_CHG,4))
	  ALLOCATE (AI_AR_CHG(N_CHG,ND))
	  ALLOCATE (dlnAI_AR_CHG_dlnT(N_CHG,ND))
	  ALLOCATE (COOL_CHG(N_CHG,ND))
	  ALLOCATE (LEV_CHG(N_CHG,4))
	  ALLOCATE (EQ_CHG(N_CHG,4))
	  ALLOCATE (EQION_CHG(N_CHG,4))
	  ALLOCATE (CHG_REACTION_AVAILABLE(N_CHG))
	  INITIALIZE_ARRAYS=.TRUE.
	END IF
C
	IF(INITIALIZE_ARRAYS)THEN
	  LEV_CHG(:,:)=0.0D0
	  EQ_CHG(:,:)=0.0D0
	  EQION_CHG(:,:)=0.0D0
	  Z_CHG(:,:)=0.0D0
	  INITIALIZE_ARRAYS=.FALSE.
!
! We must set AI_AR_CHG to unity, as we multiply it by data for each
! species in the charge exchange reaction. AI_AR_CHG allows the inverse
! reaction rate to be computed, dlnAI_AR_CHG_dlnT its variation with
! temperature.
!
	  AI_AR_CHG(:,:)=1.0D0
	  dlnAI_AR_CHG_dlnT(:,:)=0.0D0
	  COOL_CHG(:,:)=0.0D0
	END IF
!
! Check that CHARGE reactions do not refer to individual j states.
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
! Now determine whether the present species is in the CHARGE exchange reaction
! list.
!
	DO J=1,N_CHG
	  LEVEL_SET=.FALSE.
	  G_CHG_VEC(1:ND)=0.0D0
	  dG_CHG_VEC(1:ND)=0.0D0
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
	            I_S=F_TO_S(I_F)
	            LEV_CHG(J,K)=EQSPEC+I_S-1
	            EQ_CHG(J,K)=LEV_CHG(J,K)
	            EQION_CHG(J,K)=EQION
	            Z_CHG(J,K)=ZION-1.0D0
	            LEVEL_SET=.TRUE.
	          END IF
	          IF(I_S .EQ. F_TO_S(I_F))THEN
	            IF(K .EQ. 1 .OR. K .EQ. 4)THEN
	              T_VEC(1:ND)=HDKT*(EDGE_F(I_F)-EDGE_F(1))/T(1:ND)
	            ELSE
	              T_VEC(1:ND)=HDKT*EDGE_F(I_F)/T(1:ND)
	            END IF
	            DO L=1,ND
	              T1=EXP(T_VEC(L))
	              G_CHG_VEC(L)=G_CHG_VEC(L)+G_F(I_F)*T1
	              dG_CHG_VEC(L)=dG_CHG_VEC(L)-G_F(I_F)*T_VEC(L)*T1
	            END DO
	          ELSE 
	            WRITE(LUER,*)'Inconsistent level IDs in SET_CHG_EXCH'
		    WRITE(LUER,*)' Charge exchange reaction:',J
		    WRITE(LUER,*)' Species:',K
	            STOP
	          END IF
	        END IF
	      END DO
C
	      IF( LEVEL_SET)THEN
C
C Compute the mean energy change for each charge exchange reaction.
C If positive, energy is effectively removed from the electron thermal pool.
C NB:  G_CHG_VEC is proportional to the level population
C     dG_CHG_VEC is proportional to the population weighted by the level energy.
C
	        IF(K .EQ. 1 .OR. K .EQ. 3)THEN
	          DO L=1,ND
	            COOL_CHG(J,L)=COOL_CHG(J,L) -
	1               ABS(dG_CHG_VEC(L)/G_CHG_VEC(L))*T(L)/HDKT
	          END DO
	        ELSE
	          DO L=1,ND
	            COOL_CHG(J,L)=COOL_CHG(J,L) +
	1               ABS(dG_CHG_VEC(L)/G_CHG_VEC(L))*T(L)/HDKT
	          END DO
	        END IF
C
C If the charge reaction involves the final ionization stage of an atomic
C species, we need to set the DATA when we have passed data for the lower
C ionization stage. Because of our conventions, we need only check when
C K=2 or 3.
C
		IF(K .EQ. 2 .AND. EQSPEC+N_S .EQ. EQHYD)THEN
	          G_CHG_VEC(1:ND)=G_CHG_VEC(1:ND)*GION
	          dG_CHG_VEC(1:ND)=dG_CHG_VEC(1:ND)*GION
	          LEV_CHG(J,4)=EQHYD
	          Z_CHG(J,4)=ZION
	        END IF
	        IF( K .EQ. 3 .AND. EQSPEC+N_S .EQ. EQHYD )THEN
	          G_CHG_VEC(1:ND)=G_CHG_VEC(1:ND)*GION
	          dG_CHG_VEC(1:ND)=dG_CHG_VEC(1:ND)*GION
	          LEV_CHG(J,1)=EQHYD
	          Z_CHG(J,1)=ZION
	        END IF
C
	        IF(K .EQ. 1 .OR. K .EQ. 2)THEN
	          AI_AR_CHG(J,1:ND)=AI_AR_CHG(J,1:ND)*G_CHG_VEC(1:ND)
	          dlnAI_AR_CHG_dlnT(J,1:ND)=dlnAI_AR_CHG_dlnT(J,1:ND) + 
	1                            dG_CHG_VEC(1:ND)/G_CHG_VEC(1:ND)
	        ELSE 
	          AI_AR_CHG(J,1:ND)=AI_AR_CHG(J,1:ND)/G_CHG_VEC(1:ND)
	          dlnAI_AR_CHG_dlnT(J,1:ND)=dlnAI_AR_CHG_dlnT(J,1:ND) -
	1                            dG_CHG_VEC(1:ND)/G_CHG_VEC(1:ND)
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
C
	    END IF			!Species verification.
	  END DO			!K
!
500	  CONTINUE
	END DO		!J: Which charge reaction
!
	RETURN
	END
