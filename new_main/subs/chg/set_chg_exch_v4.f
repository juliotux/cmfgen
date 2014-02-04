!
! Subroutine to compute the following:
!	  AI_AR_CHG
!	  dlnAI_AR_CHG_dlnT
!	  COOL_CHG
! The first array is used for computing the reverse reaction rate so that
! LTE is recovered at depth. The second is used to compute its variation.
! The thiird array is used to help comput the charge exchange cooling rate.
!
!     Y(n+) + X([m-1]+)  <--> Y([n-1]+) + X(M+)
!
	SUBROUTINE SET_CHG_EXCH_V4(ID,LEVEL_NAMES,EDGE_F,G_F,F_TO_S,
	1            GION,N_F,N_S,ND,EQSPEC,EQHYD,T)
	USE CHG_EXCH_MOD_V3
	IMPLICIT NONE
!
! Altered 11-Apr-2003 : Simplified: No uses LEV_IN_ION_CHG to match levels with
!                         transition.
! Altered 18-Apr-2001 : ID inserted in call.
!                       Principal matching of levels is now done elsewehere 
!                         (SET_CHG_LEV_ID_V3).
! Altered 09-Oct-1999 : Error reporting inproved.
! Created 20-Aug-1998 : BASED on V1. Very different calls and a change in
!                          philosiphy of how super levels are managed.
!
	INTEGER ID    	!Ion identifier
	INTEGER N_S		!Number of super levels in atom
	INTEGER N_F		!Number of full levels in atom
	INTEGER ND		!Number of depth points
	INTEGER EQSPEC	!Eqaution of G.S. in BA matrix
	INTEGER EQHYD		!Sepecies equation in BA amtrix
!
	REAL*8 EDGE_F(N_F)
	REAL*8 G_F(N_F)
	INTEGER F_TO_S(N_F)
!
	CHARACTER(LEN=*) LEVEL_NAMES(N_F)
!
	REAL*8 T(ND)
	REAL*8 GION		!Statistical weight of ion
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
	INTEGER I,J,K,L
	INTEGER I_S,I_F
	INTEGER, PARAMETER :: LU=119
	CHARACTER(LEN=30) LOC_NAME
	LOGICAL LEVEL_SET
	LOGICAL, SAVE :: FIRST_TIME
	DATA FIRST_TIME/.TRUE./
!
!	WRITE(LUER,*)'Entering SET_CHG_EXCH_LEV_ID'
	IF(.NOT. DO_CHG_EXCH)RETURN
!
	IF(FIRST_TIME)THEN
	  WRITE(LU,*)' '
	  WRITE(LU,*)'Check of reverse charge exchange factors at depth 1.'
	  WRITE(LU,*)'To check, regroup according to reaction ID (but same order'
	  WRITE(LU,*)'within a reaction). The first column is the reaction number'
          WRITE(LU,*)'as given in CHG_EXCH_CHK. The 2nd column indicates the column'
	  WRITE(LU,*)'in charge exchange data list, while ID indicates the ion.'
	  WRITE(LU,*)'The last AI/AR ratio is the relevant conversion factor for a '
	  WRITE(LU,*)'given charge exchange reacton.'
	  WRITE(LU,*)' '
	  WRITE(LU,'(1X,A,2X,A,2X,A,6X,A,9X,A,8X,A)')
	1               'Reaction','Col','ID','AI/AR','G_EXP','GION'
	  FIRST_TIME=.FALSE.
	END IF
!
! We must set AI_AR_CHG to unity, as we multiply it by data for each
! species in the charge exchange reaction. AI_AR_CHG allows the inverse
! reaction rate to be computed, dlnAI_AR_CHG_dlnT its variation with
! temperature.
!
	IF(INITIALIZE_ARRAYS)THEN
	  AI_AR_CHG(:,:)=1.0D0
	  dlnAI_AR_CHG_dlnT(:,:)=0.0D0
	  COOL_CHG(:,:)=0.0D0
	  INITIALIZE_ARRAYS=.FALSE.
	  WRITE(LU,*)' '
	  WRITE(LU,'(1X,A,ES14.4)')'Temperature/10^4 K=',T(1)
	  WRITE(LU,*)' '
	END IF
!
! Redetermine the level for each ionic species. This is done so we
! can compute the reverse reaction rate which will recover LTE
! at depth.
!
	DO J=1,N_CHG
	  LEVEL_SET=.FALSE.
	  G_CHG_VEC(1:ND)=0.0D0
	  dG_CHG_VEC(1:ND)=0.0D0
	  LEVEL_SET=.FALSE.
	  DO K=1,4
	    IF(ID .EQ. ID_ION_CHG(J,K))THEN
	      DO I_F=1,N_F
	        I_S=F_TO_S(I_F)
	        IF(I_S .EQ. LEV_IN_ION_CHG(J,K))THEN
	          LEVEL_SET=.TRUE.
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
	        END IF
	      END DO
!
	      IF( LEVEL_SET)THEN
!
! Compute the mean energy change for each charge exchange reaction.
! If positive, energy is effectively removed from the electron thermal pool.
! NB:  G_CHG_VEC is proportional to the level population
!     dG_CHG_VEC is proportional to the population weighted by the level energy.
!
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
!
! If the charge reaction involves the final ionization stage of an atomic
! species, we need to set the DATA when we have passed data for the lower
! ionization stage. Because of our conventions, we need only check when
! K=2 or 3.
!
	        IF( (K .EQ. 2 .OR. K .EQ. 3) .AND. EQSPEC+N_S .EQ. EQHYD)THEN
	          G_CHG_VEC(1:ND)=G_CHG_VEC(1:ND)/GION
	          dG_CHG_VEC(1:ND)=dG_CHG_VEC(1:ND)/GION
	        END IF
!
	        IF(K .EQ. 1 .OR. K .EQ. 2)THEN
	          AI_AR_CHG(J,1:ND)=AI_AR_CHG(J,1:ND)*G_CHG_VEC(1:ND)
	          WRITE(LU,'(1X,3I6,3ES14.4)')J,K,ID_ION_CHG(J,K),AI_AR_CHG(J,1),G_CHG_VEC(1),GION
	          dlnAI_AR_CHG_dlnT(J,1:ND)=dlnAI_AR_CHG_dlnT(J,1:ND) + 
	1                            dG_CHG_VEC(1:ND)/G_CHG_VEC(1:ND)
	        ELSE 
	          AI_AR_CHG(J,1:ND)=AI_AR_CHG(J,1:ND)/G_CHG_VEC(1:ND)
	          WRITE(LU,'(1X,3I6,3ES14.4)')J,K,ID_ION_CHG(J,K),AI_AR_CHG(J,1),G_CHG_VEC(1),GION
	          dlnAI_AR_CHG_dlnT(J,1:ND)=dlnAI_AR_CHG_dlnT(J,1:ND) -
	1                            dG_CHG_VEC(1:ND)/G_CHG_VEC(1:ND)
	        END IF
	      ELSE
	        WRITE(LUER,*)'***********************************'
	        WRITE(LUER,*)' *** WARNNING in SET_CHG_EXCH_V3***'
	        WRITE(LUER,*)' Species match, but no name match'
	        WRITE(LUER,*)' Charge exchange reaction:',J
		WRITE(LUER,*)' Species:',K
		WRITE(LUER,*)' Check for level naming consistency'
	        WRITE(LUER,*)'***********************************'
	      END IF			!Level set
!
	    END IF			!Species verification.
	  END DO			!K
	END DO			        !J: Which charge reaction
!
	RETURN
	END
