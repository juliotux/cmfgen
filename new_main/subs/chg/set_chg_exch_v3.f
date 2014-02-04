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
	SUBROUTINE SET_CHG_EXCH_V3(ID,LEVEL_NAMES,EDGE_F,G_F,F_TO_S,
	1            GION,N_F,N_S,ND,EQSPEC,EQHYD,T)
	USE CHG_EXCH_MOD_V3
	IMPLICIT NONE
!
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
	CHARACTER*(*) LEVEL_NAMES(N_F)
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
	CHARACTER*(30) LOC_NAME
	LOGICAL LEVEL_SET
!
	IF(.NOT. DO_CHG_EXCH)RETURN
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
	END IF
!
! Redetermine the level for each ionic species. This is done so we
! can compute the reverse reaction rate which will recover LTE
! at depth.
!
	DO J=1,N_CHG
	  IF(CHG_REACTION_AVAILABLE(J))THEN
	    LEVEL_SET=.FALSE.
	    G_CHG_VEC(1:ND)=0.0D0
	    dG_CHG_VEC(1:ND)=0.0D0
	    DO K=1,4
	      IF(ID .EQ. ID_ION_CHG(J,K))THEN
	        I_S=0
	        DO I_F=1,N_F
	          LOC_NAME=LEVEL_NAMES(I_F)
	          L=INDEX(LOC_NAME,'[')
	          IF(L .NE. 0)LOC_NAME=LOC_NAME(1:L-1)
	          IF(LOC_NAME .EQ. LEV_NAME_CHG(J,K) .OR.
	1                 LOC_NAME .EQ. ALT_LEV_NAME_CHG(J,K))THEN
	            IF(I_S .EQ. 0)THEN
	              I_S=F_TO_S(I_F)
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
	              WRITE(LUER,*)'Inconsistent level IDs in SET_CHG_EXCH_V3'
	              WRITE(LUER,*)'Charge exchange reaction:',J
	              WRITE(LUER,*)'Species:',K
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
	1                 ABS(dG_CHG_VEC(L)/G_CHG_VEC(L))*T(L)/HDKT
	            END DO
	          ELSE
	            DO L=1,ND
	              COOL_CHG(J,L)=COOL_CHG(J,L) +
	1                 ABS(dG_CHG_VEC(L)/G_CHG_VEC(L))*T(L)/HDKT
	            END DO
	          END IF
C
C If the charge reaction involves the final ionization stage of an atomic
C species, we need to set the DATA when we have passed data for the lower
C ionization stage. Because of our conventions, we need only check when
C K=2 or 3.
C
	          IF( (K .EQ. 2 .OR. K .EQ. 3) .AND. EQSPEC+N_S .EQ. EQHYD)THEN
	            G_CHG_VEC(1:ND)=G_CHG_VEC(1:ND)*GION
	            dG_CHG_VEC(1:ND)=dG_CHG_VEC(1:ND)*GION
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
	          WRITE(LUER,*)' *** WARNNING in SET_CHG_EXCH_V3***'
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
	  END IF        !Reaction available ?
	END DO		!J: Which charge reaction
!
	RETURN
	END
