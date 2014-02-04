!
! Subroutine to increment the statistical equilibrium equations for each 
! depth point given the value of the mean intensity at each depth point.
!
! Subroutine also increments the QFV matrix that describe the  variation 
! of the SE quations with respect to RJ.
!
! Routine also increments the ionization equilibrium equations.
!
	SUBROUTINE EVALSE_QWVJ_V6(ID,
	1                    WSE,HN,HNST,NLEV,ION_LEV,
	1                    DI,DIST,N_DI,
	1                    JREC,JPHOT,NT,ND)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 05-Apr-2001: Changed to utilize STEQ_DATA_MOD.
!
! Altered 03-Sep-1997: QFV replaced by QFV_R, QFV_P. Allows us to spped up
!                        calculation of BA loop when the photioization
!                        cross-section is held fixed.
!                        NB: QFV = QFV_R*EMHNUKT - QFV_P
!
! Altered 29-Sep-95 : Version V4 (based on V3)
!                     Extensive changes to call (ordering AND number of 
!                       arguments)
!                     Now handles ionizations to excited states directly.
!                     Multiple states possible.
!                     No longer any need to pass HBST_P as computed on the
!                       fly.
!
! Created 09-Jul-93 - Based on EVALSE and QWVFGEN (combines functions
!                     of both routines).
!                     Routine now works for X-ray ionization,
!                     and ionizations to excited states.
!                     Note the call has 4 additional variables wrt
!                     EVALSE.
!
	INTEGER ID            !Ionic species identifier.
	INTEGER NLEV          !Number of atomic levels
	INTEGER N_DI          !Number of atomic levels in next ioization satge.
	INTEGER ION_LEV	!Levl ID of target in DI (i.e. the ION).
	INTEGER NT		!Total number of levels
	INTEGER ND		!Number of depth points
!
! NB --- NION is the total number of ionic species i.e. for
! HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and CV).
!
	REAL*8 HN(NLEV,ND),HNST(NLEV,ND)
	REAL*8 WSE(NLEV,ND)
	REAL*8 DI(N_DI,ND),DIST(N_DI,ND)
	REAL*8 JREC(ND)
	REAL*8 JPHOT(ND)
!
! Constants for opacity etc.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
! Local variables.
!
	INTEGER I,J
	REAL*8 NETR
!
! REV_HNST referes to the LTE population  of the level defined with respect
! to the actual destination (target) level.
!
	REAL*8 REV_HNST
!
	REAL*8 SUM_SE,SUM_VJ_R,SUM_VJ_P
	REAL*8 B_RAT
!
	IF(ION_LEV .EQ. 0)RETURN
!
! The net ionization (collisional and radaitive) to the last ionization stage
! must be zero from the sum of the previous equilibrum equations. Hence
! there is no need for a rate equation for the final species - it is 
! preserved for the abundance equation.
!
! If there only ionizations to the ground state, the net ionization
! term could be neglected from the rate equation for that level.
! However, ionizations to excited levels, and Auger ionization (to
! another ionization stage) mean the terms have to be explicitly included.
!
! TMP_HST is the LTE population relative to the target level in the ion.
! REV_HNST= HNST * B(ION_LEV)/B(1) where b is the deparure coefficient.
!
	DO J=1,ND
	  SUM_SE=0.0D0
	  SUM_VJ_R=0.0D0
	  SUM_VJ_P=0.0D0
	  B_RAT=(DIST(1,J)/DIST(ION_LEV,J))*(DI(ION_LEV,J)/DI(1,J))
!	  B_RAT=(DI(ION_LEV,J)/DIST(ION_LEV,J))*(DIST(1,J)/DI(1,J))
	  DO I=1,NLEV
	    REV_HNST=HNST(I,J)*B_RAT
	    NETR=WSE(I,J)*( REV_HNST*JREC(J)-HN(I,J)*JPHOT(J) )
	    SE(ID)%STEQ(I,J)=SE(ID)%STEQ(I,J)+NETR
	    SUM_SE=SUM_SE+NETR
	    SE(ID)%QFV_R(I,J)=SE(ID)%QFV_R(I,J)+WSE(I,J)*REV_HNST
	    SUM_VJ_R=SUM_VJ_R+WSE(I,J)*REV_HNST
	    SE(ID)%QFV_P(I,J)=SE(ID)%QFV_P(I,J)+WSE(I,J)*HN(I,J)
	    SUM_VJ_P=SUM_VJ_P+WSE(I,J)*HN(I,J)
	  END DO
!
! We use .LT. because X-rays may ionize to a level not included.
! This way we don't need an additional check.
!
	  I=SE(ID)%ION_LEV_TO_EQ_PNT(ION_LEV)
	  SE(ID)%STEQ(I,J)=SE(ID)%STEQ(I,J)-SUM_SE
	  SE(ID)%QFV_R(I,J)=SE(ID)%QFV_R(I,J)-SUM_VJ_R
	  SE(ID)%QFV_P(I,J)=SE(ID)%QFV_P(I,J)-SUM_VJ_P
!
	END DO
!
	RETURN
	END
