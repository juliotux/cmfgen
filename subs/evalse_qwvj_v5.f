C
C Subroutine to increment the statistical equilibrium equations for each 
C depth point given the value of the mean intensity at each depth point.
C
C Subroutine also increments the QFV matrix that describe the  variation 
C of the SE quations with respect to RJ.
C
C Routine also increments the ionization equilibrium equations.
C
	SUBROUTINE EVALSE_QWVJ_V5(STEQ,QFV_R,QFV_P,
	1                    WSE,HN,HNST,NLEV,NST,ION_LEV,
	1                    DI,DIST,N_DI,GS_IONEQ,
	1                    JREC,JPHOT,SPEC_EQ,NT,ND,
	1                    STEQION,QFVION_R,QFVION_P,EQUAT,NION)
	IMPLICIT NONE
C
C Altered 03-Sep-1997: QFV replaced by QFV_R, QFV_P. Allows us to spped up
C                        calculation of BA loop when the photioization
C                        cross-section is held fixed.
C                        NB: QFV = QFV_R*EMHNUKT - QFV_P
C
C Altered 29-Sep-95 : Version V4 (based on V3)
C                     Extensive changes to call (ordering AND number of 
C                       arguments)
C                     Now handles ionizations to excited states directly.
C                     Multiple states possible.
C                     No longer any need to pass HBST_P as computed on the
C                       fly.
C
C Created 09-Jul-93 - Based on EVALSE and QWVFGEN (combines functions
C                     of both routines).
C                     Routine now works for X-ray ionization,
C                     and ionizations to excited states.
C                     Note the call has 4 additional variables wrt
C                     EVALSE.
C
	INTEGER NLEV          !Number of atomic levels
	INTEGER N_DI          !Number of atomic levels in ION
	INTEGER NST		!Equation number for species
	INTEGER GS_IONEQ	!Equation number of target species (if G.S)
	INTEGER ION_LEV	!Levl ID of target in DI (i.e. the ION).
	INTEGER SPEC_EQ	!Equation number of abundance equation
	INTEGER NT		!Total number of levels
	INTEGER ND		!Number of depth points
	INTEGER EQUAT		!Equation number in ioization matrix
        INTEGER NION		!Numer of Eqns. in ionization matrix.
C
C NB --- NION is the total number of ionic species i.e. for
C HI,HII,CI,CII,CIII,CIV,CV would have NION=5 (dont count HII and CV).
C
	REAL*8 STEQ(2-NST:NT-NST+1,ND),STEQION(NION,ND)
	REAL*8 QFV_R(2-NST:NT-NST+1,ND),QFV_P(2-NST:NT-NST+1,ND)
	REAL*8 QFVION_R(NION,ND),QFVION_P(NION,ND)
	REAL*8 HN(NLEV,ND),HNST(NLEV,ND)
	REAL*8 WSE(NLEV,ND)
	REAL*8 DI(N_DI,ND),DIST(N_DI,ND)
	REAL*8 JREC(ND)
	REAL*8 JPHOT(ND)
C
C Constants for opacity etc.
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
C
C Local variables.
C
	INTEGER I,J
	REAL*8 NETR
C
C REV_HNST referes to the LTE population  of the level defined with respect
C to the actual destination (target) level.
C
	REAL*8 REV_HNST
C
	REAL*8 SUM_SE,SUM_VJ_R,SUM_VJ_P
	REAL*8 B_RAT
C
	IF(ION_LEV .EQ. 0)RETURN
C
C The net ionization (collisional and radaitive) to the last ionization stage
C must be zero from the sum of the previous equilibrum equations. Hence
C there is no need for a rate equation for the final species - it is 
C preserved for the abundance equation.
C
C If there only ionizations to the ground state, the net ionization
C term could be neglected from the rate equation for that level.
C However, ionizations to excited levels, and Auger ionization (to
C another ionization stage) mean the terms have to be explicitly included.
C
C TMP_HST is the LTE population relative to the target level in the ion.
C REV_HNST= HNST * B(ION_LEV)/B(1) where b is the deparure coefficient.
C
	DO J=1,ND
	  SUM_SE=0.0D0
	  SUM_VJ_R=0.0D0
	  SUM_VJ_P=0.0D0
	  B_RAT=(DI(ION_LEV,J)/DIST(ION_LEV,J))*(DIST(1,J)/DI(1,J))
	  DO I=1,NLEV
	    REV_HNST=HNST(I,J)*B_RAT
	    NETR=WSE(I,J)*( REV_HNST*JREC(J)-HN(I,J)*JPHOT(J) )
	    STEQ(I,J)=STEQ(I,J)+NETR
	    SUM_SE=SUM_SE+NETR
	    QFV_R(I,J)=QFV_R(I,J)+WSE(I,J)*REV_HNST
	    SUM_VJ_R=SUM_VJ_R+WSE(I,J)*REV_HNST
	    QFV_P(I,J)=QFV_P(I,J)+WSE(I,J)*HN(I,J)
	    SUM_VJ_P=SUM_VJ_P+WSE(I,J)*HN(I,J)
	  END DO
C
C We use .LT. because X-rays may ionize to a level not included.
C This way we don't need an additional check.
C
	  IF(GS_IONEQ .LT. SPEC_EQ)THEN
	    I=GS_IONEQ-NST+ION_LEV		!IONEQ-NST+1 + (ION_LEV-1)
	    STEQ(I,J)=STEQ(I,J)-SUM_SE
	    QFV_R(I,J)=QFV_R(I,J)-SUM_VJ_R
	    QFV_P(I,J)=QFV_P(I,J)-SUM_VJ_P
	  END IF
C
C NB: We do not increment the ionization equation for a species by
C     ionizations/recombinations from/to the lower ionization state.
C     This is satisfactory provided there are no transitions between
C     states differing by a charge of 2 --- such as occurs with
C     Auger ionization. In such cases the X-ray ionizations muyt be 
C     incorporated in a special way. 
C
	  IF(EQUAT .NE. 0)THEN
	    STEQION(EQUAT,J)=STEQION(EQUAT,J)+SUM_SE
	    QFVION_R(EQUAT,J)=QFVION_R(EQUAT,J)+SUM_VJ_R
	    QFVION_P(EQUAT,J)=QFVION_P(EQUAT,J)+SUM_VJ_P
	  END IF
	END DO
C
	RETURN
	END
