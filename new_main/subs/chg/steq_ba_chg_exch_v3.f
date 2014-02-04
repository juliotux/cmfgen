!
! Routine to MODIFY the the statistical equilibrium equations (STEQ) and
! variation of the S.E. Eqns. (BA) for charge exchange reactions. The
! routines RD_CHG_EXCH, and SET_CHG_EXCH for EACH species, must have been
! previously called for the present iteration.
!
	SUBROUTINE STEQ_BA_CHG_EXCH_V3(POPS,T,NT,ND,DIAG_INDX,UPDATE_BA)
	USE CHG_EXCH_MOD_V3
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered  08-Jan-2001 : Bug fix (never used/tested until now).
! Altered  12-Dec-2000 : Formula only used during fitting interval.
!                        Outside this interval, the values at the limits are used.
!                        Installed TYPE_CHG=3 option.
! Modified 20-Aug-1998
! Created  24-Jun-1998
!
	INTEGER NT
	INTEGER ND
	INTEGER DIAG_INDX
!
	REAL*8 POPS(NT,ND)
	REAL*8 T(ND)
	LOGICAL UPDATE_BA
!
	INTEGER LU_ER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER J			!Reaction index
	INTEGER L			!Depth index
	INTEGER E1,E2,E3,E4
	INTEGER L1,L2,L3,L4
	INTEGER AL1,AL2,AL3,AL4
	INTEGER BL1,BL2,BL3,BL4
	INTEGER NIV
	INTEGER II
	INTEGER ID,AID,BID
!
	REAL*8 ALPHA_REC
	REAL*8 ALPHA_ION
	REAL*8 ALPHA_1
	REAL*8 ALPHA_2
	REAL*8 ALPHA_3
	REAL*8 ALPHA_4
	REAL*8 dRATEdT
	REAL*8 dlnALPHA_RECdlnT
	REAL*8 dlnALPHA_IONdlnT
	REAL*8 FRD_R
	REAL*8 REV_R
	REAL*8 T1
	REAL*8 TVAL
!
	IF(.NOT. DO_CHG_EXCH)RETURN
!
	DO J=1,N_CHG
	  IF(CHG_REACTION_AVAILABLE(J))THEN
!                                
! Reaction has the form:             
!                       Y(n+) + X([m-1]+)   <-->   Y([n-1)+) + X(m+)
! Thus species 2 and 3 are those which represent the recombined ion.
!
! NB: We have the folowing associations:
!                                       A: 2: 4:
!                                       B: 1: 3:
!
! L1,...,L4 are used to refer to the absolute level in POPS.
!
	    L1=LEV_IN_POPS_CHG(J,1);   L2=LEV_IN_POPS_CHG(J,2)
	    L3=LEV_IN_POPS_CHG(J,3);   L4=LEV_IN_POPS_CHG(J,4)
!
! We use AID to refer to the ION state 2. AL1,...,AL4 then refer to
! to the variable identification in the SE equations associated with
! that ION state.
! 
	    AID=ID_ION_CHG(J,2)
	    AL1=SE(AID)%LNK_TO_IV(L1)
	    AL2=SE(AID)%LNK_TO_IV(L2)
	    AL3=SE(AID)%LNK_TO_IV(L3)
	    AL4=SE(AID)%LNK_TO_IV(L4)
!
! We use BID to refer to the ION state 3. BL1,...,BL4 then refer to
! to the variable identification in the SE equations associated with
! that ION state.
!
	    BID=ID_ION_CHG(J,3)
	    BL1=SE(BID)%LNK_TO_IV(L1)
	    BL2=SE(BID)%LNK_TO_IV(L2)
	    BL3=SE(BID)%LNK_TO_IV(L3)
	    BL4=SE(BID)%LNK_TO_IV(L4)
!
! Define variables to avoid complicated notation.
!
! Since E2 and E3 refer to the recombined ion, the associated SE equation
! is simple the level in the ion.
!
	    E2=LEV_IN_ION_CHG(J,2)
	    E3=LEV_IN_ION_CHG(J,3)
!
! Since E1 and E4 refer to the "next" ionization stage, we have to determined
! explicitly which equation is assoicated with them.
!
	    E4=SE(AID)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_CHG(J,4))
	    E1=SE(BID)%ION_LEV_TO_EQ_PNT(LEV_IN_ION_CHG(J,1));
	    IF(E1*E2*E3*E4 .EQ. 0)THEN
	      LU_ER=ERROR_LU()
	      WRITE(LU_ER,*)'Error in STEQ_BA_CHG_EXCH_V3'
	      WRITE(LU_ER,*)'At least ne of E1,...,E4 is zero'
	      WRITE(LU_ER,*)E1,E2,E3,E4
	      STOP
	    END IF
!
! II refers to the diagonal band in the Tridiagonal matrix BA.
!
	    II=DIAG_INDX
!
	    DO L=1,ND
	      TVAL=MAX(T(L),TLO_CHG(J))
	      TVAL=MIN(TVAL,THI_CHG(J))
!
! Evaluate cross section.
!
	      IF(TYPE_CHG(J) .EQ. 1)THEN
	         ALPHA_REC=COEF_CHG(J,1)*(TVAL**COEF_CHG(J,2))
	         dlnALPHA_RECdlnT=COEF_CHG(J,2)
	      ELSE IF(TYPE_CHG(J) .EQ. 2)THEN
	         T1=COEF_CHG(J,3)*EXP(COEF_CHG(J,4)*TVAL)
	         ALPHA_REC=COEF_CHG(J,1)*(TVAL**COEF_CHG(J,2))*
	1             (1.0D0+T1)
	         dlnALPHA_RECdlnT=COEF_CHG(J,2) +
	1              COEF_CHG(J,4)*T1*TVAL/(1.0D0+T1)
	      ELSE IF(TYPE_CHG(J) .EQ. 3)THEN
	         T1=EXP(COEF_CHG(J,3)*TVAL)
                 ALPHA_REC=COEF_CHG(J,1)*(TVAL**COEF_CHG(J,2))*T1
	         dlnALPHA_RECdlnT=COEF_CHG(J,2) + COEF_CHG(J,3)*TVAL
	      ELSE
	        LU_ER=ERROR_LU()
	        WRITE(LU_ER,*)'Error in STEQ_BA_CHG_EXCH'
	        WRITE(LU_ER,*)'Invalid type of charge exchange reaction'
		WRITE(LU_ER,*)'Type of charge reaction=',TYPE_CHG(J)
	        STOP
	      END IF
	      IF(TVAL .GT. THI_CHG(J))dlnALPHA_RECdlnT=0.0D0
	      IF(TVAL .LT. TLO_CHG(J))dlnALPHA_RECdlnT=0.0D0
!
	      ALPHA_ION=ALPHA_REC*AI_AR_CHG(J,L)
	      dlnALPHA_IONdlnT=dlnALPHA_RECdlnT+dlnAI_AR_CHG_dlnT(J,L)
!
	      FRD_R=ALPHA_REC*POPS(L1,L)*POPS(L2,L)
	      REV_R=ALPHA_ION*POPS(L3,L)*POPS(L4,L)
	      T1=FRD_R-REV_R
!
! In general, E1, E2, E3 and E4 can't be zero.
!
	      SE(AID)%STEQ(E2,L)=SE(AID)%STEQ(E2,L) - T1
	      SE(AID)%STEQ(E4,L)=SE(AID)%STEQ(E4,L) + T1
!
	      SE(BID)%STEQ(E1,L)=SE(BID)%STEQ(E1,L) - T1
	      SE(BID)%STEQ(E3,L)=SE(BID)%STEQ(E3,L) + T1
!
	      IF(UPDATE_BA)THEN
	        dRATEdt=(FRD_R*dlnALPHA_RECdlnT - 
	1                   REV_R*dlnALPHA_IONdlnT)/T(L)
	        ALPHA_1=ALPHA_REC*POPS(L1,L)
	        ALPHA_2=ALPHA_REC*POPS(L2,L)
	        ALPHA_3=ALPHA_ION*POPS(L3,L)
	        ALPHA_4=ALPHA_ION*POPS(L4,L)
!
! Reaction has the form:
!                       Y(n+) + X([m-1]+)   <-->   Y([n-1)+) + X(m+)
! Thus species 2 and 3 are those which represent the recombined ion.
!
	        ID=AID;        NIV=SE(ID)%N_IV
	        IF(AL1 .NE. 0)SE(ID)%BA(E2,AL1,II,L)=SE(ID)%BA(E2,AL1,II,L) - ALPHA_2
	        IF(AL2 .NE. 0)SE(ID)%BA(E2,AL2,II,L)=SE(ID)%BA(E2,AL2,II,L) - ALPHA_1
	        IF(AL3 .NE. 0)SE(ID)%BA(E2,AL3,II,L)=SE(ID)%BA(E2,AL3,II,L) + ALPHA_4
	        IF(AL4 .NE. 0)SE(ID)%BA(E2,AL4,II,L)=SE(ID)%BA(E2,AL4,II,L) + ALPHA_3
	                      SE(ID)%BA(E2,NIV,II,L)=SE(ID)%BA(E2,NIV,II,L) - dRATEdT
!
	        IF(AL1 .NE. 0)SE(ID)%BA(E4,AL1,II,L)=SE(ID)%BA(E4,AL1,II,L) + ALPHA_2
	        IF(AL2 .NE. 0)SE(ID)%BA(E4,AL2,II,L)=SE(ID)%BA(E4,AL2,II,L) + ALPHA_1
	        IF(AL3 .NE. 0)SE(ID)%BA(E4,AL3,II,L)=SE(ID)%BA(E4,AL3,II,L) - ALPHA_4
	        IF(AL4 .NE. 0)SE(ID)%BA(E4,AL4,II,L)=SE(ID)%BA(E4,AL4,II,L) - ALPHA_3
	                      SE(ID)%BA(E4,NIV,II,L)=SE(ID)%BA(E4,NIV,II,L) + dRATEdT
!	      
	        ID=BID;        NIV=SE(ID)%N_IV
	        IF(BL1 .NE. 0)SE(ID)%BA(E3,BL1,II,L)=SE(ID)%BA(E3,BL1,II,L) + ALPHA_2
	        IF(BL2 .NE. 0)SE(ID)%BA(E3,BL2,II,L)=SE(ID)%BA(E3,BL2,II,L) + ALPHA_1
	        IF(BL3 .NE. 0)SE(ID)%BA(E3,BL3,II,L)=SE(ID)%BA(E3,BL3,II,L) - ALPHA_4
	        IF(BL4 .NE. 0)SE(ID)%BA(E3,BL4,II,L)=SE(ID)%BA(E3,BL4,II,L) - ALPHA_3
	                       SE(ID)%BA(E3,NIV,II,L)=SE(ID)%BA(E3,NIV,II,L) + dRATEdT
!
	        IF(BL1 .NE. 0)SE(ID)%BA(E1,BL1,II,L)=SE(ID)%BA(E1,BL1,II,L) - ALPHA_2
	        IF(BL2 .NE. 0)SE(ID)%BA(E1,BL2,II,L)=SE(ID)%BA(E1,BL2,II,L) - ALPHA_1
	        IF(BL3 .NE. 0)SE(ID)%BA(E1,BL3,II,L)=SE(ID)%BA(E1,BL3,II,L) + ALPHA_4
	        IF(BL4 .NE. 0)SE(ID)%BA(E1,BL4,II,L)=SE(ID)%BA(E1,BL4,II,L) + ALPHA_3
	                      SE(ID)%BA(E1,NIV,II,L)=SE(ID)%BA(E1,NIV,II,L) - dRATEdT
!
	      END IF	!Update BA?
!
	    END DO		!L (depth)
	  END IF
	END DO	!J (which reaction)
C
C As we have call STEQ_SE(ID)%BA_CHG_EXCH we set INITIALIZE_ARRAYS to tell
C SET_CHG_EXCH that the arrays must be initialized when it is next
C called.
C
	INITIALIZE_ARRAYS=.TRUE.
C
	RETURN
	END
