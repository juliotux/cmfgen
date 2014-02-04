!
! Routine to MODIFY the the statistical equilibrium equations (STEQ) and
! variation of the S.E. Eqns. (BA) for charge exchange reactions. The
! routines RD_CHG_EXCH, and SET_CHG_EXCH for EACH species, must have been
! previously called for the present iteration.
!
	SUBROUTINE STEQ_BA_CHG_EXCH(BA,STEQ,BAION,STEQION,
	1               POPS,T,NT,ND,NUM_BNDS,NION,DIAG_INDX,
	1               UPDATE_BA)
	USE CHG_EXCH_MOD
	IMPLICIT NONE
!
! Altered  12-Dec-2000 : Formula only used during fitting interval.
!                        Outside this interval, the values at the limits are used.
!                        Installed TYPE_CHG=3 option.
! Modified 20-Aug-1998
! Created  24-Jun-1998
!
	INTEGER*4 NT
	INTEGER*4 ND
	INTEGER*4 NUM_BNDS
	INTEGER*4 NION
	INTEGER*4 DIAG_INDX
!
	REAL*8 STEQ(NT,ND)
	REAL*8 BA(NT,NT,NUM_BNDS,ND)
	REAL*8 STEQION(NION,ND)
	REAL*8 BAION(NION,NT,NUM_BNDS,ND)
!
	REAL*8 POPS(NT,ND)
	REAL*8 T(ND)
	LOGICAL UPDATE_BA
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	INTEGER*4 LU_ER
	INTEGER*4 ERROR_LU
	EXTERNAL ERROR_LU
!
! Local variables.
!
	INTEGER*4 J,L
	INTEGER*4 E1,E2,E3,E4
	INTEGER*4 L1,L2,L3,L4
	INTEGER*4 I1,I2,I3,I4
	INTEGER*4 II
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
	DO L=1,ND
!
	  DO J=1,N_CHG
	    IF(CHG_REACTION_AVAILABLE(J))THEN
	      TVAL=MAX(T(L),TLO_CHG(J))
	      TVAL=MIN(TVAL,THI_CHG(J))
!
! Define variables to avoid complicated notation.
!
	      L1=LEV_CHG(J,1);   L2=LEV_CHG(J,2)
	      L3=LEV_CHG(J,3);   L4=LEV_CHG(J,4)
	      E1=EQ_CHG(J,1);    E2=EQ_CHG(J,2)
	      E3=EQ_CHG(J,3);    E4=EQ_CHG(J,4)
	      I1=EQION_CHG(J,1); I2=EQION_CHG(J,2)
	      I3=EQION_CHG(J,3); I4=EQION_CHG(J,4)
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
! In general E2 and E3 can't be zero. This is checked in VERIFY_CHG_EXCH.
! E1 and E4 could be zero if the charge reaction involves the ION for
! the last included ionization stage.
!
	      IF(E1 .NE. 0)STEQ(E1,L)=STEQ(E1,L) - T1
	      IF(E2 .NE. 0)STEQ(E2,L)=STEQ(E2,L) - T1
	      IF(E3 .NE. 0)STEQ(E3,L)=STEQ(E3,L) + T1
	      IF(E4 .NE. 0)STEQ(E4,L)=STEQ(E4,L) + T1
!
!NB: We need to test I2 and I3 since some species (eg H) do not necessarily
!      have space set aside for an ION equation.
!
	      IF(I2 .NE. 0)STEQION(I2,L)=STEQION(I2,L) - T1
	      IF(I3 .NE. 0)STEQION(I3,L)=STEQION(I3,L) + T1
!
	      IF(UPDATE_BA)THEN
	        dRATEdt=(FRD_R*dlnALPHA_RECdlnT - 
	1                   REV_R*dlnALPHA_IONdlnT)/T(L)
	        ALPHA_1=ALPHA_REC*POPS(L1,L)
	        ALPHA_2=ALPHA_REC*POPS(L2,L)
	        ALPHA_3=ALPHA_ION*POPS(L3,L)
	        ALPHA_4=ALPHA_ION*POPS(L4,L)
	      
	        II=DIAG_INDX
	        IF(E1 .NE. 0)THEN
	          BA(E1,L1,II,L)=BA(E1,L1,II,L) - ALPHA_2
	          BA(E1,L2,II,L)=BA(E1,L2,II,L) - ALPHA_1
	          BA(E1,L3,II,L)=BA(E1,L3,II,L) + ALPHA_4
	          BA(E1,L4,II,L)=BA(E1,L4,II,L) + ALPHA_3
	          BA(E1,NT,II,L)=BA(E1,NT,II,L) - dRATEdT
	        END IF
	        IF(E2 .NE. 0)THEN
	          BA(E2,L1,II,L)=BA(E2,L1,II,L) - ALPHA_2
	          BA(E2,L2,II,L)=BA(E2,L2,II,L) - ALPHA_1
	          BA(E2,L3,II,L)=BA(E2,L3,II,L) + ALPHA_4
	          BA(E2,L4,II,L)=BA(E2,L4,II,L) + ALPHA_3
	          BA(E2,NT,II,L)=BA(E2,NT,II,L) - dRATEdT
	        END IF
!
	        IF(E3 .NE. 0)THEN
	          BA(E3,L1,II,L)=BA(E3,L1,II,L) + ALPHA_2
	          BA(E3,L2,II,L)=BA(E3,L2,II,L) + ALPHA_1
	          BA(E3,L3,II,L)=BA(E3,L3,II,L) - ALPHA_4
	          BA(E3,L4,II,L)=BA(E3,L4,II,L) - ALPHA_3
	          BA(E3,NT,II,L)=BA(E3,NT,II,L) + dRATEdT
	        END IF
	        IF(E4 .NE. 0)THEN
	          BA(E4,L1,II,L)=BA(E4,L1,II,L) + ALPHA_2
	          BA(E4,L2,II,L)=BA(E4,L2,II,L) + ALPHA_1
	          BA(E4,L3,II,L)=BA(E4,L3,II,L) - ALPHA_4
	          BA(E4,L4,II,L)=BA(E4,L4,II,L) - ALPHA_3
	          BA(E4,NT,II,L)=BA(E4,NT,II,L) + dRATEdT
	        END IF
!
	        IF(I2 .NE. 0)THEN
	          BAION(I2,L1,II,L)=BAION(I2,L1,II,L) - ALPHA_2
	          BAION(I2,L2,II,L)=BAION(I2,L2,II,L) - ALPHA_1
	          BAION(I2,L3,II,L)=BAION(I2,L3,II,L) + ALPHA_4
	          BAION(I2,L4,II,L)=BAION(I2,L4,II,L) + ALPHA_3
	          BAION(I2,NT,II,L)=BAION(I2,NT,II,L) - dRATEdT
	        END IF
	        IF(I3 .NE. 0)THEN
	          BAION(I3,L1,II,L)=BAION(I3,L1,II,L) + ALPHA_2
	          BAION(I3,L2,II,L)=BAION(I3,L2,II,L) + ALPHA_1
	          BAION(I3,L3,II,L)=BAION(I3,L3,II,L) - ALPHA_4
	          BAION(I3,L4,II,L)=BAION(I3,L4,II,L) - ALPHA_3
	          BAION(I3,NT,II,L)=BAION(I3,NT,II,L) + dRATEdT
	        END IF
!
	      END IF	!Update BA?
!
	    END IF
	  END DO	!J (which reaction)
	END DO		!L (depth)
C
C As we have call STEQ_BA_CHG_EXCH we set INITIALIZE_ARRAYS to tell
C SET_CHG_EXCH that the arrays must be initialized when it is next
C called.
C
	INITIALIZE_ARRAYS=.TRUE.
C
	RETURN
	END
