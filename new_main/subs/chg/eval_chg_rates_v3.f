!
! Routine to evaluate the charge exchange reactions rates for output to
! the PRRR files. This will allow a check that ionization equilibrium
! is satisfied.
!
! Also computed net cooling/heating rate for each charge exchange reaction.
!
	SUBROUTINE EVAL_CHG_RATES_V3(CHG_PR,CHG_RR,SPECIES,POPS,T,ND,NT)
	USE CHG_EXCH_MOD_V3
	IMPLICIT NONE
!
! Altered 12-Dec-2000 : Fixed reaction evaluation for TYPE_CHG=3
! Altered 10-Sep-2000 : THI_CHG limit not being correctly utilized.
! Altered 01-Oct-1999 : TLO_CHG and THI_CHG installed.
! Altered 20-Aug-1998
! Created 24-Jun-1998
!
	INTEGER NT
	INTEGER ND
	REAL*8 CHG_PR(ND)
	REAL*8 CHG_RR(ND)
!
	REAL*8 POPS(NT,ND)
	REAL*8 T(ND)
	CHARACTER*(*) SPECIES
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL*8, PARAMETER :: H=6.6261965D-27
!
! Local variables.
!
	INTEGER J,L
	INTEGER L1,L2,L3,L4
!
	REAL*8 ALPHA_REC
	REAL*8 ALPHA_ION
	REAL*8 FRD_R
	REAL*8 REV_R
	REAL*8 T1
	REAL*8 TVAL
!
! We zero the charge exchange PR (ionization) and (RR) recombination vectors
! as more than 1 reaction might contribute to each species. They must also
! be zero even if no charge reactions are included since they are tested in
! WRRECOMCHKB.
!
	CHG_PR(:)=0.0D0
	CHG_RR(:)=0.0D0
	IF(.NOT. DO_CHG_EXCH)RETURN
!
	DO J=1,N_CHG
	  IF(CHG_REACTION_AVAILABLE(J))THEN
!
! The way we have developed our ionization equtions, charge exchange
! reactions only effect the lowest ionization stage for each atomic species.
!
	    IF( SPECIES .EQ. SPEC_ID_CHG(J,3) .OR. 
	1       SPECIES .EQ. SPEC_ID_CHG(J,2) )THEN
!
! Define variables to avoid complicated notation.
!
	      L1=LEV_IN_POPS_CHG(J,1);   L2=LEV_IN_POPS_CHG(J,2)
	      L3=LEV_IN_POPS_CHG(J,3);   L4=LEV_IN_POPS_CHG(J,4)
!
	      DO L=1,ND
!
! Evaluate cross section. Outside the fitting range, we simply use the
! cross-sections at the ends of the fitting range. This avoids problems
! with potentially bad fitting formula that blow up (or become negative).
! The fitting formula usually are only really needed in the fiiting range
! anyway.
!
	        TVAL=MAX(T(L),TLO_CHG(J))
	        TVAL=MIN(TVAL,THI_CHG(J))
	        IF(TYPE_CHG(J) .EQ. 1)THEN
	           ALPHA_REC=COEF_CHG(J,1)*(TVAL**COEF_CHG(J,2))
	        ELSE IF(TYPE_CHG(J) .EQ. 2)THEN
	           T1=COEF_CHG(J,3)*EXP(COEF_CHG(J,4)*TVAL)
	           ALPHA_REC=COEF_CHG(J,1)*(TVAL**COEF_CHG(J,2))*(1.0D0+T1)
	        ELSE IF(TYPE_CHG(J) .EQ. 3)THEN
	           T1=EXP(COEF_CHG(J,3)*TVAL)
	           ALPHA_REC=COEF_CHG(J,1)*(TVAL**COEF_CHG(J,2))*T1
	        ELSE
	           WRITE(6,'(/,1X,A)')'Error -- invalid charge exchange reaction TYPE in SET_CHG_EXCH_V4'
	           STOP
	        END IF
!
	        ALPHA_ION=ALPHA_REC*AI_AR_CHG(J,L)
!
	        FRD_R=ALPHA_REC*POPS(L1,L)*POPS(L2,L)
	        REV_R=ALPHA_ION*POPS(L3,L)*POPS(L4,L)
	        IF(SPECIES .EQ. SPEC_ID_CHG(J,3))THEN
	          CHG_PR(L)=CHG_PR(L)+REV_R
	          CHG_RR(L)=CHG_RR(L)+FRD_R
	        ELSE
	          CHG_PR(L)=CHG_PR(L)+FRD_R
	          CHG_RR(L)=CHG_RR(L)+REV_R
	        END IF
!
! Factor of 10^15 is because the frequency (in COOL_CHG) is in units of
! 10^15 Hz. We must only to do this evaluation once for each charge
! exchange reaction. Before excution COOL_CHG contains the net change
! in energy (in units of 10^15 Hz) for the reaction. If negative, reaction
! is EXOTHERMIC. After evaluation, COOL_CHG contains the electron cooling
! rate in ergs/cm^3/sec.
!
	        IF( SPECIES .EQ. SPEC_ID_CHG(J,2))THEN
	          COOL_CHG(J,L)=COOL_CHG(J,L)*(FRD_R-REV_R)*H*1.0D+15
	        END IF
!
	      END DO			!L (depth)
	    END IF
	  END IF
	END DO				!J (which reaction)
!
	RETURN
	END
