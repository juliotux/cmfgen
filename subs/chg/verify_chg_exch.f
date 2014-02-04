!
! Routine checks which charge exchange reactions can be utilized (i.e. those
! with all species present), and whether the species in the charge exchange 
! reactions have the correct ordering.
!
	SUBROUTINE VERIFY_CHG_EXCH()
	USE CHG_EXCH_MOD
	IMPLICIT NONE
!
! Altered 09-Oct-1999 : Do not check Z's for unavailable reactions.
! Altered 20-Aug-1998
! Created  20-Jun-1998
!
	INTEGER*4 J
!
	DO J=1,N_CHG
!
	  CHG_REACTION_AVAILABLE(J)=.TRUE.
	  IF( (LEV_CHG(J,1) .EQ. 0) .OR.
	1      (LEV_CHG(J,2) .EQ. 0) .OR.
	1      (LEV_CHG(J,3) .EQ. 0) .OR.
	1      (LEV_CHG(J,4) .EQ. 0) )THEN
	    CHG_REACTION_AVAILABLE(J)=.FALSE.
	    WRITE(LUER,*)'Charge exchange reaction unavailable. No. ',J
	  ELSE
	    IF(Z_CHG(J,1) .NE. Z_CHG(J,3)+1)THEN
	      WRITE(LUER,*)'Invalid charge exchange: Reaction # ',J
	      STOP
	    END IF
	    IF(Z_CHG(J,2) .NE. Z_CHG(J,4)-1)THEN
	      WRITE(LUER,*)'Invalid charge exchange: Reaction # ',J
	      STOP
	    END IF
	    IF(EQ_CHG(J,2) .EQ. 0)THEN
	      WRITE(LUER,*)'Error in charge exchange reaction:',J
	      WRITE(LUER,*)'Species 2 must have an associated equation'
	      STOP
	    END IF
	    IF(EQ_CHG(J,3) .EQ. 0)THEN
	      WRITE(LUER,*)'Error in charge exchange reaction:',J
	      WRITE(LUER,*)'Species 3 must have an associated equation'
	      STOP
	    END IF
	  END IF		!Invalid reaction
	END DO
!
	RETURN
	END
