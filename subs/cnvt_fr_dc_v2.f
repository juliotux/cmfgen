!
! Routine to convert from Log(Departure coefficients) to POPULATIONS.
! In addition, the SPECIES population (==SUM) is incremented.
! The ION contribution is only added in if the next ionization
! state is not present.
!
! The ion population for the previous ionization stage is also set.
!
	SUBROUTINE CNVT_FR_DC_V2(C2,LOG_C2LTE,DC2,NC2,DCI,SUM,ND,FIRST,CIII_PRES)
	IMPLICIT NONE
!
! Created: 19-Nov-2010 : Based on CNVT_FR_DC
!
	LOGICAL FIRST,CIII_PRES
	INTEGER NC2,ND
	REAL*8 C2(NC2,ND)
	REAL*8 LOG_C2LTE(NC2,ND)
	REAL*8 DC2(ND)
	REAL*8 DCI(ND)
	REAL*8 SUM(ND)
!
	INTEGER I,J
!
	IF(FIRST)THEN
	  DO J=1,ND
	    SUM(J)=0.0D0
	  END DO
	END IF
!
	DO J=1,ND
	  DO I=1,NC2
	    C2(I,J)=EXP(C2(I,J)+LOG_C2LTE(I,J))
	    SUM(J)=SUM(J)+C2(I,J)
	  END DO
	END DO
!
	IF(.NOT. CIII_PRES)THEN
	  DO J=1,ND
	    SUM(J)=SUM(J)+DC2(J)
	  END DO
	END IF
!
! Set ion population for next lower ionization stage,
!
	DO J=1,ND
	  DCI(J)=C2(1,J)
	END DO
!
!	FLUSH(UNIT=170)
!	WRITE(170,'(A)')'NEW'
!	WRITE(170,'(5ES16.6)')C2(1:NC2,1)	
!	WRITE(170,'(5ES16.6)')C2(1:NC2,ND)	
!	FLUSH(UNIT=170)
!
	FIRST=.FALSE.
	RETURN
	END
