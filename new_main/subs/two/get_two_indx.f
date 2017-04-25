!
! Returns index value J such that FREQ lies between NU(J) and NU(J+1).
!
	INTEGER FUNCTION GET_TWO_INDX(ML,FREQ,NU,NW)
	IMPLICIT NONE
!                              
! Created 08-July-2015 : Based on GET_INDX_DP
!
	INTEGER ML
	INTEGER NW
	REAL*8 FREQ
	REAL*8 NU(NW)
!
	INTEGER J,K
	INTEGER ILOW,IHIGH
!
! This routine assumes NU monotonically decreases with frequency.
!
	IF(FREQ .GE. NU(1))THEN
	  GET_TWO_INDX=1
	END IF
	IF(FREQ .LE. NU(NW))THEN
	  GET_TWO_INDX=NW-1
	  RETURN
	END IF
!
	IF(ML .EQ. 0)THEN
	  ILOW=1
	  IHIGH=NW
	ELSE
	  ILOW=ML
	  IHIGH=ML+1
	  K=1
	  DO WHILE(FREQ .LT. NU(IHIGH))
	    ILOW=IHIGH
	    IHIGH=IHIGH+K
	    K=K*2
	  END DO
	  DO WHILE(FREQ .GT. NU(ILOW))
	    IHIGH=ILOW
	    ILOW=ILOW-K
	    K=K*2
	  END DO
	END IF
!
! Here LOW and HIGH refer to the size of the INDEX - not R which
! gets smaller as the index decreases.
!
	DO WHILE( (IHIGH-ILOW) .GT. 1)
	  J=(ILOW+IHIGH)/2
	  IF(FREQ .LT. NU(J))THEN
	    ILOW=J
	  ELSE
	    IHIGH=J
	  END IF
	END DO
	GET_TWO_INDX=ILOW
!
	RETURN
	END
