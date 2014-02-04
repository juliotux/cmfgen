!
! Returns index value J such that RVAL lies between R(J) and
! R(J+1).
!
	INTEGER FUNCTION GET_INDX_DP(RVAL,R,NW)
	IMPLICIT NONE
!                              
! Altered 02-May-1990 - Bug fix: If on two consecutive calls, RVAL was the
!                       same but the R array different, the returned index
!                       was incorrect (referring to the first array). No
!                       longer use RVAL_SAV, but directly check R array.
!                       Routine now correctly returns 1 or NW-1 if RVAL is 
!                       outside R range. Previously only J and not GET_INDX_DP
!                       was being set. (Minor correction on 8/May/91).
!
! Altered 20-Nov-1990 - Call to SEED_RITE changed.
! Altered 12-Nov-1990 = RVAL check altered. Routine nolonger STOPS program.
! Altered 23-Apr-1990 - Check that RVAL is inside range given by R included.
! Created 11-Apr-1990
!
	INTEGER NW
	REAL*8 R(NW),RVAL
!
	INTEGER J,ILOW,IHIGH,ILOW_SAV
	SAVE ILOW_SAV
	DATA ILOW_SAV/1/
!
!
! NW may have changed bewteen this and previous call, hence ILOW_SAV may
! lie outside the valid range (1:NW).
!
	IF(ILOW_SAV .GE. NW)ILOW_SAV=NW-1	
!
! Have two cases to find RVAL for, depending on whether the array R
! increases or decreases with its index.
!
	IF(R(1) .GT. R(NW))THEN
!
! Need t compare RVAl with R array (and not RSAV) in case value is
! same, but array is different.
!
	  IF( (RVAL .GE.  R(ILOW_SAV+1)) .AND. 
	1       (RVAL .LE. R(ILOW_SAV)) )THEN
            GET_INDX_DP=ILOW_SAV
	    RETURN
	  END IF
	  IF(RVAL .GT. R(1))THEN
	    WRITE(2,*)'Error in GET_INDX_DP - RVAL too large.'
	    WRITE(2,70)'R(1)=',R(1),'RVAL=',RVAL
70	    FORMAT(1X,1P,(3X,A,E12.4))
	    RVAL=R(1)
	    GET_INDX_DP=1
	    RETURN
	  END IF
	  IF(RVAL .LT. R(NW))THEN
	    WRITE(2,*)'Error in GET_INDX_DP - RVAL too small.'
	    WRITE(2,70)'R(NW)=',R(NW),'RVAL=',RVAL
	    RVAL=R(NW)
	    GET_INDX_DP=NW-1
	    RETURN
	  END IF
!
! Here LOW and HIGH refer to the size of the INDEX - not R which
! gets smaller as the index decreases.
!
	  ILOW=1
	  IHIGH=NW
	  DO WHILE( (IHIGH-ILOW) .GT. 1)
	    J=(ILOW+IHIGH)/2
	    IF(RVAL .LT. R(J))THEN
	      ILOW=J
	    ELSE
	      IHIGH=J
	    END IF
	  END DO
	  ILOW_SAV=ILOW
	  GET_INDX_DP=ILOW
	  RETURN
	ELSE
!
! Need t compare RVAl with R array (and not RSAV) in case value is
! same, but array is different.
!
	  IF( (RVAL .LE.  R(ILOW_SAV+1)) .AND. 
	1       (RVAL .GE. R(ILOW_SAV)) )THEN
            GET_INDX_DP=ILOW_SAV
	    RETURN
	  END IF
	  IF(RVAL .LT. R(1))THEN
	    WRITE(2,*)'Error in GET_INDX_DP - RVAL too small.'
	    WRITE(2,70)'R(1)=',R(1),'RVAL=',RVAL
	    RVAL=R(1)
	    GET_INDX_DP=1
	    RETURN
	  END IF
	  IF(RVAL .GT. R(NW))THEN
	    WRITE(2,*)'Warning in GET_INDX_DP - RVAL too large.'
	    WRITE(2,70)'R(NW)=',R(NW),'RVAL=',RVAL
	    RVAL=R(NW)
	    GET_INDX_DP=NW-1
	    RETURN
	  END IF
!
	  ILOW=1
	  IHIGH=NW
	  DO WHILE( (IHIGH-ILOW) .GT. 1)
	    J=(ILOW+IHIGH)/2
	    IF(RVAL .LT. R(J))THEN
	      IHIGH=J
	    ELSE
	      ILOW=J
	    END IF
	  END DO
	  ILOW_SAV=ILOW
	  GET_INDX_DP=ILOW
	  RETURN
	END IF
!
	END FUNCTION GET_INDX_DP
