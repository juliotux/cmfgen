C
C Subroutine to perform liner interpolations for RDPOP.
C The edges of the new mesh must lie on or inside the edges of
C the old mesh. R may be monotonically increasing, or decreasing.
C
C Improved error handling (5-OCT-1983)
C
C Altered 11-Apr-1989 - Implicit NONE INSTALLED. Main call changed to
C                       LIN_INTERP, and LINPOP and LININT entry points 
C                       installed. All three routines are identical
C                       and were all included so save editing a lot of
C                       routines that call LINPOP and LINIT. LIN_INTERP
C                       is the preferred call.
C
	SUBROUTINE LIN_INTERP(R,V,ND,U,W,NIN)
	IMPLICIT NONE
	INTEGER NIN,ND
	REAL*8 R(ND),V(ND),U(NIN),W(NIN)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local varaibles.
C
	REAL*8 T1
	INTEGER J,L
C
	ENTRY LININT(R,V,ND,U,W,NIN)
	ENTRY LINPOP(R,V,ND,U,W,NIN)
C
	IF((R(2)-R(1)) .GT. 0.0D0 )THEN
	  IF(R(ND) .GT. U(NIN) .OR. R(1) .LT. U(1))THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in LIN_INTERP/LINPOP - values outside range'
	    WRITE(J,*)R(ND),'should be LE',U(NIN)
	    WRITE(J,*)R(1),'should be GE ',U(1)
	    WRITE(J,*)'ND=',ND,'NIN=',NIN
	    STOP
	  END IF
	  L=2
	  DO J=1,ND
50	    IF(R(J) .LE. U(L))THEN
	      T1=(U(L)-R(J))/(U(L)-U(L-1))
	      V(J)=T1*W(L-1)+(1.0D0-T1)*W(L)
	    ELSE
	      L=L+1
	      GOTO 50
	    END IF
	  END DO
	ELSE
C
C**********************************************************************
C**********************************************************************
C
	  IF(R(1) .GT. U(1) .OR. R(ND) .LT. U(NIN))THEN
	    J=ERROR_LU()
	    WRITE(J,*)'Error in LIN_INTERP/LINPOP - values outside range'
	    WRITE(J,*)R(ND),'should be GE',U(NIN)
	    WRITE(J,*)R(1),'should be LE ',U(1)
	    WRITE(J,*)'ND=',ND,'NIN=',NIN
	    STOP
	  END IF
	  L=NIN-1		!-1 to stay in array bounds
	  DO J=ND,1,-1
500	    IF(R(J) .LE. U(L))THEN
	      T1=(U(L+1)-R(J))/(U(L+1)-U(L))
	      V(J)=T1*W(L)+(1.0D0-T1)*W(L+1)
	    ELSE
	      L=L-1
	      GOTO 500
	    END IF
	  END DO
C
	END IF
C
	RETURN
	END
