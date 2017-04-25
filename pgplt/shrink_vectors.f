!
! Subroutine to shrink number of data points in plotting arrays.
! This procedure cannot be undone.
!
! The error section has not been tested.
!
	SUBROUTINE SHRINK_VECTORS(CUT_ACCURACY)
	USE MOD_CURVE_DATA
	IMPLICIT NONE
!
! Created 7-Mar-2015
!
	REAL*4 CUT_ACCURACY
	REAL*4, ALLOCATABLE :: NEW_XVEC(:)
	REAL*4, ALLOCATABLE :: NEW_DATA(:)
	INTEGER NX
	INTEGER IP
	INTEGER I,J
!
	NX=MAXVAL(NPTS(1:NPLTS))
	ALLOCATE (NEW_XVEC(1:NX))
	ALLOCATE (NEW_DATA(1:NX))
!
	DO IP=1,NPLTS
	  CALL CUT_POINTS_FROM_PLOT(NEW_XVEC,NEW_DATA,NX,
	1        CD(IP)%XVEC,CD(IP)%DATA,NPTS(IP),CUT_ACCURACY)
	  WRITE(6,'(A,I2,2(A,I7))')'Curve ',IP,':  N(old)=',NPTS(IP),':  N(new)=',NX
!
! Although I am overwriting, this should work since I >= J (not tested).
! The use of an exact equality is okay, since CUT_POINTS_FROM_PLOT only removes
! data points -- it does not create new X values.
!
	  IF(ERR(IP))THEN
	    I=1
	    DO J=1,NX
	      IF(NEW_XVEC(J) .EQ. CD(IP)%XVEC(I))THEN
	        CD(IP)%EMIN(J)=CD(IP)%EMIN(I)
	        CD(IP)%EMAX(J)=CD(IP)%EMAX(I)
	      ELSE
	        I=I+1
	      END IF
	    END DO
	  END IF
	  CD(IP)%XVEC(1:NX)=NEW_XVEC(1:NX)
	  CD(IP)%DATA(1:NX)=NEW_DATA(1:NX)
	  NPTS(IP)=NX
	END DO
!
	DEALLOCATE (NEW_XVEC,NEW_DATA)
!
	RETURN
	END
