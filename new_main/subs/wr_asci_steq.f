!
! Routine to write a 2D Matrix out as a "MATRIX" . A maximum of ten
! numbers are written across the page.
!
! Altered 13-Dec-1989 - Implicit none installed. I index written out.
!
	SUBROUTINE WR_ASCI_STEQ(NION,ND,MES,LU)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
	INTEGER ND,NION,LU
	CHARACTER*(*) MES
!
	INTEGER MS,MF,ML,I,J,ID
!
	WRITE(LU,'(/,1X,A)')MES
!
	MS=1
	DO ML=0,ND-1,10
	  MF=ML+10
	  IF(MF .GT. ND)MF=ND
	  WRITE(LU,'(/)')
	  DO ID=1,NION
	    IF(SE(ID)%XzV_PRES)THEN
	      WRITE(LU,110)1,(SE(ID)%STEQ(1,J),J=MS,MF)
	      DO I=2,SE(ID)%N_SE
	        WRITE(LU,120)I,(SE(ID)%STEQ(I,J),J=MS,MF)
	      END DO
	    END IF
	  END DO
	  MS=MS+10
	END DO
!
110	FORMAT(1X,I4,'** ',1P,10E12.4)
120	FORMAT(1X,I4,' * ',1P,10E12.4)
	RETURN
	END
