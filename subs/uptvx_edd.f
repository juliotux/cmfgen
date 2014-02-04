C
C Routine to update the matrices which describe the V flux
C equation at each frequency. Subroutine assumes there is no
C T (temperature) variation (i.e Doppler profile is fixed.)
C Routine to be used with EDDLINE only.
C (Different sign convention for HU and HL compared with
C FORMSOL.
C
C Could used UPDATE_TVX if provided a zero TXOLD array.
	SUBROUTINE UPTVX_EDD(TVX,TX,HU,HL,HS,
	1                     RHS_dHdCHI,DEPTH_PHI,NI,NM)
	IMPLICIT NONE
C
C Altered 16-Sep-1994 - Bug Fix. I dimension if first DO block was going to ND
C                       instead of ND-1. Bug will have no effect since the 3
C                       H vectors at I=ND are (should be) zero.
C                       
C Created  6-Jun-1989 - Based on UPDATE_TVX which in turn was based on UPVNOT.
C
	INTEGER NI,NM
	REAL*8 TVX(NI-1,NI,NM)
	REAL*8 TX(NI,NI,NM)
	REAL*8 HU(NI),HL(NI),HS(NI)
	REAL*8 RHS_dHdCHI(NI-1,NI)
	REAL*8 DEPTH_PHI(NI)
C
C Local varaiables.
C
	INTEGER I,J,K
C
	DO K=1,NM
	  DO J=1,NI
	    DO I=1,NI-1
	      TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1               + HS(I)*TVX(I,J,K)
	    END DO
	  END DO
	END DO
C
	DO J=1,NI
	  DO I=1,NI-1
	    TVX(I,J,1)=TVX(I,J,1)+RHS_dHdCHI(I,J)*DEPTH_PHI(J)
	  END DO
	END DO
C
	IF(NM .EQ. 4)THEN
	  DO J=1,NI
	    DO I=1,NI-1
	      TVX(I,J,3)=TVX(I,J,3)+RHS_dHdCHI(I,J)
	    END DO
	  END DO
	END IF
C
	RETURN
	END
