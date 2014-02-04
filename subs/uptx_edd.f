C
C Routine to compute the matrices which depend on the integration
C at the previous frequency. It is assumed that TA,TB and TC
C (the tridiagonal matrix to be inverted) have already been
C modified by a call to THOMAS2D. The sign of VB etc is for
C use with EDDLINE.
C
C We could use UPDATE_TX by setting UA(I)=UC(I)=0.0, but having a separate
C routine will save some time, and in addition we don't require the
C extra array TXOLD.
C
	SUBROUTINE UPTX_EDD(TX,TVX,KI,TA,TB,TC,U,VB,VC,
	1                    ML_NE_ONE,NI,NM)
	IMPLICIT NONE
C
C Altered 28-May-1996 - Call to DP_ZERO removed.
C Created 06-Jun-1989 - Based on UPDATE_TX which was originally
C                       based on UPDATEU.
C
	INTEGER NI,NM
	REAL*8 TX(NI,NI,NM)
	REAL*8 TVX(NI-1,NI,NM)
	REAL*8 KI(NI,NI,NM)
	REAL*8 TA(NI),TB(NI),TC(NI)
	REAL*8 U(NI),VB(NI),VC(NI)
	LOGICAL ML_NE_ONE
C
C Local varables.
C
	INTEGER I,J,K
C
	IF( .NOT. ML_NE_ONE)TX(:,:,:)=0.0D0
C
	DO K=1,NM
C
C Check to see if the matrices need to be modified by the
C calculations at the preceeding frequency.
C
	  IF(ML_NE_ONE)THEN
	    DO J=1,NI
	      TX(1,J,K)=U(1)*TX(1,J,K) + VC(1)*TVX(1,J,K)
	      DO I=2,NI-1
 	        TX(I,J,K)=U(I)*TX(I,J,K)
	1               + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	      END DO
 	      TX(NI,J,K)= U(NI)*TX(NI,J,K) + VB(NI)*TVX(NI-1,J,K)
	    END DO
	  END IF
C
	  DO J=1,NI
	    DO I=1,NI
	      TX(I,J,K)=TX(I,J,K)+KI(I,J,K)
	    END DO
	  END DO
C
C Solve the simultaneous equations.
C
	  CALL SIMPTH(TA,TB,TC,TX(1,1,K),NI,NI)
	END DO
C
	RETURN
	END
