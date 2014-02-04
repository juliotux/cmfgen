C
C Routine to compute the net X-ray recombination rate via K shell ionization,
C and the net cooling rate.
C
C NB: IN outer regions, X-ray ionizationshould provide a net heating.
C
	SUBROUTINE X_RRR_COOL_V4(NET_X_RR,X_BFCR,WSE_X_A,WCR_X_A,
	1                     HN_A,HNST_A,N_A,HN_B,HNST_B,N_B,
	1                     JREC,JPHOT,JREC_CR,JPHOT_CR,ML,ND,FLAG)
	IMPLICIT NONE
C
C Altered 17-Sep-1997 : JREC,JPHOT now place RJ so that we can handle a
C                          constant ohotoiozation cross-section.
C                          As call very different, changed to v4.
C Altered 27-Oct-1995 : Changed to _V3
C                       EDGE variables deleted from call.
C                       WCR_X_A inserted
C Altered 06-Mar-1995 : Dimensioniong of WSE changed to (N,ND) from (N,NCF).
C                        _V2 append to name.
C Created 18-Jul-1994 : Based on PRRRCOOL
C
	INTEGER N_A,N_B,ND,ML
	REAL*8 NET_X_RR(ND),X_BFCR(ND)
	REAL*8 HN_A(N_A,ND),HNST_A(N_A,ND)
	REAL*8 HN_B(N_B,ND),HNST_B(N_B,ND)
	REAL*8 WSE_X_A(N_A,ND)
	REAL*8 WCR_X_A(N_A,ND)
	REAL*8 JREC(ND)			!In (2hv^3/c^2 + J) EXP(-hv/kT)/v  dv
	REAL*8 JPHOT(ND)		!In J/v dv
	REAL*8 JREC_CR(ND)		!In (2hv^3/c^2 + J) EXP(-hv/kT)   dv
	REAL*8 JPHOT_CR(ND)		!In J dv
C
	INTEGER I,J
	REAL*8 A1
	LOGICAL FLAG
C
C If ML=1 and flag is set, then zero all arrays. FLAG should be
C false for incrementing rates due to ionizations/recombinations to
C the 2p level.
C
	IF(ML .EQ. 1 .AND. FLAG)THEN
	  DO J=1,ND
	    NET_X_RR(J)=0.0D0
	    X_BFCR(J)=0.0D0
	  END DO
	END IF
C
C IF WSE_X_A(I,1) is zero, then it is zero for all depths since the weights
C are depth independent.
C
	DO I=1,N_A
	  IF(WSE_X_A(I,1) .NE. 0)THEN
	    DO J=1,ND
	      A1=HNST_A(I,J)*HNST_B(1,J)/HN_B(1,J)
	      NET_X_RR(J)=NET_X_RR(J) + 
	1       WSE_X_A(I,J)*( A1*JREC(J)-HN_A(I,J)*JPHOT(J) )
	      X_BFCR(J)=X_BFCR(J)   +   
	1       A1*( WCR_X_A(I,J)*JREC_CR(J) + WSE_X_A(I,J)*JREC(J) ) -
	1       HN_A(I,J)*( WCR_X_A(I,J)*JPHOT_CR(J) + WSE_X_A(I,J)*JPHOT(J) )
	    END DO
	  END IF
	END DO
C
	RETURN
	END
