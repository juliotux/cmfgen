C
C Routine to compute the net X-ray recombination rate via K shell ionization,
C and the net cooling rate.
C
C NB: IN outer regions, X-ray ionizationshould provide a net heating.
C
	SUBROUTINE X_RRR_COOL_V6(NET_X_RR,X_BFCR,WSE_X_A,WCR_X_A,
	1                     HN_A,LOG_HNST_A,N_A,HN_B,LOG_HNST_B,N_B,
	1                     JREC,JPHOT,JREC_CR,JPHOT_CR,INIT_ARRAYS,ND,FLAG)
	IMPLICIT NONE
!
! Altered 05-Apr-2011 : Adapted from X_RRR_COOL_V5 (29-Nov-2010).
!                         LOG_HNST_A replaces HNST_A in call.
!                         Changed to facilitate a larger dynamic range in the LTE populations.
! Altered 10-Sep-2003 : Bug fix. I had the wrong J associated with the quadrature weights
!                         W in the expression for X_BFCR.
! Altered 12-Oct-2003 : If JREC is zero, the recombination term is not computed.
!                         This is to avoid floating overflows at low temperatures.
!                         In practice, the X-ray recombination term will be effectively zero
!                           at such temperatures, and hence can be neglected.
C Altered 11-Jun-2002 : Bug fix: Factor of H was dropped from X_BFCR.
C Altered 14-May-2001 : Bug fixed. Arrays were not being initialized
C                         correctly. Using continuum bands, ML may not be
C                         one on first call. Replaced ML by INIT_ARRAYS.
C                         Changed to V5.
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
	INTEGER N_A,N_B,ND
	REAL*8 NET_X_RR(ND),X_BFCR(ND)
	REAL*8 HN_A(N_A,ND),LOG_HNST_A(N_A,ND)
	REAL*8 HN_B(N_B,ND),LOG_HNST_B(N_B,ND)
	REAL*8 WSE_X_A(N_A,ND)
	REAL*8 WCR_X_A(N_A,ND)
	REAL*8 JREC(ND)			!In (2hv^3/c^2 + J) EXP(-hv/kT)/v  dv
	REAL*8 JPHOT(ND)		!In J/v dv
	REAL*8 JREC_CR(ND)		!In (2hv^3/c^2 + J) EXP(-hv/kT)   dv
	REAL*8 JPHOT_CR(ND)		!In J dv
	LOGICAL INIT_ARRAYS
C
	INTEGER I,J
	REAL*8 A1
	REAL*8 H
	LOGICAL FLAG
C
C As we multiply by hv
C
        H=6.6261965D-12 
C
C If INIT_ARRAYS and FLAG is set, then zero all arrays. FLAG should be
C false for incrementing rates due to ionizations/recombinations to
C the 2p level.
C
	IF(INIT_ARRAYS .AND. FLAG)THEN
	  DO J=1,ND
	    NET_X_RR(J)=0.0D0
	    X_BFCR(J)=0.0D0
	  END DO
	END IF
C
C IF WSE_X_A(I,1) is zero, then it is zero for all depths since the weights
C are depth independent.
C
! The check on A1 prevents overflow.
!
	DO I=1,N_A
	  IF(WSE_X_A(I,1) .NE. 0)THEN
	    DO J=1,ND
	      A1=LOG_HNST_A(I,J)+LOG_HNST_B(1,J)-LOG(HN_B(1,J))
	      IF(JREC(J) .GT. 0 .AND. A1 .LE. 280.0D0)THEN
	         A1=EXP(A1)
	      ELSE
	         A1=0.0D0
	      END IF
	      NET_X_RR(J)=NET_X_RR(J) + 
	1       WSE_X_A(I,J)*( A1*JREC(J)-HN_A(I,J)*JPHOT(J) )
	      X_BFCR(J)=X_BFCR(J)   +   
	1       (  A1*( WCR_X_A(I,J)*JREC(J) + WSE_X_A(I,J)*JREC_CR(J) ) -
	1       HN_A(I,J)*( WCR_X_A(I,J)*JPHOT(J) + WSE_X_A(I,J)*JPHOT_CR(J) )  )*H
	    END DO
	  END IF
	END DO
C
	RETURN
	END
