C
C Subroutine to add BA_PAR to to the full BA matrix. At the same time
C BA_PAR is zeroed.
C
C Routine may be used to update both BAION and BA.
C To update BA, pass NT for NION.
C
C This routine should not be called on very frequency. Rather it should
C be called every 50 or so frequencies. In this way the BA and BAION matrices
C should suffer less cancelation effects due due to the addition of large
C positive and negative terms.
C
C Utilizing fact that consecutive frequency terms should be correlated, and
C hence similar in size. While some minor cancellation, should be much less
C then adding to full BA matrix in which terms have arbitrary size.
C
	SUBROUTINE ADD_PAR_TO_FULL(BA,BA_PAR,NION,NT,NUM_BNDS,ND)
	IMPLICIT NONE
C
C Created:   28-Feb-1995
C
	INTEGER NION,NT,NUM_BNDS,ND
	REAL*8 BA(NION,NT,NUM_BNDS,ND)
	REAL*8 BA_PAR(NION,NT,ND)
C
	INTEGER I,J,L,K
C
	K=(NUM_BNDS+1)/2
	DO L=1,ND
	  DO J=1,NT
	    DO I=1,NION
	      BA(I,J,K,L)=BA(I,J,K,L)+BA_PAR(I,J,L)
	      BA_PAR(I,J,L)=0.0D0
	    END DO
	  END DO
	END DO
C
	RETURN
	END
