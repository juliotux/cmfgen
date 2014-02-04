C
C Subroutine to increment the variation matrix BAION due to the variation of J.
C
C BAION_PAR is used for the diagonal variation only. It is updated on each call
C rather than BAION to improve numerical stability. BAION_PAR should contain
C terms of similar size. BA_PAR will need to be added to BA after every
C approximately every 50 frequencies. In this way the BA and BAION matrices
C should suffer less cancelation effects due to the addition of large positive
C and negative terms.
C
C Utilizing the fact that consecutive frequency terms should be correlated, and
C hence similar in size. While some minor cancellation, should be much less
C then adding to full BA matrix in which terms have arbitrary size.
C
	SUBROUTINE BAION_UPDATE_V4(BAION,BAION_PAR,QFVION_R,QFVION_P,
	1               VJ,T,POPS,RJ,NU,FQW,NEW_CONT,FINAL_FREQ,
	1               dJ_CHK_FAC,NION,NT,NUM_BNDS,ND,DST,DEND)
	USE BA_J_DATA_MOD_V4
	IMPLICIT NONE
C
C Altered: 01-Feb-1997 :: dJ_CHK_FAC put in call (changed from V2 to V3)
C                           Replaces parameter RMAX_FAC.
C                           The smaller dJ_CHK_FAC, the more accurate the
C                           computation of the BAION matrix.
C                           dj_CHK_FAC is normally around 1.0D-04.
C                           Larger values give less accuracy for BAION,
C                           but allow faster computation.
C                         
C Altered: 16-Aug-1996 :: COMP_VEC installed to improve vectorization.
C                           Improvement will depend on how many times
C                           innermost loop is executed.
C
C Altered   24-May-1996  -- Cleaning (IONE, RMAX_FAC etc, no DABS)
C Altered:   2-May-1995  -- DST,DEND installed (Now V2)
C Created:  28-Feb-1995
C
	INTEGER NION,NT,NUM_BNDS,ND,DST,DEND
	REAL*8 BAION(NION,NT,NUM_BNDS,ND)
	REAL*8 BAION_PAR(NION,NT,ND)
	REAL*8 QFVION_R(NION,ND)
	REAL*8 QFVION_P(NION,ND)
	REAL*8 VJ(NT,NUM_BNDS,ND)
	REAL*8 POPS(NT,ND)
	REAL*8 RJ(ND)
	REAL*8 T(ND)
	REAL*8 NU
	REAL*8 FQW
	REAL*8 dJ_CHK_FAC
C
	LOGICAL NEW_CONT
	LOGICAL FINAL_FREQ
C
C Constants for opacity etc.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ

	REAL*8 COMP_VEC(NT)
C
	REAL*8 T1,T2
	INTEGER I,J,K,L,LS
	INTEGER DIAG_INDX,BNDST,BNDEND
C
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
C
	DIAG_INDX=(NUM_BNDS+1)/2
C
C NB: L refers to the depth of the appropriate ion/recom. equation.
C     K refers to the variable depth.
C     J refers to the variable.
C     I refers to the ion/recom equation.
C
C LS is used to refer to the variable depth in POPS (which is dimensioned
C (NT,ND)  ---  not with NUM_BNDS.
C
 	CALL TUNE(IONE,'BAION_UP')
	IF(NEW_CONT .AND. FINAL_FREQ)THEN
	  T2=FQW/NU
	  DO L=DST,DEND
	    T1=T2*EXP(-HDKT*NU/T(L))
	    QFVION_R(:,L)=T1*QFVION_R(:,L)-T2*QFVION_P(:,L)
	  END DO
	  DO L=DST,DEND
	    BNDST=MAX( 1+DIAG_INDX-L, 1 )
	    BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	    DO K=BNDST,BNDEND
	      LS=L+K-DIAG_INDX
	      COMP_VEC(1:NT)=ABS(VJ(1:NT,K,L)*POPS(1:NT,LS))-RJ(L)*dJ_CHK_FAC
	      IF(K .EQ. DIAG_INDX)THEN
   	        DO J=1,NT	  	  	  		!Variable
                  IF( COMP_VEC(J) .GE. 0 )THEN
	            DO  I=1,NION
	              BAION_PAR(I,J,L)=BAION_PAR(I,J,L) +
	1                  QFVION_R(I,L)*VJ(J,K,L)
	            END DO
	          END IF
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	
                  IF( COMP_VEC(J) .GE. 0 )THEN
	            DO  I=1,NION
	              BAION(I,J,K,L)=BAION(I,J,K,L) +
	1                    QFVION_R(I,L)*VJ(J,K,L)
	            END DO
	          END IF
	        END DO
	      END IF
	    END DO
	  END DO
	ELSE IF(FINAL_FREQ)THEN
C
C NB: We use VJ_P to compute COMP_VEC as this is defined as Int[ (VJ/v) dv].
C            RJ_SUM is defined in the same way (i..e., Int[ (J/v) dv]
C
	  DO L=DST,DEND
	    BNDST=MAX( 1+DIAG_INDX-L, 1 )
	    BNDEND=MIN( ND+DIAG_INDX-L, NUM_BNDS )
	    DO K=BNDST,BNDEND
	      LS=L+K-DIAG_INDX
	      COMP_VEC(:)=ABS(VJ_P(:,K,L)*POPS(:,LS))-RJ_SUM(L)*dJ_CHK_FAC
	      IF(K .EQ. DIAG_INDX)THEN
   	        DO J=1,NT	  	  	  		!Variable
                  IF( COMP_VEC(J) .GE. 0 )THEN
	            DO  I=1,NION
	              BAION_PAR(I,J,L)=BAION_PAR(I,J,L) +
	1               (QFVION_R(I,L)*VJ_R(J,K,L)-QFVION_P(I,L)*VJ_P(J,K,L))
	            END DO
	          END IF
	        END DO
	      ELSE
   	        DO  J=1,NT	  	  	
                  IF( COMP_VEC(J) .GE. 0 )THEN
	            DO  I=1,NION
	              BAION(I,J,K,L)=BAION(I,J,K,L) +
	1               (QFVION_R(I,L)*VJ_R(J,K,L)-QFVION_P(I,L)*VJ_P(J,K,L))
	            END DO
	          END IF
	        END DO
	      END IF
	    END DO
	  END DO
	END IF
	CALL TUNE(ITWO,'BAION_UP')
C
	RETURN
	END
