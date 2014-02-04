!
! Routine to increment the charge conservation equation, and the
! variation charge equation. This was originally done in
! steqheii.
!
	SUBROUTINE STEQNE_V4(ED,NT,DIAG_INDX,ND,COMPUTE_BA,DST,DEND)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 14-Mar-2001 :  COMPUTE_BA installed (changed to V3)
!                        BA no longer modifed if COMPUTE_BA is FALSE.
!                        LNK_F_TO_IV installed.
! Altered 26-May-1996 :  ERROR_LU installed.
! Altered 05-Oct-1989 :  DST,DEND variable installed to avoid reading in
!                          enitre BA matrix for each ion.
! Created 15-Feb-1989
!
	INTEGER NT
	INTEGER DIAG_INDX
	INTEGER ND
	INTEGER DST,DEND
!
	REAL*8 ED(ND)
	LOGICAL COMPUTE_BA
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
! Local varaiables.
!
	INTEGER K,M,JJ
!
	DO K=DST,DEND
	  STEQ_ED(K)=STEQ_ED(K)-ED(K)
	END DO
!
	IF(COMPUTE_BA)THEN
 	  DO K=DST,DEND
	    BA_ED(NT-1,DIAG_INDX,K)=BA_ED(NT-1,DIAG_INDX,K)-1.0D0
	  END DO
	END IF
!
	RETURN
	END
