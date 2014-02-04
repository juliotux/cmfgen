C
C Routine to alter the Statistical Equilibrium equations so that
C a particular population is held fixed. DST AND DEND are used to
C minimize the reading of the BA matrix.
C
	SUBROUTINE FIXPOP_T(NT,ND,DIAG_INDX,DST,DEND,DESC)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
	!
	INTEGER NT
	INTEGER ND
	INTEGER DIAG_INDX
	INTEGER DST
	INTEGER DEND
	CHARACTER*(*) DESC
C
C Varaibles to allow information to be output regarding the number
C of levels and depths where a population was held fixed.
C
	REAL*8 T1
	INTEGER, SAVE, ALLOCATABLE :: CNT(:)
!
	INTEGER LUER
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	BA_T(:,:,DST:DEND)=0.0D0
	BA_T(NT,DIAG_INDX,DST:DEND)=1.0D0
	STEQ_T(DST:DEND)=0.0D0
!
	IF(ALLOCATED(CNT)) THEN
	  IF( SIZE(CNT) .NE. NT)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Iconsistent dynamic allocation of CNT in FIXPOP'
	    STOP
	  END IF
	ELSE
	  ALLOCATE(CNT(NT))
	END IF
C
	IF(DST .EQ. 1)CNT(NT)=0	  !First time in routine this pass.
	CNT(NT)=CNT(NT)+(DEND-DST+1)
C
	IF(DEND .EQ. ND .AND. CNT(NT) .NE. 0)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,100)CNT(NT)
100	  FORMAT(1X,' T held fixed at ',I4,' depths')
	END IF
C
	RETURN
	END
