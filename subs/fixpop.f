C
C Routine to alter the Statistical Equilibrium equations so that
C a particular population is held fixed. DST AND DEND are used to
C minimize the reading of the BA matrix.
C
	SUBROUTINE FIXPOP(BA,STEQ,NT,NUM_BNDS,ND,DST,DEND,
	1                     EQSPEC,NSPEC,FIX_NSPEC,DESC,
	1                     POP,POPVEC,SPEC_PRES,FIX_IMPURITY)
	IMPLICIT NONE
C
C Altered 24-May-1996 --- CNT is now dynamically allocated.
C Altered 15-Jan-1991 --- CNT put in SAVE statement for CRAY compatibility.
C Altered 29-Aug-1991 --- Extensive changes. Call changed: FIX_IMPURIYTY
C                         option installed.
C Created 18-Oct-1989.
C
	LOGICAL FIX_IMPURITY,SPEC_PRES
	INTEGER NT,ND,NUM_BNDS,DST,DEND
	INTEGER EQSPEC,NSPEC,FIX_NSPEC
	REAL*8 BA(NT,NT,NUM_BNDS,ND),STEQ(NT,ND)
	REAL*8 POP(NSPEC,ND),POPVEC(ND)
	CHARACTER*(*) DESC	
C
C Local variables.
C
	INTEGER I,J,K,L,NDIAG,FIX_N,ERROR_LU,LUER
	LOGICAL LOC_IMP
	EXTERNAL ERROR_LU
C
C Varaibles to allow information to be output regarding the number
C of levels and depths where a population was held fixed.
C
	REAL*8 T1
	INTEGER, SAVE, ALLOCATABLE :: CNT(:)
 
C FIX_NSPEC takes priority in determining the number of levels
C to be fixed. This is necessary to fix T, for example. For this
C case, FIX_IMPURITY should be false (as no POPVEC).
C
	IF(.NOT. SPEC_PRES)RETURN
	IF(FIX_NSPEC .EQ. 0 .AND. .NOT. FIX_IMPURITY)RETURN
	IF(FIX_NSPEC .NE. 0)THEN
	  FIX_N=MIN( ABS(FIX_NSPEC),NSPEC )
	  LOC_IMP=.FALSE.
	ELSE
	  FIX_N=NSPEC
	  LOC_IMP=.TRUE.
	END IF
C
	IF(ALLOCATED(CNT)) THEN
	  IF( SIZE(CNT) .NE. NT)THEN
	    I=ERROR_LU()
	    WRITE(I,*)'Iconsistent dynamic allocation of CNT in FIXPOP'
	    STOP
	  END IF
	ELSE
	  ALLOCATE(CNT(NT))
	END IF
C
	IF(DST .EQ. 1)CNT(EQSPEC)=0	!First time in routine this pass.
C
	LUER=ERROR_LU()
	NDIAG=(NUM_BNDS+1)/2
	DO L=DST,DEND
C
C Determine whether this depth is to be held fixed, and if so
C update depth counter.
C
	  IF(LOC_IMP)THEN
	    T1=0
	    DO J=1,NSPEC
	      T1=T1+POP(J,L)
	    END DO
	    IF( T1/POPVEC(L) .GT. 1.0D-15 )GOTO 500
	  END IF
	  CNT(EQSPEC)=CNT(EQSPEC)+1
C	
	  DO K=1,NUM_BNDS
	    DO J=1,NT
	      DO I=EQSPEC,EQSPEC+FIX_N-1
	        BA(I,J,K,L)=0.0D0
	      END DO
	    END DO
	  END DO
C
	  IF(NUM_BNDS .EQ. ND)THEN
	    K=L
	  ELSE
	    K=NDIAG
	  END IF
	  DO I=EQSPEC,EQSPEC+FIX_N-1
	    BA(I,I,K,L)=1.0D0
	    STEQ(I,L)=0.0D0
	  END DO
C
500	  CONTINUE
	END DO
C
	IF(DEND .EQ. ND .AND. CNT(EQSPEC) .NE. 0)THEN
	  WRITE(LUER,100)FIX_N,NSPEC,DESC,CNT(EQSPEC)
100	  FORMAT(1X,I3,' levels of',I4,' fixed for ',A,' at',I4,' depths')
	END IF
C
	RETURN
	END
