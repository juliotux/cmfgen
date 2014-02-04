C
C
C Compute the new population levels and T using a block matrix,
C Gausse-Siedel iterative technique . The matrix is preconditioned
C by zeroing the block immediately below the diagonal.
C
	SUBROUTINE NEWGSIT(BA,STEQ,TBA,FQ,PIVOT,NV,ND,
	1                    MSOL,REPA,MATELIM)
	IMPLICIT NONE
C
C Altered 24-May-1996 - DOUBLE PRECISION declaration removed.
C                       ERROR_lu installed.
C                       LIMIT on NV (Through WXX) removed.
C
C Altered 04-Mar-1988 - MATELIM variable installed. Indicates the number
C                         of sub-matrices below the diagonal to be eliminated.
C Created  1/Dec/1987 - Based on NAGGSIT
C
	INTEGER NV,ND,MATELIM
	REAL*8 BA(NV,NV,ND,ND),STEQ(NV,ND),PIVOT(NV,ND)
	REAL*8 FQ(NV,ND),TBA(NV,NV,ND),REPA
	LOGICAL MSOL,TEST
C
	INTEGER I,J,K,L,LIMIT,IT
	REAL*8 RELAX
	REAL*8 WXX(NV)
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	LUER=ERROR_LU()
	MSOL=.TRUE.
C
C Prepare large matrix for iteration with the "MATELIM" blocks below the
C diagonal block eliminated. PIVOT and TBA are work arrays.
C
	CALL MODIFYBA(BA,STEQ,PIVOT,TBA,NV,ND,MATELIM,MSOL)
	IF(.NOT. MSOL)THEN
C
C Note that in this case can still continue with other technique, even though
C matrix is in a different form. Ideally would read in BA again.
C
	  WRITE(LUER,*)'Dam - NEWGSIT method doesnt work either'
	  RETURN
	END IF
C
C Initialize solution vector
C
	FQ(:,:)=0.0D0
C
C Main iteration loop    !!!!
C
 
	LIMIT=NV*ND/3
	IF(LIMIT .LT. 50)LIMIT=50
	DO 1000 IT=1,LIMIT
	  RELAX=REPA
	  IF(IT .EQ. 1)RELAX=1.0D0
	  TEST=.TRUE.			!CONVERGENCE TEST PARAMETER
C
	  DO 900 I=ND,1,-1		!EQUATION DEPTH
C
C Evaluate the steq euations at depth I interms of the unknown
C perturbations at depth I and using the estimates for the
C perturbations at the other depths.
C
	  DO J=1,NV
	    WXX(J)=STEQ(J,I)
	  END DO
C
	  DO L=1,ND
	    IF(L .EQ. I)GOTO 1500
	    DO K=1,NV
	      DO J=1,NV
	        WXX(J)=WXX(J)-BA(J,K,L,I)*FQ(K,L)
	      END DO
	    END DO
1500	    CONTINUE
	  END DO
C
C Solve for the perturbations at depth I. Because the diagonal elements
C of the sub-matrix are zero below the diagonal, and one on the diagonal,
C this is a trivial exercise.
C
	  DO J=NV,1,-1
	    DO K=J+1,NV
	      WXX(J)=WXX(J)-BA(J,K,I,I)*FQ(K,I)
	    END DO
C
C Check if desired accuracy has been obtained or whether solution
 
	    IF(ABS(WXX(J)-FQ(J,I)) .GT. ABS(FQ(J,I)/1000.0D0)
	1)    TEST=.FALSE.
	    IF(ABS(FQ(J,I)) .GT. 1.0D+10)THEN
	      MSOL=.FALSE.
	      WRITE(LUER,*)'NEWGSIT iteration blowing up.'
	      CALL WR2D(FQ,NV,ND,
	1          'Iterated solution array - blowing up',39)
	      RETURN
	    END IF
	    FQ(J,I)=FQ(J,I)+RELAX*(WXX(J)-FQ(J,I))
	  END DO
C
C End depth loop
C
900	CONTINUE
	IF(TEST)GOTO 1100
C
C End iteration loop.
C
	IF(IT .EQ. 1 .OR. IT .EQ. 20)CALL WR2D(FQ,NV,ND,
	1   'Iterated solution array',39)
1000	CONTINUE
C
C not enough iterations !!!!
C
	MSOL=.FALSE.
	WRITE(LUER,'(A)')' Not sufficient iterations'
	RETURN
C
C End iteration loop   !!!!
C
1100	CONTINUE
C
C Store results in STEQ array.
C
	DO I=1,ND
	  DO J=1,NV
	    STEQ(J,I)=FQ(J,I)
	  END DO
	END DO
C
	WRITE(LUER,'(A)')' Number of iterations required was:-'
	WRITE(LUER,*)IT
C
	RETURN
	END
