	SUBROUTINE SOLVEBA_V6(STEQ,POPS,
	1              DIAG_INDX,NT,NION,NUM_BNDS,ND,
	1              MAXCH,METH_SOL,SUCCESS,SCALE_OPT,
	1              CHANGE_LIM,LAMBDA_IT,FIX_IMPURITY)
	IMPLICIT NONE
!
	INTEGER NT
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER DIAG_INDX
!
	REAL*8 STEQ(NT,ND)
	REAL*8 POPS(NT,ND)
	REAL*8 MAXCH,CHANGE_LIM
	CHARACTER*(*) METH_SOL,SCALE_OPT
	LOGICAL SUCCESS,LAMBDA_IT
	LOGICAL FIX_IMPURITY
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
!
!
! Local variables.
!
	REAL*8 SCALE,MINSCALE,T1,T2,T3
	REAL*8 INCREASE,DECREASE
	REAL*8 BIG_LIM,LIT_LIM
	INTEGER I,J,K,L,JJ,IINC,IDEC
C
	LUER=ERROR_LU()
C
	CALL SET_CASE_UP(SCALE_OPT,IONE,IZERO)
	IF(SCALE_OPT(1:5) .NE. 'LOCAL' .AND. SCALE_OPT(1:4) .NE. 'NONE'
	1   .AND. SCALE_OPT(1:6) .NE. 'GLOBAL'
	1   .AND. SCALE_OPT(1:5) .NE. 'MAJOR')THEN
	  SCALE_OPT='MAJOR'
	  WRITE(LUER,*)'Warning - Invalid scale option in SOLVEBA',
	1          ' MAJOR scaling assumed'
	END IF
!
! Solve for the perturbations. ND: The scaling of BA is now done within
! CMF_BLKBAND.
!
	IF(METH_SOL(1:4) .EQ. 'DIAG' .OR. METH_SOL(1:3) .EQ. 'TRI')THEN
	  CALL TUNE(IONE,'BLKBAND')
!
! Perform the solution for depth L. The IONE refers to DIAG_INDX, NUM_BNDS,
! and ND respectively.
!
	    CALL CMF_BLKBAND_V2(STEQ,POPS,METH_SOL,SUCCESS,
	1              DIAG_INDX,NT,NION,NUM_BNDS,ND,FIX_IMPURITY)
	    IF(.NOT. SUCCESS)THEN
	      WRITE(LUER,*)'Error in CMF_BLKBAND_V2 - shutting code down'
	      STOP
	    END IF
!
	  CALL TUNE(ITWO,'BLKBAND')
!
	ELSE
	  WRITE(LUER,*)'Error - invalid solution method in SOLVEBA'
	  WRITE(LUER,*)'Solution method is ',METH_SOL
	  STOP
	END IF
!
C 
C
C Adjust scale parameter so that the biggest decrease in any variable
C (except T) is a factor of 20, the biggest increase is a factor of CHANGE_LIM.
C T is limited to a maximum change of 20% . Three options are used.
C In case I (LOCAL), the changes are scaled local according to the largest
C single change. In case II (NONE) no scaling is done - the changes for EACH
C variable are limited to the values given above.  Case III (MAJOR) is a
C combination of LOCAL and NONE scaling. We use NONE scaling for species
C whose population is significantly less then the Electron density. In case IV
C (GLOBAL), the biggest change at any depth is used to scale all changes.
C This option is probaly obsolete.
C
	BIG_LIM=(CHANGE_LIM-1.0D0)/CHANGE_LIM
	LIT_LIM=-CHANGE_LIM
	MINSCALE=1.0D0
	IF(SCALE_OPT(1:5) .EQ. 'LOCAL')THEN
	  DO I=1,ND
	    T1=BIG_LIM			!Prevents division by zero and insures
	    T2=LIT_LIM 			!SCALE=1 if small changes.
	    DO J=1,NT-1
	      T1=MAX(T1,STEQ(J,I))   	!Note + means decrease
	      T2=MIN(T2,STEQ(J,I))      !Note - means increase
	    END DO
	    SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2 )
C
C Limit the change in T to a maximum of 20%.
C
	    T1=0.2D0
	    T3=MAX( T1,ABS(STEQ(NT,I)) )
	    SCALE=MIN( T1/T3,SCALE )
C
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0D0-STEQ(J,I)*SCALE)
	    END DO
	    MINSCALE=MIN(SCALE,MINSCALE)
	  END DO
	  WRITE(LUER,'(A,1PE12.4)')
	1  ' The local minimum value of scale is:',MINSCALE
	ELSE IF(SCALE_OPT(1:4) .EQ. 'NONE')THEN
	  DO I=1,ND
	    DO J=1,NT-1
	      IF(STEQ(J,I) .GT. BIG_LIM)THEN
	        POPS(J,I)=POPS(J,I)*(1.0D0-BIG_LIM)
	        MINSCALE=MIN( BIG_LIM/STEQ(J,I),MINSCALE )
	      ELSE IF(STEQ(J,I) .LT. LIT_LIM)THEN
	        POPS(J,I)=POPS(J,I)*(1.0D0-LIT_LIM)
	        MINSCALE=MIN( LIT_LIM/STEQ(J,I),MINSCALE )
	      ELSE
	        POPS(J,I)=POPS(J,I)*(1.0D0-STEQ(J,I))
	      END IF
	    END DO
C
C Limit T to a 20% change if it is a variable.
C
	    SCALE=0.2D0/MAX( 0.2D0,ABS(STEQ(NT,I)) )
	    POPS(NT,I)=POPS(NT,I)*(1.0D0-STEQ(NT,I)*SCALE)
	    MINSCALE=MIN(SCALE,MINSCALE)
	  END DO
C
	  WRITE(LUER,'(A,1PE12.4)')
	1  ' The minimum value of scale for all species is:',MINSCALE
	ELSE IF(SCALE_OPT(1:5) .EQ. 'MAJOR')THEN
	  DO I=1,ND
	    T1=BIG_LIM			!Prevents division by zero and insures
	    T2=LIT_LIM 			!SCALE=1 if small changes.
	    DO J=1,NT-1
	      IF(POPS(J,I) .GT. 1.0E-10*POPS(NT-1,I))THEN
	        T1=MAX(T1,STEQ(J,I))   		!Note + means decrease
	        T2=MIN(T2,STEQ(J,I))            !Note - means increase
	      END IF
	    END DO
	    SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2 )
C
C Limit the change in T to a maximum of 20%.
C
	    T3=MAX( 0.2D0,ABS(STEQ(NT,I)) )
	    SCALE=MIN( 0.2D0/T3,SCALE )
	    MINSCALE=MIN(SCALE,MINSCALE)
C
	    DO J=1,NT
	      T1=STEQ(J,I)*SCALE
	      IF(T1 .GT. BIG_LIM)T1=BIG_LIM
	      IF(T1 .LT. LIT_LIM)T1=LIT_LIM
	      POPS(J,I)=POPS(J,I)*(1.0D0-T1)
	    END DO
	  END DO
	  WRITE(LUER,'(A,1PE12.4)')
	1   ' The minimum value of scale for Major species is:',MINSCALE
C
	ELSE			!Global Scaling !
	  T1=BIG_LIM		!Prevents division by zero and insures
	  T2=LIT_LIM   		!SCALE=1 if small changes.
	  DO I=1,ND
	    DO J=1,NT-1
	      T1=MAX(STEQ(J,I),T1)			!Note + means decrease
	      T2=MIN(STEQ(J,I),T2) 			!Note - means increase
	    END DO
	  END DO
	  SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2)
C
C Limit the change in T to a maximum of 20%.
C
	  T3=MAX( 0.2D0,ABS(STEQ(NT,I)) )
	  SCALE=MIN( 0.2D0/T3,SCALE )
	  WRITE(LUER,'(A,1PE12.4)')' The value of scale is:',SCALE
C
C Update the population levels (and the temperature) .
C
	  DO I=1,ND
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0D0-STEQ(J,I)*SCALE)
	    END DO
	  END DO
	END IF
C
C Check whether the desired accuracy in the 'population parameters'
C has been obtained:
C
C   STEQ is +ve when the populations are to decrease.
C   STEQ is -ve when the populations are to increase.
C
	DECREASE=0.0D0
	IDEC=0
	INCREASE=0.0D0
	IINC=0
	DO I=1,ND
	  DO J=1,NT
	    IF(STEQ(J,I) .GT. DECREASE)THEN
	      DECREASE=STEQ(J,I)
	      IDEC=I
	    ELSE IF(STEQ(J,I) .LT. INCREASE)THEN
	      INCREASE=STEQ(J,I)
	      IINC=I
	    END IF
	  END DO
	END DO
C
	DECREASE=100.0D0*DECREASE
	INCREASE=100.0D0*INCREASE
	IF(LAMBDA_IT)THEN
 	  WRITE(LUER,9000)IINC,ABS(INCREASE),'  (LAMBDA)'
	  WRITE(LUER,9200)IDEC,DECREASE,'  (LAMBDA)'
	ELSE
 	  WRITE(LUER,9000)IINC,ABS(INCREASE)
	  WRITE(LUER,9200)IDEC,DECREASE
	END IF
9000	FORMAT(' Maximum % increase at depth ',I3,' is',1PE10.2: A)
9200	FORMAT(' Maximum % decrease at depth ',I3,' is',1PE10.2: A)
C
C Convert decrease into more useful form. If DECREASE=100, the
C new population value is zero and change is infinite. i.e we
C convert decrease to form 100(old-new)/new. Currently in the form
C 100(old-new)/old [same as INCREASE]. We allow for big decreases -
C limit MAXCH to 1.0D+07 - This value does not halt program execution.
C
	MAXCH=-INCREASE
	IF(DECREASE .LT. 99.999D0)THEN
	  DECREASE=100.0D0*DECREASE/(100.0D0-DECREASE)
	ELSE
	  DECREASE=1.0D+07
	END IF
	MAXCH=MAX(MAXCH,DECREASE)
C
	RETURN
	END
