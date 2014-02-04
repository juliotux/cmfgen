	SUBROUTINE SOLVEBA_V8(STEQ,POPS,
	1              DIAG_INDX,NT,NION,NUM_BNDS,ND,
	1              MAXCH,METH_SOL,SUCCESS,SCALE_OPT,CHANGE_LIM,T_MIN,
	1              BA_COMPUTED,WR_BA_INV,WR_PRT_INV,LAMBDA_IT)
	IMPLICIT NONE
!
! Altered 23-May-2006 : Extra column output to CORRECTION_SUM.
! Altered 19-May-2005 : T_MIN inserted into call.
! Altered 11-Feb-2004 : Bug fixed with calculation of SCALE in 'GLOBAL' scale option.  
! Altered 30-Mar-2003 : Output to CORRECTION_SUM file included.
! Altered 28-Jan-2002 : Changed to V& as BA_COMPUTED and WR_BA_INV now passed in the CALL.
!                       CMF_BLKBAND_V3 now used.
!
	INTEGER NT
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER ND
	INTEGER DIAG_INDX
!
	REAL*8 STEQ(NT,ND)
	REAL*8 POPS(NT,ND)
	REAL*8 MAXCH
	REAL*8 CHANGE_LIM
	REAL*8 T_MIN
	CHARACTER*(*) METH_SOL
	CHARACTER*(*) SCALE_OPT
	LOGICAL SUCCESS
	LOGICAL BA_COMPUTED 
	LOGICAL WR_BA_INV
	LOGICAL WR_PRT_INV
	LOGICAL LAMBDA_IT
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
	INTEGER I,J,K,IINC,IDEC,IOS,LU_SUM
	INTEGER COUNT(7)
	LOGICAL LOC_WR_BA_INV
!
	LUER=ERROR_LU()
!
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
	    LOC_WR_BA_INV=WR_BA_INV
	    IF(LAMBDA_IT)LOC_WR_BA_INV=.FALSE.
	    CALL CMF_BLKBAND_V3(STEQ,POPS,METH_SOL,SUCCESS,
	1              DIAG_INDX,NT,NION,NUM_BNDS,ND,
	1              BA_COMPUTED,LOC_WR_BA_INV,WR_PRT_INV)
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
! 
!
! Adjust scale parameter so that the biggest decrease in any variable
! (except T) is a factor of 20, the biggest increase is a factor of CHANGE_LIM.
! T is limited to a maximum change of 20% . Three options are used.
! In case I (LOCAL), the changes are scaled local according to the largest
! single change. In case II (NONE) no scaling is done - the changes for EACH
! variable are limited to the values given above.  Case III (MAJOR) is a
! combination of LOCAL and NONE scaling. We use NONE scaling for species
! whose population is significantly less then the Electron density. In case IV
! (GLOBAL), the biggest change at any depth is used to scale all changes.
! This option is probaly obsolete.
!
	IF(CHANGE_LIM .LE. 1.0D0)THEN
          WRITE(LUER,'(A,1PE12.4)')' Error in SOLVE_BA_V7'
          WRITE(LUER,'(A,1PE12.4)')' Maximum change for normal iteration must be > 1.'
	  STOP
	END IF
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
!
! Limit the change in T to a maximum of 20%, and ensure T > T_MIN.
!
	    T1=0.2D0
	    T3=MAX( T1,ABS(STEQ(NT,I)) )
	    SCALE=MIN( T1/T3,SCALE )
	    IF(STEQ(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0D0-STEQ(NT,I)*SCALE) .LT. T_MIN)THEN
	      SCALE=(1.0D0-T_MIN/POPS(NT,I))/STEQ(NT,I)
	    END IF
!
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
	    IF(POPS(NT,I) .LT. T_MIN)POPS(NT,I)=T_MIN
!
! Limit T to a 20% change if it is a variable.
!
	    SCALE=0.2D0/MAX( 0.2D0,ABS(STEQ(NT,I)) )
	    POPS(NT,I)=POPS(NT,I)*(1.0D0-STEQ(NT,I)*SCALE)
	    MINSCALE=MIN(SCALE,MINSCALE)
	  END DO
!
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
!
! Limit the change in T to a maximum of 20%, and ensure T > T_MIN.
!
	    T3=MAX( 0.2D0,ABS(STEQ(NT,I)) )
	    SCALE=MIN( 0.2D0/T3,SCALE )
	    MINSCALE=MIN(SCALE,MINSCALE)
	    IF(STEQ(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0D0-STEQ(NT,I)*SCALE) .LT. T_MIN)THEN
	        SCALE=(1.0D0-T_MIN/POPS(NT,I))/STEQ(NT,I)
	    END IF
	    IF(SCALE .GT. 1.0D0)SCALE=1.0D0		!i.e. will not force T to T_MIN
!
	    DO J=1,NT
	      T1=STEQ(J,I)*SCALE
	      IF(T1 .GT. BIG_LIM)T1=BIG_LIM
	      IF(T1 .LT. LIT_LIM)T1=LIT_LIM
	      POPS(J,I)=POPS(J,I)*(1.0D0-T1)
	    END DO
	  END DO
	  WRITE(LUER,'(A,1PE12.4)')
	1   ' The minimum value of scale for Major species is:',MINSCALE
!
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
	  DO I=1,ND
	    IF(STEQ(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
	1                      POPS(NT,I)*(1.0D0-STEQ(NT,I)*SCALE) .LT. T_MIN)THEN
	        SCALE=(1.0D0-T_MIN/POPS(NT,I))/STEQ(NT,I)
	    END IF
	  END DO
!
! Limit the change in T to a maximum of 20%.
!
	  DO I=1,ND
	    T3=MAX( 0.2D0,ABS(STEQ(NT,I)) )
	  END DO
	  SCALE=MIN( 0.2D0/T3,SCALE )
	  WRITE(LUER,'(A,1PE12.4)')' The value of scale is:',SCALE
!
! Update the population levels (and the temperature) .
!
	  DO I=1,ND
	    DO J=1,NT
	      POPS(J,I)=POPS(J,I)*(1.0D0-STEQ(J,I)*SCALE)
	    END DO
	  END DO
	END IF
!
! Determine maximum corrections to the 'population parameters', and output
! summary file:
!
	LU_SUM=7
	CALL GEN_ASCI_OPEN(LU_SUM,'CORRECTION_SUM','UNKNOWN',' ',' ',IZERO,IOS)
	WRITE(LU_SUM,'(A)')' '
	WRITE(LU_SUM,'(A)')' Summary of changes at each depth '
	WRITE(LU_SUM,'(A,I6)')' NT=',NT
	WRITE(LU_SUM,'(A)')' '
	WRITE(LU_SUM,'(3X,A,7(1X,A))')'Depth',' 100.0%','  10.0%','   1.0%',
	1        '   0.1%','  0.01%',' 0.001%','0.0001%'
!
!   STEQ is +ve when the populations are to decrease.
!   STEQ is -ve when the populations are to increase.
!
	DECREASE=0.0D0
	IDEC=0
	INCREASE=0.0D0
	IINC=0
	DO I=1,ND
	  COUNT(:)=0
	  DO J=1,NT
	    IF(STEQ(J,I) .GT. DECREASE)THEN
	      DECREASE=STEQ(J,I)
	      IDEC=I
	    ELSE IF(STEQ(J,I) .LT. INCREASE)THEN
	      INCREASE=STEQ(J,I)
	      IINC=I
	    END IF
	    T1=1.0D+06*ABS(STEQ(J,I))+1.0D-31		!to ensure non-zero.
	    K=LOG10(T1)+1; K=MIN(7,K)
	    IF(K > 0)COUNT(1:K)=COUNT(1:K)+1
	  END DO
	  WRITE(LU_SUM,'(8I8)')I,COUNT(7:1:-1)
	END DO
	CLOSE(LU_SUM)
!
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
!
! Convert decrease into more useful form. If DECREASE=100, the
! new population value is zero and change is infinite. i.e we
! convert decrease to form 100(old-new)/new. Currently in the form
! 100(old-new)/old [same as INCREASE]. We allow for big decreases -
! limit MAXCH to 1.0D+07 - This value does not halt program execution.
!
	MAXCH=-INCREASE
	IF(DECREASE .LT. 99.999D0)THEN
	  DECREASE=100.0D0*DECREASE/(100.0D0-DECREASE)
	ELSE
	  DECREASE=1.0D+07
	END IF
	MAXCH=MAX(MAXCH,DECREASE)
!
	RETURN
	END
