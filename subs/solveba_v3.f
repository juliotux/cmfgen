C
C Altered 24-JUL-82
C Altered 16-APR_85 - Call changed and nag routines called in subrouines.
C Altered 26-May-86 - Call changed. Local scaling of the solution is now an
C                     option.
C Altered 22-Jan-88 - NEWGSIT installed. Only used if other option not used.
C
C Population version - Altered to allow constant T
C Altered 04-Apr-88 - NEWGSIT installed. MATELIM variable installed.
C                     MATELIM=0 corresponds to normal block matrix
C                     iterative technique. Call to SOLVEBA changed.
C                     Implicit None installed.
c
C
C Altered 7-Feb-88  Blocktri installed, TBA made NT*NT*ND.
C                   NEWGSIT was commented out. This is only a
C                   temporary measure.
C
C Altered 14-Feb-88 BABLKBAND installed. Character string now choses
C                   method of solution. Still operates on large
C                   BA matrix.
C
C Created 16-Feb-1989 - Based on NAGSOLUT. Name change. Call changed.
C                       BA can now be of dimension BA(n,n,NUM_BNDS,nd).
C                       BA is no longer zeroed inside this routine.
C
C Altered 05-Sep-1989 - INCREASE and DECREASE variables installed. Now
C                       indicate at which depth biggest change(s) occur.
C
C Altered 19-Dec-1989 - LCLSCALE changed to a character option.
C                       'NONE' scaling option introduced.
C Altered 28-Dec-1989 - CHANGE_LIM variable intoduced. This
C                       varaible can now be changed externally,
C                       although correction to T is still limited to
C                       20%.
C
C Altered 28-Mar-1990 - METH_SOL checked to insure valid option.
C
C Altered 16-Aug-1990 - J,K loops switched in NORM section.
C
C Altered 10-Jun-1991 - MAJOR option installed. Some cleaning and tidying
C                       of option section.
C Altered 26-May-1996 - String continuations across lines removed.
C                       ERROR_LU installed.
C                       IZERO, IONE, ITWO installed for calls.
C                       0.2D0 replced bt T1 in call to MAX, MIN.
C Altered 12-Dec-1997 - WM and PIVOT removed from call, and dynamic dimensioning
!                          niw used for these variables.
! Altered 16-Jun-2000 - Change to V2: LAMBDA_IT pased in call so that it can be 
!                       output with the change to OUTGEN.
!
	SUBROUTINE SOLVEBA_V3(BA,STEQ,TBA,POPS,
	1   DIAG_IND,NT,NUM_BNDS,ND,
	1   REPA,MATELIM,MAXCH,METH_SOL,
	1   SUCCESS,SCALE_OPT,CHANGE_LIM,LAMBDA_IT)
	IMPLICIT NONE
C
	INTEGER NT,ND,MATELIM,DIAG_IND,NUM_BNDS
	REAL*8 BA(NT,NT,NUM_BNDS,ND),STEQ(NT,ND)
	REAL*8 TBA(NT,NT,ND)				!Working array
	REAL*8 POPS(NT,ND),REPA,MAXCH,CHANGE_LIM
	CHARACTER*(*) METH_SOL,SCALE_OPT
	LOGICAL SUCCESS,LAMBDA_IT
C
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
C
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: IONE=1
	INTEGER, PARAMETER :: ITWO=2
C
C Work arrays.
C
	REAL*8 WM(NT*ND)
	REAL*8 PIVOT(NT*ND)	
C
C Local variables.
C
	REAL*8 SCALE,MINSCALE,T1,T2,T3,INCREASE,DECREASE
	REAL*8 BIG_LIM,LIT_LIM
	INTEGER I,J,K,L,IST,IEND,NEWI,IINC,IDEC
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
C
C Normalize the equations to improve solution stability.
C
	CALL TUNE(IONE,'NORM')
	DO L=1,ND
	  IF(NUM_BNDS .EQ. ND)THEN
	    IST=1
	    IEND=1
	  ELSE
	    IST=MAX(DIAG_IND+1-L,1)
	    IEND=MIN(ND+DIAG_IND-L,NUM_BNDS)
	  END IF
	  DO I=IST,IEND
	    IF(ND .EQ. NUM_BNDS)THEN
	      NEWI=I
	    ELSE
	      NEWI=L+I-DIAG_IND
	    END IF
	    DO J=1,NT
	      DO K=1,NT
	        BA(K,J,I,L)=BA(K,J,I,L)*POPS(J,NEWI)
	      END DO
	    END DO
	  END DO
	END DO
	CALL TUNE(ITWO,'NORM')
C
C Solve for the perturations using either the iterative
C technique or GAUSSIAN elimination.
C
	IF(METH_SOL(1:4) .EQ. 'DIAG' .OR.
	1       METH_SOL(1:3) .EQ. 'TRI' .OR.
	1       METH_SOL(1:3) .EQ. 'PEN' )THEN
	  CALL TUNE(IONE,'BLKBAND')
	  CALL BLKBAND(BA,STEQ,TBA(1,1,1),TBA(1,1,3),METH_SOL,
	1      SUCCESS,DIAG_IND,NT,NUM_BNDS,ND)	
	  CALL TUNE(ITWO,'BLKBAND')
	  IF(.NOT. SUCCESS)THEN
	    WRITE(LUER,*)'Error in BLOCKTRI - shutting code down'
	    STOP
	  END IF
	ELSE IF(METH_SOL(1:4) .EQ. 'GSIT')THEN
	  IF(ND .NE. NUM_BNDS)THEN
	    WRITE(LUER,*)'Error - Can''t call NEWGSIT since ',
	1        'ND .ne. NUM_BNDS in SOLVEBA'
	    STOP
	  END IF
	  CALL TUNE(IONE,'NEWGSIT')
	  WRITE(LUER,*)' NEWGSIT used for solution'
	  CALL NEWGSIT(BA,STEQ,TBA,WM,PIVOT,NT,ND,SUCCESS,REPA,MATELIM)
	  CALL TUNE(ITWO,'NEWGSIT')
	  IF(.NOT. SUCCESS)METH_SOL='MIN'
	ELSE IF(METH_SOL(1:3) .EQ. 'MIN')THEN
	  IF(ND .NE. NUM_BNDS)THEN
	    WRITE(LUER,*)'Error - Can''t call MINSOL since ',
	1        'ND .ne. NUM_BNDS in SOLVEBA'
	    STOP
	  END IF
	  CALL TUNE(IONE,'NAGMINSOL')
	  WRITE(LUER,*)' NAGMINSOL used for solution'
	  CALL NAGMINSOL(BA,STEQ,TBA,WM,NT,ND)
	  CALL TUNE(ITWO,'NAGMINSOL')
	ELSE
	  WRITE(LUER,*)'Error - invalid solution method in SOLVEBA'
	  WRITE(LUER,*)'Solution method is ',METH_SOL
	  STOP
	END IF
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
