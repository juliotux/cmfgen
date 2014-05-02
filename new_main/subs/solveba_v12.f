	SUBROUTINE SOLVEBA_V12(STEQ,POPS,POP_ATOM,DIAG_INDX,NT,NION,NUM_BNDS,ND,
	1              MAXCH,METH_SOL,SUCCESS,SCALE_OPT,LAM_SCALE_OPT,CHANGE_LIM,T_MIN,
	1              BA_COMPUTED,WR_BA_INV,WR_PRT_INV,LAMBDA_IT,
	1              MAIN_COUNTER,SET_POPS_D2_EQ_D1)
	IMPLICIT NONE
!
! Created 26-Mar-2012 : Changed to V11 and added POP_ATOM to call.
! Altered 12-MAr-2013 : We can now set the populations at depth 2 equal to those at depth 1 for 
!                          a FULL iteration (i.e., non=LAMBDA). This is designed to overcome a convergence
!                          issue with SN models. This should only be done when depth 2 is very close to 
!                          depth 1.
! Altered 08-Nov-2011 : Updated maximum correction output with info on whether BA was computed.
! Altered 15-Nov-2010 : FIDDLE_POP_CORRECTIONS added to assist convergence. This allows
!                         -ver relaxation to be used at some depths.
! Altered 23-Apr-2010 : DO_LEVEL_CHECK, MAX_INC_VEC, & MAX_DEC_VEC variables installed.
!                         Return 10th maximum value, rather than biggest correction.
!                         Done to help with a few non-converging levels.
! Altered 23-Nov-2007 : Changed to V9
!                       LAM_SCALE_OPT inserted in call.
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
	INTEGER MAIN_COUNTER
!
	REAL*8 STEQ(NT,ND)
	REAL*8 POPS(NT,ND)
	REAL*8 POP_ATOM(ND)
	REAL*8 MAXCH
	REAL*8 CHANGE_LIM
	REAL*8 T_MIN
	CHARACTER*(*) METH_SOL
	CHARACTER*(*) SCALE_OPT
	CHARACTER*(*) LAM_SCALE_OPT
	LOGICAL SUCCESS
	LOGICAL BA_COMPUTED 
	LOGICAL WR_BA_INV
	LOGICAL WR_PRT_INV
	LOGICAL LAMBDA_IT
	LOGICAL SET_POPS_D2_EQ_D1
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
	INTEGER, PARAMETER :: NV=10
	REAL*8 MAX_INC_VEC(NV)
	REAL*8 MAX_DEC_VEC(NV)
!
	REAL*8 SCALE,MINSCALE,T1,T2,T3
	REAL*8 INCREASE,DECREASE
	REAL*8 BIG_LIM,LIT_LIM
	INTEGER I,J,K,IINC,IDEC,IOS,LU_SUM
	INTEGER COUNT(7)
	LOGICAL LOC_WR_BA_INV
	LOGICAL DO_LEVEL_CHK
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
! This is done to prevent corrections at depth 2 setting the nature of future iterations.
!
	IF(SET_POPS_D2_EQ_D1 .AND. .NOT. LAMBDA_IT)THEN
	  STEQ(:,2)=0.0D0
	END IF
!
! 
!
! Lambda iterations (for Ne fixed) should generally have corrections < unity.
! Due to instabilities, large -ve corrections can sometimes arrise. This
! limits these corrections, and potentilly wil help facilitate convergence.
! 
	IF(LAMBDA_IT .AND. LAM_SCALE_OPT(1:5) .EQ. 'LIMIT')THEN
	  COUNT(1)=0
	  DO I=1,ND
	    DO J=1,NT-1
	      IF(STEQ(J,I) .GT. 1.1)THEN
	        STEQ(J,I)=0.999D0
	        COUNT(1)=COUNT(1)+1
	      END IF
	    END DO
	  END DO
	  IF(COUNT(1) .NE. 0)
	1   WRITE(LUER,*)'Warning -- using  LIMIT option for LAMBDA iteration in SOLVEBA_V9.f'
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
 	  WRITE(LUER,9000)IINC,ABS(INCREASE),'  (LAMBDA)',
	1          '          --- iteration ',MAIN_COUNTER
	  WRITE(LUER,9200)IDEC,DECREASE,'  (LAMBDA)',
	1          '          --- iteration ',MAIN_COUNTER
	ELSE IF(BA_COMPUTED)THEN
 	  WRITE(LUER,9000)IINC,ABS(INCREASE),'  (BA computed)',
	1          '     --- iteration ',MAIN_COUNTER
	  WRITE(LUER,9200)IDEC,DECREASE,'  (BA computed)',
	1          '     --- iteration ',MAIN_COUNTER
	ELSE
 	  WRITE(LUER,9000)IINC,ABS(INCREASE),'  (BA not computed)',
	1          ' --- iteration ',MAIN_COUNTER
	  WRITE(LUER,9200)IDEC,DECREASE,'  (BA not computed)',
	1          ' --- iteration ',MAIN_COUNTER
	END IF
9000	FORMAT(' Maximum % increase at depth ',I4,' is',ES10.2,A,A,I4)
9200	FORMAT(' Maximum % decrease at depth ',I4,' is',ES10.2,A,A,I4)
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
!Determine the largest set of changes:
!We determine the largest NV increases, and the largest NV decreases.
!
	MAX_INC_VEC=0.0D0
	MAX_DEC_VEC=0.0D0
	DO J=1,ND
	  DO I=1,NT
	    T1=STEQ(I,J)
	    IF(T1 .LE. 0.0D0)THEN
	      IF(T1 .LT. MAX_INC_VEC(NV))THEN
	        DO K=NV-1,1,-1
	          IF(T1 .GT. MAX_INC_VEC(K))THEN
	            MAX_INC_VEC(K+2:NV)=MAX_INC_VEC(K+1:NV-1)
	            MAX_INC_VEC(K+1)=T1
	            EXIT
	          END IF
	          IF(K .EQ. 1)THEN
	            MAX_INC_VEC(2:NV)=MAX_INC_VEC(1:NV-1)
	            MAX_INC_VEC(1)=T1
	          END IF
	        END DO
	      END IF
	    ELSE IF(STEQ(I,J) .GT. 0.0D0)THEN
	      IF(T1 .GT. MAX_DEC_VEC(NV))THEN
	        DO K=NV-1,1,-1
	          IF(T1 .LT. MAX_DEC_VEC(K))THEN
	            MAX_DEC_VEC(K+2:NV)=MAX_DEC_VEC(K+1:NV-1)
	            MAX_DEC_VEC(K+1)=T1
	            EXIT
	          END IF
	          IF(K .EQ. 1)THEN
	            MAX_DEC_VEC(2:NV)=MAX_DEC_VEC(1:NV-1)
	            MAX_DEC_VEC(1)=T1
	          END IF
	        END DO
	      END IF
	    END IF
	  END DO
	END DO
!
! DO_LEVEL_CHECL allows us to check for oscillating variables, or a few non-coverging levels.
! If the fast majority of the corrections (in this case all bar 2*NV) are less than 10%,
! we set DO_LEVEL_CHK as TRUE. When DO_LEVEL_CHECK is true we will do the following:
!     (a) Ignore the levels with corrections > 10% when computing scale (MAJOR only)
!     (b) Scale the corrections > 10% down by a factor of 3 (-ve relaxation).
!
	DO_LEVEL_CHK=.FALSE.
	IF(MAX_INC_VEC(NV) .GT. -0.1D0 .AND. MAX_DEC_VEC(NV) .LT. 0.1D0)DO_LEVEL_CHK=.TRUE.
!
! Rather than return the largest change, return the max/min of the NV correction.
!
	IF(DO_LEVEL_CHK)THEN
	  MAXCH=MAX( ABS(MAX_INC_VEC(NV)), MAX_DEC_VEC(NV)/(1.0D0-MIN(MAX_DEC_VEC(NV),0.9999D0)))
	  MAXCH=100.0D0*MAXCH
	END IF
	WRITE(6,*)'Maximm changes as returned by SOLVEBA_V9 is',MAXCH
	WRITE(6,'(A,10ES10.2)')' DEC_VEC: ',MAX_DEC_VEC(1:10)
	WRITE(6,'(A,10ES10.2)')' INC_VEC: ',MAX_INC_VEC(1:10)
!
!
!**********************************************************************************
!**********************************************************************************
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
!**********************************************************************************
!**********************************************************************************

	IF(CHANGE_LIM .LE. 1.0D0)THEN
          WRITE(LUER,'(A,1PE12.4)')' Error in SOLVEBA_V9'
          WRITE(LUER,'(A,1PE12.4)')' Maximum change for normal iteration must be > 1.'
	  STOP
	END IF
	BIG_LIM=(CHANGE_LIM-1.0D0)/CHANGE_LIM
!	LIT_LIM=-CHANGE_LIM
	LIT_LIM=1.0D0-CHANGE_LIM
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
!
!
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
!
!
	ELSE IF(SCALE_OPT(1:5) .EQ. 'MAJOR')THEN
	  CALL FIDDLE_POP_CORRECTIONS(POPS,STEQ,T_MIN,CHANGE_LIM,SCALE_OPT,LAMBDA_IT,LU_SUM,NT,ND)
!
!	  DO I=1,ND
!	    T1=BIG_LIM			!Prevents division by zero and insures
!	    T2=LIT_LIM 			!SCALE=1 if small changes.
!	    DO J=1,NT-1
!	      IF(POPS(J,I) .GT. 1.0D-10*POPS(NT-1,I))THEN
!	        T1=MAX(T1,STEQ(J,I))   		!Note + means decrease
!	        T2=MIN(T2,STEQ(J,I))            !Note - means increase
!	      END IF
!	    END DO
!	    IF(DO_LEVEL_CHK)THEN
!	      SCALE=1.0D0			!As majority of corrections < 10%
!	    ELSE
!	      SCALE=MIN( BIG_LIM/T1, LIT_LIM/T2 )
!	    END IF
!
! Limit the change in T to a maximum of 20%, and ensure T > T_MIN.
!
!	    T3=MAX( 0.2D0,ABS(STEQ(NT,I)) )
!	    SCALE=MIN( 0.2D0/T3,SCALE )
!	    MINSCALE=MIN(SCALE,MINSCALE)
!	    IF(STEQ(NT,I) .NE. 0 .AND. POPS(NT,I) .GT. T_MIN .AND.
!	1                      POPS(NT,I)*(1.0D0-STEQ(NT,I)*SCALE) .LT. T_MIN)THEN
!	        SCALE=(1.0D0-T_MIN/POPS(NT,I))/STEQ(NT,I)
!	    END IF
!	    IF(SCALE .GT. 1.0D0)SCALE=1.0D0		!i.e. will not force T to T_MIN
!!
!	    DO J=1,NT
!	      T1=STEQ(J,I)*SCALE
!	      IF(T1 .GT. BIG_LIM)T1=BIG_LIM
!	      IF(T1 .LT. LIT_LIM)T1=LIT_LIM
!	      IF(DO_LEVEL_CHK .AND. ABS(T1) .GT. 0.1D0)T1=0.3D0*T1
!!	      IF(I .EQ. 28 .OR. I .EQ. 29)T1=0.3D0*T1
!	      POPS(J,I)=POPS(J,I)*(1.0D0-T1)
!	   END DO
!	  END DO
!	  WRITE(LUER,'(A,1PE12.4)')
!	1   ' The minimum value of scale for Major species is:',MINSCALE
!
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
	IF(SET_POPS_D2_EQ_D1 .AND. .NOT. LAMBDA_IT)THEN
	  POPS(:,2)=POPS(:,1)*POP_ATOM(2)/POP_ATOM(1)
	END IF
!
	RETURN
	END
