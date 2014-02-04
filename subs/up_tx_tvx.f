C
C Routine to compute the matrices TX and TVX which depend on the integration
C at the previous frequency.
C
C TX(, , )  describes the variation of J
C TVX(, , ) describes the variation of RSQH
C
C It is assumed that TA,TB and TC (the tridiagonal matrix to be inverted) have
C already been modified by a call to THOMAS.
C
C KI is a 2 dimensional matrix with
C
C       KI( , ,1) variation of transfer equation w.r.t. CHI
C       KI( , ,2) variation of transfer equation w.r.t. ETA
C
C Upon entry TX( , ,1)  and TX( , ,2) should be zero, as these now reflect
C change in J with respect to CHI and ETA at the current frequency.
C The other matrices are used to reflect changes in J at an earlier frequency.
C
C In general K=1 denotes dCHI and K=2 denotes dETA.
C
	SUBROUTINE UP_TX_TVX(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                       OLD_TX,ND,NM_TX,NM_KI,
	1                       INIT,DO_THIS_TX_MATRIX)
	IMPLICIT NONE
C
C Altered 28-May-1996 : Calls to DP_ZERO removed.
C                       (, at end of subroutine specification deleted).
C
C Created 10-Mar-1995.   Based on UPTX_J_EDD_V3
C                        EPS_A,EP_B etc passed in call.
C                        Argument odering altered.
C                        TX and TVX modified in the same routine, since
C                          TVX may depend on TX at the previus frequency.
C
	INTEGER ND,NM_TX,NM_KI
	REAL*8 TX(ND,ND,NM_TX)
	REAL*8 TVX(ND-1,ND,NM_TX)
	REAL*8 KI(ND,ND,NM_KI)
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 PSIPREV_MOD(ND),VB(ND),VC(ND)
	REAL*8 HU(ND),HL(ND),HS(ND)
	REAL*8 RHS_dHdCHI(ND-1,ND)
C
C NB: _A denotes that EPS(I) multiples RJ(I)
C     _B denotes that EPS(I) multiples RJ(I+1)
C
	REAL*8 EPS_A(ND),EPS_B(ND)
	REAL*8 EPS_PREV_A(ND),EPS_PREV_B(ND)
	LOGICAL INIT,DO_THIS_TX_MATRIX(NM_TX)
C
C Work Array.
C
	REAL*8 OLD_TX(ND,ND)
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local varables.
C
	INTEGER I,J,K
	LOGICAL USE_EPS
C
C Determine whether we are using the G eddington factor to describe N
C (in terms of H) or whether N is also being described interms of J.
C
C If only G is being used, a faster section of code is excuted.
C Need to check both EPS_A and EPS_PREV_A because only 1 frequency may be using
C RSQN_ON_RSQJ. NO need to check _B, since differ by constant (factor) only.
C
	USE_EPS=.FALSE.			!G Eddington factor only
	DO I=1,ND
	  IF(EPS_A(I) .NE. 0.0D0 .OR. EPS_PREV_A(I) .NE. 0.0D0)USE_EPS=.TRUE.
	END DO
C
	IF(NM_TX .LT. 2 .OR. NM_KI .LT. 2)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Invalid NM in UP_TX_TVX'
	  WRITE(I,*)'NM_TX=',NM_TX
	  WRITE(I,*)'NM_KI=',NM_KI
	  STOP
	END IF
C
C INIT will be true for the very first frequency. We initialize all storage
C locations, even those not in use.
C
	IF(INIT)THEN
	  TX(:,:,:)=0.0D0       !ND,ND,NM_TX
	  TVX(:,:,:)=0.0D0      !(ND-1),ND*NM_TX
	  OLD_TX(:,:)=0.0D0     !ND,ND
	END IF
C
C Now modify the matrices, operating on each matrix (labeled by K) separately.
C We use OLD_TX to store TX( , ,K) at the previous frequency. Only necessary
C when N is being (at least partially) specified in terms of J.
C
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(I,J,OLD_TX)
	DO K=1,NM_TX
	  IF(DO_THIS_TX_MATRIX(K))THEN
	    IF(USE_EPS .AND. .NOT. INIT)THEN
	      DO J=1,ND
	        DO I=1,ND
	          OLD_TX(I,J)=TX(I,J,K)
	        END DO
	      END DO
	      DO J=1,ND
	        TX(1,J,K)=PSIPREV_MOD(1)*OLD_TX(1,J) + VC(1)*TVX(1,J,K)
	        DO I=2,ND-1
 	          TX(I,J,K)=PSIPREV_MOD(I)*OLD_TX(I,J)
	1                 + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	1                 + ( EPS_PREV_B(I)*OLD_TX(I+1,J)
	1                       - EPS_PREV_A(I-1)*OLD_TX(I-1,J) )
	        END DO
 	        TX(ND,J,K)= PSIPREV_MOD(ND)*OLD_TX(ND,J) +
	1                       VB(ND)*TVX(ND-1,J,K)
	      END DO
	    ELSE IF(.NOT. INIT)THEN
	      DO J=1,ND
	        TX(1,J,K)=PSIPREV_MOD(1)*TX(1,J,K) + VC(1)*TVX(1,J,K)
	        DO I=2,ND-1
 	          TX(I,J,K)=PSIPREV_MOD(I)*TX(I,J,K)
	1                 + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	        END DO
 	        TX(ND,J,K)= PSIPREV_MOD(ND)*TX(ND,J,K) + VB(ND)*TVX(ND-1,J,K)
	      END DO
	    END IF
C
	    IF(K .EQ. 1 .OR. K .EQ. 2)THEN
	      DO J=1,ND
	        DO I=1,ND
	          TX(I,J,K)=TX(I,J,K)+KI(I,J,K)
	        END DO
	      END DO
	    END IF
C
C Solve the simultaneous equations.
C
	    CALL SIMPTH(TA,TB,TC,TX(1,1,K),ND,ND)
C
	    IF(USE_EPS)THEN
	      DO J=1,ND
	        DO I=1,ND-1
	           TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1               + HS(I)*TVX(I,J,K) +
	1               (EPS_PREV_A(I)*OLD_TX(I,J)-EPS_A(I)*TX(I,J,K)) +
	1               (EPS_PREV_B(I)*OLD_TX(I+1,J)-EPS_B(I)*TX(I+1,J,K))
	        END DO
	      END DO
	    ELSE
	      DO J=1,ND
	        DO I=1,ND-1
	           TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1               + HS(I)*TVX(I,J,K)
	        END DO
	      END DO
	    END IF
C
	    IF(K .EQ. 1)THEN
	      DO J=1,ND
	        DO I=1,ND-1
	          TVX(I,J,K)=TVX(I,J,K)+RHS_dHdCHI(I,J)
	        END DO
	      END DO
	    END IF
C
	  END IF	!DO_THIS_MATRIX
	END DO		!K
!$OMP END PARALLEL DO
 
	RETURN
	END
