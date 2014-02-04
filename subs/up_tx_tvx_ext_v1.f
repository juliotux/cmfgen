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
C KI is a 3 dimensional matrix with
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
C This version asumes we are computing J at NDEXT points based on primary
C data at ND nodes.
C
C This routine can superced UP_TX_TVX by passing with ND=NDEXT
C
	SUBROUTINE UP_TX_TVX_EXT_V1(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       EPS_A,EPS_B,EPS_PREV_A,EPS_PREV_B,
	1                       OLD_TX,NDEXT,ND,NM_TX,NM_KI,
	1                       INIT,DO_THIS_TX_MATRIX)
	IMPLICIT NONE
C
C Created : 15-May-1987  Based on UP_TX_TVX
C
	INTEGER NDEXT,ND,NM_TX,NM_KI
	REAL*8 TX(NDEXT,ND,NM_TX)
	REAL*8 TVX(NDEXT-1,ND,NM_TX)
	REAL*8 KI(NDEXT,ND,NM_KI)
	REAL*8 TA(NDEXT),TB(NDEXT),TC(NDEXT)
	REAL*8 PSIPREV_MOD(NDEXT),VB(NDEXT),VC(NDEXT)
	REAL*8 HU(NDEXT),HL(NDEXT),HS(NDEXT)
	REAL*8 RHS_dHdCHI(NDEXT-1,ND)
C
C NB: _A denotes that EPS(I) multiples RJ(I)
C     _B denotes that EPS(I) multiples RJ(I+1)
C
	REAL*8 EPS_A(NDEXT),EPS_B(NDEXT)
	REAL*8 EPS_PREV_A(NDEXT),EPS_PREV_B(NDEXT)
	LOGICAL INIT,DO_THIS_TX_MATRIX(NM_TX)
C
C Work Array.
C
	REAL*8 OLD_TX(NDEXT,ND)
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
	DO I=1,NDEXT
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
	  TX(:,:,:)=0.0D0       !NDEXT,ND,NM_TX
	  TVX(:,:,:)=0.0D0      !(NDEXT-1),ND*NM_TX
	  OLD_TX(:,:)=0.0D0     !NDEXT,ND
	END IF
C
C Now modify the matrices, operating on each matrix (labeled by K) separately.
C We use OLD_TX to store TX( , ,K) at the previous frequency. Only necessary
C when N is being (at least partially) specified in terms of J.
C
	DO K=1,NM_TX
	  IF(DO_THIS_TX_MATRIX(K))THEN
	    IF(USE_EPS .AND. .NOT. INIT)THEN
	      DO J=1,ND
	        DO I=1,NDEXT
	          OLD_TX(I,J)=TX(I,J,K)
	        END DO
	      END DO
	      DO J=1,ND
	        TX(1,J,K)=PSIPREV_MOD(1)*OLD_TX(1,J) + VC(1)*TVX(1,J,K)
	        DO I=2,NDEXT-1
 	          TX(I,J,K)=PSIPREV_MOD(I)*OLD_TX(I,J)
	1                 + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	1                 + ( EPS_PREV_B(I)*OLD_TX(I+1,J)
	1                       - EPS_PREV_A(I-1)*OLD_TX(I-1,J) )
	        END DO
 	        TX(NDEXT,J,K)= PSIPREV_MOD(NDEXT)*OLD_TX(NDEXT,J) +
	1                       VB(NDEXT)*TVX(NDEXT-1,J,K)
	      END DO
	    ELSE IF(.NOT. INIT)THEN
	      DO J=1,ND
	        TX(1,J,K)=PSIPREV_MOD(1)*TX(1,J,K) + VC(1)*TVX(1,J,K)
	        DO I=2,NDEXT-1
 	          TX(I,J,K)=PSIPREV_MOD(I)*TX(I,J,K)
	1                 + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	        END DO
 	        TX(NDEXT,J,K)= PSIPREV_MOD(NDEXT)*TX(NDEXT,J,K)
	1                 + VB(NDEXT)*TVX(NDEXT-1,J,K)
	      END DO
	    END IF
C
	    IF(K .EQ. 1 .OR. K .EQ. 2)THEN
	      DO J=1,ND
	        DO I=1,NDEXT
	          TX(I,J,K)=TX(I,J,K)+KI(I,J,K)
	        END DO
	      END DO
	    END IF
C
C Solve the simultaneous equations.
C
	    CALL SIMPTH(TA,TB,TC,TX(1,1,K),NDEXT,ND)
C
	    IF(USE_EPS)THEN
	      DO J=1,ND
	        DO I=1,NDEXT-1
	           TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1               + HS(I)*TVX(I,J,K) +
	1               (EPS_PREV_A(I)*OLD_TX(I,J)-EPS_A(I)*TX(I,J,K)) +
	1               (EPS_PREV_B(I)*OLD_TX(I+1,J)-EPS_B(I)*TX(I+1,J,K))
	        END DO
	      END DO
	    ELSE
	      DO J=1,ND
	        DO I=1,NDEXT-1
	           TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1               + HS(I)*TVX(I,J,K)
	        END DO
	      END DO
	    END IF
C
	    IF(K .EQ. 1)THEN
	      DO J=1,ND
	        DO I=1,NDEXT-1
	          TVX(I,J,K)=TVX(I,J,K)+RHS_dHdCHI(I,J)
	        END DO
	      END DO
	    END IF
C
	  END IF	!DO_THIS_MATRIX
	END DO		!K
 
	RETURN
	END
