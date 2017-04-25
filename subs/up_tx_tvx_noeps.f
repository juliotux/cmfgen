!
! Routine to compute the matrices TX and TVX which depend on the integration
! at the previous frequency.
!
! TX(, , )  describes the variation of J
! TVX(, , ) describes the variation of RSQH
!
! It is assumed that TA,TB and TC (the tridiagonal matrix to be inverted) have
! already been modified by a call to THOMAS.
!
! KI is a 2 dimensional matrix with
!
!       KI( , ,1) variation of transfer equation w.r.t. CHI
!       KI( , ,2) variation of transfer equation w.r.t. ETA
!
! Upon entry TX( , ,1)  and TX( , ,2) should be zero, as these now reflect
! change in J with respect to CHI and ETA at the current frequency.
! The other matrices are used to reflect changes in J at an earlier frequency.
!
! In general K=1 denotes dCHI and K=2 denotes dETA.
!
	SUBROUTINE UP_TX_TVX_NOEPS(TX,TVX,KI,TA,TB,TC,PSIPREV_MOD,
	1                       VB,VC,HU,HL,HS,RHS_dHdCHI,
	1                       OLD_TX,ND,NM_TX,NM_KI,
	1                       INIT,DO_THIS_TX_MATRIX)
	IMPLICIT NONE
!
! Altered 02-Oct-2016 : Based on UP_TX_TVX
!                       EPS variables deleted
!                       Clean, and THOMAS solution explicitly included
!                       in the code.
!
! Altered 28-May-1996 : Calls to DP_ZERO removed.
!                       (, at end of subroutine specification deleted).
!
! Created 10-Mar-1995.   Based on UPTX_J_EDD_V3
!                        EPS_A,EP_B etc passed in call.
!                        Argument odering altered.
!                        TX and TVX modified in the same routine, since
!                          TVX may depend on TX at the previus frequency.
!
	INTEGER ND,NM_TX,NM_KI
	REAL*8 TX(ND,ND,NM_TX)
	REAL*8 TVX(ND-1,ND,NM_TX)
	REAL*8 KI(ND,ND,NM_KI)
	REAL*8 TA(ND),TB(ND),TC(ND)
	REAL*8 PSIPREV_MOD(ND),VB(ND),VC(ND)
	REAL*8 HU(ND),HL(ND),HS(ND)
	REAL*8 RHS_dHdCHI(ND-1,ND)
!
	LOGICAL INIT,DO_THIS_TX_MATRIX(NM_TX)
!
! Work Array.
!
	REAL*8 OLD_TX(ND,ND)
!
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
! Local varables.
!
	INTEGER I,J,K
!
	IF(NM_TX .LT. 2 .OR. NM_KI .LT. 2)THEN
	  I=ERROR_LU()
	  WRITE(I,*)'Invalid NM in UP_TX_TVX'
	  WRITE(I,*)'NM_TX=',NM_TX
	  WRITE(I,*)'NM_KI=',NM_KI
	  STOP
	END IF
!
! INIT will be true for the very first frequency. We initialize all storage
! locations, even those not in use.
!
	IF(INIT)THEN
	  TX(:,:,:)=0.0D0       !ND,ND,NM_TX
	  TVX(:,:,:)=0.0D0      !(ND-1),ND*NM_TX
	  OLD_TX(:,:)=0.0D0     !ND,ND
	END IF
!
! Now modify the matrices, operating on each matrix (labeled by K) separately.
! We use OLD_TX to store TX( , ,K) at the previous frequency. Only necessary
! when N is being (at least partially) specified in terms of J.
!
!$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(I,J)
	DO K=1,NM_TX
	  IF(DO_THIS_TX_MATRIX(K))THEN
	    IF(.NOT. INIT)THEN
	      DO J=1,ND
	        TX(1,J,K)=PSIPREV_MOD(1)*TX(1,J,K) + VC(1)*TVX(1,J,K)
	        DO I=2,ND-1
 	          TX(I,J,K)=PSIPREV_MOD(I)*TX(I,J,K)
	1                 + VB(I)*TVX(I-1,J,K) + VC(I)*TVX(I,J,K)
	        END DO
 	        TX(ND,J,K)= PSIPREV_MOD(ND)*TX(ND,J,K) + VB(ND)*TVX(ND-1,J,K)
	      END DO
	    END IF
!
	    IF(K .EQ. 1 .OR. K .EQ. 2)THEN
	      DO J=1,ND
	        DO I=1,ND
	          TX(I,J,K)=TX(I,J,K)+KI(I,J,K)
	        END DO
	      END DO
	    END IF
!
! Solve the simultaneous equations.
!
!
!	    CALL SIMPTH(TA,TB,TC,TX(1,1,K),ND,ND)
!
! Forward elimination.
!
	    DO J=1,ND
!
	      TX(1,J,K)=TX(1,J,K)*TB(1)
              DO I=2,ND
                TX(I,J,K)=(TX(I,J,K)-TA(I)*TX(I-1,J,K))*TB(I)
              END DO
!
! Perform the back substitution. NB D(N1,J)=D(N1,J) is first step.
!
              DO I=ND-1,1,-1
                TX(I,J,K)=TX(I,J,K)+TC(I)*TX(I+1,J,K)
              END DO
!
	      DO I=1,ND-1
	        TVX(I,J,K)= HU(I)*TX(I+1,J,K) - HL(I)*TX(I,J,K)
	1            + HS(I)*TVX(I,J,K)
	      END DO
!
	      IF(K .EQ. 1)THEN
	        DO I=1,ND-1
	          TVX(I,J,K)=TVX(I,J,K)+RHS_dHdCHI(I,J)
	        END DO
	      END IF
!
	    END DO
	  END IF	!DO_THIS_MATRIX
	END DO		!K
!$OMP END PARALLEL DO
! 
	RETURN
	END
