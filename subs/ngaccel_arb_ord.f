!
! Subroutine to perform an NG acceleration on estimates obtained
! using an operator which has linear convergence. The method
! is described in detail by Auer (p101, Numerical Radiative Transfer).
!
! IF the logical variable WEIGHT is true,  weighting inversely proportional
! to the value is used.
!
	SUBROUTINE NGACCEL_ARB_ORD(RJ,PREVRJ,ND,NORD,USE_WEIGHTING)
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER NORD
	REAL*8 RJ(ND)
	REAL*8 PREVRJ(0:NORD+1,ND)
	LOGICAL USE_WEIGHTING
!
	REAL*8 CMAT(NORD,NORD)
	REAL*8 RHS(NORD)
	REAL*8 WEIGHT(ND)
	REAL*8 T1,T2
	INTEGER I,J,K
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	WEIGHT=1.0D0
	DO I=1,ND
          IF(USE_WEIGHTING .AND. PREVRJ(0,I) .NE. 0.0D0)WEIGHT(I)=1.0D0/PREVRJ(0,I)
	END DO
!
	RHS=0.0D0
	DO I=1,ND
	  T1=PREVRJ(0,I)-PREVRJ(1,I)
	  DO J=1,NORD
	    RHS(J)=RHS(J)+WEIGHT(I)*T1*(T1-(PREVRJ(J,I)-PREVRJ(J+1,I)))
	  END DO
	END DO
!
	CMAT=0.0D0
	DO I=1,ND
	  T1=PREVRJ(0,I)-PREVRJ(1,I)
	  DO K=1,NORD
	    T2=(PREVRJ(K,I)-PREVRJ(K+1,I)-T1)
	    DO J=1,NORD
	       CMAT(J,K)=CMAT(J,K)+WEIGHT(I)*T2*(PREVRJ(J,I)-PREVRJ(J+1,I)-T1)
	    END DO
	  END DO
	END DO
!
	CALL SIMQ(CMAT,RHS,NORD,K)
	IF( K .NE. 0 )THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Warning - Singular determinant in NGACCEL'
	  RJ(1:ND)=PREVRJ(0,1:ND)
	  RETURN
	END IF
!
	DO I=1,ND
	  RJ(I)=PREVRJ(0,I)
	  DO J=1,NORD
	    RJ(I)=RJ(I)+RHS(J)*(PREVRJ(J,I)-PREVRJ(0,I) )
	  END DO
	END DO
!
	RETURN
	END
