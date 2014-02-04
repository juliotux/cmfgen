!
! Subroutine to compute the value of the statistical equilibrium
! equations and the variation of the statistical equilibrium matrix for
! the advection terms: - GRAD.(v ni).
!
	SUBROUTINE STEQ_ADVEC_V3(ID,HN_S,DI_S,ADVEC_RR,N_S,R,VEL,
	1                 NUM_BNDS,ND,LAMBDA_ITERATION,COMPUTE_BA,LINEAR)
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 13-May-2004 - Included ADVEC_RR into call, and changed to V3.
! Altered 12-May-2004 - Can now use linear differencing, which will have greater
!                          convergence stability, although potentially not as accurate.
!                          Changed V2,
! Altered 05-May-2004 - Bug fix when updating BA for k=ND: MP1=M-1 not ND-1 
! Created 20-Jan-2003
!
	INTEGER ID
	INTEGER N_S
	INTEGER NUM_BNDS
	INTEGER ND
!
	REAL*8 HN_S(N_S,ND)
	REAL*8 DI_S(ND)
	REAL*8 R(ND)
	REAL*8 VEL(ND)
	REAL*8 ADVEC_RR(ND)
!
	LOGICAL LAMBDA_ITERATION
	LOGICAL COMPUTE_BA
	LOGICAL LINEAR
!
! Local variables.
!
	INTEGER K,KP1			!Depth index
	INTEGER M,MP1			!Band index
	INTEGER I			!Variable index
	INTEGER EQION
!
	REAL*8 T1,T2
	REAL*8 DERIV_FAC
	REAL*8 DERIV_CONST
	REAL*8 LOG_DERIV
!
        EQION=N_S+1			!Ion equation
	M=(NUM_BNDS/2)+1
!
! We use backward linear differencing. This should be more stable than linear
! differencing in the log-log plane.
!
! NB: Warning: We could sum ADVEC_RR over I, since we ate refering the ionizations
!     to XzV, not the ion. Different if photoioizations to multiple levels. Hoever,
!     numerial proplems could occur since HN(1,K) migh be >> DI(K).
!
! NB: The factor of 1.0D-05 arises from V/R.
!
	WRITE(6,*)'**************************'
	WRITE(6,*)'**************************'
	WRITE(6,*)'WARNNG IN STEQ_ADVEC_V3 --- reduced advection TERMS'
	WRITE(6,*)'**************************'
	WRITE(6,*)'**************************'
	DERIV_FAC=1.0D-05                         !Should be 1.0D-05
!
	IF(LINEAR)THEN
	  DO K=1,ND
	    ADVEC_RR(K)=0.0D0
	    KP1=K+1
	    IF(K .EQ. ND)KP1=ND-1
	    DERIV_CONST=DERIV_FAC/R(K)/R(K)/(R(K)-R(KP1))
	    DO I=1,N_S			!Which S.E. equationa
	      T1=R(K)*R(K)*VEL(K)*HN_S(I,K)
	      T2=R(KP1)*R(KP1)*VEL(KP1)*HN_S(I,KP1)
	      SE(ID)%STEQ(I,K)=SE(ID)%STEQ(I,K) - DERIV_CONST*(T1-T2)
	    END DO
	    T1=R(K)*R(K)*VEL(K)*DI_S(K)
	    T2=R(KP1)*R(KP1)*VEL(KP1)*DI_S(KP1)
	    SE(ID)%STEQ(EQION,K)=SE(ID)%STEQ(EQION,K)-DERIV_CONST*(T1-T2)
	    ADVEC_RR(K)=ADVEC_RR(K) + DERIV_CONST*(T1-T2)
	  END DO
!
	  IF(COMPUTE_BA)THEN
	    DO K=1,ND
	      KP1=K+1
	      MP1=M+1
	      IF(K .EQ. ND)THEN
	        KP1=ND-1
	        MP1=M-1
	      END IF
!
	      DERIV_CONST=DERIV_FAC/(R(K)-R(KP1))
	      T1= DERIV_CONST*VEL(K)
	      T2= DERIV_CONST*R(KP1)*R(KP1)*VEL(KP1)/R(K)/R(K)
	      DO I=1,N_S			!Which S.E. equation
	        SE(ID)%BA(I,I,M,K)=SE(ID)%BA(I,I,M,K)-T1
	        IF(NUM_BNDS .GE. 3 .AND. .NOT. LAMBDA_ITERATION)THEN
	          SE(ID)%BA(I,I,MP1,K)=SE(ID)%BA(I,I,MP1,K)+T2
	        END IF
	      END DO
!
! Update ion equation variation.
!
	      SE(ID)%BA(EQION,EQION,M,K)=SE(ID)%BA(EQION,EQION,M,K)-T1
	      IF(NUM_BNDS .GE. 3 .AND. .NOT. LAMBDA_ITERATION)THEN
	        SE(ID)%BA(EQION,EQION,MP1,K)=SE(ID)%BA(EQION,EQION,MP1,K)+T2
	      END IF
!
	    END DO		!loop over depth
	  END IF
	ELSE
!
! We use backward differencing in the LOG-LOG plane. This was found to be more stable
! for adiabatic cooling was included.
!
	  DO K=1,ND
	    ADVEC_RR(K)=0.0D0
	    KP1=K+1
	    IF(K .EQ. ND)KP1=ND-1
	    DERIV_CONST=DERIV_FAC*VEL(K)/R(K)
	    DO I=1,N_S			!Which S.E. equation
	      LOG_DERIV=LOG( (R(K)/R(KP1))**2 * (VEL(K)/VEL(KP1)) * (HN_S(I,K)/HN_S(I,KP1)) ) / LOG(R(K)/R(KP1))
	      SE(ID)%STEQ(I,K)=SE(ID)%STEQ(I,K)-DERIV_CONST*LOG_DERIV*HN_S(I,K)
	    END DO
	    LOG_DERIV=LOG( (R(K)/R(KP1))**2 * (VEL(K)/VEL(KP1)) * (DI_S(K)/DI_S(KP1)) ) / LOG(R(K)/R(KP1))
	    SE(ID)%STEQ(EQION,K)=SE(ID)%STEQ(EQION,K)-DERIV_CONST*LOG_DERIV*DI_S(K)
	    ADVEC_RR(K)=ADVEC_RR(K)+DERIV_CONST*LOG_DERIV*DI_S(K)
	  END DO
!
	  IF(COMPUTE_BA)THEN
	    DO K=1,ND
	      KP1=K+1
	      MP1=M+1
	      IF(K .EQ. ND)THEN
	        KP1=ND-1
	        MP1=M-1
	      END IF
!
	      DERIV_CONST=DERIV_FAC*VEL(K)/R(K)
	      DO I=1,N_S			!Which S.E. equation
	        T1= (R(K)/R(KP1))**2 * (VEL(K)/VEL(KP1)) * (HN_S(I,K)/HN_S(I,KP1))
	        T2=LOG(R(K)/R(KP1))
	        LOG_DERIV=LOG(T1)/T2
	        SE(ID)%BA(I,I,M,K)=SE(ID)%BA(I,I,M,K)-DERIV_CONST*LOG_DERIV
	        SE(ID)%BA(I,I,M,K)=SE(ID)%BA(I,I,M,K)-DERIV_CONST/T2
	        IF(NUM_BNDS .GE. 3 .AND. .NOT. LAMBDA_ITERATION)THEN
	          SE(ID)%BA(I,I,MP1,K)=SE(ID)%BA(I,I,MP1,K)+DERIV_CONST*HN_S(I,K)/T2/HN_S(I,KP1)
	        END IF
	      END DO
!
! Update ion equation variation.
!
	      T1= (R(K)/R(KP1))**2 * (VEL(K)/VEL(KP1)) * (DI_S(K)/DI_S(KP1))
	      T2=LOG(R(K)/R(KP1))
	      LOG_DERIV=LOG(T1)/T2
	      SE(ID)%BA(EQION,EQION,M,K)=SE(ID)%BA(EQION,EQION,M,K)-DERIV_CONST*LOG_DERIV
	      SE(ID)%BA(EQION,EQION,M,K)=SE(ID)%BA(EQION,EQION,M,K)-DERIV_CONST/T2
	      IF(NUM_BNDS .GE. 3 .AND. .NOT. LAMBDA_ITERATION)THEN
	        SE(ID)%BA(EQION,EQION,MP1,K)=SE(ID)%BA(EQION,EQION,MP1,K)+DERIV_CONST*DI_S(K)/T2/DI_S(KP1)
	      END IF
!
	    END DO		!loop over depth
	  END IF
	END IF			!linear or in loglog plane
!
	RETURN
	END
