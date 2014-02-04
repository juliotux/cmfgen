!
! Subroutine to compute the value of the statistical equilibrium
! equations and the variation of the statistical equilibrium matrix for
! the advection terms: - GRAD.(v ni).
!
	SUBROUTINE STEQ_ADVEC_V4(RELAXATION_PARAMETER,LINEAR,
	1                 NUM_BNDS,ND,
	1                 INCL_ADVECTION,LAMBDA_ITERATION,COMPUTE_BA)
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered 18-May-2004 - Major rewrite: Assume linear only, and modifiy computation of
!                          ion terms for ionizaton/recombination equations.
! Created 20-Jan-2003
!
	REAL*8 RELAXATION_PARAMETER
!
	INTEGER NUM_BNDS
	INTEGER ND
!
	LOGICAL LAMBDA_ITERATION
	LOGICAL COMPUTE_BA
	LOGICAL LINEAR
	LOGICAL INCL_ADVECTION 
!
! Local variables.
!
	REAL*8 SUM(NUM_IONS,ND)
	REAL*8 T1,T2
	REAL*8 DERIV_CONST
	REAL*8 UNIT_CONST
!
	INTEGER K,KP1			!Depth index
	INTEGER M,MP1			!Band index
	INTEGER I			!Variable index
	INTEGER J			!Variable index
	INTEGER ID
	INTEGER ID_STRT,ID_END
	INTEGER ISPEC
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
	LOGICAL DO_OFF_DIAG
!
	M=(NUM_BNDS/2)+1
	DO_OFF_DIAG=.FALSE.
	IF(NUM_BNDS .GE. 3 .AND. .NOT. LAMBDA_ITERATION)DO_OFF_DIAG=.TRUE.
!
	SUM(:,:)=0.0D0
	DO ID=1,NUM_IONS
	  SE(ID)%STEQ_ADV(:)=0.0D0
	  SE(ID)%STRT_ADV_ID(:)=0.0D0
	  SE(ID)%END_ADV_ID(:)=0.0D0
	END DO
	BA_ADV_TERM(:,:)=0.0D0
	IF(.NOT. INCL_ADVECTION)RETURN
!
	IF(.NOT. LINEAR)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in STEQ_ADVEC_V4'
	  WRITE(LUER,*)'Only linear option for advection terms is curently installed.'
	  LINEAR=.TRUE.
	END IF
!
! The relaxation factor should be < 1, and is used to adjust the importance of the
! advection terms. It should be 1 for the final model. It should only be used to
! help converge a model in which advection terms are very important.
!
! NB: The factor of 1.0D-05 arises from V/R.
!
	UNIT_CONST=1.0D-05*RELAXATION_PARAMETER
!
! We use backward linear differencing. This should be more stable than linear
! differencing in the log-log plane.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    DO K=1,ND
	      KP1=K+1
	      IF(K .EQ. ND)KP1=ND-1
	      DERIV_CONST=UNIT_CONST/R(K)/R(K)/(R(K)-R(KP1))
	      DO I=1,ATM(ID)%NXzV                                !Which S.E. equationa
	        T1=R(K)*R(K)*V(K)*ATM(ID)%XzV(I,K)
	        T2=R(KP1)*R(KP1)*V(KP1)*ATM(ID)%XzV(I,KP1)
	        SE(ID)%STEQ(I,K)=SE(ID)%STEQ(I,K) - DERIV_CONST*(T1-T2)
	        SUM(ID,K)=SUM(ID,K)+DERIV_CONST*(T1-T2)
	      END DO
	      IF(ID .EQ. SPECIES_END_ID(ISPEC)-1)THEN
	        T1=R(K)*R(K)*V(K)*ATM(ID)%DXzV(K)
	        T2=R(KP1)*R(KP1)*V(KP1)*ATM(ID)%DXzV(KP1)
	        SUM(ID+1,K)=SUM(ID+1,K)+DERIV_CONST*(T1-T2)
	      END IF
	    END DO
	  END DO
	END DO
!
! We now have to compute the advection terms for the ion equations. This is complicated
! because of stability issues. Want to use the smallest terms possible.
! We have this option since Sum[all j] Xj = X (X=species population).
!
! Let Xi refer to the total population of ionization stage i. Then
! 
!      Sum[j=1,i]      {vdXj/dr} = RRi  or
!
!      Sum[j=i+1,..]  {-vdXj/dr} = RRi
!
! where RRI refers to the net recombination (phot+col) to ionization stage i.
! We choose the equation for which Sum Xj is smallest.
!
	DO ISPEC=1,NUM_SPECIES
	  ID_STRT=SPECIES_BEG_ID(ISPEC)
	  ID_END=SPECIES_END_ID(ISPEC)
	  DO K=1,ND
	    DO ID=ID_STRT,ID_END-1
	      T1=0.0D0
	      DO I=ID_STRT,ID
	        DO J=1,ATM(I)%NXzV
	          T1=T1+ATM(I)%XzV(J,K)
	        END DO
	      END DO
	      T2=0.0D0
	      DO I=ID+1,ID_END-1
	        DO J=1,ATM(I)%NXzV
	          T2=T2+ATM(I)%XzV(J,K)
	        END DO
	      END DO
	      T2=T2+ATM(ID_END-1)%DXzV(K)
	      IF(T1 .GT. T2)THEN
	        DO I=ID+1,ID_END
		  SE(ID)%STEQ_ADV(K)=SE(ID)%STEQ_ADV(K)+SUM(I,K)
	        END DO
	        SE(ID)%STRT_ADV_ID(K)=ID+1
	        SE(ID)%END_ADV_ID(K)=ID_END-1          !Don't count ion.
	      ELSE
	        DO I=ID_STRT,ID
                  SE(ID)%STEQ_ADV(K)=SE(ID)%STEQ_ADV(K)-SUM(I,K)
	        END DO
	        SE(ID)%STRT_ADV_ID(K)=ID_STRT
	        SE(ID)%END_ADV_ID(K)=ID
	      END IF
	    END DO
	  END DO
	END DO
!
! We define the diagonal & off diagonal terms so that we add or subtract both of them.
!
	DO K=1,ND-1
	  BA_ADV_TERM(M,K)=UNIT_CONST*V(K)/(R(K)-R(K+1))
	  IF(DO_OFF_DIAG)BA_ADV_TERM(M+1,K)=-UNIT_CONST*V(K+1)*R(K+1)*R(K+1)/(R(K)-R(K+1))/R(K)/R(K)
	END DO
	BA_ADV_TERM(M,ND)=UNIT_CONST*V(ND)/(R(ND)-R(ND-1))
	IF(DO_OFF_DIAG)BA_ADV_TERM(M-1,K)=-UNIT_CONST*V(ND-1)*R(ND-1)*R(ND-1)/(R(ND)-R(ND-1))/R(ND)/R(ND)
!
	IF(COMPUTE_BA)THEN
	  DO ISPEC=1,NUM_SPECIES
	    ID_STRT=SPECIES_BEG_ID(ISPEC)
	    ID_END=SPECIES_END_ID(ISPEC)
	    DO ID=ID_STRT,ID_END-1
	      DO K=1,ND
	        KP1=K+1
	        MP1=M+1
	        IF(K .EQ. ND)THEN
	          KP1=ND-1
	          MP1=M-1
	        END IF
	        DERIV_CONST=UNIT_CONST/(R(K)-R(KP1))
	        T1= DERIV_CONST*V(K)
	        T2= DERIV_CONST*R(KP1)*R(KP1)*V(KP1)/R(K)/R(K)
	        DO I=1,ATM(ID)%NXzV			!Which S.E. equation
	          SE(ID)%BA(I,I,M,K)=SE(ID)%BA(I,I,M,K)-T1
	          IF(DO_OFF_DIAG)SE(ID)%BA(I,I,MP1,K)=SE(ID)%BA(I,I,MP1,K)+T2
	        END DO
	      END DO		!loop over depth
	    END DO
	  END DO
	END IF
!
	RETURN
	END
