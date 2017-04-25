!
! Simple routine to force the components of the variation matrix (C_MAT) to assume
! that we are dealing with LTE populations. In this fudge, we assume the departure
! coefficients, T and Ne are the independent variables. It assumed that we have LTE
! populations, and that the departure coeficient will remain at unit. Doing an initial
! lambda iteration (with LTE_MODEL set to true) will force the populations to be LTE.
! NB: More complicated approaches could be implmented.
!
	SUBROUTINE ADJUST_CMAT_TO_LTE(C_MAT,STEQ_VEC,DIAG_MAT,DEPTH_INDX,NT)
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
!
! Finalized: 17-Feb-2016
!
	INTEGER NT
	INTEGER DEPTH_INDX
	REAL*8 C_MAT(NT,NT)
	REAL*8 STEQ_VEC(NT)
	LOGICAL DIAG_MAT
!
	INTEGER ISPEC
	INTEGER ID
	INTEGER I
	INTEGER EQ
	INTEGER ION_EQ
!
	IF(DIAG_MAT)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,ATM(ID)%NXzV
	        EQ=SE(ID)%EQ_IN_BA(I)
	        C_MAT(EQ,1:NT)=0.0D0
	      END DO
	    END DO
	  END DO
	  STEQ_VEC(1:NT-1)=0.0D0                !As forcing LTE and number conservation.
	ELSE
	  C_MAT(1:NT-1,1:NT)=0.0D0
	END IF
!
	IF(DIAG_MAT)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      ION_EQ=SE(ID)%EQ_IN_BA(1)+ATM(ID)%NXzV
	      DO I=1,ATM(ID)%NXzV
	        EQ=SE(ID)%EQ_IN_BA(I)
	        C_MAT(EQ,EQ)=1.0D0
	        C_MAT(EQ,ION_EQ)=-ATM(ID)%XzVLTE(I,DEPTH_INDX)/ATM(ID)%DXzV(DEPTH_INDX)
	        C_MAT(EQ,NT-1)=-ATM(ID)%XzVLTE(I,DEPTH_INDX)/ED(DEPTH_INDX)
	        C_MAT(EQ,NT)=-ATM(ID)%XzVLTE(I,DEPTH_INDX)*ATM(ID)%dlnXzVLTE_dlnT(I,DEPTH_INDX)/T(DEPTH_INDX)
	      END DO
	    END DO
	  END DO
	END IF
!
	RETURN
	END
