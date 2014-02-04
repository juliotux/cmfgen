!
! Routine to evaluate the Penning contribution to the de-excitation rate of
! the highly meta-stable state, 2s_3Se of HeI.
!
!    HeI(1s_2s_3Se) + H(1s_2Se)  --->  HeI(1s2_1Se) +  H+ + e-
!
	SUBROUTINE DO_PENNING_ION(HDKT,COMPUTE_BA,DIAG_INDX,ND)
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	IMPLICIT NONE
!
! Altered: 09-Jan-2012: Bug fixed with linearization.
!
	INTEGER ND
	INTEGER DIAG_INDX
	REAL*8 HDKT
	LOGICAL COMPUTE_BA
!
	REAL*8, PARAMETER :: ALPHA=7.5D-10
	REAL*8, PARAMETER :: CI=2.07078D-22
!
	REAL*8 RATE
	REAL*8 REV_RATE
	REAL*8 dREV_RATEdT
	REAL*8 REV_RATE_COEF
	REAL*8 T1,T2
!
	INTEGER ID1		!Ion ID for H
	INTEGER ID2		!Ion ID for HeI
	INTEGER ION_EQ_H	!Equation for H+ in H SEEs.
	INTEGER IV_ED_HI        !Equation for Ne in H SEEs
	INTEGER IV_ED_HeI	!Equation for Ne in HeI SEEs
	INTEGER IV_T_HI         !Equation for T in H SEEs
	INTEGER IV_T_HeI	!Equation for T in HeI SEEs
	INTEGER IV
	INTEGER L
	INTEGER IQ		!Used as shorthand notation for ION_EQ_H
	INTEGER IED		!Used as shorthand notation for ED equations (ion dependent)
	INTEGER IT              !Used as shorthand for T equation
	INTEGER II		!Used as shorthand notationm for DIAG_INDX
!
	DO ID1=1,NUM_IONS
	  IF(ION_ID(ID1) .EQ. 'HI')THEN
	    DO ID2=1,NUM_IONS
	      IF(ION_ID(ID2) .EQ. 'HeI')THEN
!
	        ION_EQ_H=ATM(ID1)%NXzV+1
	        DO L=1,ND
	          RATE=ALPHA*ATM(ID1)%XzV_F(1,L)*ATM(ID2)%XzV_F(2,L)
	          T1=HDKT*(ATM(ID1)%EDGEXzV_F(1)+ATM(ID2)%EDGEXzV_F(2)-ATM(ID2)%EDGEXzV_F(1))/T(L)
	          REV_RATE_COEF=ALPHA*6.0D0*CI*EXP(T1)
	          REV_RATE=REV_RATE_COEF*ED(L)*ATM(ID1)%DXzV(L)*ATM(ID2)%XzV_F(1,L)
!
	          T2=RATE-REV_RATE
	          SE(ID1)%STEQ(1,L)=SE(ID1)%STEQ(1,L)-T2
	          SE(ID1)%STEQ(ION_EQ_H,L)=SE(ID1)%STEQ(ION_EQ_H,L)+T2
!
	          SE(ID2)%STEQ(2,L)=SE(ID2)%STEQ(2,L)-T2		!HeI(2s 3S) state.
	          SE(ID2)%STEQ(1,L)=SE(ID2)%STEQ(1,L)+T2		!HeI ground term.
	        END DO
!
	        IF(COMPUTE_BA)THEN
	          IV_ED_HI=SE(ID1)%N_IV-1
	          IV_T_HI=SE(ID1)%N_IV
	          IV_ED_HeI=SE(ID2)%N_IV-1
	          IV_T_HeI=SE(ID2)%N_IV
	          II=DIAG_INDX
	          DO L=1,ND
	            RATE=ALPHA*ATM(ID1)%XzV_F(1,L)*ATM(ID2)%XzV_F(2,L)
	            T1=HDKT*(ATM(ID1)%EDGEXzV_F(1)+ATM(ID2)%EDGEXzV_F(2)-ATM(ID2)%EDGEXzV_F(1))/T(L)
	            REV_RATE_COEF=ALPHA*6.0D0*CI*EXP(T1)
	            REV_RATE=REV_RATE_COEF*ED(L)*ATM(ID1)%DXzV(L)*ATM(ID2)%XzV_F(1,L)
	            dREV_RATEdT=-REV_RATE*T1/T(L)
!
	            IV=SE(ID1)%LNK_TO_IV(ATM(ID2)%EQXzV+1)
	            SE(ID1)%BA(1,1,II,L) =SE(ID1)%BA(1,1,II,L) -ALPHA*ATM(ID2)%XzV_F(2,L)
	            SE(ID1)%BA(1,IV,II,L)=SE(ID1)%BA(1,IV,II,L)-ALPHA*ATM(ID1)%XzV_F(1,L)
!
	            IQ=ION_EQ_H
	            IED=IV_ED_HI
	            IT=IV_T_HI
	            IV=SE(ID1)%LNK_TO_IV(ATM(ID2)%EQXzV)
	            SE(ID1)%BA(1,IV,II,L) =SE(ID1)%BA(1,IV,II,L) +REV_RATE_COEF*ED(L)*ATM(ID1)%DXzV(L)
	            SE(ID1)%BA(1,IQ,II,L) =SE(ID1)%BA(1,IQ,II,L) +REV_RATE_COEF*ED(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID1)%BA(1,IED,II,L)=SE(ID1)%BA(1,IED,II,L)+REV_RATE_COEF*ATM(ID1)%DXzV(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID1)%BA(1,IT,II,L) =SE(ID1)%BA(1,IT,II,L) +dREV_RATEdT
!
	            SE(ID1)%BA(IQ,1,II,L)   =SE(ID1)%BA(IQ,1,II,L)   +ALPHA*ATM(ID2)%XzV_F(2,L)
	            SE(ID1)%BA(IQ,IV+1,II,L)=SE(ID1)%BA(IQ,IV+1,II,L)+ALPHA*ATM(ID1)%XzV_F(1,L)
!
	            SE(ID1)%BA(IQ,IV,II,L) =SE(ID1)%BA(IQ,IV,II,L) -REV_RATE_COEF*ED(L)*ATM(ID1)%DXzV(L)
	            SE(ID1)%BA(IQ,IQ,II,L) =SE(ID1)%BA(IQ,IQ,II,L) -REV_RATE_COEF*ED(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID1)%BA(IQ,IED,II,L)=SE(ID1)%BA(IQ,IED,II,L)-REV_RATE_COEF*ATM(ID1)%DXzV(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID1)%BA(IQ,IT,II,L) =SE(ID1)%BA(IQ,IT,II,L) -dREV_RATEdT
!
! Now do HeI equations.
!
	            IV=SE(ID2)%LNK_TO_IV(ATM(ID1)%EQXzV)
	            SE(ID2)%BA(2,2,II,L) =SE(ID2)%BA(2,2,II,L) -ALPHA*ATM(ID1)%XzV_F(1,L)
	            SE(ID2)%BA(2,IV,II,L)=SE(ID2)%BA(2,IV,II,L)-ALPHA*ATM(ID2)%XzV_F(2,L)
	            SE(ID2)%BA(1,2,II,L) =SE(ID2)%BA(1,2,II,L) +ALPHA*ATM(ID1)%XzV_F(1,L)
	            SE(ID2)%BA(1,IV,II,L)=SE(ID2)%BA(1,IV,II,L)+ALPHA*ATM(ID2)%XzV_F(2,L)
!
	            IV=SE(ID2)%LNK_TO_IV(ATM(ID1)%EQXzV+ATM(ID1)%NXzV-1)
	            IED=IV_ED_HeI
	            IT=IV_T_HeI
	            SE(ID2)%BA(2,1,II,L)  =SE(ID2)%BA(2,1,II,L)  +REV_RATE_COEF*ED(L)*ATM(ID1)%DXzV(L)
	            SE(ID2)%BA(2,IV,II,L) =SE(ID2)%BA(2,IV,II,L) +REV_RATE_COEF*ED(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID2)%BA(2,IED,II,L)=SE(ID2)%BA(2,IED,II,L)+REV_RATE_COEF*ATM(ID1)%DXzV(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID2)%BA(2,IT,II,L) =SE(ID2)%BA(1,IT,II,L) +dREV_RATEdT
!
	            SE(ID2)%BA(1,1,II,L)  =SE(ID2)%BA(1,1,II,L)  -REV_RATE_COEF*ED(L)*ATM(ID1)%DXzV(L)
	            SE(ID2)%BA(1,IV,II,L) =SE(ID2)%BA(1,IV,II,L) -REV_RATE_COEF*ED(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID2)%BA(1,IED,II,L)=SE(ID2)%BA(1,IED,II,L)-REV_RATE_COEF*ATM(ID1)%DXzV(L)*ATM(ID2)%XzV_F(1,L)
	            SE(ID2)%BA(1,IT,II,L) =SE(ID2)%BA(1,IT,II,L) -dREV_RATEdT
!
	          END DO
	        END IF
	        EXIT
	      END IF
	    END DO
	    EXIT
	  END IF
	END DO
!
	RETURN
	END
