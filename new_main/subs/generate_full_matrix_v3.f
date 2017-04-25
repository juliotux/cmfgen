!
! Routine designed to convert the SE%STEQ and SE%BA structures into the
! full STEQ and BA matrices that are used to find the corrections.
! The conversion is done one depth at a time. At each depth, the DIAGONAL
! matrix should be passed first.
!
	SUBROUTINE GENERATE_FULL_MATRIX_V3(C_MAT,STEQ_VEC,POPS,
	1                REPLACE,ZERO_STEQ,
	1                NT,ND,NION,NUM_BNDS,BAND_INDX,DIAG_INDX,DEPTH_INDX,
	1                FIRST_MATRIX,LAST_MATRIX,USE_PASSED_REP)
	USE MOD_CMFGEN
	USE STEQ_DATA_MOD
	USE CONTROL_VARIABLE_MOD, ONLY : LTE_MODEL
	IMPLICIT NONE
!
! Altered 13-Mar-2014 : Issues with crude electron-energy balance equation when
!                          non-thermal ionization was included.
! Altered 30-Jan-2002 : Changed to V2
!                          REPLACE and ZERO_STEQ now passed in call.
! Altered 10-Sep-2001 : Using DIAG_INDX as logical variable in an IF statement
! Created 05-Apr-2001
!
!******************************************************************************
!
	INTEGER NT
	INTEGER ND
	INTEGER NION
	INTEGER NUM_BNDS
	INTEGER BAND_INDX 
	INTEGER DIAG_INDX
	INTEGER DEPTH_INDX
!
        REAL*8 POPS(NT,ND)
        REAL*8 C_MAT(NT,NT)
	REAL*8 STEQ_VEC(NT)
	LOGICAL ZERO_STEQ(NT)		!Which vectors to be zeroed.
	LOGICAL REPLACE(NION)		!Which equations to be replaced
	LOGICAL FIRST_MATRIX		!Fist call to GENERATE_FULL_MATRIX.
	LOGICAL LAST_MATRIX		!Last call to   "       "    "
!
! When true we replace those equations indicated by the passed logical vector REPLACE.
! We always use the passed logical vector REPLACE if not updating the diagonal band.
!
	LOGICAL USE_PASSED_REP
!
! Local matrices and vectors to facilitate construction of C_MAT and STEQ_VEC.
!
! We use ?_ION to create the ionization balance equations for each
! ionization stage.
!
	REAL*8 C_ION(NION,NT)
	REAL*8 STEQ_ION(NION)
!
! We use ?_NC for the Number conservation equation for each species.
! We don't replace them directly into C_MAT and STEQ_VEC for ease of
! programming.
!
	REAL*8 C_NC(NUM_SPECIES,NT)
	REAL*8 STEQ_NC(NUM_SPECIES)
!
	REAL*8 G_SUM(NT)
	REAL*8 EDGE_SUM(NT)
	REAL*8 SUM,T1
	INTEGER NS
!
! FAC is used as a scale factor to determine at what depth the ground-state
! equilibrium equation is replaced by the ionization equation.
!
!	REAL*8, PARAMETER :: FAC=1.0D+05
	REAL*8, PARAMETER :: FAC=1.0D+02
!
	LOGICAL DIAG_BAND
!
! For consistency with the old version of CMFGEN, we only replace the
! ground-state equations over a continuous set of depths. REPLACE
! is used to indicate whether the current depth is to be replaced
! (as determined from the DIAGONAL band). Replace now passed.
!
	INTEGER, SAVE, ALLOCATABLE ::  REP_CNT(:)
!
	INTEGER, SAVE :: LST_DEPTH_INDX=0
!
	INTEGER I,J,K,L,JJ,N
	INTEGER EQ
	INTEGER ID
	INTEGER ISPEC
!
	INTEGER ERROR_LU,WARNING_LU
	INTEGER LUER,LUWARN
	EXTERNAL ERROR_LU,WARNING_LU
!
! To save typing.
!
	LUER=ERROR_LU()
	LUWARN=WARNING_LU()
        K=DEPTH_INDX
	DIAG_BAND=.FALSE.
	IF(DIAG_INDX .EQ. BAND_INDX)DIAG_BAND=.TRUE.
	IF(DEPTH_INDX .NE. LST_DEPTH_INDX)THEN
	  LST_DEPTH_INDX=DEPTH_INDX
	  IF(.NOT. DIAG_BAND)THEN
	    WRITE(LUER,*)'Error in GENERATE_FULL_MATRIX'
	    WRITE(LUER,*)'Diagonal band must be passed first at each depth'
	    WRITE(LUER,*)'DEPTH_INDX=',DEPTH_INDX
	    STOP
	  END IF
	END IF
!
	C_MAT(:,:)=0.0D0
	C_ION(:,:)=0.0D0
	C_NC(:,:)=0.0D0
!
	IF(DIAG_BAND)THEN
	  STEQ_VEC(:)=0.0D0
	  STEQ_ION(:)=0.0D0
	  STEQ_NC(:)=0.0D0
!
!	  IF(FIRST_MATRIX)THEN
!            I=SE(4)%LNK_TO_IV(1101)
!	     WRITE(99,*)I
!            WRITE(99,*)SE(1)%BA(1,I,DIAG_INDX,1),SE(1)%BA(2,I,DIAG_INDX,1)
!            WRITE(99,*)SE(3)%BA(1,I,DIAG_INDX,1),SE(3)%BA(2,I,DIAG_INDX,1)
!            WRITE(99,*)SE(4)%BA(1,I,DIAG_INDX,1),SE(4)%BA(2,I,DIAG_INDX,1)
!	  END IF
!
	END IF
!
	IF( .NOT. ALLOCATED(REP_CNT))THEN
	  ALLOCATE (REP_CNT(NION)); REP_CNT(:)=0
	END IF
!
! Map the small BA array onto the full BA array
!
! NB: C_ION(1,:) refers to to the ionization/recombination equation
!                for ion 1 (e.g. CI in the carbon sequence). It is
!                dN(CI)/dt. Since the equation (i.e. BA(I,:,:,:) with 
!                I > ATMD(ID)%NXzV refers to dN/dt for the recombining 
!                level (i.e. C2) we need a - sign when we evaluate 
!                C_ION and STEQ_ION.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	    DO J=1,SE(ID)%N_IV
	      JJ=SE(ID)%LNK_TO_F(J)
	      DO I=1,SE(ID)%N_SE-1
	        EQ=SE(ID)%EQ_IN_BA(I)
                C_MAT(EQ,JJ)=C_MAT(EQ,JJ)+SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	      DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	        C_ION(ID,JJ)=C_ION(ID,JJ)-SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	      C_NC(ISPEC,JJ)=C_NC(ISPEC,JJ)+SE(ID)%BA(SE(ID)%N_SE,J,BAND_INDX,K)
	    END DO
	  END DO
	END DO
!
! Now do STEQ
!
	IF(DIAG_BAND)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      DO I=1,SE(ID)%N_SE-1
	        EQ=SE(ID)%EQ_IN_BA(I)
                STEQ_VEC(EQ)=STEQ_VEC(EQ)+SE(ID)%STEQ(I,K)
	      END DO
	      DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	        STEQ_ION(ID)=STEQ_ION(ID)-SE(ID)%STEQ(I,K)
	      END DO
	      STEQ_NC(ISPEC)=STEQ_NC(ISPEC)+SE(ID)%STEQ(SE(ID)%N_SE,K)
	    END DO
	  END DO
	  STEQ_VEC(NT-1)=STEQ_ED(K)
	  STEQ_VEC(NT)=STEQ_T(K)
	END IF
!
! Update the ionization equations for when X-rays are included.
! The following is photoionizations/recombinations which change z by 2.
! Only other process allowed are /\z=1. The corrections to the ionization 
! equations follow from simple algebraic manipulations.
!
	DO ISPEC=1,NUM_SPECIES
	  DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-2
            I=SE(ID)%XRAY_EQ
	    IF(I .NE. 0)THEN
	      DO J=1,SE(ID)%N_IV
	        JJ=SE(ID)%LNK_TO_F(J)
	        C_ION(ID+1,JJ)=C_ION(ID+1,JJ)-SE(ID)%BA(I,J,BAND_INDX,K)
	      END DO
	      IF(DIAG_BAND)STEQ_ION(ID+1)=STEQ_ION(ID+1)-SE(ID)%STEQ(I,K)
	    END IF
	  END DO
	END DO
!
! Allow for advection terms in the ionization equations. 
!
	IF(DIAG_BAND)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      STEQ_ION(ID)=STEQ_ION(ID)+SE(ID)%STEQ_ADV(K)
	    END DO
	  END DO
	END IF
!
! Allow for advection terms in the ionization equations. Since we are using linear derivatives,
! the terms, at each depth, are identical. We only need to worry about the sign of the terms.
!
	IF(BA_ADV_TERM(DIAG_INDX,1) .NE. 0)THEN
	  DO ISPEC=1,NUM_SPECIES
	    DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	      IF(SE(ID)%STRT_ADV_ID(K) .EQ. SPECIES_BEG_ID(ISPEC))THEN
	        DO L=SE(ID)%STRT_ADV_ID(K),ID
	          DO J=1,ATM(L)%NXzV
	            JJ=SE(L)%LNK_TO_F(J)
	            C_ION(ID,JJ)=C_ION(ID,JJ)-BA_ADV_TERM(BAND_INDX,K)
	          END DO
	        END DO
	      ELSE
	        DO L=ID+1,SE(ID)%END_ADV_ID(K)
	          DO J=1,ATM(L)%NXzV
	            JJ=SE(L)%LNK_TO_F(J)
	            C_ION(ID,JJ)=C_ION(ID,JJ)+BA_ADV_TERM(BAND_INDX,K)
	          END DO
	        END DO
!
!Including DxZV for last ionization state
!
	        L=SE(ID)%END_ADV_ID(K)
	        J=ATM(L)%NXzV
	        JJ=SE(L)%LNK_TO_F(J)+1
	        C_ION(ID,JJ)=C_ION(ID,JJ)+BA_ADV_TERM(BAND_INDX,K)
	      END IF
	    END DO
	  END DO
	END IF
!
! We now replace any equations that need replacing.
!
! Insert the number conservation equation for each species.
!
	DO ISPEC=1,NUM_SPECIES
	  IF(SPECIES_PRES(ISPEC))THEN
	    I=EQ_SPECIES(ISPEC)
	    C_MAT(I,:)=C_NC(ISPEC,:)
	  END IF
	END DO
!
	DO ISPEC=1,NUM_SPECIES
	  IF(DIAG_BAND .AND. SPECIES_PRES(ISPEC))THEN
	    I=EQ_SPECIES(ISPEC)
	    STEQ_VEC(I)=STEQ_NC(ISPEC)
	  END IF
	END DO
!
! Charge conservation and radiative equilibrium equations.
!
	C_MAT(NT-1,:)=BA_ED(:,BAND_INDX,DEPTH_INDX)
	C_MAT(NT,:)  =BA_T(:,BAND_INDX,DEPTH_INDX)
!
! Crude section to create ELECTRON cooling equation.
! Must be done after all equations are done, but before GS equation is replaced.
!
! The following was omitted 13-Mar-2014. Problem with non-thermal ionization. Needs modifcation.
!
	IF(DEPTH_INDX .LE. -12)THEN
          DO ISPEC=1,NUM_SPECIES
	    IF(SPECIES_PRES(ISPEC))THEN
	      DO ID=SPECIES_BEG_ID(ISPEC),SPECIES_END_ID(ISPEC)-1
	        NS=ATM(ID)%NXzV
                G_SUM(1:NS)=0.0D0
                EDGE_SUM(1:NS)=0.0D0
                DO J=1,ATM(ID)%NXzV_F
                  I=ATM(ID)%F_TO_S_XzV(J)
                  G_SUM(I)=G_SUM(I)+ATM(ID)%GXzV_F(J)
                  EDGE_SUM(I)=EDGE_SUM(I)+ATM(ID)%EDGEXzV_F(J)*ATM(ID)%GXzV_F(J)
                END DO
!
	        IF(DIAG_BAND)THEN
	          SUM=0.0D0
                  DO I=1,ATM(ID)%NXzV
                    T1=EDGE_SUM(I)/G_SUM(I)
	            EQ=SE(ID)%EQ_IN_BA(I)
	            SUM=SUM+5.27296D-03*T1*STEQ_VEC(EQ)
	          END DO
	          STEQ_VEC(NT)=STEQ_VEC(NT)+SUM
	        END IF
!
	        DO J=1,NT
	          SUM=0.0D0
                  DO I=1,ATM(ID)%NXzV
                    T1=EDGE_SUM(I)/G_SUM(I)
	            EQ=SE(ID)%EQ_IN_BA(I)
	            SUM=SUM+5.27296D-03*T1*C_MAT(EQ,J)
	          END DO
	          C_MAT(NT,J)=C_MAT(NT,J)+SUM
	        END DO
	      END DO
	    END IF
	  END DO
	END IF
!
! Check if we wish to use the ionization conservation equation.
! We only use the diagonal band to do this. We use REPLACE
! array so that we only use the conservation equation over a
! sequence of optical depths.
!
	IF(DIAG_BAND .AND. DEPTH_INDX .EQ. 1)REP_CNT(:)=0
        IF(.NOT. USE_PASSED_REP .AND. DIAG_BAND)THEN
	  DO ID=1,NION
	    REPLACE(ID)=.FALSE.
	    IF(ATM(ID)%XzV_PRES)THEN
	      EQ=ATM(ID)%EQXzV
	      N=ATM(ID)%NXzV
              IF( ABS(C_MAT(EQ,EQ))*ATM(ID)%XzV(1,K) .GT.
	1          FAC*ABS(C_MAT(EQ,EQ+N))*ATM(ID)%DXzV(K) )THEN
	        REPLACE(ID)=.TRUE.
	        REP_CNT(ID)=REP_CNT(ID)+1
	      END IF
	    END IF
	  END DO
	END IF
!
! In all cases, we replace the ground state equation.
! We need the check on XzV_PRES since we initially set
! all the REPLACE to true.
!
	DO ID=1,NION
          IF(REPLACE(ID) .AND. ATM(ID)%XzV_PRES)THEN
	    C_MAT(ATM(ID)%EQXzV,:)=C_ION(ID,:)
	    IF(DIAG_BAND)STEQ_VEC(ATM(ID)%EQXzV)=STEQ_ION(ID)
	  END IF
	END DO
!
	IF(DIAG_BAND .AND. DEPTH_INDX .EQ. ND)THEN
	  WRITE(LUWARN,'(/,/,1X,A,/)')' Equation selection in generate_full_matrix_v3.f' 
	  DO ID=1,NION
	    IF(REP_CNT(ID) .GT. 0)THEN
              WRITE(LUWARN,'(1X,A,T9,A,I3,A)')
	1         TRIM(ION_ID(ID)),' g.s. eq. replaced by ionization eq. at ',
	1         REP_CNT(ID),' depths.'
	      END IF
	  END DO
	  WRITE(LUWARN,'(A)')' '
	END IF
!
!	IF(K .EQ. 1 .AND. DIAG_BAND)THEN
!	  CALL WR2D(STEQ_VEC,NT,1,'STEQ_VEC_D1',94)
!	  OPEN(UNIT=96,FILE='BA_ASCI_D1',STATUS='UNKNOWN')
!	    CALL WR2D(C_MAT,NT,NT,'C_MAT_VEC_D1',96)
!	  CLOSE(UNIT=96)
!	END IF
!	IF(K .EQ. 61 .AND. DIAG_BAND)THEN
!	  CALL WR2D(STEQ_VEC,NT,1,'STEQ_VEC_D10',94)
!	  OPEN(UNIT=96,FILE='BA_ASCI_D61',STATUS='UNKNOWN',POSITION='APPEND')
!	    CALL WR2D(C_MAT,NT,NT,'C_MAT_VEC_D10',96)
!	  CLOSE(UNIT=96)
!	END IF
!
! Fix any populations.
!
	IF(USE_PASSED_REP)THEN
	  DO I=1,NT
	     IF(ZERO_STEQ(I))STEQ_VEC(I)=0.0D0
	  END DO
	ELSE
	  CALL FIXPOP_IN_BA_V3(C_MAT,STEQ_VEC,ZERO_STEQ,
	1       NT,ND,NION,DIAG_BAND,K,
	1       FIRST_MATRIX,LAST_MATRIX)
	END IF
!
	IF(LTE_MODEL)THEN
	  CALL ADJUST_CMAT_TO_LTE(C_MAT,STEQ_VEC,DIAG_BAND,DEPTH_INDX,NT)
	END IF
!
! Scale the BA matrix so that we solve for the fractional corrections to
! the populations. This seems to yield better solutions for large matrices,
! especially when we subsequently condition the matrices.
!                
	L=DEPTH_INDX-DIAG_INDX+BAND_INDX
	DO J=1,NT
	  DO I=1,NT
	    C_MAT(I,J)=C_MAT(I,J)*POPS(J,L)
	  END DO
	END DO
!
	IF(K .EQ. 1 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=96,FILE='BA_ASCI_N_D1',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',96)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',96)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',96)
	  CLOSE(UNIT=96)
	END IF
	IF(K .EQ. 6 .AND. DIAG_BAND .AND. K .LE. ND)THEN
	  OPEN(UNIT=96,FILE='BA_ASCI_N_D6',STATUS='UNKNOWN')
	    CALL WR2D_MA(POPS(1,K),NT,1,'POPS_D1',96)
	    CALL WR2D_MA(STEQ_VEC,NT,1,'STEQ_VEC_D1',96)
	    CALL WR2D_MA(C_MAT,NT,NT,'C_MAT_D1',96)
	  CLOSE(UNIT=96)
	END IF
!
	RETURN
	END
