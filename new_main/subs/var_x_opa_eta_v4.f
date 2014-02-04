!
! Subroutine to compute the contribution to the opacity AND emissivity
! by K shell ionization. The K (& L) shell cross-sections are
! assumed to be independent of the level of the valence electron. In practice,
! ionizations will generally be determined by the population of the ground
! configuration.
!
	SUBROUTINE VAR_X_OPA_ETA_V4(VCHI,VETA,
	1              HN_A,HNST_A,dlnHNST_AdlnT,N_A,
	1              HN_B,HNST_B,dlnHNST_BdlnT,N_B,
	1              ED,DI,T,IMP_VAR,
	1              EQ_A,EQION,AT_NO,Z_A,
	1              NU,EMHNUKT,NT,ND,LST_DEPTH_ONLY)
	IMPLICIT NONE
!
! Altered 19-May-2002 : Changed to version V4
!                       LST_DEPTH_ONLY inserted to save time in DTDR computation.
!                       Rewritten to sum over population variable in order to save
!                           time.
!
! Altered 07-May-2001 : Inserted IMP_VAR vector in call. Only compute
!                         variation if IMP_VAR(?)=.TRUE.
! Altered 26-Oct-1995 : dlnHNST_... passed in call instead of edge.
!                         Made version 2.
! Altered 22-Jul-1994 : Extensive modifications and testing.
! Created 20-Jul-1993
!
	EXTERNAL XCROSS_V2
!
	INTEGER ND,NT,N_A,N_B,EQ_A,EQION
	REAL*8 VCHI(NT,ND),VETA(NT,ND)
	REAL*8 HN_A(N_A,ND),HNST_A(N_A,ND),dlnHNST_AdlnT(N_A,ND)
	REAL*8 HN_B(N_B,ND),HNST_B(N_B,ND),dlnHNST_BdlnT(N_B,ND)
	REAL*8 ED(ND),DI(ND),T(ND),EMHNUKT(ND),NU
	REAL*8 AT_NO,Z_A
	LOGICAL IMP_VAR(NT)
	LOGICAL LST_DEPTH_ONLY
!
! Functions called.
!
	REAL*8 XCROSS_V2
	INTEGER ERROR_LU
!
! Constants for opacity etc.
!
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
!
	REAL*8 LTE_POP_SUM(ND)
	REAL*8 dLTE_SUM_VEC(ND)
	REAL*8 LTE_POP
	REAL*8 dLTE_SUM
!
! Local constants.
!
	INTEGER I,J,LEV
	INTEGER J_D_ST
	REAL*8 ALPHA,NO_ELEC
	REAL*8 TCHI1,TETA2,TETA3
!
	INTEGER, PARAMETER :: IZERO=0
	LOGICAL, PARAMETER :: L_FALSE=.FALSE.
	LOGICAL ERROR_OUTPUT
	DATA ERROR_OUTPUT/.FALSE./
!
	NO_ELEC=AT_NO-Z_A+1
!
! We only include the variation terms if the ION and the level are regared as
! important variables. T and ED are always regarded as important, and are therefore
! not checked.
!                            
! Add in BOUND-FREE contributions (if any). XCROSS_V2 must return 0
! if frequency is to low to cause ionizations. The cross-section
! is assumed to be independent of the level of the valence electron.
!
	IF( .NOT. IMP_VAR(EQION) )RETURN
	ALPHA=XCROSS_V2(NU,AT_NO,NO_ELEC,IZERO,IZERO,L_FALSE,L_FALSE)
	IF(ALPHA .LE. 0.0D0)RETURN
	TETA2=ALPHA*TWOHCSQ*(NU**3)
!
! Set initial depth location, in case we are just computing for the
! outermost depth point.
!
	J_D_ST=1
	IF(LST_DEPTH_ONLY)J_D_ST=ND
!
	IF(NO_ELEC .GT. 3)THEN
!
! LTE_POP_SUM is the sum over all levels.
! dLTE_SUM_VE is the of dHNST_AdlnT.
!
	  LTE_POP_SUM(1:ND)=0.0D0
	  dLTE_SUM_VEC(1:ND)=0.0D0
!
!$OMP PARALLEL DO PRIVATE(I,LEV,J)
	  DO I=1,N_A
	    LEV=EQ_A+I-1
            IF( IMP_VAR(LEV) )THEN
	      DO J=J_D_ST,ND
	        VCHI(LEV,J)=VCHI(LEV,J)+ALPHA
	        LTE_POP_SUM(J)=LTE_POP_SUM(J)+HNST_A(I,J)
	        dLTE_SUM_VEC(J)=dLTE_SUM_VEC(J)+HNST_A(I,J)*dlnHNST_AdlnT(I,J)
	      END DO
	    END IF
	  END DO
!
! If LTE_POP_SUM is zero, we have no IMPORTANT variables.
!
	  IF(LTE_POP_SUM(ND) .EQ. 0)RETURN
!
	  DO J=J_D_ST,ND
!
	    LTE_POP=(EMHNUKT(J)*LTE_POP_SUM(J))*HNST_B(1,J)/HN_B(1,J)
!
! Convert to dlnHNST_AdlnT. The factor LTE_POP_SUM in inclued
! in LTE_POP later.
!
	    dLTE_SUM=dLTE_SUM_VEC(J)/LTE_POP_SUM(J)
!
	    TCHI1=ALPHA*LTE_POP
	    VCHI(EQION,J)=VCHI(EQION,J)-TCHI1/DI(J)
	    VCHI(NT-1,J)=VCHI(NT-1,J)-2.0D0*TCHI1/ED(J)
	    VCHI(NT,J)=VCHI(NT,J)-TCHI1*(HDKT*NU/T(J)+dLTE_SUM+dlnHNST_BdlnT(1,J))/T(J)
!
	    TETA3=TETA2*LTE_POP
	    VETA(EQION,J)=VETA(EQION,J)+TETA3/DI(J)
	    VETA(NT-1,J)=VETA(NT-1,J)+2.0D0*TETA3/ED(J)
	    VETA(NT,J)=VETA(NT,J)+TETA3*(HDKT*NU/T(J)+dLTE_SUM+dlnHNST_BdlnT(1,J))/T(J)
	  END DO
	ELSE		!Variation for K shell ionization of Li ions.
	  IF(.NOT. ERROR_OUTPUT)THEN
	    WRITE(ERROR_LU(),*)'**************************************'
	    WRITE(ERROR_LU(),*)'General K shell ionization for Lithium'//
	1              'ions is not yet treated in VAR_X_OPA_ETA'
	  END IF
	  ERROR_OUTPUT=.TRUE.
	END IF
!
	RETURN
	END
