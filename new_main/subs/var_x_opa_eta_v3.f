!
! Subroutine to compute the contribution to the opacity AND emissivity
! by K shell ionization.
!
	SUBROUTINE VAR_X_OPA_ETA_V3(VCHI,VETA,
	1              HN_A,HNST_A,dlnHNST_AdlnT,N_A,
	1              HN_B,HNST_B,dlnHNST_BdlnT,N_B,
	1              ED,DI,T,IMP_VAR,
	1              EQ_A,EQION,AT_NO,Z_A,
	1              NU,EMHNUKT,NT,ND)
	IMPLICIT NONE
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
! Local constants.
!
	INTEGER I,J,LEV
	REAL*8 ALPHA,LTE_POP,NO_ELEC
	REAL*8 TCHI1,TETA2,TETA3
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
	IF(NO_ELEC .GT. 3)THEN
	  DO I=1,N_A
	    LEV=EQ_A+I-1
            IF( IMP_VAR(LEV) )THEN
	      DO J=1,ND
	        LTE_POP=HNST_A(I,J)*HNST_B(1,J)/HN_B(1,J)
	        VCHI(LEV,J)=VCHI(LEV,J)+ALPHA
	        TCHI1=ALPHA*LTE_POP*EMHNUKT(J)
	        VCHI(EQION,J)=VCHI(EQION,J)-TCHI1/DI(J)
	        VCHI(NT-1,J)=VCHI(NT-1,J)-2.0*TCHI1/ED(J)
	        VCHI(NT,J)=VCHI(NT,J)-TCHI1*
	1          (HDKT*NU/T(J)+dlnHNST_AdlnT(I,J)+dlnHNST_BdlnT(1,J))/T(J)
!
	        TETA3=TETA2*LTE_POP*EMHNUKT(J)
	        VETA(EQION,J)=VETA(EQION,J)+TETA3/DI(J)
	        VETA(NT-1,J)=VETA(NT-1,J)+2.0*TETA3/ED(J)
	        VETA(NT,J)=VETA(NT,J)+TETA3*
	1        (HDKT*NU/T(J)+dlnHNST_AdlnT(I,J)+dlnHNST_BdlnT(1,J))/T(J)
	      END DO
	    END IF
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
