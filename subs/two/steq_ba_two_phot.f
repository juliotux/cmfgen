!
! Routine to increment the Statistical Equilibriaum Equations and
! the Variation equation.
!
	SUBROUTINE STEQ_BA_TWO_PHOT_RATE(STEQ,BA,POPS,NT,ND,NUM_BNDS,
	1             DIAG_INDX,UPDATE_BA,LU_OUT,WRITE_RATES)
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Created 26-Jun-1998
!
	INTEGER*4 NT,ND
	INTEGER*4 NUM_BNDS,DIAG_INDX
!
	REAL*8 STEQ(NT,ND)			!Statistical Equilibrium equations
	REAL*8 POPS(NT,ND)
	REAL*8 BA(NT,NT,NUM_BNDS,ND)		!Variation of STEQ
!
	LOGICAL UPDATE_BA
!
	INTEGER*4 LU_OUT	!Unit for summary of 2-photon rates.
	LOGICAL WRITE_RATES	!Indicates that rates should be written.
!
	INTEGER*4 LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER*4, PARAMETER :: IZERO=0
	INTEGER*4 L,J,II,IOS
	INTEGER*4 NL,NUP
	REAL*8 TA(ND)
	REAL*8 T1
!
	DO L=1,ND
	  DO J=1,N_TWO
	    IF(TWO_PHOT_AVAILABLE(J))THEN
	      NL=LOW_LEV_TWO(J)
	      NUP=UP_LEV_TWO(J)
	      T1=POPS(NUP,L)*DOWN_RATE_TWO(L,J)-POPS(NL,L)*UP_RATE_TWO(L,J)
	      STEQ(NL,L)=STEQ(NL,L)+T1
	      STEQ(NUP,L)=STEQ(NUP,L)-T1
	    END IF
	  END DO
	END DO
C
	IF(UPDATE_BA)THEN
	  II=DIAG_INDX
	  DO L=1,ND
	    DO J=1,N_TWO
	      IF(TWO_PHOT_AVAILABLE(J))THEN
	        NL=LOW_LEV_TWO(J)
	        NUP=UP_LEV_TWO(J)
	        BA(NL,NL,II,L)=BA(NL,NL,II,L)-UP_RATE_TWO(L,J)
	        BA(NUP,NL,II,L)=BA(NUP,NL,II,L)+UP_RATE_TWO(L,J)
	        BA(NL,NUP,II,L)=BA(NL,NUP,II,L)+DOWN_RATE_TWO(L,J)
	        BA(NUP,NUP,II,L)=BA(NUP,NUP,II,L)-DOWN_RATE_TWO(L,J)
	      END IF
	    END DO
	  END DO
	END IF
!
	IF(WRITE_RATES)THEN
	  CALL GEN_ASCI_OPEN(LU_OUT,'TWO_PHOT_SUM','UNKNOWN',
	1                                        ' ',' ',IZERO,IOS)
	  IF(IOS .NE. 0)THEN
	     LUER=ERROR_LU()
	     WRITE(LUER,*)'Unable to open TWO_PHOT_SUM in',
	1                   ' STEQ_BA_TWO_PHOT'
	     INITIALIZE_TWO=.TRUE.
	     RETURN
	  END IF
!
	  DO J=1, N_TWO
	    NL=LOW_LEV_TWO(J)
	    NUP=UP_LEV_TWO(J)
	    WRITE(LU_OUT,*)' '
	    WRITE(LU_OUT,*)' '
	    IF(TWO_PHOT_AVAILABLE(J))THEN
	      WRITE(LU_OUT,'(X,I4,3X,L1,A,3X,A,3X,A)')
	1             J,.TRUE.,TRIM(SPEC_ID_TWO(J)),
	1             TRIM(LOW_NAME_TWO(J)),TRIM(UP_NAME_TWO(J))
	      WRITE(LU_OUT,*)'Down rate:'
!
	      CALL WRITV(DOWN_RATE_TWO(1,J),ND,'Down rate',LU_OUT)
	      CALL WRITV(UP_RATE_TWO(1,J),ND,'Up rate',LU_OUT)
!
	      TA(1:ND)=DOWN_RATE_TWO(1:ND,J)*POPS(NUP,1:ND)
	      CALL WRITV(TA,ND,'Total # of down transitions',LU_OUT)
	      TA(1:ND)=UP_RATE_TWO(1:ND,J)*POPS(NL,1:ND)
	      CALL WRITV(TA,ND,'Total # of upward transitions',LU_OUT)
	    ELSE
	      WRITE(LU_OUT,'(X,I4,3X,L1,A,3X,A,3X,A)')
	1             J,.FALSE.,TRIM(SPEC_ID_TWO(J)),
	1             TRIM(LOW_NAME_TWO(J)),TRIM(UP_NAME_TWO(J))
	    END IF
	  END DO
	END IF
!
! Indicates that STEQ and BA arrays have been updated. Thus ths rates
! can be zeroed. Done in SET_TWO_PHOT.
!
	INITIALIZE_TWO=.TRUE.
!
	RETURN
	END
