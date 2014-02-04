!
! Routine to increment the Statistical Equilibriaum Equations and
! the Variation equation.
!
	SUBROUTINE STEQ_BA_TWO_PHOT_RATE_V3(POPS,NT,ND,
	1             DIAG_INDX,UPDATE_BA,LU_OUT,WRITE_RATES)
	USE STEQ_DATA_MOD 
	USE TWO_PHOT_MOD
	IMPLICIT NONE
!
! Created 26-Jun-1998
!
	INTEGER NT
	INTEGER ND
	INTEGER DIAG_INDX
!
	REAL*8 POPS(NT,ND)
!
	LOGICAL UPDATE_BA
!
	INTEGER LU_OUT	!Unit for summary of 2-photon rates.
	LOGICAL WRITE_RATES	!Indicates that rates should be written.
!
	INTEGER LUER,ERROR_LU
	EXTERNAL ERROR_LU
!
	INTEGER, PARAMETER :: IZERO=0
	INTEGER L,J,II,IOS,ID
	INTEGER NL,NUP
	INTEGER ION_NL,ION_NUP
	REAL*8 TA(ND)
	REAL*8 T1
!
	DO L=1,ND
	  DO J=1,N_TWO
	    IF(TWO_PHOT_AVAILABLE(J))THEN
	            ID=ION_ID_TWO(J)
	        ION_NL=ION_LOW_LEV_TWO(J)
 	       ION_NUP=ION_UP_LEV_TWO(J)
	            NL=LOW_LEV_TWO(J)
	           NUP=UP_LEV_TWO(J)
!
	      T1=POPS(NUP,L)*DOWN_RATE_TWO(L,J)-POPS(NL,L)*UP_RATE_TWO(L,J)
	      SE(ID)%STEQ(ION_NL,L)=SE(ID)%STEQ(ION_NL,L)+T1
	      SE(ID)%STEQ(ION_NUP,L)=SE(ID)%STEQ(ION_NUP,L)-T1
	    END IF
	  END DO
	END DO
C
	IF(UPDATE_BA)THEN
	  II=DIAG_INDX
	  DO L=1,ND
	    DO J=1,N_TWO
	      IF(TWO_PHOT_AVAILABLE(J))THEN
	        ID=ION_ID_TWO(J)
	        NL=LOW_LEV_TWO(J)
	        NUP=UP_LEV_TWO(J)
	        ION_NL=ION_LOW_LEV_TWO(J)
	        ION_NUP=ION_UP_LEV_TWO(J)
		SE(ID)%BA(ION_NL,ION_NL,II,L)=SE(ID)%BA(ION_NL,ION_NL,II,L)-UP_RATE_TWO(L,J)
	        SE(ID)%BA(ION_NUP,ION_NL,II,L)=SE(ID)%BA(ION_NUP,ION_NL,II,L)+UP_RATE_TWO(L,J)
	        SE(ID)%BA(ION_NL,ION_NUP,II,L)=SE(ID)%BA(ION_NL,ION_NUP,II,L)+DOWN_RATE_TWO(L,J)
	        SE(ID)%BA(ION_NUP,ION_NUP,II,L)=SE(ID)%BA(ION_NUP,ION_NUP,II,L)-DOWN_RATE_TWO(L,J)
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
	      WRITE(LU_OUT,'(1X,I4,3X,L1,A,3X,A,3X,A)')
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
	      WRITE(LU_OUT,'(1X,I4,3X,L1,A,3X,A,3X,A)')
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
