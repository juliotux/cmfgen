!
! Routine to place X-ray bound-free edge frequencies into a vector.
!
	SUBROUTINE SET_X_FREQ_V2(FREQ,NCF,NCF_MAX,
	1                   MAX_CONT_FREQ,ZCORE,ZION,
	1                   NI_PRES,N2_PRES)
	USE XRAY_DATA_MOD
	IMPLICIT NONE
!
! Altered 20-Mar-2014: Bug fix -- WARN_LU and ERROR_LU were not necessarily defined
!                         after the first entry.
! Altered 24-Jan-2014: WARNING_LU installed
! Altered 22-Dec-2004: Error message improved.
! Created 23-Oct-2000: Bsed on SET_X_FRQ
!
	INTEGER NCF,NCF_MAX
	REAL*8 MAX_CONT_FREQ		!Units 10^15 Hz
	REAL*8 FREQ(NCF_MAX)
	REAL*8 ZCORE,ZION
	LOGICAL NI_PRES,N2_PRES 
!
	INTEGER ERROR_LU,LUER,WARNING_LU,LUWARN
	EXTERNAL ERROR_LU,WARNING_LU
!
	REAL*8 CONV_FAC
	INTEGER I,J,K,LOOP_PQN_MAX
	INTEGER IZ,NE
	LOGICAL, SAVE :: FIRST=.TRUE.
	LOGICAL, SAVE :: FIRST_PQN=.TRUE.
!
	IF( .NOT. NI_PRES)RETURN
	IF( .NOT. N2_PRES)RETURN
!
	CONV_FAC=0.24191				!ev to 10^15 Hz
	IZ=ZCORE
	NE=ZCORE-ZION+1
!
	IF(IZ .GT. 30)THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error in XCROSS_V2'
	  WRITE(LUER,*)'No X-ray data available for species'//
	1                     ' with Atomic No.  > 30'
	  STOP
	END IF
!
! We only wish to return photoionization edges for the core
! (i.e. non-valence) electrons. This section of the code should
! be consistent with that in XCROSS_V2,
!
	DO J=0,ANG_MAX_X
	  LOOP_PQN_MAX=3
	  IF(J .GT. 1 .AND. NE .LT. 20)LOOP_PQN_MAX=2
	  IF(NE .LT. 18)LOOP_PQN_MAX=2
	  IF(NE .LE. 10)LOOP_PQN_MAX=1
	  IF(NE .LE. 2)RETURN
	  DO I=1,LOOP_PQN_MAX
	    IF( SIG_0_X(IZ,NE,I,J) .NE. 0)THEN
              IF(NCF+1 .GT. NCF_MAX)THEN
	        LUER=ERROR_LU()
	        WRITE(LUER,*)'Error in SET_X_FREQ --- NCF_MAX too small'
	        STOP
	      END IF
	      NCF=NCF+1
	      FREQ(NCF)=E_THRESH_X(IZ,NE,I,J)*CONV_FAC
!
! If IZ > 20, we only output error to OUTGEN when PQN > 1. For most of
! CMFGEN's applications, don't need to consider X-ray iozation from
! n=1 when IZ > 20.
!
              IF(FREQ(NCF) .GT. MAX_CONT_FREQ)THEN
	        K=0
	        LUWARN=WARNING_LU()
	        LUER=ERROR_LU()
	        IF(IZ .GT. 20 .AND. I .EQ. 1)THEN
	          IF(FIRST_PQN)THEN 
	            K=LUWARN
	            FIRST_PQN=.FALSE.
	          END IF
	        ELSE IF(FIRST)THEN
	          K=LUER
	          FIRST=.FALSE.
	        END IF
	        IF(K .NE. 0)THEN
	          WRITE(K,*)' '
	          WRITE(K,*)'*************** Warning -- Warning -- Warning ****************'
	          WRITE(K,*)'Max. cont. freq may be too small in in SET_X_FREQ_V2'
	          WRITE(K,*)'Max. cont. should generally be set to 1000 when X-rays present'
	          WRITE(K,*)'Need to allow for ionization from inner shells. Generally can'
	          WRITE(K,*)'ignore ionization from n=1 (=PQN) state of iron group elements'
	          WRITE(K,*)'since these can also ionize from n=2 sate. A list of effected'
	          WRITE(K,*)'ionization routes follows:'
	          WRITE(K,*)
	        END IF
	        IF(IZ .GT. 20 .AND. I .EQ. 1)THEN
	          WRITE(LUWARN,'(1X,4(A4,I2,3X),3X,A,ES10.2)')' IZ=',IZ,' NE=',NE,'PQN=',I,
	1                     'ANG=',J,'Edge Freq(10^15 Hz)=',FREQ(NCF)
                ELSE
	          WRITE(LUER,'(1X,4(A4,I2,3X),3X,A,ES10.2)')' IZ=',IZ,' NE=',NE,'PQN=',I,
	1                     'ANG=',J,'Edge Freq(10^15 Hz)=',FREQ(NCF)
	        END IF
	        NCF=NCF-1
	      END IF
	    END IF
	  END DO
	END DO
!
	RETURN
	END
