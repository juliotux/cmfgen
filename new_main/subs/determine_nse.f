!
! Subroutine to determine the number of statistical equilibrium equations
! for each species. 
!
      SUBROUTINE DETERMINE_NSE(NION,XRAYS)
      USE MOD_CMFGEN
      USE STEQ_DATA_MOD
      USE CHG_EXCH_MOD_V3
      IMPLICIT NONE
!
! Altered 16-May-2003 : Small bug fix for so that charge exchange reactions that
!                          are not available are corectly handled. Only important
!                          if unvailable if species present but level name 
!                          doesn't match.
! Created 17-Mar-2001
!
      INTEGER NION            !Total number of ioization stages.
      LOGICAL XRAYS
!
! Local variables
!
      INTEGER ID
      INTEGER I
      INTEGER J
      INTEGER K
      INTEGER MK
      INTEGER L
      LOGICAL IONIZING_LEVEL_FOUND
!
      INTEGER, ALLOCATABLE :: ION_LEV_PNT(:)
!
      WRITE(93,*)' '
      WRITE(93,*)'Summary of photoiization routes'
      WRITE(93,*)'Eqn. in Spec BA ',' Lev in ion ',' Equation #'
      WRITE(93,*)' '
!
      DO ID=1,NION
!
	IF(ATM(ID)%XzV_PRES)THEN
	  ALLOCATE (ION_LEV_PNT(ATM(ID)%NXzV+30))
	  ION_LEV_PNT(:)=0
!
          SE(ID)%XzV_PRES=.TRUE.
!
! Set the total number of equations for this ION. This is
! equal to:
!             NXzV  (i.e. one for each level in the atom)
!            +1     (number conservation equation)
!            +NI    (number of states atom can ionize to).
! NI requires some thought, since it is set by the 
! number of photoionization routes, X-ray ionization,
! and charge exchange reactions.
!
          SE(ID)%N_SE=ATM(ID)%NXzV+ATM(ID)%N_XzV_PHOT
	  DO J=1,ATM(ID)%N_XzV_PHOT
	    ION_LEV_PNT(ATM(ID)%NXzV+J)=ATM(ID)%XzV_ION_LEV_ID(J)
	    WRITE(93,'(1X,A,1X,I2,A)')TRIM(ION_ID(ID)),ATM(ID)%XzV_ION_LEV_ID(J),
	1                         ' photoionization route'
	  END DO
!
! Now do the charge exchange rections. Since charge reactions have the form
!
! We only need to consider species 2 and 3.
!
          IF(DO_CHG_EXCH)THEN
	    DO K=2,3
	      IF(K .EQ. 2)MK=4
	      IF(K .EQ. 3)MK=1
	      DO J=1,N_CHG
	        IF(CHG_REACTION_AVAILABLE(J) .AND. ID_ION_CHG(J,K) .EQ. ID)THEN
	          IONIZING_LEVEL_FOUND=.FALSE.
	          DO L=ATM(ID)%NXzV+1,SE(ID)%N_SE
	             IF(ION_LEV_PNT(L) .EQ. LEV_IN_ION_CHG(J,MK))THEN
	               IONIZING_LEVEL_FOUND=.TRUE.
	               EXIT
	             END IF
	          END DO
	          IF(.NOT. IONIZING_LEVEL_FOUND)THEN
	             SE(ID)%N_SE=SE(ID)%N_SE+1
	             L=SE(ID)%N_SE
	          END IF
	          ION_LEV_PNT(L)=LEV_IN_ION_CHG(J,MK)
                END IF
	      END DO
	    END DO
          END IF
!
! ION_LEV_PNT refers to the level in the next ionization.
! For Xrays we want the ground state, 2 stages up.
!
	   IF(XRAYS .AND. ATM(ID+1)%XzV_PRES)THEN
	      SE(ID)%N_SE=SE(ID)%N_SE+1 
	      SE(ID)%XRAY_EQ=SE(ID)%N_SE
	      ION_LEV_PNT(SE(ID)%N_SE)=ATM(ID+1)%NXzV+1
	   ELSE
	      SE(ID)%XRAY_EQ=0
	   END IF
!
! Finally add in the number conservation equation, which is always the
! last equation.
!
	   SE(ID)%N_SE=SE(ID)%N_SE+1 
!
	   ALLOCATE (SE(ID)%EQ_TO_ION_LEV_PNT(SE(ID)%N_SE))
	   SE(ID)%EQ_TO_ION_LEV_PNT(1:SE(ID)%N_SE)=ION_LEV_PNT(1:SE(ID)%N_SE)
	   K=MAXVAL(ION_LEV_PNT)
	   ALLOCATE (SE(ID)%ION_LEV_TO_EQ_PNT(K))
	   SE(ID)%ION_LEV_TO_EQ_PNT(:)=0
	   DO I=ATM(ID)%NXzV+1,SE(ID)%N_SE-1
	     SE(ID)%ION_LEV_TO_EQ_PNT(ION_LEV_PNT(I))=I
	     WRITE(93,*)SE(ID)%ION_LEV_TO_EQ_PNT(ION_LEV_PNT(I)),ION_LEV_PNT(I),I
	   END DO
	   DEALLOCATE (ION_LEV_PNT)
!
	END IF		!ionization stage is present
      END DO		!loop over ID
!
      RETURN
      END    
