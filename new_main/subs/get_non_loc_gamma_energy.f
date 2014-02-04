!
! Subroutine to read in non-local rdiactive energy deposition.
! The data should be computed by a separate calculation, and is
! read in from:
!                current_nonlocal_decay_energy.dat
!
	SUBROUTINE GET_NON_LOCAL_GAMMA_ENERGY(V,ND,LU)
	USE CONTROL_VARIABLE_MOD
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
	INTEGER ND
	INTEGER LU
	REAL*8 V(ND)
!
! Variables for reading in the non-local energy deposition from outside file
!
	INTEGER NDTMP
	REAL*8, ALLOCATABLE :: VTMP(:),EDEPTMP(:),EDEPNEW(:)
	REAL*8 T1
	INTEGER LUER
	INTEGER I
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
        IF (GAMRAY_TRANS .EQ. 'LOCAL') THEN
!
        ELSE IF (GAMRAY_TRANS.EQ.'ABS_TRANS') THEN
	   CALL DO_GAM_ABS_APPROX(RADIOACTIVE_DECAY_ENERGY,ND)
!
	ELSE IF (GAMRAY_TRANS.EQ.'NONLOCAL') THEN
!
	   CLOSE(LU)
	   OPEN(LU,FILE='current_nonlocal_decay_energy.dat',STATUS='UNKNOWN')
!
	   READ(LU,*) NDTMP
	   READ(LU,*) T1
	   IF (T1 .NE. SN_AGE_DAYS) THEN
	      WRITE(LUER,'(A80)') 'New time for Gamma-ray transport calculation is incompatible'
	      WRITE(LUER,'(A50)') 'Monte Carlo new time [days]',T1
	      WRITE(LUER,'(A50)') 'CMFGEN new time [days]',SN_AGE_DAYS
	      STOP
	   ENDIF
!
	   ALLOCATE (VTMP(NDTMP),EDEPTMP(NDTMP),EDEPNEW(ND))
	   DO I=1,NDTMP
	      READ(LU,*) VTMP(I),EDEPTMP(I)
	   ENDDO
	   CLOSE(LU)
!
	   VTMP(NDTMP) = V(ND)
	   CALL LIN_INTERP(V,EDEPNEW,ND,VTMP,EDEPTMP,NDTMP)
!
	   OPEN(LU,FILE='check_edep.dat',STATUS='UNKNOWN')
	   WRITE(LU,'(I5,A50)') ND,' !Number of depth points'
	   WRITE(LU,'(12X,A,5X,A,5X,A)')'Velocity','    Local Decay','Non-local Decay'
	   WRITE(LU,'(12X,A,5X,A,5X,A)')'  km/s  ','    ergs/cm^3/s','   ergs/cm^3/s'
	   DO I=1,ND
	      WRITE(LU,'(3ES20.8)')V(I),RADIOACTIVE_DECAY_ENERGY(I),EDEPNEW(I)
	   ENDDO
	   WRITE(LU,'(/,/,I5,A50)') NDTMP,' !Number of depth points in MC computation'
	   DO I=1,NDTMP
	      WRITE(LU,'(3ES20.8)') VTMP(I),EDEPTMP(I)
	   ENDDO
!
	   CLOSE(LU)
!
! We overwrite the edep computed before assuming local deposition
!
	   RADIOACTIVE_DECAY_ENERGY = EDEPNEW
	   DEALLOCATE (VTMP,EDEPTMP,EDEPNEW)
!
	ELSE
	   WRITE(LUER,*)'Error in GET_NON_LOCAL_GAMMA_ENERGY'
	   WRITE(LUER,*)'Unrecognized gamma-ray transport option'
	   WRITE(LUER,*)GAMRAY_TRANS
	   STOP
	ENDIF
!
	RETURN
	END
