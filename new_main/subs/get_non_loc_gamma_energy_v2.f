!
! Subroutine to read in non-local rdiactive energy deposition.
! The data should be computed by a separate calculation, and is
! read in from:
!                current_nonlocal_decay_energy.dat
!
	SUBROUTINE GET_NON_LOCAL_GAMMA_ENERGY_V2(R,V,ND,LU)
	USE CONTROL_VARIABLE_MOD
	USE NUC_ISO_MOD
	IMPLICIT NONE
!
! Altered 29-Jan-2014 - Call DO_GAM_ABS_APPROX_V2
!                         Compute emitted and absorbed decau luminosities.
! Altered 05-Jan-2014 - R added to call.
!                         Now output energy emitted & absorbed to check_edep.
!                         Changed to V2.
!
	INTEGER ND
	INTEGER LU
	REAL*8 R(ND)
	REAL*8 V(ND)
!
! Variables for reading in the non-local energy deposition from outside file
!
	INTEGER NDTMP
	REAL*8, ALLOCATABLE :: VTMP(:),EDEPTMP(:),EDEPNEW(:)
	REAL*8 T1,T2
	REAL*8 CONV_CONST
	REAL*8 WRK(ND)
	INTEGER LUER
	INTEGER I
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
!
	LUER=ERROR_LU()
        IF (GAMRAY_TRANS .EQ. 'LOCAL') THEN
!
        ELSE IF (GAMRAY_TRANS.EQ.'ABS_TRANS') THEN
!
! Before the call, WRK must contain the locally EMITTED energy. After the call,
! RADIOACTIVE_DECAY_ENERGY will contain the locally ABSORBED energy.
! 
	   WRK=RADIOACTIVE_DECAY_ENERGY
	   CALL DO_GAM_ABS_APPROX_V2(RADIOACTIVE_DECAY_ENERGY,WRK,KINETIC_DECAY_ENERGY,ND)
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
!
!
! The conversion factor:
!   (a) 4pi r^2 dr --> 4pi .10^30 (as R is in units of 10^10 cm).
!   (b) The extra factor of 4 arises as ATAN(1.0D0) is pi/4.
!   (c) /Lsun to convert to solar luminosities.
!
	   CONV_CONST=16.0D0*ATAN(1.0D0)*1.0D+30/3.826D+33
	   DO I=1,ND
	     WRK(I)=RADIOACTIVE_DECAY_ENERGY(I)*R(I)*R(I)
	   END DO
	   CALL LUM_FROM_ETA(WRK,R,ND)
	   T2=SUM(WRK)*CONV_CONST
	   WRITE(LU,'(A)')'!'
	   WRITE(LU,'(A,ES13.5,A)')'!            Radioactive energy emitted is: ',T2,' Lsun'
	   DO I=1,ND
	     WRK(I)=EDEPNEW(I)*R(I)*R(I)
	   END DO
	   CALL LUM_FROM_ETA(WRK,R,ND)
	   T1=SUM(WRK)*CONV_CONST
	   WRITE(LU,'(A,ES13.5,A)')'!            Radioactive energy absorbed is:',T1,' Lsun'
	   WRITE(LU,'(A,ES13.5,A)')'!Fraction of radioactive energy absorbed is:',T1/T2
	   WRITE(LU,'(A)')'!'
!
! Switched ordering of EDEPNEW and RADIOACTIVE_DECAY_ENERGY to be consistent
! with that in 'current_nonlocal_decay_energy.dat'
!
	   WRITE(LU,'(A,12X,A,5X,A,5X,A)')'!','Velocity','Non-Local Decay','    local Decay'
	   WRITE(LU,'(A,12X,A,5X,A,5X,A)')'!','  km/s  ','    ergs/cm^3/s','   ergs/cm^3/s'
	   DO I=1,ND
	      WRITE(LU,'(3ES20.8)')V(I),EDEPNEW(I),RADIOACTIVE_DECAY_ENERGY(I)
	   ENDDO
!
	   WRITE(LU,'(/,/,I5,A50)') NDTMP,' !Number of depth points in MC computation'
	   DO I=1,NDTMP
	      WRITE(LU,'(3ES20.8)') VTMP(I),EDEPTMP(I)
	   ENDDO
	   WRITE(LU,'(A)')'!'
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
