!
! Routine to convolve a STARK profile with a DOPPLER profile. The Doppler
! profile is characterized by a Doppler parameter DLAM_TURB. For the
! Butler & Lemmke stark tables, this should NOT include the thermal
! contribution, since this has already been included.
!
	SUBROUTINE CONV_STRK_V2(PRO,PROF_LAM,NF,
	1                 STARK,DWS,NWS,
	1                 SYM_STARK,TRAP_QUAD,NORMALIZE_PROFILE,
	1                 DLAM_THERM,DLAM_TURB,WAVE,
	1                 ELOG,TLOG,DIAGNOSTICS)
	IMPLICIT NONE
!
! Altered 01-Sep-2007: ELOG & TLOG inserted for diagnostic purposes (changed to V2).
! Altered 24-Sep-2003:  For symmetric profiles, normalization was being done twice.
!                         This didn't effect the output, provided NORMALIZE_PROFILE
!                         was true (default for CMF_FLUX).
! Altered  5-Sep-2003:  Routine was not working correctly for large turbulent
!                         velocities (> 200 km/s). Routine may be somewhat
!                         slower.
!
	INTEGER NF
!
	REAL*8 PRO(NF)		!Should contain LOG10 of the STARK profile. 
	REAL*8 PROF_LAM(NF)	!Offset from line center in Angstroms.
!
	INTEGER NWS
	REAL*8 STARK(NWS)
	REAL*8 DWS(NWS)
	REAL*8 ELOG
	REAL*8 TLOG
	LOGICAL SYM_STARK
	LOGICAL TRAP_QUAD
	LOGICAL NORMALIZE_PROFILE
	LOGICAL DIAGNOSTICS
!
! The thermal velocity is used as an indicator of the structure in the
! passed STARK profile. If the STARK profile has NOT been convolved with
! the thermal distribution, it should characterize the structure
! variation in the RAW STARK profile.
!
! DLAM_TURB is the width due to pure TURBULENCE. If the RAW stark
! profile has NOT already been convolved with a Doppler profile, 
! the Doppler contributions should be included.
!
	REAL*8 DLAM_THERM		!Thermal velocity
	REAL*8 DLAM_TURB		!Turbulent velocity
	REAL*8 WAVE			!Central wavelength of line in Ang.
!                
! Local variables
!
	REAL*8 MAX_DWS,MIN_DWS
	REAL*8 DLAM
	REAL*8 MAX_INT_LAM
	REAL*8 MIN_INT_LAM
	REAL*8 SLOPE_RHS,SLOPE_LHS
	REAL*8 VAL_RHS,VAL_LHS
	REAL*8 T1
	REAL*8 SQRT_PI
	REAL*8 LAM_VAL
	REAL*8 STARK_VAL
!
	REAL*8, ALLOCATABLE :: DLAM_INT(:)
	REAL*8, ALLOCATABLE :: STARK_INT(:)
	SAVE DLAM_INT,STARK_INT
	REAL*8, SAVE :: WAVE_SAVE=0.0D0
!
	INTEGER NI,NG
	INTEGER I,J,K,L,IOS
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	CHARACTER(LEN=80) STRING
	LOGICAL SIMP_QUAD
!
	IF(DWS(NWS) .LT. DWS(1))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error on CONV_STRK_V2'
	  WRITE(LUER,*)'DWS should increase with pixel'
	  WRITE(LUER,*)WAVE,DWS(1),DWS(NWS)
	  STOP
	END IF
	IF(PROF_LAM(NF) .LT. PROF_LAM(1))THEN
	  LUER=ERROR_LU()
	  WRITE(LUER,*)'Error on CONV_STRK_V2'
	  WRITE(LUER,*)'PROF_LAM should increase with pixel'
	  WRITE(LUER,*)WAVE,PROF_LAM(1),PROF_LAM(NWS)
	  STOP
	END IF
!
! We can either use Simpson's rule, of the Trapazoidal rule, for the
! convolution quadrature.
!
	SIMP_QUAD=.NOT. TRAP_QUAD
	MAX_DWS=MAXVAL(DWS)
	MIN_DWS=MINVAL(DWS)
!
! We choose the minimum so as to resolve both profiles accurately.
!
	DLAM=MIN(DLAM_THERM,DLAM_TURB)/5.0D0
!
	PRO(1:NF)=0.0D0
	SQRT_PI=1.772453851D0
!
! Determine profile limits over which we perform the interpolation. These
! are set by the smaller of the tabulated range, and the range required
! to determine the requested Stark profile.
!
	MAX_INT_LAM=MIN( MAX_DWS, PROF_LAM(NF)+5.2*DLAM_TURB )
	IF(SYM_STARK)THEN
	  MIN_INT_LAM=0.0D0
	ELSE
	  MIN_INT_LAM=MAX( MIN_DWS, PROF_LAM(1)-5.2*DLAM_TURB )
	END IF
	NI=(MAX_INT_LAM-MIN_INT_LAM)/DLAM+1
!
! Allocate working vectors, if needed.
!
!	IF(ALLOCATED(DLAM_INT) .AND. .LT. NI)THEN
	IF(ALLOCATED(DLAM_INT))THEN
	  DEALLOCATE(DLAM_INT)
	  DEALLOCATE(STARK_INT)
	END IF
	IF(.NOT. ALLOCATED(DLAM_INT))THEN
	  ALLOCATE(DLAM_INT(NI),STAT=IOS)
	  IF(IOS .EQ. 0)ALLOCATE(STARK_INT(NI),STAT=IOS)
	  IF(IOS .NE. 0)THEN
	    LUER=ERROR_LU()
	    WRITE(LUER,*)'Error on CONV_STRK_V2'
	    WRITE(LUER,*)'Unable to allocate DLAM_INT and STARK_INT'
	    STOP
	  END IF
	END IF
!
! Tabulated the wavelength array, equally spaced in wavelength.
!
	DO I=1,NI
	  DLAM_INT(I)=MIN_INT_LAM+(I-1)*DLAM
	END DO
	DLAM_INT(NI)=MAX_INT_LAM
	CALL MON_INTERP_FAST(STARK_INT,NI,1,DLAM_INT,NI,STARK,NWS,DWS,NWS)
!
! Since LOG10 of profile is tabulated.
!
	STARK_INT(1:NI)=10**(STARK_INT(1:NI))
!
! Generate values for extrapolation beyond tabulated profile limits.
!
	VAL_RHS=10.0D0**STARK(NWS)
	SLOPE_RHS=(STARK(NWS)-STARK(NWS-1)) / LOG10(DWS(NWS)/DWS(NWS-1)) 
	IF(.NOT. SYM_STARK)THEN
	  VAL_LHS=10.0D0**STARK(1)
	  SLOPE_LHS=(STARK(1)-STARK(2)) / LOG10(DWS(1)/DWS(2))
	END IF
!
! NG is the number of points (on each side) used to evaluate the
! Gaussian profile for the convolution.
!
	NG=5.0D0*DLAM_TURB/DLAM
!
! Can now perform the convolution. The method is dependent on whether the 
! tabulated profile is symmetric. NB: Unfortunately the wavelengths
! at which the profile is to be determined are generally not symmetric
! about line center.
!
! The normalization of the profile is done after the full profile is 
! computed.
!
	IF(SYM_STARK)THEN
!
! Perform the convolution quadrature. We first find the location of the current
! frequency in interpolation array. If convolution bandpass extends beyond
! tabulated range, we extrapolate profile in log-log plane using a power
! law. Constants common to the quadrature are included separately.
!
	  DO J=1,NF
	    K=NINT((ABS(PROF_LAM(J))-DLAM_INT(1))/DLAM)+1
	    DO I=K-NG,K+NG
	      LAM_VAL=(I-1)*DLAM
	      L=I
	      IF(I .LE. 0)L=1-I
	      IF(L .LE. NI)THEN
	        STARK_VAL=STARK_INT(L)
	      ELSE
	        STARK_VAL=VAL_RHS*(LAM_VAL/DWS(NWS))**SLOPE_RHS
	      END IF
	      IF(SIMP_QUAD .AND. MOD(I-K,2) .EQ. 0)STARK_VAL=STARK_VAL*2.0D0
	      PRO(J)=PRO(J)+STARK_VAL*EXP( -( (LAM_VAL-ABS(PROF_LAM(J)))/DLAM_TURB )**2 )
	    END DO
	  END DO
!
	ELSE 			!Not symmetric.
!
! Can now perform the quadrature. We first find the location of the current
! frequency in interpolation array. If convolution bandpass extends beyond
! tabulated range, we extrapolate profile in log-log plane using a power
! law.
!
	  DO J=1,NF
	    K=NINT((PROF_LAM(J)-DLAM_INT(1))/DLAM)+1
	    DO I=K-NG,K+NG
	      LAM_VAL=DLAM_INT(1)+(I-1)*DLAM
	      IF(I .LT. 1)THEN
	        STARK_VAL=VAL_LHS*(LAM_VAL/DWS(1))**SLOPE_LHS
	      ELSE IF(I .GT. NI)THEN
	        STARK_VAL=VAL_RHS*(LAM_VAL/DWS(NWS))**SLOPE_RHS
	      ELSE
	        STARK_VAL=STARK_INT(I)
	      END IF
	      IF(SIMP_QUAD .AND. MOD(I-K,2) .EQ. 0)STARK_VAL=STARK_VAL*2.0D0
	      PRO(J)=PRO(J)+STARK_VAL*EXP( -( (LAM_VAL-PROF_LAM(J))/DLAM_TURB )**2 )
	    END DO
	  END DO
	END IF		!Symmetric / anti-symmetric profile
!
! Perform the quadrature normalization.
!
	T1=DLAM/SQRT_PI/DLAM_TURB
	IF(SIMP_QUAD)T1=T1/1.5D0
	PRO(1:NF)=PRO(1:NF)*T1
!
! Normalize the profile to have unit area if desired.
!
	IF(NORMALIZE_PROFILE)THEN
	  T1=0.0D0
	  DO I=1,NF-1
	    T1=T1+( 1.0D0/(PROF_LAM(I)+WAVE)-1.0D0/(PROF_LAM(I+1)+WAVE) )*
	1         (PRO(I+1)+PRO(I))
	  END DO
	  T1=T1*0.5D0*2.99794D+18
	  PRO(1:NF)=PRO(1:NF)/T1
	  IF(ABS(T1-1.0D0) .GT. 0.3)THEN
            LUER=ERROR_LU()
	    IF(WAVE .NE. WAVE_SAVE)THEN
	      BACKSPACE(LUER,IOSTAT=IOS)
	      IF(IOS .EQ. 0)THEN
	        READ(LUER,'(A)')STRING
	        BACKSPACE(LUER)
	        WRITE(LUER,'(/,A)')TRIM(STRING)
	      END IF
	      WRITE(LUER,*)'Possible error in CONV_STRK_V2'
	      WRITE(LUER,*)'Profile normalization constant differs from 1 by more than 30%'
	      IF(WAVE .GT. 100.0D0 .AND. WAVE .LT. 1.0D+05)THEN
	        WRITE(LUER,'(3(A,F15.8,A,3X))')' Wave=',WAVE,'A','Lam_ST(1)=',PROF_LAM(1),'A',
	1                      'Lam_END(NF)=',PROF_LAM(NF),'A'
	      ELSE
	        WRITE(LUER,'(3(A,ES16.8,A,3X))')' Wave=',WAVE,'A','Lam_ST(1)=',PROF_LAM(1),'A',
	1                      'Lam_END(NF)=',PROF_LAM(NF),'A'
	      END IF
	      WAVE_SAVE=WAVE
	      WRITE(LUER,'(5A14)')'Norm Const.','    Log(Ne)',
	1                      '     Log(T)',' dLAM_THERM',' dLAM_TURB'
	    END IF
	    WRITE(LUER,'(5ES14.4)')T1,ELOG,TLOG,DLAM_THERM,DLAM_TURB
	  END IF
	END IF
!
	IF(DIAGNOSTICS)THEN
	  DLAM_INT(1:NI)=DLAM_INT(1:NI)+WAVE
	  DLAM_INT(1:NI)=2.998D+18/DLAM_INT(1:NI)
	  T1=0.0D0
	  DO I=1,NI-1
	    T1=T1+0.5D0*(DLAM_INT(I)-DLAM_INT(I+1))*(STARK_INT(I+1)+STARK_INT(I))
	  END DO
	  IF(SYM_STARK)T1=T1*2.0D0
	  WRITE(6,*)'The area under the INTERP profile is',T1
	  DLAM_INT(1:NI)=2.998D+18/DLAM_INT(1:NI)
	  DLAM_INT(1:NI)=DLAM_INT(1:NI)-WAVE
	  STARK_INT(1:NI)=DLOG10(STARK_INT(1:NI))
	  CALL DP_CURVE(NI,DLAM_INT,STARK_INT)
	END IF
!
	RETURN
	END
