!
! Routine to compute limits on the intrinsic line profile:
!
!RETURNED:
!        VEC_STRT_FREQ - Frequency at which intrinsic line absorption profile
!                          is assumed to start. Profile limits are assumed
!                          to be symmetric about line center.
!        VEC_VDOP_MIN  - Minimum Doppler width for line
!	
!INPUT:
!        ED_IN     - Electron density (/cm^3) (Vector, length ND)
!        TEMP_IN   - Temperature (10^4 K)     (Vector, length ND)
!        VTURB_IN  - Turbulent velocity in km/s (function of depth - length ND).
!        CHIL      - Line opacity             (Vector, length ND)
!
!        PROF_TYPE - Indicate type of intrinsic line profile; Profile types
!                      curently implemented are:
!                                               DOPPLER
!                                               VOIGT
!                                               STARK (Hydrogenic species only)
!
!        AMASS_IN  - Atomic mass for Doppler profile in amu's.
!        Z_IN      - Ion charge.
!        NL        - Lower level of transition.
!        NUP       - Upper level of transition.
!        ND        - Nuber of (Ne,T) values profile is to be computed for.
!        GAM_RAD   - Natural (radiative) broadening parameter for Voigt profile
!        CAM_COL   - Collisional broadening parameter for Voigt profile
!                          (quadratic stark effect).
!
!	LIMIT_SET_BY_OPACITY - Determines the STRT_FREQ of the intrinsic
!                                 line absorption profile based on the
!                                 ratio of line to electron scattering opacity.
!       DOP_LIMIT  - Truncate the Doppler profile when the ratio of the LINE 
!                      opacity to the Electron scattering opacity is DOP_LIMIT.
!                      Should be of order 10^{-3} or less.
!                      The profile CANNOT be truncated insiede 3.5 Doppler 
!                      widths from line center
!       VOIGT_LIMIT  - Truncate the VOIGT profile when the ratio of the LINE 
!                      opacity to the Electron scattering opacity is DOP_LIMIT.
!                      Should be of order 10^{-3} or less.
!                      The profile CANNOT be truncated inside 3.5 Doppler 
!                      widths from line center
!       V_PROF_LIMIT - In km/s. Used to determine limits for profiles (excluding
!                         Doppler and Voigt).
!        MAX_PROF_ED - Maximum electron density used to compute Stark profile.
!                         Not used at present, but may be need later.
!
	SUBROUTINE SET_PROF_LIMITS_V4(VEC_STRT_FREQ,VEC_VDOP_MIN,
	1             CHIL,ED_IN,TEMP_IN,VTURB_IN,ND,
	1             PROF_TYPE,PROF_LIST_LOCATION,
	1             NU_ZERO,NL,NUP,SPECIES_IN,AMASS_IN,Z_IN,
	1             GAM_RAD,GAM_COL,TDOP,AMASS_DOP,VTURB,
	1             DOP_LIMIT,VOIGT_LIMIT,V_PROF_LIMIT,MAX_PROF_ED,
	1             LIMIT_SET_BY_OPACITY)
	USE MOD_STRK_LIST
!
! Altered 18-May-2015 : Bug fix: When WAVE was outside valid range, GET_INDX_DP was returning an error 
!                          to fort.2 (not OUTGEN) but was also returning a valid index. As a result, some
!                          lines (in this case belonging to HeI) would use the wrong profile.
!                          LOC_GAM_COL is no longer set to zer0. Uses passed value (GAM_COL).
!                          Removed limit on whenVOIGT profiles used -- may need to be updated.  
! Altered 14-May-2015 : Limit Ne to MAX_PROF_ED.
! Altered 31-Jan-2014 : Added V_PROF_LIMIT & MAX_PROF_ED to call. Changed to V4.
! Altered  6-Jan-2014 : V3 was added to repositry 6_jan-2014.
!                           (taken from cur_cm_25jun13 development version).
! Altered 01-May-2003 : Minor bug fix LO_GAM_COL was being incorrectly set.
! Altered 27-Nov-2001 : Bug fix. Incorrect frequency limits for DOP_FIX option.
! Altered 06-Apr-2000 : LINE_TO_CONT_RATIO was zero at one depth, causing a divide
!                         by zero.
!
	IMPLICIT NONE
	INTEGER ND
	INTEGER NL
	INTEGER NUP
!
	REAL*8 VEC_STRT_FREQ
	REAL*8 VEC_VDOP_MIN
	REAL*8 DOP_LIMIT
	REAL*8 VOIGT_LIMIT
!
	REAL*8 CHIL(ND) 
	REAL*8 ED_IN(ND)
	REAL*8 TEMP_IN(ND)
	REAL*8 VTURB_IN(ND)
!
	CHARACTER(LEN=*) SPECIES_IN
	REAL*8 AMASS_IN
	REAL*8 Z_IN
	REAL*8 NU_ZERO
	REAL*8 GAM_RAD
	REAL*8 GAM_COL
	REAL*8 TDOP
	REAL*8 AMASS_DOP
	REAL*8 VTURB
	REAL*8 V_PROF_LIMIT
	REAL*8 MAX_PROF_ED
	INTEGER PROF_LIST_LOCATION
	CHARACTER*(*) PROF_TYPE
!
	LOGICAL LIMIT_SET_BY_OPACITY
!
	REAL*8 SPEED_OF_LIGHT
	INTEGER GET_INDX_DP
	EXTERNAL SPEED_OF_LIGHT,GET_INDX_DP
!
! Local variables
!
	INTEGER I,J
	REAL*8 ESEC(ND)
	REAL*8 ED_MOD(ND)
	REAL*8 TMP_VEC(ND)
	REAL*8 NU_DOP(ND)
	REAL*8 LINE_TO_CONT_RATIO(ND)
	REAL*8 dNU
	REAL*8 C_KMS
	REAL*8 WAVE
	REAL*8 LOC_GAM_RAD
	REAL*8 LOC_GAM_COL
	REAL*8 PROF_LINE_CENTER
	REAL*8 T1,T2
	REAL*8 LOCAL_V_PROF_LIMIT
!
	INTEGER, PARAMETER :: NUM_DOP=6
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()		!km/s
	LOC_GAM_RAD=GAM_RAD
	LOC_GAM_COL=GAM_COL
	PROF_LIST_LOCATION=0
	LOCAL_V_PROF_LIMIT=V_PROF_LIMIT
	WAVE=0.01D0*C_KMS/NU_ZERO
!
! The option assumes fixed width Doppler profiles, and recovers exactly the
! same option as was installed in CMFGEN prior to the installations of variable 
! Doppler widths.
!
	IF(PROF_TYPE .EQ. 'DOP_FIX')THEN
	  VEC_VDOP_MIN=12.85D0*SQRT( TDOP/AMASS_DOP + (VTURB/12.85D0)**2 )
          NU_DOP(1:ND)=NU_ZERO*VEC_VDOP_MIN/C_KMS				!kms
          VEC_STRT_FREQ=NU_ZERO+NUM_DOP*NU_DOP(1)
	  RETURN
	END IF
!
	IF(PROF_TYPE .EQ. 'DOP_SPEC')THEN
	  VEC_VDOP_MIN=12.85D0*SQRT( TDOP/AMASS_IN + (VTURB/12.85D0)**2 )
          NU_DOP(1:ND)=NU_ZERO*VEC_VDOP_MIN/C_KMS				!kms
          VEC_STRT_FREQ=NU_ZERO+NUM_DOP*NU_DOP(1)
	  RETURN
	END IF
!
	IF(PROF_TYPE(1:4) .EQ. 'LIST')THEN
	  WAVE=0.01D0*C_KMS/NU_ZERO
	  IF(WAVE .LT. LST_WAVE(1))THEN
	    I=1
	  ELSE IF(WAVE .GT. LST_WAVE(N_LST))THEN
	    I=N_LST
	  ELSE
	    I=GET_INDX_DP(WAVE,LST_WAVE,N_LST)
	  END IF
	  DO J=MAX(1,I-2),MIN(I+2,N_LST)
	    IF( TRIM(SPECIES_IN) .EQ. TRIM(LST_SPECIES(J)) .AND.
	1          ABS((WAVE-LST_WAVE(J))/WAVE) .LT. 1.0D-05)THEN
	      PROF_TYPE=LST_TYPE(J)
	      PROF_LIST_LOCATION=J
	      IF(LST_V_PROF_LIMIT(J) .NE. 0)LOCAL_V_PROF_LIMIT=LST_V_PROF_LIMIT(J)
	      IF(PROF_TYPE .EQ. 'VOIGT')THEN
	        LOC_GAM_RAD=LST_GAM_RAD(J)
	        LOC_GAM_COL=LST_GAM_COL(J)
	      END IF
	      EXIT		!As found profile.
	    END IF
	  END DO
	END IF
!
! For high HI and HeII line profiles we automatically chooses the HZ_STARK option
! if the pofile type has not already been set.
!
	IF(PROF_TYPE(1:4) .EQ. 'LIST')THEN
	  IF(SPECIES_IN .EQ. 'HI' .OR. SPECIES_IN .EQ. 'He2')PROF_TYPE='HZ_STARK'
	END IF
!
! We assume all lines are dominated by the Doppler profile at line center.
! We ensure LINE_TO_CONT ratio is not zero, to prevent division by zero.
!                    
        TMP_VEC(1:ND)=( VTURB_IN(1:ND)/12.85D0  )**2
	ESEC(1:ND)=6.65D-15*ED_IN(1:ND)
	ED_MOD=ED_IN
	DO I=1,ND
	  ED_MOD(I)=MIN(ED_IN(I),MAX_PROF_ED)
	END DO
!                       
        T1=1.0D-15/1.77245385095516D0         !1.0D-15/SQRT(PI)
	VEC_VDOP_MIN=1.0D+50			!Very large number 
        DO I=1,ND
          NU_DOP(I)=12.85D0*SQRT( TEMP_IN(I)/AMASS_IN + TMP_VEC(I) )	!kms
	  VEC_VDOP_MIN=MIN(VEC_VDOP_MIN,NU_DOP(I))			!kms
          NU_DOP(I)=NU_DOP(I)*NU_ZERO/C_KMS				!10^15 Hz
          PROF_LINE_CENTER=T1/NU_DOP(I)
	   LINE_TO_CONT_RATIO(I)=ABS(CHIL(I))*PROF_LINE_CENTER/ESEC(I)
	   IF(LINE_TO_CONT_RATIO(I) .EQ. 0)LINE_TO_CONT_RATIO(I)=1.0D-50
!	   IF(LINE_TO_CONT_RATIO(I) .GT. 1.0E+04 .AND. PROF_TYPE .EQ. 'LIST_VGT')THEN
	   IF(PROF_TYPE .EQ. 'LIST_VGT')THEN
	     LOC_GAM_RAD=GAM_RAD
	     LOC_GAM_COL=1.55D+04*(ABS(GAM_COL)**0.667D0)		!GAM_COL
	     PROF_TYPE='VOIGT'
	   END IF 
	END DO
!
! If the line profile has not been set, we choose DOPPLER as the default.
!
	IF(PROF_TYPE(1:4) .EQ. 'LIST')THEN
	  PROF_TYPE='DOPPLER'
	END IF
!
	VEC_STRT_FREQ=0.0D0
	IF(PROF_TYPE .EQ. 'DOPPLER' .OR. PROF_TYPE .EQ. 'VOIGT')THEN
	  IF(LIMIT_SET_BY_OPACITY)THEN
	    DO I=1,ND
              T1=-LOG(DOP_LIMIT/LINE_TO_CONT_RATIO(I))
	      IF(T1 .GT. 0)T1=SQRT(T1)
	      T1=MAX(3.5D0,T1)
              VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+T1*NU_DOP(I))
	    END DO
	  ELSE
	    DO I=1,ND
              VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+NUM_DOP*NU_DOP(I))
	    END DO
	  END IF
	  IF(PROF_TYPE .EQ. 'DOPPLER')RETURN
	END IF
!
! NB: 6.7005D-09=1.0D-15*SQRT(1.0D+15 /4  / PI**1.5)
! The factor of 1.0D-15 arises since dNU has to be in units of 10^15 Hz.
! The factor of 1.0D+15 arisis since NU_DOP is in units of 10^15 Hz.
!
	dNU=0.0D0
	IF(PROF_TYPE .EQ. 'VOIGT')THEN
	  IF(LIMIT_SET_BY_OPACITY)THEN
	    DO I=1,ND            
	      T2=LOC_GAM_RAD+LOC_GAM_COL*ED_MOD(I)
              T1=6.7005D-09*SQRT(T2*NU_DOP(I)*LINE_TO_CONT_RATIO(I)/VOIGT_LIMIT)
	      dNU=MAX(T1,dNU)
	    END DO
	    dNU=MIN(dNU,NU_ZERO*5.0D+03/C_KMS)
	    VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+dNU)
	  ELSE
!
! Truncates profile when profile is down to 1.0E-06 of the value at line
! center.
!
	    DO I=1,ND
	      T2=LOC_GAM_RAD+LOC_GAM_COL*ED_MOD(I)
              T1=6.7005D-09*SQRT(1.0D+06*T2*NU_DOP(I))
	      dNU=MAX(T1,dNU)
	    END DO
	    dNU=MIN(dNU,NU_ZERO*5.0D+03/C_KMS)
	    VEC_STRT_FREQ=MAX(VEC_STRT_FREQ,NU_ZERO+dNU)
	  END IF
	  RETURN
	END IF
!
! For all other profile types: CRUDE.
!
	VEC_STRT_FREQ=NU_ZERO*(1.0D0+LOCAL_V_PROF_LIMIT/C_KMS)
!
!	WRITE(6,'(ES12.4,2I6,3X,A10,2F6.2)')0.2998D+04/NU_ZERO,NL,NUP,SPECIES_IN,AMASS_IN,Z_IN
!	T1=(VEC_STRT_FREQ/NU_ZERO-1.0D0)*C_KMS
!	WRITE(177,'(X,A10,2E14.5)')PROF_TYPE,NU_ZERO,T1
!
	RETURN
	END
