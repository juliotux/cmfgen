      SUBROUTINE STRK_BS_HHE(PROF,FREQ,NF,
	1             ED,T,VTURB,ND,
	1             NU_ZERO,AMASS,TRAP_QUAD)
C
C
C routine to obtain stark profile for a tabulated transition
C         as a function of electron density, temperature and lambda.
C input:
C     ED(ND)      : electron density in cgs 
C     T(ND)       : temperature in units of 1e4 K 
C     DXLAM (NLAM): displacement in Angstroms from line center
C                  
C output: loge of stark profile is returned in PROF(NLAM,ND)
C
C  Altered 5-Sep-03 : Calculation of DLAM_THERM is now correct.     
C  Altered Feb-99  LEMKE'S HYD PROFILES ALSO INCLUDED
C  created Nov-95  Based in the interpolation loop of routine
C                  utable from Keith. Changes have been introduced
C                  so that routine does not return PROF=0 when electron
C                  density is below ranges in table (MID variable)
C                  and extrapolation is performed. Before electron density
C                  and temperature where forced to be within ranges of table.
C                  
C            
C
	USE STRK_MOD_HHE
	IMPLICIT NONE

	INTEGER ND
	INTEGER NF
!
	REAL*8 PROF(ND,NF)
	REAL*8 FREQ(NF)
	REAL*8 ED(ND)
	REAL*8 T(ND)
	REAL*8 VTURB(ND)
!
	REAL*8 NU_ZERO
	LOGICAL TRAP_QUAD
!
	REAL*8 WAVE
	REAL*8 AMASS
!
	REAL*8 TLOG
	REAL*8 ELOG
	REAL*8 ST1,ST2
	REAL*8 C_KMS
!
	REAL*8, ALLOCATABLE :: STARK(:)
	REAL*8, ALLOCATABLE :: TMP_DWS(:)
!
	REAL*8 DLAM(NF)
	REAL*8 LOC_PROF(NF)
!
	REAL*8 T_WGT
	REAL*8 ED_WGT
!
	REAL*8 DLAM_THERM
	REAL*8 DLAM_TURB
!               
	INTEGER I,ID
	INTEGER T_INDX
	INTEGER ED_INDX
	INTEGER L0,L1,L2,L3,ML
	LOGICAL DIAGNOSTICS
!
	DIAGNOSTICS=.FALSE.
	C_KMS=2.99794D+05
	WAVE=0.01D0*C_KMS/NU_ZERO
	DLAM(1:NF)=0.01D0*C_KMS/FREQ(1:NF)-WAVE
!
	IF(ALLOCATED(STARK))THEN
	   DEALLOCATE(STARK)
	   DEALLOCATE(TMP_DWS)
	END IF
	ALLOCATE(STARK(STKTA_NWS))
	ALLOCATE(TMP_DWS(STKTA_NWS))
!
! Obtain profile at each depth.
!
	DO ID = 1,ND
	  TLOG = LOG10(T(ID))+4.0    ! T was in 1.e4 K
	  TLOG = MAX(STKTA_TS(1),MIN(STKTA_TS(STKTA_NTS),TLOG))
!
! Find T interval
!
          DO I = 2,STKTA_NTS
            T_INDX=I-1
            IF (TLOG .LT. STKTA_TS(I))EXIT
	  END DO
          T_WGT =(TLOG-STKTA_TS(T_INDX))/(STKTA_TS(T_INDX+1)-STKTA_TS(T_INDX))
!
! Find electron density interval
!
          ELOG = LOG10(ED(ID))                                    
          ELOG = MAX(STKTA_ES(1),MIN(STKTA_ES(STKTA_NES),ELOG))
          DO I = 2,STKTA_NES                                     
            ED_INDX = I-1                                           
            IF (ELOG .LT. STKTA_ES(I))EXIT
	  END DO
          ED_WGT = (ELOG-STKTA_ES(ED_INDX))/
	1                 (STKTA_ES(ED_INDX+1)-STKTA_ES(ED_INDX))
!                           
! Obtain the interpolated profile. NB: The porifles aren't tabulated
! on a fine enough grid. Thus interpolated profiles can be inaccurate.
! As a result of the inaccuracy, profiles need to be normalized, by
! up to a factor of 1.2
!
! Test value: T=3.1, Log Ne=13.3 (obviously not table value)
!
	  L0 = STKTA_NWS* (T_INDX-1+STKTA_NTS*(ED_INDX-1))
	  L1 = STKTA_NWS* (T_INDX+STKTA_NTS*(ED_INDX-1))
	  L2 = STKTA_NWS* (T_INDX-1+STKTA_NTS*ED_INDX)
	  L3 = STKTA_NWS* (T_INDX+STKTA_NTS*ED_INDX)
	  DO ML=1,STKTA_NWS
	    ST1=T_WGT*(STKTA_PS(ML+L1)-STKTA_PS(ML+L0))+STKTA_PS(ML+L0)
	    ST2=T_WGT*(STKTA_PS(ML+L3)-STKTA_PS(ML+L2))+STKTA_PS(ML+L2)
	    STARK(ML)=ED_WGT*(ST2-ST1)+ST1
	  END DO
!
	  DLAM_THERM=WAVE*12.86D0*SQRT(T(ID)/AMASS)/C_KMS
	  DLAM_TURB=WAVE*VTURB(ID)/C_KMS
!
! The LEMKE profiles are not tabulated to 0 offset. We set the
! first frequency to zero, to allow the convolution.
!                            
	  IF(STKTA_WSCA)THEN
	   TMP_DWS=STKTA_DWS*1.25D-9*(10**(ELOG*2.0D0/3.0D0))
	   IF(STKTA_QHALF .AND. TMP_DWS(1) .NE. 0)TMP_DWS(1)=0.0D0
	  ELSE
	   TMP_DWS=STKTA_DWS
	  END IF
          CALL CONV_STRK_V2(LOC_PROF,DLAM,NF,
	1          STARK,TMP_DWS,STKTA_NWS,
	1          STKTA_QHALF,TRAP_QUAD,.TRUE.,
	1          DLAM_THERM,DLAM_TURB,WAVE,
	1          ELOG,TLOG,DIAGNOSTICS)
!
	  PROF(ID,1:NF)=LOC_PROF(1:NF)
!
	END DO
!
	RETURN
	END
