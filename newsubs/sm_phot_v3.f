!============================================================================
!
	SUBROUTINE SM_PHOT_V3(NU,CROSS,NCROSS,NSM_MAX,
	1                        SIG_GAU_KMS,FRAC_SIG_GAU,
	1                        CUT_ACCURACY,ABOVE)
!
! Altered: 8-Oct-2008 : Bug fixed 
!                       Ensure cross section falls off at least as 1/nu
!                         in extrapolation region.
!
!============================================================================
!
! Routine to smooth an arbitrary photionization cross-section. Cross-section
! is smoothed with a Gaussian whose width is equal to SIG_GAU times the
! current smoothed frequency. Thus cross-section area is not conserved
! exactly. Allows wider freqeuncy spacing as the frequency is increased.
!
! On entry
!         NU            Should contain the frequency in units of nu_zero
!         CROSS         The cross section.
!         NCROSS        Number of cross sections.
!         SIG_GAU       The sigma of the gaussian smoothing routine.
!         DEL_NU        Frequency spacing for smoothed cross-section.
!                         (should be fraction of SIG_GAU).
!         CUT_ACCURACY  Accuracy to limit omission of data points.
!                          Assumes inear interpolation.
!         ABOVE         Smooth from above threshold only?
!
! On exit
!        NU_SM       Contains the evenly spaced frequencies at which
!                    the cross-section is tabulated.
!        CROSS_SM    Smoothed Cross-Section
!        NSM         Number of frequency points in smooth cross-section.
!
! Altered 24-May-2005 : Fixed to handle case where mirroring of cross-section
!                        about NU=1 gives negative frequencies.
! Created 19-Apr-2004 : Based on SM_PHOT_V2
!                       Designed to be incorporated into CMFGEN
!                       Some cleaning.
!			See SM_PHOT_V2 for earlier changes.
!
	IMPLICIT NONE
	INTEGER NCROSS
!
	INTEGER NSM_MAX
	REAL*8 NU(NSM_MAX)
	REAL*8 CROSS(NSM_MAX)
!
	REAL*8 SIG_GAU_KMS
	REAL*8 FRAC_SIG_GAU
	REAL*8 CUT_ACCURACY
	LOGICAL ABOVE
!
! Work arrays.
!
	REAL*8 NU_SM(NSM_MAX)
	REAL*8 CROSS_SM(NSM_MAX)
!
	LOGICAL FLAG
	EXTERNAL EQUAL
	LOGICAL EQUAL
!
! Local variables
!
	INTEGER, PARAMETER :: NZ=10000
	INTEGER, PARAMETER :: NLOC=200000
	REAL*8 Z(NZ)
	REAL*8 NU_FINE(NLOC)
	REAL*8 CROSS_FINE(NLOC)
	INTEGER NFINE
!
	REAL*8 T1,T2,DIFF,SUM,SIG_NEW
	REAL*8 LAST_FREQ
	REAL*8 DNU
	REAL*8 SIG_GAU
	INTEGER I,J,K,ML
	INTEGER IST,IEND,IBEG,N_INS
	INTEGER NSM
!
	INTEGER NU_INDX_ST
	REAL*8 NU_ST
!
	REAL*8 MAX_DEL
	REAL*8 PTS_PER_SIG
!
	PTS_PER_SIG=4.0D0			!Used for the convolution.
	SIG_GAU=SIG_GAU_KMS/2.998D+05		!C does not need to be accurate
	DNU=FRAC_SIG_GAU*SIG_GAU
!
! If ABOVE is TRUE, we ignore cross-section (except for point immediately
! proceeding threshold) when smoothing cross-section.  We reflect the cross 
! section about the threshold frequency. This ensures that we conserve area.
! If ABOVE is FALSE, we smooth acrossa the threshold. 
!
	IF(ABOVE)THEN
!
	  I=1
	  DO WHILE(NU(I) .LT. 1.0D0)
	    I=I+1
	  END DO
	  IBEG=I
!
! Find range of Gaussian.
!
	  DO WHILE(NU(I) .LT. 1.0D0+6.0D0*SIG_GAU)
	    I=I+1
	  END DO
	  N_INS=I-IBEG+1
!
! How we proceed depends on whether NU=1 is present.
! If there is a large step size in NU near NU=1, we may get -ve frequencies
! This will not effect the answer, since we are just integrating the
! cross-section. For reasonable a reasonable SIG_GAU (< 10,000 km/s) we
! will only integrate in the positive frequency range anyway.
!
	  IF(NU(IBEG) .NE. 1.0D0)THEN
	    IF(IBEG .EQ. 1)THEN
	      T2=0.0D0
	    ELSE
	      T1=(NU(IBEG-1)-1.0D0)/(NU(IBEG-1)-NU(IBEG))
	      T2=(1.0D0-T1)*CROSS(IBEG-1)+T1*CROSS(IBEG)
	    END IF
	    N_INS=N_INS+1
	    IF(N_INS .LT. IBEG-1)N_INS=IBEG-1
	    DO J=NCROSS,IBEG,-1
	      K=J+N_INS-IBEG+1
	      NU(K)=NU(J)
	      CROSS(K)=CROSS(J)
	    END DO
	    NU(N_INS)=1.0D0
	    CROSS(N_INS)=T2
	    DO J=N_INS-1,1,-1
	      K=2*N_INS-J
	      NU(J)=2.0D0-NU(K)
	      CROSS(J)=CROSS(K)
	    END DO
	    NCROSS=NCROSS-IBEG+1 + N_INS
	  ELSE
	    IF(N_INS .LT. IBEG-1)N_INS=IBEG-1
	    DO J=NCROSS,IBEG,-1
	      K=J+N_INS-IBEG+1
	      NU(K)=NU(J)
	      CROSS(K)=CROSS(J)
	    END DO
	    DO J=N_INS,1,-1
	      K=2*N_INS-J+2
	      NU(J)=2.0D0-NU(K)
	      CROSS(J)=CROSS(K)
	    END DO
	    NCROSS=NCROSS-IBEG+1+ N_INS
	  END IF
	ELSE
!
! Usually we will have enough points below threshold. For simplicty,
! we assume a constant cross-section below these points.
!
	 IF(NU(1) .GT. 1.0D0-6.0D0*SIG_GAU)THEN
	   DO J=NCROSS,1,-1
	     K=J+1
	     NU(K)=NU(J)
	     CROSS(K)=CROSS(J)
	   END DO
	   NU(1)=1.0D0-6.0D0*SIG_GAU
	   CROSS(1)=CROSS(2)
	   NCROSS=NCROSS+1
	  END IF
	END IF
!
! This test should be done earlier.
!
	IF(NCROSS+1 .GT. NSM_MAX)THEN
	  WRITE(6,*)'Error in SM_PHOT_V3'
	  WRITE(6,*)'NSM_MAX is too small'
	  STOP
	END IF
!
! Extend cross-sections at high frequencies. We ensure cross-section
! falls off as least as 1/nu.
!
	NCROSS=NCROSS+1
	NU(NCROSS)=NU(NCROSS-1)*(1.0D0+6.0D0*SIG_GAU)
	T1=LOG(NU(NCROSS-1)/NU(NCROSS-2))
	T2=LOG(CROSS(NCROSS-1)/CROSS(NCROSS-2))/T1
	IF(T2 .GT. -1.0D0)T2=-1.0D0
	CROSS(NCROSS)=CROSS(NCROSS-1)*(NU(NCROSS)/NU(NCROSS-1))**T2
!
!	CALL DP_CURVE(NCROSS,NU,CROSS)
!
! Interpolate cross-section onto a sufficiently fine grid so that the Gaussian
! correctly samples the cross-section data. At present we assume linear 
! interpolation. This is satisfactory given the accuracy of the cross-sections.
!
! NB: The Gaussian is assumed to have constant width in velocity space.
!     For this reason SIG_GAU is mutipled by NU(ML). For frequencies < 1,
!     we use SIG_GAU. This prevents DIFF from becoming arbitrarily small if
!     NU should approach zero.
!
	I=1
	NU_FINE(1)=NU(1)
	CROSS_FINE(1)=CROSS(1)
	DO ML=2,NCROSS
	  MAX_DEL=SIG_GAU*NU(ML)/PTS_PER_SIG
	  IF(NU(ML) .LE. 1.0D0)MAX_DEL=SIG_GAU/PTS_PER_SIG
	  DIFF=NU(ML)-NU(ML-1)
	  IF( DIFF .GT. MAX_DEL)THEN
            J=DIFF/MAX_DEL+1
	    T2=DIFF/J
	    DO K=1,J-1
	      IF(I+K .GT. NLOC)GOTO 999
	      NU_FINE(I+K)=NU_FINE(I)+K*T2
	      T1=(NU_FINE(I+K)-NU(ML-1))/(NU(ML)-NU(ML-1))
	      CROSS_FINE(I+K)=(1.0D0-T1)*CROSS(ML-1)+T1*CROSS(ML)
	    END DO
	    I=I+J-1
	  END IF
	  I=I+1
	  IF(I .GT. NLOC)GOTO 999
	  NU_FINE(I)=NU(ML)
	  CROSS_FINE(I)=CROSS(ML)
	END DO
	NFINE=I
!
!	CALL DP_CURVE(NFINE,NU_FINE,CROSS_FINE)
!
! Define the smooth mesh. DNU is the frequency spacing for the 
! smoothed cross-section output.
!
	FLAG=.TRUE.
	NSM=1
	NU_SM(1)=1.0D0
	LAST_FREQ=NU(NCROSS)
	DO WHILE(FLAG)
	  IF(NSM+1 .GT. NSM_MAX)THEN
	    WRITE(6,10)NSM,NSM_MAX
 10         FORMAT(' Error NSM_MAX TOO small : NSM =',I6,': NSM_MAX =',I6)
	    WRITE(6,*)'DNU=',DNU
	    STOP
	  END IF
	  NU_SM(NSM+1)=NU_SM(1)*(1.0D0+DNU)**NSM
	  IF(NU_SM(NSM+1) .GT. LAST_FREQ)EXIT
	  NSM=NSM+1
	END DO
C 
C Now determine smooth cross section for each (smooth) frequency.
C
C We use Z to generate for gaussian over the desired freqency band.
C We also integrate Z to get the proper normailization.
C
C IST  is index for ML at which left side of Gaussian begins.
C IEND is index for ML at which right side of Gaussian ends.
C
	IST=1
	IEND=2
	DO ML=1,NSM
	  SIG_NEW=SIG_GAU*NU_SM(ML)
	  DO WHILE (NU_FINE(IST) .LT. NU_SM(ML)-5.0D0*SIG_NEW)
	    IST=IST+1
	  END DO
	  DO WHILE ( (NU_FINE(IEND)-NU_SM(ML)) .LT. 5.0D0*SIG_NEW)
	    IF(IEND .EQ. NFINE)GOTO 100
	    IEND=IEND+1
	  END DO
100	CONTINUE
!
! Determine the smoothing Gaussian.
!
	  IF(IEND-IST+1 .GT. NZ)THEN
	    WRITE(6,*)'Error in SM_PHOT_V3: NZ too small'
	    STOP
	  END IF
	  DO I=IST,IEND
	    J=I-IST+1
            Z(J)=EXP( -0.5D0*( (NU_SM(ML)-NU_FINE(I))/SIG_NEW )**2 )
	  END DO
!
! Perform the integration for the current frequency.
!
	  CROSS_SM(ML)=0.0D0
	  SUM=0.0D0
	  DO I=IST,IEND-1
	    J=I-IST+1
	    CROSS_SM(ML)=CROSS_SM(ML)+ (NU_FINE(I+1)-NU_FINE(I))*
	1         (CROSS_FINE(I)*Z(J)+CROSS_FINE(I+1)*Z(J+1))
	    SUM=SUM + (NU_FINE(I+1)-NU_FINE(I))*(Z(J)+Z(J+1))
	  END DO
	  IF(SUM .NE. 0.0D0)THEN
	    CROSS_SM(ML)=CROSS_SM(ML)/SUM
	  ELSE
	    WRITE(6,*)'Zero SUM in SM_PHOT_V3'
	    WRITE(6,*)IST,IEND,ML,NSM
	    WRITE(6,*)'NU_SM(ML)=',NU_SM(ML)
	    WRITE(6,*)'NU_FINE(IST)=',NU_FINE(IST)
	    WRITE(6,*)'NU_FINE(IEND)=',NU_FINE(IEND)
	    WRITE(6,*)'SIG_GAU=',SIG_GAU
	    DO J=1,NCROSS
	      WRITE(30,'(I6,2ES16.6)')J,NU(J),CROSS(J) 
	    END DO
	    STOP
	  END IF
	END DO				!Frequency loop
!
! Omit pints that are uneessary to maintain an accuracy of CUT_ACCURACY in the
! cross-section assuming linear interpolation.
!
        CALL CUT_POINTS_V3(NU,CROSS,NCROSS,NU_SM,CROSS_SM,NSM,CUT_ACCURACY)
!
	RETURN
!
999	WRITE(6,*)'Error --- NLOC too small in SM_PHOT_V3'
	WRITE(6,*)'NU(1)=',NU(1)
	WRITE(6,*)'NU(NCROSS)=',NU(NCROSS)
	WRITE(6,*)'DNU=',DNU
!
	STOP
	END
