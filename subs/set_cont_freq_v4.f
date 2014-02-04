!
! Routine to generate frequencies for radiative transfer program given a set of 
! bound-free edges, and maximum and minimum frequencies. The frequency spacing is 
! controlled by passable parameters.
!
! Every bound-free edged is defined by two frequencies, except where two bound-free 
! edges coincide to 1 part in 10^11. Generally points are inserted so that Simpson's 
! rule can be used, although in some regions where the bound-free edges are close the
! trapezoidal rule will be used.
!
	SUBROUTINE SET_CONT_FREQ_V4(NEW_FREQ,FREQ,INDX,
	1                        SMALL_RAT,BIG_AMP,DNU_MAX,MAX_FREQ,MIN_FREQ,
	1                        dV_LEV,AMP_DIS,MIN_FREQ_LEV_DIS,
	1                        DELV_CONT,DELV_XRAY,NU_END_XRAY,N,NCF,NCF_MAX,LUOUT)
	IMPLICIT NONE
!
! Altered 18-Aug-2004 : Changed way we insert points after a bound-free edge.
!                       Initial spacing is DNU_MAX. The spacing increases by
!                         0.5*DNU_MAX on each additional quadrature point. This should
!                         work since spacing is fine enough when hv/kT is small, and 
!                         integral is dominated by first quadrature point when hv/kT
!                         is very large. There is no specific maximum set, but in 
!                         practice this will be set by DELV_CONT. The meaning of
!                         BIG_AMP has changed --- typical value should be 0.5.
!                         Spacing increases by BIG_AMP*DNU_MAX on successive nodes.
!
!                       We no longer insert extra points to use Simpson's rule. We
!                         no longer use Simpson's rule since ether are so many edges,
!                         and we use the same quadrature weights FOR ALL integrals.
!                         As a consequence spacing is now DNU_MAX near a bound-free 
!                         edge, not DNU_MAX/2.
!                         [Changed to V4 so that keep old routine] 
!
! Altered 26-May-1996 : ERROR_LU installed.
! Cleaned
! Created 29-Mar-1990 (Based on GEN_FREQ).
!
	INTEGER N,NCF,NCF_MAX,LUOUT,INDX(NCF_MAX)
	REAL*8 FREQ(NCF_MAX),NEW_FREQ(NCF_MAX)
!
	REAL*8 MAX_FREQ		!Maximum continuum frequency
	REAL*8 MIN_FREQ		!Minimum continuum frequency
	REAL*8 DNU_MAX		! Twice the maximum frequency spacing near/above
				!  bound-free edge: i.e. dNU < 0.5* DNU_MAX
	REAL*8 BIG_AMP  	!Amplification. dNU increases by a factor
				! BIG_AMP as we move away from the b.f. edge.
				! for frequencies above SWITCH_FREQ.
   	REAL*8 SMALL_RAT	!Used to define frequency spacing for
				! frequencies less than SWITCH_FREQ.
				! dNU/NU=SMALL_RAT-1
!
! Parameters for installing extra frequencies near bound-free edges
! (low frequency side) to allow for level dissolution.
!
	REAL*8 dV_LEV			!Spacing near b-f edge.
	REAL*8 AMP_DIS			!Amplification factor for dNU as we
					!  move to smaller frequencies.
	REAL*8 MIN_FREQ_LEV_DIS		!Indicates that the extra frequencies
					! should only be installed for
					! frequencies above MIN_FREQ_LEV_DIS.
!
	REAL*8 DELV_CONT
	REAL*8 DELV_XRAY
	REAL*8 NU_END_XRAY
!
        REAL*8, PARAMETER :: RZERO=0.0D0
        REAL*8, PARAMETER :: RHALF=0.5D0
        REAL*8, PARAMETER :: RONE=1.0D0
        REAL*8, PARAMETER :: RTWO=2.0D0
!
	INTEGER ERROR_LU,LUER
	EXTERNAL ERROR_LU
!
	REAL*8 T1,T2,dV_NEW,dV
	REAL*8 C_KMS
!
	INTEGER I,J,K,L,ML,NPTS,K_BEG
	LOGICAL NUMER,EQUAL
	REAL*8 FAC,EQUAL_FAC
	REAL*8 UP,LOW,DELF,DIFF,RAT,SWITCH_FREQ
!
	LUER=ERROR_LU()
	C_KMS=2.998D+05			!Doesn't have to be accurate
!
! Sort frequencies into numerical order. New freq is used as a
! work array.
!
	N=N+2
	FREQ(N-1)=MIN_FREQ
	FREQ(N)=MAX_FREQ
	NUMER=.TRUE.
	CALL INDEXX(N,FREQ,INDX,NUMER)
	CALL SORTDP(N,FREQ,INDX,NEW_FREQ)
!
! Not equal as MIN_FREQ, and MAX_FREQ have been inserted in FREQ array.
!
	IF(MIN_FREQ .NE. FREQ(1))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V4 - MIN_FREQ too big'
	  WRITE(LUER,*)'MIN_FREQ =',MIN_FREQ
	  WRITE(LUER,*)'FREQ(1) =',FREQ(1)
	  STOP
	END IF
	IF(MAX_FREQ .NE. FREQ(N))THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V4 - MAX_FREQ too small'
	  WRITE(LUER,*)'MAX_FREQ =',MAX_FREQ
	  WRITE(LUER,*)'FREQ(N) =',FREQ(N)
	  STOP
	END IF
!
	EQUAL_FAC=1.0D-11
	FAC = RONE + EQUAL_FAC
!
! SMAL_FAC is the ratio used to set the frequency spacing for
! small frequencies (frequencies less than 1.0 approximately).
! Two points are inserted so that Simpson's rule can be used.
!
	IF(SMALL_RAT .LE. RONE .OR. SMALL_RAT .GT. 3)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V4 - SMALL_RAT Outside ',
	1           'valid range ',SMALL_RAT
	  STOP
	END IF
!
! DNU_MAX is the maximum frequency spacing for frequencies adjacent
! to a bound-free edge. For large v, we require constant spacing since
! error in integral is a function of {del v }.
!
	IF(DNU_MAX .LE. RZERO .OR. DNU_MAX .GT. RTWO)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V4 - DNU_MAX Outside ',
	1           'valid range ',DNU_MAX
	  STOP
	END IF
!
! G is used to amplify the spacing for large frequencies. Close
! to the edge, the spacing is less than DNU_MAX, but for every
! two points, the spacing increases by a factor of G.
!
	IF(BIG_AMP .LE. RZERO .OR. BIG_AMP .GT. 2)THEN
	  WRITE(LUER,*)'Error in SET_CONT_FREQ_V4 - BIG_AMP ',
	1                 'Outside valid range ',BIG_AMP
	  STOP
	END IF
!
	SWITCH_FREQ=DNU_MAX/(SMALL_RAT-RONE)
!
! Now begin inserting extra points into frequency array.
!
	K=1
	NEW_FREQ(1)=MIN_FREQ
	DO I=1,N-1
	  IF( EQUAL(FREQ(I),FREQ(I+1),EQUAL_FAC) )THEN
	    FREQ(I+1)=FREQ(I)
	  ELSE
	    LOW=FREQ(I)
	    UP=FREQ(I+1)
	    IF(I .NE. 1)THEN
	      LOW=LOW*FAC
	      K=K+1
	      NEW_FREQ(K)=LOW
	    END IF
	    IF(I .NE. N-1)UP=UP/FAC
	    RAT=UP/LOW
	    DIFF=UP-LOW
!
! The insertion of frequencies is treated differently depending on whether the
! EDGE_FREQ (stored in FREQ) is greater than, or less than, SWITCH_FREQ.
!
! If EDGE_FREQ < SWITCH_FREQ  we insert frequencies so that FREQ(K+1)/FREQ(K)
! < SMALL_RAT.
!
! If EDGE_FREQ > SWITCH_FREQ  we insert frequencies so that FREQ(K+1)-FREQ(K)
! > DNU_MAX at the  bound-free edge. The spacing is ramped up by a factor
! G as we move to higher frequencies away from the bound-free edge.
!
! Note the spacing is actually a factor of 2 better than implied by DNU_MAX
! and SMALL_RAT because we insert points midway so that Simpson's rule can
! be used.
!
	    IF(LOW .LT. SWITCH_FREQ)THEN
	      IF(RAT .GT. SMALL_RAT)THEN
	        NPTS=LOG(RAT)/LOG(SMALL_RAT) + 1
	        DELF=(RAT)**(RONE/NPTS)
	        DO J=1,NPTS-1
	          K=K+1
	          NEW_FREQ(K)=NEW_FREQ(K-1)*DELF
	        END DO
	        K=K+1
	        NEW_FREQ(K)=UP
	      ELSE
	        K=K+1
	        NEW_FREQ(K)=UP
	      END IF
	    ELSE IF(DIFF .GT. DNU_MAX)THEN
	      T1=BIG_AMP*DNU_MAX
	      T2=2*DNU_MAX-T1
	      NPTS=0.5D0*(-T2+SQRT(T2*T2+8.0D0*T1*DIFF))/T1
	      DELF=DNU_MAX
	      DO J=1,NPTS-1
	        K=K+1
	        NEW_FREQ(K)=NEW_FREQ(K-1)+DELF
	        DELF=DELF+T1
	      END DO
	      K=K+1
	      NEW_FREQ(K)=UP
	    ELSE
	      K=K+1
	      NEW_FREQ(K)=UP
	    END IF			!< switch_freq if
	  END IF			!Equal if
!
! Now we will insert a few extra points before the bound-free edge.
! This will allow for blending of levels due to level dissolution near the
! bound-free edge.
!
! This is only done for frequencies above MIN_FREQ_LEV_DIS.
!
	  IF(NEW_FREQ(K) .GT. MIN_FREQ_LEV_DIS)THEN
	    K_BEG=K-1
500	    CONTINUE
	    dV=2.998D+05*(NEW_FREQ(K)-NEW_FREQ(K_BEG))/NEW_FREQ(K)
	    IF(dV .GT. dV_LEV .AND. I .LT . N-1)THEN
	      T1= LOG( dV/dV_LEV*(AMP_DIS-RONE) + RONE) / LOG(AMP_DIS)
	      J=NINT(T1)-1
	      IF(K_BEG .NE. 1)THEN
	        IF( dV_LEV*AMP_DIS**(J+1) .LT. 2.998D+05*
	1               (NEW_FREQ(K_BEG)-NEW_FREQ(K_BEG-1))/NEW_FREQ(K))THEN
	          K_BEG=K_BEG-1
	          GOTO 500
	        END IF
	      END IF
	      dV_NEW=dV*(AMP_DIS-RONE)/(AMP_DIS**(J+1) -RONE)
	      UP=NEW_FREQ(K)
	      NEW_FREQ(K_BEG+J+1)=UP
	      T1=0.0D0
	      DO L=J,1,-1
	        T1=T1+dV_NEW*(AMP_DIS**(J-L))/2.998D+05
	        NEW_FREQ(K_BEG+L)=UP*(RONE-T1)
	      END DO
	      K=K_BEG+J+1
	    END IF	
	  END IF
!
	  IF(K .GT. NCF_MAX)THEN
	    WRITE(LUER,*)'Error NCF too small in SET_CONT_FREQ_V4'
	    WRITE(LUER,*)N,NCF_MAX,SMALL_RAT,BIG_AMP,DNU_MAX
	    OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	      DO J=1,NCF
	        WRITE(LUOUT,100)NEW_FREQ(J)
	      END DO
	    CLOSE(LUOUT)
	    STOP
	  END IF
!
	END DO				!I loop
!
	NCF=K
!
! Now sort frequencies into numerical decreasing order. Then check
! frequency array is monotonic.
!
	DO I=1,NCF/2
	  DELF=NEW_FREQ(I)
	  NEW_FREQ(I)=NEW_FREQ(NCF-I+1)
	  NEW_FREQ(NCF-I+1)=DELF
	END DO
!
! Add in extra frequencies to ensure adequate sampling of continuum.
! NU_EVAL is used as a temporary array.
!
	FREQ(1:NCF)=NEW_FREQ(1:NCF)	
        K=1
        DO ML=2,NCF
          T1=C_KMS*(FREQ(ML-1)-FREQ(ML))/FREQ(ML)
          IF(T1 .GT. 1.25D0*DELV_XRAY .AND. FREQ(ML) .GT. NU_END_XRAY)THEN
            J=T1/DELV_XRAY/1.2D0
            DO L=1,J
              K=K+1
              IF(K .GT. NCF_MAX)EXIT
              NEW_FREQ(K)=FREQ(ML-1)-L*(FREQ(ML-1)-FREQ(ML))/(J+1.0D0)
            END DO
          ELSE IF(T1 .GT. 1.25D0*DELV_CONT .AND. FREQ(ML) .LE. NU_END_XRAY)THEN
            J=T1/DELV_CONT/1.2D0
            DO L=1,J
              K=K+1
              IF(K .GT. NCF_MAX)EXIT
              NEW_FREQ(K)=FREQ(ML-1)-L*(FREQ(ML-1)-FREQ(ML))/(J+1.0D0)
            END DO
          END IF
          K=K+1
          IF(K .GT. NCF_MAX)THEN
            WRITE(LUER,*)'Error NCF_MAX is too small in SET_CONT_FREQ_V4'
            WRITE(LUER,*)'NCF_MAX=',NCF_MAX
            STOP
          END IF
          NEW_FREQ(K)=FREQ(ML)
        END DO
        NCF=K
!
	OPEN(UNIT=LUOUT,FILE='CFDAT_OUT',STATUS='UNKNOWN')
	  DO I=2,NCF
	    WRITE(LUOUT,100)NEW_FREQ(I),2.998D+05*(NEW_FREQ(I+1)/NEW_FREQ(I)-1.0D0)
	  END DO
100	  FORMAT(1X,F22.16,2X,ES14.4)
	CLOSE(LUOUT)
!
	DO I=1,NCF-1
	  IF( NEW_FREQ(I) .LE. NEW_FREQ(I+1) )THEN
	    WRITE(LUER,*)'Error in SET_CONT_FREQ_V4 - frequency array not',
	1             ' monotonic'
	    STOP
	  END IF
	END DO
!
! Ensure FREQ array is zeroed (as probably will use OBSF).
!
	DO I=1,NCF
	  FREQ(I)=0.0D0
	END DO
!
	RETURN
	END
