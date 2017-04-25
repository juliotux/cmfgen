!
! Subroutine to compute the CMF line frequencies. The spacing of the line 
! frequencies is determined by passed parameters which specify the
! the resonance zone extent, spacing in Doppler widths etc.
!
! The CMF line frequencies are merged with the previously adopted continuum
! frequency set. Unnecessary frequencies are eliminated. Edge frequencies are
! retained.
!
! NU_STRT_LINE and NU_CONT must be ordered monotonically from highest to
! lowest frequency.
!
	SUBROUTINE INS_LINE_V6(
	1		FREQ,LINES_THIS_FREQ,NFREQ,NFREQ_MAX,
	1		NU_LINE,NU_STRT_LINE,VEC_MIN_VDOP,TRANS_TYPE,
	1               LINE_ST_INDX,LINE_END_INDX,N_LINES,
	1		NU_CONT,NCF,FRAC_DOP,VINF,dV_CMF_PROF,
	1               dV_CMF_WING,ES_WING_EXT,R_CMF_WING_EXT,
	1               INCLUDE_LINE_CENTERS)
	IMPLICIT NONE
!
! Altered 06-Jul-2015: dNU was not being compared to dV_CMF_PRO (GT instead of LT).
! Altered 23-Jan-2014: Cleaned and fixed issues with frequencies at bound-free edges.
! Created 03-Jan-2014: Essentially the same as INS_LINE_V5 except specing in resonance zone
!                         is actually FRAC_DOP rather than 05*FRAC_DOP. INS_LINE_V5 was left
!                         unchanged so that CMF_FLUX did not change.
!
	INTEGER NCF,NFREQ_MAX,N_LINES
	INTEGER NFREQ				!Returned
!
! Vectors returned by subroutine:
!
! Line+continuum frequencies
!
	REAL*8 FREQ(NFREQ_MAX)			!Continuum frequencies
	INTEGER LINES_THIS_FREQ(NFREQ_MAX)	!Indicates that this frequency 
						!  has line contributions,
	INTEGER LINE_ST_INDX(N_LINES)		!Start index for the line 
						!  in the NEW frequency array.
	INTEGER LINE_END_INDX(N_LINES)		!End index for the line 
						! in the NEW frequency array.
!
! Passed vectors.
!
	REAL*8 NU_CONT(NCF)			!Continuum frequencies
	REAL*8 NU_LINE(N_LINES)			!Line frequencies
	REAL*8 NU_STRT_LINE(N_LINES)		!Start frequency of resoance zone.
	REAL*8 VEC_MIN_VDOP(N_LINES)		!Minimum doppler velocity for line.
	CHARACTER*(*) TRANS_TYPE(N_LINES)
!
! Passed constants:
	REAL*8 VINF		!Terminal velocity of wind.
	REAL*8 FRAC_DOP		!Indicates dNU across line in Doppler widths.
	REAL*8 dV_CMF_PROF	!Indicate spacing in profile but outside
                                !  resonance zone (in km/s).
	REAL*8 dV_CMF_WING	!Indicate spacing in wings (i.e. outside 
				!  intrinsic profile) (in km/s).
!
! R_CMF_WING_EXT indicates how far profile should extend beyond red edge
! of RESONANCE zone. This will normally be greater than 2VINF as electron
! scattering by cold electrons broadens the line to the red side. Value
! is set on the basis of cold electrons.
!
! ES_WING_EXT is useful when have non-coherent electron scattering.
! Used for both blue and red sides of the line profile.
!
	REAL*8 ES_WING_EXT
	REAL*8 R_CMF_WING_EXT
!
	LOGICAL INCLUDE_LINE_CENTERS
!
! 
	REAL*8 NU_END_LINE(N_LINES)
!
! Local variables.
!
	REAL*8 dNU_on_NU	!Actual spacing used across intrinsic line 
				!  profile given by dNU =NU*dNU_on_NU
!
	REAL*8 ES_BLUE_WING_EXT		!In km/s
	REAL*8 ES_RED_WING_EXT
	REAL*8 CUR_RED_PROF_EXT		!10^15 Hz
	REAL*8 EDGE_SEP_FAC
	REAL*8 MIN_FREQ_RAT
!
	INTEGER INDX		!Current frequency index.
	INTEGER LN_INDX		!Current line whose frequencies we are installing.
	INTEGER NUM_RES_LINES
	INTEGER LOCAL_N_LINES
!
	INTEGER ML		!Continuum frequency index
	INTEGER I,J,K		!Miscellaneous loop variables.
	INTEGER LU_ER
	INTEGER LU
	REAL*8 C_KMS
	REAL*8 dNU
	REAL*8 dNU_NEXT
	REAL*8 T1
!
	LOGICAL EDGE_FREQ(NCF)
	LOGICAL VERBOSE
!
! External functions
!
	INTEGER ERROR_LU
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	LU_ER=ERROR_LU()
	CALL GET_VERBOSE_INFO(VERBOSE)
!
! Check validity of some of the parameters.
!
	IF(R_CMF_WING_EXT .LT. 2.0D0 .OR. R_CMF_WING_EXT .GT. 20.0D0)THEN
	  WRITE(LU_ER,*)'Invalid value for R_CMF_WING_EXT in INS_LINE_V6'
	  WRITE(LU_ER,*)'R_CMF_WING_EXT=',R_CMF_WING_EXT
	  STOP
	END IF
	IF(ES_WING_EXT .LT. 0.0D0 .OR. ES_WING_EXT .GT. 20000.0D0)THEN
	  WRITE(LU_ER,*)'Invalid value for ES_WING_EXT in INS_LINE_V6'
	  WRITE(LU_ER,*)'ES_WING_EXT=',ES_WING_EXT
	  STOP
	END IF
	IF(dV_CMF_PROF .LT. 1.0D0 .OR. dV_CMF_PROF .GT. 1.0D+05)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_PROF in INS_LINE_V6'
	  WRITE(LU_ER,*)'dV_CMF_PROF=',dV_CMF_PROF
	  STOP
     	END IF
	IF(dV_CMF_WING .LT. 1.0D0 .OR. dV_CMF_WING .GT. 1.0D+05)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_WING in INS_LINE_V6'
	  WRITE(LU_ER,*)'dV_CMF_WING=',dV_CMF_WING
	  STOP
	END IF
!
	T1=1.0D0-(6.0D0*MAXVAL(VEC_MIN_VDOP)+3.0D0*VINF)/C_KMS
	T1=1.0D0/(1.0D0+(6.0D0*MAXVAL(VEC_MIN_VDOP)+3.0D0*VINF)/C_KMS)
	T1=1.0D0/(1.0D0+(6.0D0*MAXVAL(VEC_MIN_VDOP)+3.0D0*MAX(VINF,3000.0D0))/C_KMS)
	IF(T1 .LT. 0.5D0)T1=0.5D0
	LOCAL_N_LINES=N_LINES
	DO WHILE(NU_CONT(NCF) .GT. NU_LINE(LOCAL_N_LINES)*T1)
	  LOCAL_N_LINES=LOCAL_N_LINES-1
	END DO
	IF(N_LINES .NE. LOCAL_N_LINES)THEN
	  WRITE(LU_ER,'(/,A)')' Warning from INS_LINE_V6'
	  WRITE(LU_ER,'(1X,I5,A,A)')N_LINES-LOCAL_N_LINES,
	1        ' weak lines in extreme ',
	1        'IR will be ignored as outside continuum range.'
	  WRITE(LU_ER,*)'    Min(Nu_CONT)=',NU_CONT(NCF)
	  WRITE(LU_ER,*)'    Min(Nu_LINE)=',NU_LINE(N_LINES)
	END IF
!
! Check that frequencies are monontonically decreaing.
!
	DO I=1,NCF-1
	  IF(NU_CONT(I) .LT. NU_CONT(I+1))THEN
	    WRITE(LU_ER,*)'Error in INS_LINE_V6'
	    WRITE(LU_ER,*)'Continuum frequencies not monotonically decreasing'
	    STOP
	  END IF
	END DO
!
	DO I=1,N_LINES-1
	  IF(NU_STRT_LINE(I) .LT. NU_STRT_LINE(I+1))THEN
	    WRITE(LU_ER,*)'Error in INS_LINE_V6:'
	    WRITE(LU_ER,*)'Start line frequencies not monotonically decreasing'
	    STOP
	  END IF
	END DO
	DO I=1,N_LINES
	  IF(NU_STRT_LINE(I) .LE. NU_LINE(I) .AND. TRANS_TYPE(I) .EQ. 'BLANK')THEN
	    WRITE(LU_ER,*)'Error in INS_LINE_V6:'
	    WRITE(LU_ER,*)'Inconsistent line and start line frequencies'
	    WRITE(LU_ER,*)I-1,NU_LINE(I-1),NU_STRT_LINE(I-1)
	    WRITE(LU_ER,*)I,NU_LINE(I),NU_STRT_LINE(I)
	    WRITE(LU_ER,*)I+1,NU_LINE(I+1),NU_STRT_LINE(I+1)
	    STOP
	  END IF
	END DO
!
! We assume that the resonance zone is symmetrical about line center.
!
	NU_END_LINE(:)=NU_LINE(:)-(NU_STRT_LINE(:)-NU_LINE(:))
!
! 
!
! We assume hat both lines and continuum are ordered from highest to
! lowest frequencies.
!
	dNU_on_NU=FRAC_DOP*MINVAL(VEC_MIN_VDOP)/C_KMS
!
        WRITE(6,*)'In INS_LINE_V6 FRAC_DOP is',FRAC_DOP
        WRITE(6,*)'In INS_LINE_V6 VDop is',MINVAL(VEC_MIN_VDOP)
!
! To avoid numerical instabilities in the iteration procedure when solving 
! for the corrections we ensure that the frequencies bracketing a bound-free 
! edge are EDGE_SEP_FAC*FRAC_DOP Doppler widths appart. We adjust the lower 
! frequency to ensure this. MIN_FREQ_RAT is the minimum ratio allowed between
! successive frequencies.
!
	EDGE_SEP_FAC=0.1D0
	MIN_FREQ_RAT=1.0D0+EDGE_SEP_FAC*dNU_on_NU
!
! Define edges of the e.s blue and red wings. The form of ES_RED_WING_EXTENT
! ensures that it is always positive, even when VINF is clse to c.
!
	ES_BLUE_WING_EXT=1.0D0+ES_WING_EXT/C_KMS		!v/v(o)
	ES_RED_WING_EXT=1.0D0/( 1.0D0+(ES_WING_EXT+R_CMF_WING_EXT*VINF)/C_KMS)
!
! Determine continuum frequencies bracketing bound-free edges. We keep
! these in our final continuum list. To avoid numerical instabilities
! in the iteration procedure when solving for the corrections we ensure
! that the frequencies bracketing a bound-free edge are EDGE_SEP_FAC*FRAC_DOP
! Doppler widths appart. We adjust the lower frequency to ensure this.
!
	EDGE_FREQ(1:NCF)=.FALSE.
	I=2
	DO WHILE (I .LT. NCF)
	  IF( ABS(NU_CONT(I-1)/NU_CONT(I)-1.0D0) .LT. 1.0D-08)THEN
	    EDGE_FREQ(I-1)=.TRUE.
	      EDGE_FREQ(I)=.TRUE.
	    I=I+2
	  ELSE
	    I=I+1
	  END IF
	END DO
!
! Find the first line that is to be included as a blanketed line.
!
	LN_INDX=1
	DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1          TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	  LN_INDX=LN_INDX+1
	END DO
!
! If no lines to be included, set FREQ directly to continuum frequencies
! and return.
!
	IF(VERBOSE)THEN
	  CALL GET_LU(LU,'LU in INS_LINE_V6')
	  OPEN(UNIT=LU,FILE='FREQUENCY_CONT_SEP',STATUS='UNKNOWN')
            DO ML=2,NCF
              WRITE(LU,'(X,I6,3X,3ES15.6,3X,L1)')ML,NU_CONT(ML),0.01D0*C_KMS/NU_CONT(ML),
	1             C_KMS*(NU_CONT(ML)-NU_CONT(ML-1))/NU_CONT(ML),EDGE_FREQ(ML)
	    END DO
	  CLOSE(LU)
	END IF
!
	IF(LN_INDX .GT. N_LINES)THEN
	  FREQ(1:NCF)=NU_CONT(1:NCF)
	  NFREQ=NCF
	  RETURN
	END IF
! 
!
! Compute the frequency grid. We use a combination of CONTINUUM and LINE
! frequencies. Continuum frequencies are inserted simultaneously with the
! line frequencies (rather than after the line grid is created) to avoid 
! the need for extra temporary storage for LINES_THIS_FREQ, and to avoid 
! having to alter LINE_ST_INDX etc.
!
! All edge frequencies are included.
!
	ML=1
	INDX=1
	FREQ(1)=NU_CONT(1)
	CUR_RED_PROF_EXT=10.0D0*NU_CONT(1)
!
	LINES_THIS_FREQ(:)=0
	LINE_ST_INDX(:)=0
	LINE_END_INDX(:)=0
!
	dNU_NEXT=0
	DO WHILE (ML .LE. NCF)
!
! Check that we can include the next frequency.
!
	  IF(INDX+1 .GT. NFREQ_MAX)GOTO 9999
!
          IF( NU_CONT(ML) .GE. FREQ(INDX))THEN
	     IF(NU_CONT(ML)-dNU_ON_NU .LE. FREQ(INDX) .AND. EDGE_FREQ(ML))THEN
	       dNU_NEXT=FREQ(INDX)*(1.0D0-1.0D0/MIN_FREQ_RAT)
	       INDX=INDX+1
	       FREQ(INDX)=FREQ(INDX-1)-dNU_NEXT
	     ELSE 
	     END IF
	     ML=ML+1			!Ignore as past this freq.
!
! If continuum frequency is within 0.2*FRAC_DOP doppler widths of the last set 
! frequency, C there is no need to use it, unless it is a bound-free edge 
! frequency.
!
          ELSE IF( NU_CONT(ML) .GE. FREQ(INDX)*(1.0D0-dNU_ON_NU)
	1                       .AND. .NOT. EDGE_FREQ(ML) )THEN
	     ML=ML+1			!Use current set frequency.
!
          ELSE IF( NU_CONT(ML)*MIN_FREQ_RAT .GE. FREQ(INDX) .AND. EDGE_FREQ(ML))THEN
	     dNU_NEXT=FREQ(INDX)*(1.0D0-1.0D0/MIN_FREQ_RAT)
	     ML=ML+1
	    INDX=INDX+1
	    FREQ(INDX)=FREQ(INDX-1)-dNU_NEXT
!
	  ELSE IF( LN_INDX .GT. LOCAL_N_LINES)THEN
!
	    dNU_NEXT=FREQ(INDX)-NU_CONT(ML)
	    IF(INDX .GT. 1)THEN
	      T1=2.0D0*(FREQ(INDX-1)-FREQ(INDX))
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    END IF
	    IF(NU_CONT(ML) .GT. CUR_RED_PROF_EXT)THEN
	      T1=FREQ(INDX)*dV_CMF_PROF/C_KMS
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    ELSE IF(NU_CONT(ML) .GT. CUR_RED_PROF_EXT*ES_RED_WING_EXT*
	1                    (1.0D0+2.0D0*VINF/C_KMS) )THEN
	      T1=FREQ(INDX)*dV_CMF_WING/C_KMS
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    END IF
	    INDX=INDX+1
	    FREQ(INDX)=FREQ(INDX-1)-dNU_NEXT
!
	  ELSE
!
	    dNU_NEXT=FREQ(INDX)-NU_CONT(ML)
!
! Spacing for electron-scattering wings.
!
	    IF(NU_CONT(ML) .LT. NU_STRT_LINE(LN_INDX)*ES_BLUE_WING_EXT .AND.
	1       NU_CONT(ML) .GT. NU_LINE(LN_INDX)*ES_RED_WING_EXT)THEN
	      T1=FREQ(INDX)*dV_CMF_WING/C_KMS
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    END IF
!
	    IF(FREQ(INDX)-dNU_NEXT .LT. CUR_RED_PROF_EXT)THEN
	      T1=FREQ(INDX)*dV_CMF_PROF/C_KMS
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    END IF
!
! If we are in a resonance zone, spacing will be even smaller. NUM_RES_ZONES
! is used to represent the total number of lines for which the next
! frequency is reward of the blue edge of the resonance zone.
!
	    NUM_RES_LINES=0
	    K=LN_INDX
	    DO WHILE( (FREQ(INDX)-dNU_NEXT) .LE. NU_STRT_LINE(K) )
	      IF(TRANS_TYPE(K)(1:3) .EQ. 'BLA')THEN
	        dNU=FREQ(INDX)*FRAC_DOP*VEC_MIN_VDOP(K)/C_KMS
	        dNU_NEXT=MIN(dNU,dNU_NEXT)
	        LINES_THIS_FREQ(INDX+1)=LINES_THIS_FREQ(INDX+1)+1
	      END IF
	      NUM_RES_LINES=NUM_RES_LINES+1
	      K=K+1
	      IF(K .GT. LOCAL_N_LINES)GOTO 100
	    END DO
100	    CONTINUE
!
! If requested we include the center of each line. We only do this provided
! its separation from the previous frequency is less than 10% of the typical
! frequency spacing around line center (FRAC_DOP*VEC_MIN_VDOP).
!
! As the lines are ordered by start frequency (and not line center) we
! check all NUM_RES_LINES.
!
	    IF(INCLUDE_LINE_CENTERS)THEN
	      DO K=LN_INDX,LN_INDX+NUM_RES_LINES-1
	         IF(TRANS_TYPE(K)(1:3) .EQ. 'BLA')THEN
	           dNU=FREQ(INDX)-NU_LINE(K)
	           IF(dNU .GT. FRAC_DOP*VEC_MIN_VDOP(K)*FREQ(INDX)*0.1D0)THEN
	             dNU_NEXT=MIN(dNU_NEXT,dNU)
	           END IF
	         END IF
	      END DO
	    END IF
!
	    INDX=INDX+1
	    FREQ(INDX)=FREQ(INDX-1)-dNU_NEXT
!
! Set location of resonance zone. Because the reosnance zones can have
! different widths, we do it in the following statements.
!
	    DO K=LN_INDX,LN_INDX+NUM_RES_LINES-1
	      IF( TRANS_TYPE(K)(1:3) .EQ. 'BLA')THEN
	        IF( FREQ(INDX) .LE. NU_END_LINE(K) .AND. 
	1             LINE_END_INDX(K) .EQ. 0)LINE_END_INDX(K)=INDX
		IF( FREQ(INDX) .LE. NU_STRT_LINE(K) .AND. 
	1               LINE_ST_INDX(K) .EQ. 0)LINE_ST_INDX(K)=INDX
	      END IF
	   END DO
!
! Increase the line index (LN_INDX) if we have gone outside the resonance zone,
! and determine the current extent of the red edge of the profile,
! allowing for velocity shifts.
!
	    DO WHILE( LN_INDX .LE. N_LINES .AND. FREQ(INDX) .LE. 
	1         NU_END_LINE(LN_INDX) )
	      IF(TRANS_TYPE(LN_INDX)(1:3) .EQ. 'BLA')THEN
	        T1=NU_END_LINE(LN_INDX)/(1.0D0+2.0D0*VINF/C_KMS)
	        CUR_RED_PROF_EXT=MIN(CUR_RED_PROF_EXT,T1)
	      END IF
	      LN_INDX=LN_INDX+1
	    END DO
!
! Find the next line that is to be included as a blanketed line.
!
	    DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1          TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	       LN_INDX=LN_INDX+1
	    END DO
!	 
	  END IF
	END DO
!
! Set the number of frequencies.
!
	NFREQ=INDX
!
! Test for monotocity of frequencies, and determine MINIMUM frequency
! spacing in velocity space.
!
	T1=10000.0D0
	DO ML=1,NFREQ-1
	  T1=MIN( T1 , C_KMS*(FREQ(ML)-FREQ(ML+1))/FREQ(ML) )
	  IF(FREQ(ML) .LE. FREQ(ML+1))THEN
	    WRITE(LU_ER,*)' Invalid frequency grid computed in INS_LINE_V6'
	    WRITE(LU_ER,*)' ML=',ML
	    DO I=MAX(1,ML-20),MIN(NCF,ML+20)
	      WRITE(LU_ER,*)I,FREQ(I)
	    END DO
	    STOP
	  END IF
	END DO
	WRITE(LU_ER,'(1X,A,1PE9.2,A)')
	1          'Minimum frequency spacing is:',T1,'km/s'
!
! Test that all lines treated in blanketing mode have LINE_ST_INDX and
! LINE_END_INDX defined.
!
	DO I=1,LOCAL_N_LINES
	  IF(TRANS_TYPE(I) .EQ. 'BLANK')THEN
	    IF(LINE_ST_INDX(I) .EQ. 0 .OR. LINE_END_INDX(I) .EQ. 0 .OR.
	1        LINE_END_INDX(I) .LE. LINE_ST_INDX(I) .OR.
	1        LINE_ST_INDX(MAX(I-1,1)) .GT. LINE_ST_INDX(I))THEN
	      WRITE(LU_ER,*)' Invalid LINE_ST_INDX or LINE_END_INDX in INS_LINE_V6'
	      WRITE(LU_ER,*)'             INDX=',I
		      WRITE(LU_ER,*)'    LOCAL_N_LINES=',LOCAL_N_LINES
		      WRITE(LU_ER,*)'          N_LINES=',N_LINES
	      WRITE(LU_ER,*)'            NFREQ=',NFREQ
	      WRITE(LU_ER,*)'        NFREQ_MAX=',NFREQ_MAX
	      WRITE(LU_ER,*)'      FREQ(ML_ST)=',FREQ(LINE_ST_INDX(I))
	      WRITE(LU_ER,*)'          NU_LINE=',NU_LINE(I)
	      WRITE(LU_ER,*)'     NU_STRT_LINE=',NU_STRT_LINE(I)
	      WRITE(LU_ER,*)'LINE_ST_INDX(I-1)=',LINE_ST_INDX(MAX(I-1,1))
	      WRITE(LU_ER,*)'  LINE_ST_INDX(I)=',LINE_ST_INDX(I)
	      WRITE(LU_ER,*)' LINE_END_INDX(I)=',LINE_END_INDX(I)
	      STOP
	    END IF
	  END IF
	END DO
!
	IF(VERBOSE)THEN
	  OPEN(UNIT=LU,FILE='FREQUENCY_GRID',STATUS='UNKNOWN')
	  DO I=2,NFREQ
	    WRITE(63,'(1X,I6,1P,2E12.4)')I,FREQ(I),3.0D+05*(FREQ(I-1)-FREQ(I))/FREQ(I-1)
	  END DO
	  CLOSE(UNIT=LU)
	END IF
!
	IF(VERBOSE)THEN
	  OPEN(UNIT=LU,FILE='LINE_FREQUENCY_DATA',STATUS='UNKNOWN')
	    DO I=1,LOCAL_N_LINES
	      IF(TRANS_TYPE(I) .EQ. 'BLANK')THEN
	        WRITE(77,'(3I8,1P,4E15.5)')
	1          I,LINE_ST_INDX(I),LINE_END_INDX(I),
	1          NU_LINE(I),0.01D0*C_KMS/NU_LINE(I),
	1          C_KMS*( FREQ(LINE_ST_INDX(I))-NU_LINE(I))/NU_LINE(I),
	1          C_KMS*( NU_LINE(I)-FREQ(LINE_END_INDX(I)) )/NU_LINE(I)
	      END IF
	    END DO
	  CLOSE(LU)
	END IF
!
	RETURN
!
9999	CONTINUE
	LU_ER=ERROR_LU()
	WRITE(LU_ER,*)'Error --- insufficient frequencies to store'//
	1               ' both line and continuum frequencies'
	WRITE(LU_ER,*)'ML(cont. index)= ',ML
	WRITE(LU_ER,*)'   NCF(# cont.)= ',NCF
	WRITE(LU_ER,*)'        LN_INDX= ',LN_INDX
	WRITE(LU_ER,*)'          NLINE= ',N_LINES
	WRITE(LU_ER,*)'           INDX= ',INDX
	WRITE(LU_ER,*)'       NFREQ_MAX= ',NFREQ_MAX
!
	STOP
	END
