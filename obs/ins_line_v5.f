C
C Subroutine to compute the CMF line frequencies. The spacing of the line 
C frequencies is determined by passed parameters which specify the
C the resonance zone extent, spacing in Doppler widths etc.
C
C The CMF line frequencies are merged with the previously adopted continuum
C frequency set. Unnecessary frequencies are eliminated. Edge frequencies are
C retained.
C
C NU_STRT_LINE and NU_CONT must be ordered monotonically from highest to
C lowest frequency.
c
	SUBROUTINE INS_LINE_V5(
	1		FREQ,LINES_THIS_FREQ,NFREQ,NFREQ_MAX,
	1		NU_LINE,NU_STRT_LINE,VEC_MIN_VDOP,TRANS_TYPE,
	1               LINE_ST_INDX,LINE_END_INDX,N_LINES,
	1		NU_CONT,NCF,FRAC_DOP,VINF,dV_CMF_PROF,
	1               dV_CMF_WING,ES_WING_EXT,R_CMF_WING_EXT,
	1               INCLUDE_LINE_CENTERS)
	IMPLICIT NONE
C
C Altered: 07-Jul-2011: Enhanced error checking and messages.
C Altered: 25-Apr-2010: Changed check of whether lines outside of continuum range.
C                         Previously cut all lines for VINF --> C.
C Altered: 18-Jul-2008: Changed terms of the form 1-a*Vinf/C_kms to 1/(1-a*Vinf/C_kms).
C                         This is to prevent problems when Vinf is close to c. 
C Altered: 13-Apr-2003: Bug fix: dNU replaced by dNU_NEXT in if statement (line 287)
C Created: 22-Dec-1998: Complete rewrite of INS_LINE_V5.
C                       Written to handle lines with different intrinsic
C                         profile extents.
C
	INTEGER NCF,NFREQ_MAX,N_LINES
	INTEGER NFREQ				!Returned
C
C Vectors returned by subroutine:
C
C Line+continuum frequencies
	REAL*8 FREQ(NFREQ_MAX)			!Continuum frequencies
	INTEGER LINES_THIS_FREQ(NFREQ_MAX)	!Indicates that this frequency 
						!  has line contributions,
C
	INTEGER LINE_ST_INDX(N_LINES)		!Start index for the line 
						!  in the NEW frequency array.
	INTEGER LINE_END_INDX(N_LINES)	!End index for the line 
						! in the NEW frequency array.
C
C Passed vectors.
C
	REAL*8 NU_CONT(NCF)		!Continuum frequencies
	REAL*8 NU_LINE(N_LINES)		!Line frequencies
	REAL*8 NU_STRT_LINE(N_LINES)	!Start frequency of resoance zone.
	REAL*8 VEC_MIN_VDOP(N_LINES)	!Minimum doppler velocity for line.
	CHARACTER*(*) TRANS_TYPE(N_LINES)
C
C Passed constants:
	REAL*8 VINF		!Terminal velocity of wind.
	REAL*8 FRAC_DOP		!Indicates dNU across line in Doppler widths.
	REAL*8 dV_CMF_PROF	!Indicate spacing in profile but outside
                                !  resonance zone (in km/s).
	REAL*8 dV_CMF_WING	!Indicate spacing in wings (i.e. outside 
				!  intrinsic profile) (in km/s).
C
C R_CMF_WING_EXT indicates how far profile should extend beyond red edge
C of RESONANCE zone. This will normally be greater than 2VINF as electron
C scattering by cold electrons broadens the line to the red side. Value
C is set on the basis of cold electrons.
C
C ES_WING_EXT is useful when have non-coherent electron scattering.
C Used for both blue and red sides of the line profile.
C
	REAL*8 ES_WING_EXT
	REAL*8 R_CMF_WING_EXT
C
	LOGICAL INCLUDE_LINE_CENTERS
C
C 
	REAL*8 NU_END_LINE(N_LINES)
C
C Local variables.
C
	REAL*8 dNU_on_NU	!Actual spacing used across intrinsic line 
				!  profile given by dNU =NU*dNU_on_NU
C
	REAL*8 ES_BLUE_WING_EXT		!In km/s
	REAL*8 ES_RED_WING_EXT
	REAL*8 CUR_RED_PROF_EXT		!10^15 Hz
	REAL*8 EDGE_SEP_FAC
	REAL*8 MIN_FREQ_RAT
C
	INTEGER INDX		!Current frequency index.
	INTEGER LN_INDX	!Current line whose frequencies we are 
				!   installing.
	INTEGER NUM_RES_LINES
	INTEGER LOCAL_N_LINES
C
	INTEGER ML		!Continuum frequency index
	INTEGER I,J,K		!Miscellaneous loop variables.
	INTEGER LU_ER
	REAL*8 C_KMS
	REAL*8 dNU
	REAL*8 dNU_NEXT
	REAL*8 T1
C
	LOGICAL EDGE_FREQ(NCF)
C
C External functions
C
	INTEGER ERROR_LU
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
C
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	LU_ER=ERROR_LU()
C
C Check validity of some of the parameters.
C
	IF(R_CMF_WING_EXT .LT. 2.0D0 .OR. R_CMF_WING_EXT .GT. 20.0D0)THEN
	  WRITE(LU_ER,*)'Invalid value for R_CMF_WING_EXT in INS_LINE'
	  WRITE(LU_ER,*)'R_CMF_WING_EXT=',R_CMF_WING_EXT
	  STOP
	END IF
	IF(ES_WING_EXT .LT. 0.0D0 .OR. ES_WING_EXT .GT. 20000.0D0)THEN
	  WRITE(LU_ER,*)'Invalid value for ES_WING_EXT in INS_LINE'
	  WRITE(LU_ER,*)'ES_WING_EXT=',ES_WING_EXT
	  STOP
	END IF
	IF(dV_CMF_PROF .LT. 1.0D0 .OR. dV_CMF_PROF .GT. 1.0D+05)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_PROF in INS_LINE'
	  WRITE(LU_ER,*)'dV_CMF_PROF=',dV_CMF_PROF
	  STOP
     	END IF
	IF(dV_CMF_WING .LT. 1.0D0 .OR. dV_CMF_WING .GT. 1.0D+05)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_WING in INS_LINE'
	  WRITE(LU_ER,*)'dV_CMF_WING=',dV_CMF_WING
	  STOP
	END IF
C
	T1=1.0D0-(6.0D0*MAXVAL(VEC_MIN_VDOP)+3.0D0*VINF)/C_KMS
	T1=1.0D0/(1.0D0+(6.0D0*MAXVAL(VEC_MIN_VDOP)+3.0D0*VINF)/C_KMS)
	IF(T1 .LT. 0.5D0)T1=0.5D0
	LOCAL_N_LINES=N_LINES
	DO WHILE(NU_CONT(NCF) .GT. NU_LINE(LOCAL_N_LINES)*T1)
	  LOCAL_N_LINES=LOCAL_N_LINES-1
	END DO
	IF(N_LINES .NE. LOCAL_N_LINES)THEN
	  WRITE(LU_ER,*)'Warning from INS_LINE_V5'
	  WRITE(LU_ER,'(1X,I5,A,A)')N_LINES-LOCAL_N_LINES,
	1        ' weak lines in extreme ',
	1        'IR will be ignored as outside continuum range.'
	  WRITE(LU_ER,*)'Min(Nu_CONT)=',NU_CONT(NCF)
	  WRITE(LU_ER,*)'Min(Nu_LINE)=',NU_LINE(N_LINES)
	END IF
C
C Check that frequencies are monontonically decreaing.
C
	DO I=1,NCF-1
	  IF(NU_CONT(I) .LT. NU_CONT(I+1))THEN
	    WRITE(LU_ER,*)'Error in INS_LINE'
	    WRITE(LU_ER,*)'Continuum frequencies not monotonically decreasing'
	    STOP
	  END IF
	END DO
C
	DO I=1,N_LINES-1
	  IF(NU_STRT_LINE(I) .LT. NU_STRT_LINE(I+1))THEN
	    WRITE(LU_ER,*)'Error in INS_LINE:'
	    WRITE(LU_ER,*)'Start line frequencies not monotonically decreasing'
	    STOP
	  END IF
	END DO
	DO I=1,N_LINES
	  IF(NU_STRT_LINE(I) .LE. NU_LINE(I) .AND. TRANS_TYPE(I) .EQ. 'BLANK')THEN
	    WRITE(LU_ER,*)'Error in INS_LINE:'
	    WRITE(LU_ER,*)'Inconsistent line and start line frequencies'
	    WRITE(LU_ER,*)I-1,NU_LINE(I-1),NU_STRT_LINE(I-1)
	    WRITE(LU_ER,*)I,NU_LINE(I),NU_STRT_LINE(I)
	    WRITE(LU_ER,*)I+1,NU_LINE(I+1),NU_STRT_LINE(I+1)
	    STOP
	  END IF
	END DO
C
C We assume that the resonance zone is symmetrical about line center.
C
	NU_END_LINE(:)=NU_LINE(:)-(NU_STRT_LINE(:)-NU_LINE(:))
C
C 
C
C We assume hat both lines and continuum are ordered from highest to
C lowest frequencies.
C
	dNU_on_NU=FRAC_DOP*MINVAL(VEC_MIN_VDOP)/C_KMS
C
c To avoid numerical instabilities in the iteration procedure when solving 
C for the corrections we ensure that the frequencies bracketing a bound-free 
C edge are EDGE_SEP_FAC*FRAC_DOP Doppler widths appart. We adjust the lower 
C frequency to ensure this. MIN_FREQ_RAT is the minimum ratio allowed between
C successive frequencies.
C
	EDGE_SEP_FAC=0.5D0                              !was 0.1
	MIN_FREQ_RAT=1.0D0+EDGE_SEP_FAC*dNU_on_NU
C
C Define edges of the e.s blue and red wings. The form of ES_RED_WING_EXTENT
C ensures that it is always positive, even when VINF is clse to c.
C
	ES_BLUE_WING_EXT=1.0D0+ES_WING_EXT/C_KMS		!v/v(o)
	ES_RED_WING_EXT=1.0D0/( 1.0D0+(ES_WING_EXT+R_CMF_WING_EXT*VINF)/C_KMS)
C
C Determine continuum frequencies bracketing bound-free edges. We keep
C these in our final continuum list. To avoid numerical instabilities
C in the iteration procedure when solving for the corrections we ensure
C that the frequencies bracketing a bound-free edge are EDGE_SEP_FAC*FRAC_DOP
C Doppler widths appart. We adjust the lower frequency to ensure this.
C
	EDGE_FREQ(1:NCF)=.FALSE.
	I=2
	DO WHILE (I .LT. NCF)
	  IF( ABS(NU_CONT(I-1)/NU_CONT(I)-1.0D0) .LT. 1.0D-06)THEN
!	  IF( ABS(NU_CONT(I-1)/NU_CONT(I)-1.0D0) .LT. 1.0D-04)THEN
	    IF(NU_CONT(I)/MIN_FREQ_RAT .GT. MIN_FREQ_RAT*NU_CONT(I+1))THEN
	      EDGE_FREQ(I-1)=.TRUE.
	      EDGE_FREQ(I)=.TRUE.
	      NU_CONT(I)=NU_CONT(I)/MIN_FREQ_RAT
	      I=I+2
	    ELSE
	      EDGE_FREQ(I-1)=.TRUE.
	      K=1
	      DO WHILE ( NU_CONT(I)/MIN_FREQ_RAT .LT. NU_CONT(I+K))
	        K=K+1
	      END DO
	      EDGE_FREQ(I+K)=.TRUE.
	      I=I+K+1
	    END IF
	  END IF
	  I=I+1
	END DO
C
C Find the first line that is to be included as a blanketed line.
C
	LN_INDX=1
	DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1          TRANS_TYPE(MIN(LN_INDX,N_LINES))(1:3) .NE. 'BLA')
	  LN_INDX=LN_INDX+1
	END DO
C
C If no lines to be included, set FREQ directly to continuum frequencies
C and return.
C
        DO ML=2,NCF
          WRITE(19,'(X,I6,3X,3ES15.6)')ML,NU_CONT(ML),0.01D0*C_KMS/NU_CONT(ML),
	1              C_KMS*(NU_CONT(ML)-NU_CONT(ML-1))/NU_CONT(ML)
	END DO
!
	IF(LN_INDX .GT. N_LINES)THEN
	  FREQ(1:NCF)=NU_CONT(1:NCF)
	  NFREQ=NCF
!	  FREQ(1:2)=NU_CONT(1:2)
!	  J=2
!	  I=3
!	  DO WHILE(I .LE. NCF)
!	    IF( (FREQ(J)-NU_CONT(I)) .GT. 2.5D0*(FREQ(J-1)-FREQ(J)))THEN
!	      FREQ(J+1)=FREQ(J)-2.0D0*(FREQ(J-1)-FREQ(J))
!	    ELSE
!	      FREQ(J+1)=NU_CONT(I)
!	      I=I+1
!	    END IF
!	    J=J+1
!          WRITE(17,'(X,2I6,3X,5ES15.6)')I,J,NU_CONT(I),FREQ(J), 
!	1              (FREQ(J)-NU_CONT(I)),  2.5D0*(FREQ(J-1)-FREQ(J)), FREQ(J+1)
!	  END DO 
!	  NFREQ=J
!	  LU_ER=ERROR_LU()
!	  WRITE(LU_ER,*)'Warning --- no line inserted in INS_LINE_V5'
!
!         DO ML=2,NFREQ
!            WRITE(18,'(X,I6,3X,3ES15.6)')ML,FREQ(ML),0.01D0*C_KMS/FREQ(ML),
!	1                C_KMS*(FREQ(ML)-FREQ(ML-1))/FREQ(ML)
!	  END DO
!	  STOP
!
	  RETURN
	END IF
C 
C
C Compute the frequency grid. We use a combination of CONTINUUM and LINE
C frequencies. Continuum frequencies are inserted simultaneously with the
C line frequencies (rather than after the line grid is created) to avoid 
C the need for extra temporary storage for LINES_THIS_FREQ, and to avoid 
C having to alter LINE_ST_INDX etc.
C
C All edge frequencies are included.
C
	ML=1
	INDX=1
	FREQ(1)=NU_CONT(1)
	CUR_RED_PROF_EXT=10.0D0*NU_CONT(1)
C
	LINES_THIS_FREQ(:)=0
	LINE_ST_INDX(:)=0
	LINE_END_INDX(:)=0
C
	DO WHILE (ML .LE. NCF)
C
C Check that we can include the next frequency.
C
	  IF(INDX+1 .GT. NFREQ_MAX)GOTO 9999
C
          IF( NU_CONT(ML) .GE. FREQ(INDX) )THEN
	     ML=ML+1			!Ignore as past this freq.
C
C If continuum frequency is within 0.2*FRAC_DOP doppler widths of the last set 
C frequency, C there is no need to use it, unless it is a bound-free edge 
C frequency.
C
          ELSE IF( NU_CONT(ML) .GE. FREQ(INDX)/(1.0D0+0.2D0*dNU_on_NU) )THEN
!	1                       .AND. .NOT. EDGE_FREQ(ML) )THEN
	     ML=ML+1			!Use current set frequency.
C
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
C
	    dNU_NEXT=FREQ(INDX)-NU_CONT(ML)
	    IF(INDX .GT. 1)THEN
	      T1=2.0D0*(FREQ(INDX-1)-FREQ(INDX))
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    END IF
C
C Spacing for electron-scattering wings.
C
	    IF(NU_CONT(ML) .LT. NU_STRT_LINE(LN_INDX)*ES_BLUE_WING_EXT .AND.
	1       NU_CONT(ML) .GT. NU_LINE(LN_INDX)*ES_RED_WING_EXT)THEN
	      T1=FREQ(INDX)*dV_CMF_WING/C_KMS
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    END IF
C
	    IF(FREQ(INDX)-dNU_NEXT .GT. CUR_RED_PROF_EXT)THEN
	      T1=FREQ(INDX)*dV_CMF_PROF/C_KMS
	      dNU_NEXT=MIN(dNU_NEXT,T1)
	    END IF
C
C If we are in a resonance zone, spacing will be even smaller. NUM_RES_ZONES
C is used to represent the total number of lines for which the next
C frequency is reward of the blue edge of the resonance zone.
C
	    NUM_RES_LINES=0
	    K=LN_INDX
	    DO WHILE( (FREQ(INDX)-dNU_NEXT) .LE. NU_STRT_LINE(K) )
	      IF(TRANS_TYPE(K)(1:3) .EQ. 'BLA')THEN
	        dNU=0.5D0*FREQ(INDX)*FRAC_DOP*VEC_MIN_VDOP(K)*
	1             SQRT( 1.0D0+ C_KMS*ABS(1.0D0-NU_LINE(K)/FREQ(INDX))/
	1                   VEC_MIN_VDOP(K) )/C_KMS
	        dNU_NEXT=MIN(dNU,dNU_NEXT)
	        LINES_THIS_FREQ(INDX+1)=LINES_THIS_FREQ(INDX+1)+1
	      END IF
	      NUM_RES_LINES=NUM_RES_LINES+1
	      K=K+1
	      IF(K .GT. LOCAL_N_LINES)GOTO 100
	    END DO
100	    CONTINUE
C
C If requested we include the center of each line. We only do this provided
C its separation from the previous frequency is less than 10% of the typical
C frequency spacing around line center (FRAC_DOP*VEC_MIN_VDOP).
C
C As the lines are ordered by start frequency (and not line center) we
C check all NUM_RES_LINES.
C
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
C
	    INDX=INDX+1
	    FREQ(INDX)=FREQ(INDX-1)-dNU_NEXT
C
C Set location of resonace zone. Because the reosnance zones can have
C different widths, we do it in the following statements.
C
	    DO K=LN_INDX,LN_INDX+NUM_RES_LINES-1
	      IF( TRANS_TYPE(K)(1:3) .EQ. 'BLA')THEN
	        IF( FREQ(INDX) .LE. NU_END_LINE(K) .AND. 
	1             LINE_END_INDX(K) .EQ. 0)LINE_END_INDX(K)=INDX
		IF( FREQ(INDX) .LE. NU_STRT_LINE(K) .AND. 
	1               LINE_ST_INDX(K) .EQ. 0)LINE_ST_INDX(K)=INDX
	      END IF
	   END DO
C
C Increase the line index (LN_INDX) if we have gone outside the resonance zone,
C and determine the current extent of the red edge of the profile,
C allowing for velocity shifts.
C
	    DO WHILE( LN_INDX .LE. N_LINES .AND. FREQ(INDX) .LE. 
	1         NU_END_LINE(LN_INDX) )
	      IF(TRANS_TYPE(LN_INDX)(1:3) .EQ. 'BLA')THEN
	        T1=NU_END_LINE(LN_INDX)/(1.0D0+2.0D0*VINF/C_KMS)
	        CUR_RED_PROF_EXT=MIN(CUR_RED_PROF_EXT,T1)
	      END IF
	      LN_INDX=LN_INDX+1
	    END DO
C
C Find the next line that is to be included as a blanketed line.
C
	    DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1          TRANS_TYPE(MIN(LN_INDX,N_LINES))(1:3) .NE. 'BLA')
	       LN_INDX=LN_INDX+1
	    END DO
C	 
	  END IF
	END DO
C
C Set the number of frequencies.
C
	NFREQ=INDX
C
C Test for monotocity of frequencies, and determine MINIMUM frequency
C spacing in velocity space.
C
	T1=10000.0D0
	DO ML=1,NFREQ-1
	  T1=MIN( T1 , C_KMS*(FREQ(ML)-FREQ(ML+1))/FREQ(ML) )
	  IF(FREQ(ML) .LE. FREQ(ML+1))THEN
	    WRITE(LU_ER,*)' Invalid frequency grid computed in INS_LINE_V5'
	    WRITE(LU_ER,*)' ML=',ML
	    DO I=MAX(1,ML-20),MIN(NCF,ML+20)
	      WRITE(LU_ER,*)I,FREQ(I)
	    END DO
	    STOP
	  END IF
	END DO
	WRITE(LU_ER,'(1X,A,1PE9.2,A)')
	1          'Minimum frequency spacing is:',T1,'km/s'
C
C Test that all lines treated in blanketing mode have LINE_ST_INDX and
C LINE_END_INDX defined.
C
	DO I=1,LOCAL_N_LINES
	  IF(TRANS_TYPE(I) .EQ. 'BLANK')THEN
	    IF(LINE_ST_INDX(I) .EQ. 0 .OR. LINE_END_INDX(I) .EQ. 0 .OR.
	1        LINE_END_INDX(I) .LE. LINE_ST_INDX(I) .OR.
	1        LINE_ST_INDX(MAX(I-1,1)) .GT. LINE_ST_INDX(I))THEN
	      WRITE(LU_ER,*)' Invalid LINE_ST_INDX or LINE_END_INDX in INS_LINE_V5'
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
C
	DO I=2,NFREQ
	  WRITE(63,'(1X,I6,1P,2E12.4)')I,FREQ(I),3.0D+05*(FREQ(I-1)-FREQ(I))/FREQ(I-1)
	END DO
C
	DO I=1,LOCAL_N_LINES
	   IF(TRANS_TYPE(I) .EQ. 'BLANK')THEN
	     WRITE(77,'(3I8,1P,4E15.5)')
	1       I,LINE_ST_INDX(I),LINE_END_INDX(I),
	1       NU_LINE(I),0.01D0*C_KMS/NU_LINE(I),
	1      C_KMS*( FREQ(LINE_ST_INDX(I))-NU_LINE(I))/NU_LINE(I),
	1      C_KMS*( NU_LINE(I)-FREQ(LINE_END_INDX(I)) )/NU_LINE(I)
	  END IF
	END DO
C
	RETURN
C
9999	CONTINUE
	LU_ER=ERROR_LU()
	WRITE(LU_ER,*)'Error --- insufficient frequencies to store'//
	1               ' both line and continuum frequencies'
	WRITE(LU_ER,*)'ML= ',ML
	WRITE(LU_ER,*)'LN_INDX= ',LN_INDX
	WRITE(LU_ER,*)'INDX= ',INDX
	WRITE(LU_ER,*)'NFREQ_MAX= ',NFREQ_MAX
	STOP
	END
