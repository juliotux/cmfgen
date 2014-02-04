!
! Subroutine to compute the CMF line frequencies. The spacing of the line 
C frequencies is determined by passed parameters which specify the
C the resonance zone extent, spacing in Doppler widths etc.
C
C The CMF line frequencies are merged with the previously adopted continuum
C frequency set. Unnecessary frequencies are eliminated. Edge frequencies are
C retained.
C
	SUBROUTINE INS_LINE_V4(
	1		FREQ,LINES_THIS_FREQ,NFREQ,NFREQ_MAX,
	1		NU_LINE,TRANS_TYPE,
	1               LINE_ST_INDX,LINE_END_INDX,N_LINES,
	1		NU_CONT,NCF,
	1		V_DOP,FRAC_DOP,MAX_DOP,VINF,dV_CMF_PROF,
	1               dV_CMF_WING,ES_WING_EXT,R_CMF_WING_EXT)
	IMPLICIT NONE
C
C Altered 12-Dec-1997 : If no lines are to be inserted, FREQ is set to
C                         NU_CONT, and there is an immediate return.
C                         All edge frequencies [defined by 
C                             ABS( NU(ML)/MU(M+1)-1 ) < 0.000001 ]
C                         are now retained in the list.
C                         Changes based/requested by PACO.
C
C Altered 15-Aug-1997 : Bug fix: Inserting to many pints in profile zone
C                         when some lines not in BLANK mode. LST_LN_INDX
C                         variable now used.
C Altered 05-Dec-1996 : NDOP made intgere (bug fix)
C Altered 30-May-1996 : Warning about lines below minimum continuum frequency
C                         installed.
C Altered 15-May-1996 : ES_WING_EXT introduced (Now V4). May be zero if
C                         coherent scattering (in km/s).
C                       Some cleaning done to minimize frequencies which are
C                         unnecessarily close. 
C                       R_CMF_WING_EXT now refers only to the wing caused by
C                         coherent scattering.
C                        
C Altered 28-Feb-1995 : dv_CMF_WING,CMF_WING_EXT inserted (_V2)
C Altered 02-Feb-1995 : Not correctly inserting all continuum points.
C
	INTEGER*4 NCF,NFREQ_MAX,N_LINES
	INTEGER*4 NFREQ				!Returned
C
C Vectors returned by subroutine:
C
C Line+continuum frequencies
	REAL*8 FREQ(NFREQ_MAX)			!Continuum frequencies
	INTEGER*4 LINES_THIS_FREQ(NFREQ_MAX)	!Indicates that this frequency 
						!  has line contributions,
C
	INTEGER*4 LINE_ST_INDX(N_LINES)		!Start index for the line 
						!  in the NEW frequency array.
	INTEGER*4 LINE_END_INDX(N_LINES)	!End index for the line 
						! in the NEW frequency array.
C
C Passed vectors.
C
	REAL*8 NU_CONT(NCF)		!Continuum frequencies
	REAL*8 NU_LINE(N_LINES)		!Line frequencies
	CHARACTER*(*) TRANS_TYPE(N_LINES)
C
C Passed constants:
	REAL*8 VINF		!Terminal velocity of wind.
	REAL*8 V_DOP		!Doppler velocity (km/s).
	REAL*8 FRAC_DOP		!Indicates dNU across line in Doppler widths.
	REAL*8 MAX_DOP		!Half the extent of intrinsic profile
				!  in Doppler widths,
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
C Local variables.
C
	REAL*8 dNU_on_NU	!Actual spacing used across intrinsic line 
				!  profile given by dNU =NU*dNU_on_NU
	INTEGER*4 NDOP		!Number of frequencies across intrinsic 
				!  profile.
	REAL*8 RES_EXTENT	!Maximum frequency in line is NU*RES_EXTENT
C
	REAL*8 BLUE_WING_EXT	!In km/s
	REAL*8 RED_WING_EXT
	REAL*8 APP_RES_EXT	!No units.
	REAL*8 APP_ESBW_EXT
	REAL*8 EDGE_SEP_FAC
	REAL*8 MIN_FREQ_RAT
C
	INTEGER*4 INDX		!Current frequency index.
	INTEGER*4 LN_INDX	!Current line whose frequencies we are 
				!   installing.
	INTEGER*4 LST_LN_INDX	!Last line whose frequencies we 
				!   installed.
C
	INTEGER*4 ML		!Continuum frequency index
	INTEGER*4 I,J,K		!Miscellaneous loop variables.
	INTEGER*4 LU_ER
	REAL*8 C_KMS
	REAL*8 DELF
	REAL*8 MIN_FREQ
	REAL*8 SWITCH_FREQ
	REAL*8 TEMP_FREQ
	REAL*8 T1
C
	LOGICAL EDGE_FREQ(NCF)
C
C External functions
C
	INTEGER*4 ERROR_LU
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
C
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	LU_ER=ERROR_LU()
C
C Check validity of some of the parameters.
C
	IF(R_CMF_WING_EXT .LT.2 .OR. R_CMF_WING_EXT .GT. 20)THEN
	  WRITE(LU_ER,*)'Invalid value for R_CMF_WING_EXT in INS_LINE'
	  WRITE(LU_ER,*)'R_CMF_WING_EXT=',R_CMF_WING_EXT
	  STOP
	END IF
	IF(ES_WING_EXT .LT. 0 .OR. ES_WING_EXT .GT. 20000)THEN
	  WRITE(LU_ER,*)'Invalid value for ES_WING_EXT in INS_LINE'
	  WRITE(LU_ER,*)'ES_WING_EXT=',ES_WING_EXT
	  STOP
	END IF
	IF(MAX_DOP .LT.2 .OR. MAX_DOP .GT. 20)THEN
	  WRITE(LU_ER,*)'Invalid value for MAX_DOP in INS_LINE'
	  WRITE(LU_ER,*)'MAX_DOP=',MAX_DOP
	  STOP
	END IF
	IF(dV_CMF_PROF .LT. 1.0 .OR. dV_CMF_PROF .GT. 1.0E+05)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_PROF in INS_LINE'
	  WRITE(LU_ER,*)'dV_CMF_PROF=',dV_CMF_PROF
	  STOP
     	END IF
	IF(dV_CMF_WING .LT. 1.0 .OR. dV_CMF_WING .GT. 1.0E+05)THEN
	  WRITE(LU_ER,*)'Invalid value for dV_CMF_WING in INS_LINE'
	  WRITE(LU_ER,*)'dV_CMF_WING=',dV_CMF_WING
	  STOP
	END IF
	IF(NU_CONT(NCF) .GT. NU_LINE(N_LINES))THEN
	  LN_INDX=N_LINES-1
	  DO WHILE(NU_CONT(NCF) .GT. NU_LINE(LN_INDX))
	    LN_INDX=LN_INDX-1
	  END DO
	  WRITE(LU_ER,*)'Warning from INS_LINE_V4'
	  WRITE(LU_ER,'(1X,I5,A,A)')N_LINES-LN_INDX,' weak lines in ',
	1        'extreme IR will be ignored as outside continuum range.'
	  WRITE(LU_ER,*)'Min(Nu_CONT)=',NU_CONT(NCF)
	  WRITE(LU_ER,*)'Min(Nu_LINE)=',NU_LINE(N_LINES)
	END IF
C
C 
C
C We assume that both lines and continuum are ordered from highest to
C lowest frequencies.
C
	dNU_on_NU=FRAC_DOP*V_DOP/C_KMS
	NDOP=NINT(MAX_DOP/FRAC_DOP)
	RES_EXTENT=(1.0D0+dNU_on_NU)**NDOP
C
c To avoid numerical instabilities in the iteration procedure when solving 
C for the corrections we ensure that the frequencies bracketing a bound-free 
C edge are EDGE_SEP_FAC*FRAC_DOP Doppler widths appart. We adjust the lower 
C frequency to ensure this. MIN_FREQ_RAT is the minimum ratio allowed between
C successive frequencies.
C
	EDGE_SEP_FAC=0.1D0
	MIN_FREQ_RAT=1.0D0+EDGE_SEP_FAC*dNU_on_NU
C
C Define approximate edges of the resonance zone and the e.s blue wing.
C These limit getting frequencies unnecessarily close.
C
	BLUE_WING_EXT=ES_WING_EXT+(NDOP+2)*V_DOP*FRAC_DOP	!In km/s
	RED_WING_EXT=ES_WING_EXT+R_CMF_WING_EXT*VINF
C                                                  
	APP_RES_EXT=RES_EXTENT*(1+0.5D0*dNU_on_NU)	!No units
	APP_ESBW_EXT=1.0D0+(BLUE_WING_EXT+0.2D0*dV_CMF_WING)/C_KMS
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
	1          TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	  LN_INDX=LN_INDX+1
	END DO
C
C If no lines to be included, set FREQ directly to continuum frequencies
C and return.
C
	IF(LN_INDX .GT. N_LINES)THEN
	  FREQ(1:NCF)=NU_CONT(1:NCF)
	  NFREQ=NCF
	  LU_ER=ERROR_LU()
	  WRITE(LU_ER,*)'Warning --- no line inserted in INS_LINE_V4'
	  RETURN
	END IF
C 
C
C Compute the frequency grid. We use a combination of CONTINUUM and LINE
C frequencies. Continuum frequencies are inserted simultaneously with the
C line frequencies (rather than after the line grid is created) to avoid 
C the need for extra temporary storage for LINE_THIS_FREQ, and to avoid 
C having to alter LINE_ST_INDX etc.
C
C All edge frequencies are included.
C
	ML=1
	INDX=1
	FREQ(1)=NU_CONT(1)
	LINES_THIS_FREQ(1)=0
C
	DO WHILE (ML .LE. NCF)
          IF( NU_CONT(ML) .GE. FREQ(INDX) )THEN
	     ML=ML+1			!Ignore as past this freq.
C
C If continuum frequency is within 0.2 doppler widths of last set frequency,
C there is no need to use it, unless it is a bound-free edge frequency.
C
          ELSE IF( NU_CONT(ML) .GE. FREQ(INDX)/(1.0D0+0.2D0*dNU_on_NU) 
	1                       .AND. .NOT. EDGE_FREQ(ML) )THEN
	     ML=ML+1			!Use current set frequency.
C
	  ELSE IF( LN_INDX .GT. N_LINES .OR. NU_CONT(ML) .GT. 
	1          NU_LINE(LN_INDX)*MAX(APP_RES_EXT,APP_ESBW_EXT) )THEN
C
C No nearby line, so continuum point becomes next frequency point.
C
	    INDX=INDX+1
	    IF(INDX .GT. NFREQ_MAX)GOTO 9999
	    FREQ(INDX)=NU_CONT(ML)
	    LINES_THIS_FREQ(INDX)=0
	    ML=ML+1
	  ELSE
C
C Add in additional frequencies for this line.
C
C First, we allow for electron scattering on the blue side of the line profile.
C Strictly only need if scattering is non-coherent. Insertion of frequencies
C only occurs until we hit the resonance zone.
C
C Insert points at beginning (high frequency side) of electron scattering
C wing if needed.
C
	    IF(FREQ(INDX) .GT. NU_LINE(LN_INDX)*APP_ESBW_EXT)THEN
	      T1=NU_LINE(LN_INDX)*(1+BLUE_WING_EXT/C_KMS)
	      DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND. 
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=NU_CONT(ML)
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	        ML=ML+1
	      END DO
              IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=T1
	        LINES_THIS_FREQ(INDX)=0
	      END IF
	    END IF
C
C Now insert points in the blue e.s. wing, until we reach the resonance
C zone.
C
! NB: The following slab of code simple ensures that all EDGE frequencies
!     are include in revised frequency list. Code similar to this must
!     be enacted for each new frequency. T1 is the new frequency to be set.
!
!	      DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
!	        IF(EDGE_FREQ(ML))THEN
!	          INDX=INDX+1
!	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
!	          FREQ(INDX)=NU_CONT(ML)
!	          LINES_THIS_FREQ(INDX)=0
!	        END IF
!	        ML=ML+1
!	      END DO
C
	    T1=FREQ(INDX)/(1.0D0+1.1*dV_CMF_WING/C_KMS) 
	    IF(T1 .GT. NU_LINE(LN_INDX)*RES_EXTENT)THEN
	      DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND. 
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=NU_CONT(ML)
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	        ML=ML+1
	      END DO
              IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=T1
	        LINES_THIS_FREQ(INDX)=0
	      END IF
C
C Only add in extra points if resonance zone is 1.2dV_CMF_WING away.
C
	      I=INT( (FREQ(INDX)-NU_LINE(LN_INDX)*RES_EXTENT)/
	1               FREQ(INDX)*C_KMS/dV_CMF_WING - 0.1)
	      DELF=(FREQ(INDX)-NU_LINE(LN_INDX)*RES_EXTENT)/(I+1)
	      DO J=1,I
	        T1=FREQ(INDX)-DELF
	        DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND. 
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	            INDX=INDX+1
	            IF(INDX .GT. NFREQ_MAX)GOTO 9999
	            FREQ(INDX)=NU_CONT(ML)
	            LINES_THIS_FREQ(INDX)=0
	          END IF
	          ML=ML+1
	        END DO
	        IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=T1
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	      END DO
	    END IF
C
C Now begin the resonance zone.
C
	    LINE_ST_INDX(LN_INDX)=INDX+1
	    DO I=-NDOP,NDOP
	      T1=NU_LINE(LN_INDX)*(1.0D0+dNU_on_NU)**(-I)
	      DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	        IF(EDGE_FREQ(ML) .AND. 
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=NU_CONT(ML)
	          LINES_THIS_FREQ(INDX)=1
	        END IF
	        ML=ML+1
	      END DO
	      IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT )THEN
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=T1
	        LINES_THIS_FREQ(INDX)=1
	      END IF
	    END DO
	    LINE_END_INDX(LN_INDX)=INDX
C
C Ready for next line or continuum point.
C
	    LST_LN_INDX=LN_INDX
	    LN_INDX=LN_INDX+1
	    DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1                           TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	      LN_INDX=LN_INDX+1
	    END DO
C
C Check if overlapping line [Overlap of intrinsic profile only!]
C
100	    CONTINUE
	    IF(LN_INDX .LE. N_LINES)THEN
	      IF(FREQ(INDX) .LT. NU_LINE(LN_INDX)*RES_EXTENT)THEN
C
C Check how far line extends back
C
	        J=INDX
	        DO WHILE(FREQ(J) .LT. NU_LINE(LN_INDX)*RES_EXTENT)
	          LINES_THIS_FREQ(J)=LINES_THIS_FREQ(J)+1
	          J=J-1
	        END DO
	        LINE_ST_INDX(LN_INDX)=J+1
C
	        DO WHILE(FREQ(INDX) .GT. NU_LINE(LN_INDX)/RES_EXTENT)
	          T1=FREQ(INDX)/(1.0D0+dNU_on_NU)
	          DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. T1)
	          IF(EDGE_FREQ(ML) .AND. 
	1                NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	              INDX=INDX+1
	              IF(INDX .GT. NFREQ_MAX)GOTO 9999
	              FREQ(INDX)=NU_CONT(ML)
	              LINES_THIS_FREQ(INDX)=1
	            END IF
	            ML=ML+1
	          END DO
	          IF(T1 .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	            INDX=INDX+1
	            IF(INDX .GT. NFREQ_MAX)GOTO 9999
	            FREQ(INDX)=T1
	            LINES_THIS_FREQ(INDX)=1
	          END IF
	        END DO
C
	        LINE_END_INDX(LN_INDX)=INDX
	        LST_LN_INDX=LN_INDX
	        LN_INDX=LN_INDX+1
	        DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1                           TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	          LN_INDX=LN_INDX+1
	        END DO
C
C Check another line
C
	        GOTO 100
	      END IF
	    END IF
C
C Now we need to extend line to allow for the outflow velocity. We only extend
C until the blue edge of the next line.
C
C We extend the frequencies with spacing dV_CMF_PROF (in km/s) from the
C resonance zone edge for another 2Vinf/kms. Beyond this freq (called
C SWITCH_FREQ) we use until dV_VMF_WING, which is the spacing in the e.s. wings.
C
	    MIN_FREQ=NU_LINE(LST_LN_INDX)/RES_EXTENT/
	1              (1.+RED_WING_EXT/C_KMS)
C
C The V_CMF_PROF is added to allow for some bleeding.
C
	    SWITCH_FREQ=NU_LINE(LST_LN_INDX)/RES_EXTENT/
	1              (1.+(2.0D0*VINF+dV_CMF_PROF)/C_KMS)
C
C We check that the minimum frequency does not extend beyond the
C resonance zone of the next line. As we will put a frequency at
C the beginning of the resonance zone, we back of a little bit.
C
	    T1=RES_EXTENT*(1.0D0+0.3D0*MIN(dV_CMF_PROF,dV_CMF_WING)/C_KMS)
	    IF(LN_INDX .GT. N_LINES)THEN
	    ELSE IF(MIN_FREQ .LT. NU_LINE(LN_INDX)*T1)THEN
	      MIN_FREQ=NU_LINE(LN_INDX)*T1
	    END IF
	    DO WHILE(FREQ(INDX) .GT. MIN_FREQ)
	      IF( FREQ(INDX) .GT. SWITCH_FREQ)THEN
	        TEMP_FREQ=FREQ(INDX)/(1.0D0+dV_CMF_PROF/C_KMS)
	      ELSE
	        TEMP_FREQ=FREQ(INDX)/(1.0D0+dV_CMF_WING/C_KMS)
	      END IF
C
C We need to check again, since we didn't know the frequency step size.
C
	      IF(TEMP_FREQ .GT. MIN_FREQ)THEN
	        DO WHILE(ML .LT. NCF .AND. NU_CONT(ML) .GT. TEMP_FREQ)
	        IF(EDGE_FREQ(ML) .AND. 
	1              NU_CONT(ML) .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	            INDX=INDX+1
	            IF(INDX .GT. NFREQ_MAX)GOTO 9999
	            FREQ(INDX)=NU_CONT(ML)
	            LINES_THIS_FREQ(INDX)=0
	          END IF
	          ML=ML+1
	        END DO
                IF(TEMP_FREQ .LT. FREQ(INDX)/MIN_FREQ_RAT)THEN
	          INDX=INDX+1
	          IF(INDX .GT. NFREQ_MAX)GOTO 9999
	          FREQ(INDX)=TEMP_FREQ
	          LINES_THIS_FREQ(INDX)=0
	        END IF
	      ELSE
	        GOTO 200	!Exit from this section
	      END IF
	    END DO
200	    CONTINUE
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
	    WRITE(LU_ER,*)' Invalid frequency grid computed in INS_LINE_V4'
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
	DO I=2,NFREQ
	  WRITE(63,'(1X,I6,1P,2E12.4)')I,FREQ(I),
	1            3.0D+05*(FREQ(I-1)-FREQ(I))/FREQ(I-1)
	END DO
C
!	DO I=1,N_LINES
!	  J=LINE_ST_INDX(I)
!	  K=LINE_END_INDX(I)
!          WRITE(83,'(1P,3E16.6,0P,2I7,A)')NU_LINE(I),FREQ(J),FREQ(K),J,K,TRIM(TRANS_TYPE(I))
!        END DO
!
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
