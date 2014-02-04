	SUBROUTINE INS_LINE_OBS_V3(
	1		FREQ,NFREQ,NFREQ_MAX,
	1               NU_LINE,TRANS_TYPE,N_LINES,INCL_ALL_LINES,
	1		NU_MAX,NU_MIN,VINF,
	1               dV_OBS_PROF,dV_OBS_WING,dV_OBS_BIG,
	1               OBS_PRO_EXT_RAT,ES_WING_EXT,V_DOP)
	IMPLICIT NONE
C
C Altered 16-Jul-2008 : Adjusted MAX_B_EXTENT etc to accomodate velocities close to C. 
C Created 23-Nov-1998 : Based on INS_LINE_OBS_V2.
C                       Variable TRANS_TYPE and INC_ALL_LINES included in call.
C                       If INC_ALL_LINES is FALSE, only LINES treated in BLANK mode 
C                       are used to define the observers frequency grid.
C Altered 25-May-1996 : MIN(LN_INDX,N_LINES) inserted in 2 IF staements because
C                        of order in which CRAY executes the IF statements.
C Altered 17-May-1996 : Better computation of frequency grid to allow for
C                         coherent and noncoherent electron scattering, and to
C                         handle P Cygni stars.
C	                OBS_PRO_EXT_RAT,ES_WING_EXT,V_DOP parameters inserted.
C                       Now V2.
C Altered 12-Feb-2006 : Fixed bug to ensure that last frequency is NU_MIN.
C
	INTEGER NFREQ_MAX,N_LINES
	INTEGER NFREQ				!Returned
C
C Vecters returned by subroutine:
C
C Line+continuum frequencies
	REAL*8 FREQ(NFREQ_MAX)			!New observers frequencies
C
C Passed vectors.
C
	REAL*8 NU_LINE(N_LINES)		!Line frequencies
	CHARACTER*(*) TRANS_TYPE(N_LINES)
	LOGICAL      INCL_ALL_LINES
C
C Passed constants:
	REAL*8 VINF			!Terminal velocity of wind.
	REAL*8 dV_OBS_PROF		!Spcing across line profile
	REAL*8 dV_OBS_WING		!Spacing in e.s. wing
	REAL*8 dV_OBS_BIG		!Spacing between lines (i.e., in continuum)
	REAL*8 NU_MIN			!Minimum frequency
	REAL*8 NU_MAX			!Maximum frequency
	REAL*8 OBS_PRO_EXT_RAT
	REAL*8 ES_WING_EXT		!Extent of e.s. wing (thermal; in km/s)
	REAL*8 V_DOP
C
C Local variables.
C
	REAL*8 MAX_B_EXTENT		!Blue profile extent
	REAL*8 MAX_R_EXTENT		!Red profile extent
	REAL*8 MAX_BW_EXTENT		!Blue e.s. wing extent (from line core)
	REAL*8 MAX_RW_EXTENT		!Red e.s. wing etent (from line core)
C
	REAL*8 PROF_SPACING
	REAL*8 WING_SPACING
	REAL*8 BIG_SPACING
	REAL*8 T1
C
	INTEGER INDX		!Current frequency index.
	INTEGER LN_INDX		!Current line whose frequencies we are 
				!   installing.
	INTEGER LST_LN_INDX	!Index of last line whose frequencies we
				!   installed. Needed as lines computed in
                                !   SOB or CMF mode may not be included.
C
	INTEGER I,J		!Micellaneous loop variables.
	INTEGER LU_ER
	REAL*8 C_KMS,MIN_FREQ
C
C External functions
C
	INTEGER ERROR_LU
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
C
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
C
C Check parameters.
C
	LU_ER=ERROR_LU()
	IF( (ES_WING_EXT .GT. 0 .AND. ES_WING_EXT .LT. 50) .OR.
	1      ES_WING_EXT .LT. 0)THEN
	  WRITE(LU_ER,*)'Invalid ES_WING_EXT in INS_LINE_OBS_V3'
	  WRITE(LU_ER,*)'ES_WING_EXT is measured in km/s.'
	  STOP
	END IF
	IF( OBS_PRO_EXT_RAT .LT. 1)THEN
	  WRITE(LU_ER,*)'Invalid OBS_PRO_EXT_RAT in INS_LINE_OBS_V3'
	  WRITE(LU_ER,*)'OBS_PRO_EXT_RAT=',OBS_PRO_EXT_RAT
	  WRITE(LU_ER,*)'OBS_PRO_EXT_RAT must >= to unity.'
	  STOP
	END IF
C
C We assume that both lines and continuum are ordered from highest to
C lowest frequencies.
C
	MAX_B_EXTENT  = 1.0D0+(OBS_PRO_EXT_RAT*VINF+3.0D0*V_DOP)/C_KMS
	MAX_R_EXTENT  = 1.0D0-(OBS_PRO_EXT_RAT*VINF+3.0D0*V_DOP)/C_KMS
C
C BW and RW refer to the extent of the electron scattering wings. They are 
C defined from line center.
C                               
C When Vinf is small, extent of red wing will be primarily determined by
C thermal redistribution effects, and hence ES_WING_EXT is important. When
C Vinf is large (>> Velec) the "coherent" scattering will dominate, and
C the extent of the red wing will be determined by Vinf.
C
	MAX_RW_EXTENT = 1.0D0-(ES_WING_EXT+4.0D0*VINF)/C_KMS
	MAX_BW_EXTENT = 1.0D0+(ES_WING_EXT+VINF)/C_KMS
C
C Ensures BW extent is bigger than profile extent.
C
	MAX_BW_EXTENT = MAX(MAX_BW_EXTENT,MAX_B_EXTENT+2.0D0*dV_OBS_BIG/C_KMS)
C
	MAX_R_EXTENT=MAX(0.2D0,MAX_R_EXTENT)
	MAX_B_EXTENT=MIN(1.8D0,MAX_B_EXTENT)
	MAX_RW_EXTENT=MAX(0.19D0,MAX_RW_EXTENT)
	MAX_BW_EXTENT=MIN(1.81D0,MAX_BW_EXTENT)
C
	PROF_SPACING = 1.0D0-dV_OBS_PROF/C_KMS
	WING_SPACING = 1.0D0-dV_OBS_WING/C_KMS
	BIG_SPACING  = 1.0D0-dV_OBS_BIG/C_KMS
C                                                
	LN_INDX=1
	IF(.NOT. INCL_ALL_LINES)THEN
	  DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1                       TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	    LN_INDX=LN_INDX+1
	  END DO
	END IF
C
	INDX=1         
	FREQ(1)=NU_MAX
	DO WHILE( FREQ(INDX) .GT. NU_MIN/MAX_RW_EXTENT )
	  IF(LN_INDX .GT. N_LINES .OR.
	1      FREQ(INDX)*BIG_SPACING .GT. NU_LINE(MIN(LN_INDX,N_LINES))*
	1         MAX_BW_EXTENT*(1.0D0+0.1D0*dV_OBS_WING/C_KMS) )THEN
	     INDX=INDX+1
	     IF(INDX .GT. NFREQ_MAX)GOTO 9999
	     FREQ(INDX)=FREQ(INDX-1)*BIG_SPACING
	     IF(FREQ(INDX) .LT. NU_MIN)GOTO 1000
	  ELSE
	     IF(FREQ(INDX) .GT. NU_LINE(LN_INDX)*MAX_BW_EXTENT)THEN
	       INDX=INDX+1             
	       IF(INDX .GT. NFREQ_MAX)GOTO 9999
	       FREQ(INDX)=NU_LINE(LN_INDX)*MAX_BW_EXTENT
	       IF(FREQ(INDX) .LT. NU_MIN)GOTO 1000
	     END IF
C
	     IF(FREQ(INDX) .GT. NU_LINE(LN_INDX)*MAX_B_EXTENT)THEN
	       T1=FREQ(INDX)-NU_LINE(LN_INDX)*MAX_B_EXTENT
	       I=C_KMS*T1/(FREQ(INDX)+0.5*T1)/dV_OBS_WING
	       T1=T1/(I+1)
	       DO J=1,I
	         INDX=INDX+1
	         IF(INDX .GT. NFREQ_MAX)GOTO 9999
	         FREQ(INDX)=FREQ(INDX-1)-T1
	       END DO
	       INDX=INDX+1
	       IF(INDX .GT. NFREQ_MAX)GOTO 9999
	       FREQ(INDX)=NU_LINE(LN_INDX)*MAX_B_EXTENT
	       IF(FREQ(INDX) .LT. NU_MIN)GOTO 1000
	     END IF
C
C Now do the spacing across the main line profile.
C	 
	     DO WHILE (FREQ(INDX)*PROF_SPACING
	1         .GT. NU_LINE(LN_INDX)*MAX_R_EXTENT)
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=FREQ(INDX-1)*PROF_SPACING
	        IF(FREQ(INDX) .LT. NU_MIN)GOTO 1000
	     END DO
C
C As we're still not to the red edge, we do one more.
C
	     INDX=INDX+1
	     IF(INDX .GT. NFREQ_MAX)GOTO 9999
	     FREQ(INDX)=FREQ(INDX-1)*PROF_SPACING
C
	     LST_LN_INDX=LN_INDX
	     LN_INDX=LN_INDX+1
	     IF(.NOT. INCL_ALL_LINES)THEN
	       DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1                         TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	        LN_INDX=LN_INDX+1
	       END DO
	     END IF
C
C Check if overlapping line [Overlap of intrinsic profile only!]
C [NB. Inserting red wing frequencies will be equally useful as blue wing
C   frequencies on the next line., hence compare with MAX_B_EXTENT]
C
	    DO WHILE (LN_INDX .LE. N_LINES .AND.
	1                  FREQ(INDX)*WING_SPACING .LT. 
	1                  NU_LINE(MIN(LN_INDX,N_LINES))*MAX_B_EXTENT)
	      DO WHILE (FREQ(INDX)*PROF_SPACING .GT. 
	1                   NU_LINE(LN_INDX)*MAX_R_EXTENT)
	         INDX=INDX+1
	         IF(INDX .GT. NFREQ_MAX)GOTO 9999
	         FREQ(INDX)=FREQ(INDX-1)*PROF_SPACING
	         IF(FREQ(INDX) .LT. NU_MIN)GOTO 1000
	      END DO
C
C As we're still not to the red edge, we do one more.
C
	      IF(FREQ(INDX) .GT. NU_LINE(LN_INDX)*MAX_R_EXTENT)THEN
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=FREQ(INDX-1)*PROF_SPACING
	        IF(FREQ(INDX) .LT. NU_MIN)GOTO 1000
	      END IF
	      LST_LN_INDX=LN_INDX
	      LN_INDX=LN_INDX+1
	      IF(.NOT. INCL_ALL_LINES)THEN
	         DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1                         TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	         LN_INDX=LN_INDX+1
	         END DO
	      END IF
	    END DO
C
C Now we need to extend line to allow for electron scattering. We only extend
C until the blue (profile) edge of the next line.
C
	    IF(LN_INDX .GT. N_LINES)THEN
	      MIN_FREQ=NU_MIN
	    ELSE
	      MIN_FREQ=NU_LINE(LN_INDX)*MAX_B_EXTENT
	    END IF
	    MIN_FREQ=MAX(NU_LINE(LST_LN_INDX)*MAX_RW_EXTENT,MIN_FREQ)
	    DO WHILE (FREQ(INDX)*WING_SPACING .GT. MIN_FREQ)
	      INDX=INDX+1
	      IF(INDX .GT. NFREQ_MAX)GOTO 9999
	      FREQ(INDX)=FREQ(INDX-1)*WING_SPACING
	      IF(FREQ(INDX) .LT. NU_MIN)GOTO 1000
	    END DO
	  END IF
	END DO
C
C Set the number of frequencies.
C
1000	CONTINUE
	IF(FREQ(INDX-1) .GT. NU_MIN)FREQ(INDX)=NU_MIN
	NFREQ=INDX
C
C Test for monotocity of frequencies.
C
	DO J=1,NFREQ-1
	  IF(FREQ(J) .LE. FREQ(J+1))THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)' Invalid frequency grid computed in INS_OBS_LINE'
	    WRITE(LU_ER,*)' J=',J
	    DO I=MAX(1,J-20),MIN(NFREQ,J+20)
	      WRITE(LU_ER,*)I,FREQ(I)
	    END DO
	    STOP
	  END IF
	END DO
C
	RETURN
C
9999	CONTINUE
	LU_ER=ERROR_LU()
	WRITE(LU_ER,*)'Error in INS_LINE_OBS --- insufficient frequencies '//
	1               ' to store both line and continuum frequencies'
	WRITE(LU_ER,*)'LN_INDX= ',LN_INDX
	WRITE(LU_ER,*)'INDX= ',INDX
	WRITE(LU_ER,*)'NFREQ_MAX= ',NFREQ_MAX
	STOP
	END
