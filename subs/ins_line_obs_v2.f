	SUBROUTINE INS_LINE_OBS_V2(
	1		FREQ,NFREQ,NFREQ_MAX,NU_LINE,N_LINES,
	1		NU_MAX,NU_MIN,VINF,
	1               dV_OBS_PROF,dV_OBS_WING,dV_OBS_BIG,
	1               OBS_PRO_EXT_RAT,ES_WING_EXT,V_DOP)
	IMPLICIT NONE
C
C Altered 25-May-1996 : MIN(LN_INDX,N_LINES) inserted in 2 IF staements because
C                        of order in which CRAY executes the IF statements.
C Altered 17-May-1996 : Better computation of frequency grid to allow for
C                         coherent and noncoherent electron scattering, and to
C                         handle P Cygni stars.
C	                OBS_PRO_EXT_RAT,ES_WING_EXT,V_DOP parameters inserted.
C                       Now V2.
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
C
C Passed constants:
	REAL*8 VINF		!Terminal velocity of wind.
	REAL*8 dV_OBS_PROF
	REAL*8 dV_OBS_WING
	REAL*8 dV_OBS_BIG
	REAL*8 NU_MIN
	REAL*8 NU_MAX
	REAL*8 OBS_PRO_EXT_RAT
	REAL*8 ES_WING_EXT
	REAL*8 V_DOP
C
C Local variables.
C
	REAL*8 MAX_B_EXTENT
	REAL*8 MAX_R_EXTENT
	REAL*8 MAX_BW_EXTENT
	REAL*8 MAX_RW_EXTENT
C
	REAL*8 PROF_SPACING
	REAL*8 WING_SPACING
	REAL*8 BIG_SPACING
	REAL*8 T1
C
	INTEGER INDX		!Current frequency index.
	INTEGER LN_INDX	!Current line whose frequencies we are 
				!   installing.
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
	  WRITE(LU_ER,*)'Invalid ES_WING_EXT in INS_LINE_OBS_V2'
	  WRITE(LU_ER,*)'ES_WING_EXT is measured in km/s.'
	  STOP
	END IF
	IF( OBS_PRO_EXT_RAT .LT. 1)THEN
	  WRITE(LU_ER,*)'Invalid OBS_PRO_EXT_RAT in INS_LINE_OBS_V2'
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
	MAX_BW_EXTENT = MAX(MAX_BW_EXTENT,MAX_B_EXTENT+2.0D0*dV_OBS_BIG)
C
	PROF_SPACING = 1.0D0-dV_OBS_PROF/C_KMS
	WING_SPACING = 1.0D0-dV_OBS_WING/C_KMS
	BIG_SPACING  = 1.0D0-dV_OBS_BIG/C_KMS
C                                                
	INDX=1         
	FREQ(1)=NU_MAX
	LN_INDX=1
	DO WHILE( FREQ(INDX) .GT. NU_MIN/MAX_RW_EXTENT )
	  IF(LN_INDX .GT. N_LINES .OR.
	1      FREQ(INDX)*BIG_SPACING .GT. NU_LINE(MIN(LN_INDX,N_LINES))*
	1         MAX_BW_EXTENT*(1.0D0+0.1D0*dV_OBS_WING/C_KMS) )THEN
	     INDX=INDX+1
	     IF(INDX .GT. NFREQ_MAX)GOTO 9999
	     FREQ(INDX)=FREQ(INDX-1)*BIG_SPACING
	  ELSE
	     IF(FREQ(INDX) .GT. NU_LINE(LN_INDX)*MAX_BW_EXTENT)THEN
	       INDX=INDX+1             
	       IF(INDX .GT. NFREQ_MAX)GOTO 9999
	       FREQ(INDX)=NU_LINE(LN_INDX)*MAX_BW_EXTENT
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
	     END IF
C
C Now do the spacing across the main line profile.
C	 
	     DO WHILE (FREQ(INDX)*PROF_SPACING
	1         .GT. NU_LINE(LN_INDX)*MAX_R_EXTENT)
	        INDX=INDX+1
	        IF(INDX .GT. NFREQ_MAX)GOTO 9999
	        FREQ(INDX)=FREQ(INDX-1)*PROF_SPACING
	     END DO
C
C As we're still not to the red edge, we do one more.
C
	     INDX=INDX+1
	     IF(INDX .GT. NFREQ_MAX)GOTO 9999
	     FREQ(INDX)=FREQ(INDX-1)*PROF_SPACING
C
	     LN_INDX=LN_INDX+1
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
	      END DO
C
C As we're still not to the red edge, we do one more.
C
	      INDX=INDX+1
	      IF(INDX .GT. NFREQ_MAX)GOTO 9999
	      FREQ(INDX)=FREQ(INDX-1)*PROF_SPACING
	      LN_INDX=LN_INDX+1
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
	    MIN_FREQ=MAX(NU_LINE(LN_INDX-1)*MAX_RW_EXTENT,MIN_FREQ)
	    DO WHILE (FREQ(INDX)*WING_SPACING .GT. MIN_FREQ)
	      INDX=INDX+1
	      IF(INDX .GT. NFREQ_MAX)GOTO 9999
	      FREQ(INDX)=FREQ(INDX-1)*WING_SPACING
	    END DO
	  END IF
	END DO
C
C Set the number of frequencies.
C
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
