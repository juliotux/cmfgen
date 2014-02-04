	SUBROUTINE INS_LINE_OBS_V4(
	1		FREQ,NFREQ,NFREQ_MAX,
	1               NU_LINE,NU_STRT_LINE,VEC_VMIN_VDOP,TRANS_TYPE,
	1               N_LINES,INCL_ALL_LINES,
	1		NU_MAX,NU_MIN,VINF,
	1               FRAC_DOP_OBS,dV_OBS_PROF,dV_OBS_WING,dV_OBS_BIG,
	1               OBS_PRO_EXT_RAT,ES_WING_EXT,V_DOP)
	IMPLICIT NONE
!
! Altered 02-Jul-2000 : Complete rewrite. Changed from V3 to V4.
!                       Routine now correctly handles the case where lines are
!                       are ordered by the start frequency, rather than the
!                       central frequency.
! Altered 23-Nov-1998 : Based on INS_LINE_OBS_V2.
!                       Variable TRANS_TYPE and INC_ALL_LINES included in call.
!                       If INC_ALL_LINES is FALSE, only LINES treated in BLANK mode 
!                       are used to define the observers frequency grid.
! Altered 25-May-1996 : MIN(LN_INDX,N_LINES) inserted in 2 IF staements because
!                        of order in which CRAY executes the IF statements.
! Altered 17-May-1996 : Better computation of frequency grid to allow for
!                         coherent and noncoherent electron scattering, and to
!                         handle P Cygni stars.
!	                OBS_PRO_EXT_RAT,ES_WING_EXT,V_DOP parameters inserted.
!                       Now V2.
!
	INTEGER*4 NFREQ_MAX,N_LINES
	INTEGER*4 NFREQ				!Returned
!
! Vecters returned by subroutine:
!
! Line+continuum frequencies

	REAL*8 FREQ(NFREQ_MAX)			!New observers frequencies
!
! Passed vectors.
!
	REAL*8 NU_LINE(N_LINES)		!Line frequencies
	REAL*8 NU_STRT_LINE(N_LINES)	!Frequencies of start of res. zone.
	REAL*8 VEC_VMIN_VDOP(N_LINES)	!Minium Doppler velocity in km/s.
!
	CHARACTER*(*) TRANS_TYPE(N_LINES)
	LOGICAL      INCL_ALL_LINES
!
! Passed constants:
	REAL*8 VINF		!Terminal velocity of wind.
	REAL*8 dV_OBS_PROF
	REAL*8 dV_OBS_WING
	REAL*8 dV_OBS_BIG
	REAL*8 NU_MIN
	REAL*8 NU_MAX
	REAL*8 OBS_PRO_EXT_RAT
	REAL*8 ES_WING_EXT
	REAL*8 V_DOP
	REAL*8 FRAC_DOP_OBS
!
! Local variables.
!
	REAL*8 MAX_B_EXTENT		!Maximum blueward extent of line profile
	REAL*8 MAX_R_EXTENT		!Maximum  redward extent of line profile
	REAL*8 MAX_BW_EXTENT		!Maximum blueward extent of e.s. wings
	REAL*8 MAX_RW_EXTENT		!Maximum  redward extent of e.s. wings
!
	REAL*8 PROF_SPACING
	REAL*8 WING_SPACING
	REAL*8 BIG_SPACING
	REAL*8 T1,T2
	REAL*8 dNU
	REAL*8 NU_END_LINE
!
	INTEGER*4 INDX		!Current frequency index.
	INTEGER*4 LN_INDX	!Current line whose frequencies we are 
				!   installing.
!
	INTEGER*4 I,J,K		!Micellaneous loop variables.
	INTEGER*4 LU_ER
	REAL*8 C_KMS
!
! External functions
!
	INTEGER*4 ERROR_LU
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL ERROR_LU,SPEED_OF_LIGHT
!
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
!
! Check parameters.
!
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
!
! We assume that both lines and continuum are ordered from highest to
! lowest frequencies.
!
	MAX_B_EXTENT  = 1.0D0+(OBS_PRO_EXT_RAT*VINF+3.0D0*V_DOP)/C_KMS
	MAX_R_EXTENT  = 1.0D0-(OBS_PRO_EXT_RAT*VINF+3.0D0*V_DOP)/C_KMS
!
! BW and RW refer to the extent of the electron scattering wings. They are 
! defined from line center.
!                               
! When Vinf is small, extent of red wing will be primarily determined by
! thermal redistribution effects, and hence ES_WING_EXT is important. When
! Vinf is large (>> Velec) the "coherent" scattering will dominate, and
! the extent of the red wing will be determined by Vinf.
!
	MAX_RW_EXTENT = 1.0D0-(ES_WING_EXT+4.0D0*VINF)/C_KMS
	MAX_BW_EXTENT = 1.0D0+(ES_WING_EXT+VINF)/C_KMS
!
! Ensures BW extent is bigger than profile extent.
!
	MAX_BW_EXTENT = MAX(MAX_BW_EXTENT,MAX_B_EXTENT+2.0D0*dV_OBS_BIG/C_KMS)
!
! Spacing in km/s across various parts of the frequency spectrum. NB: These
! meanings have changed from version V3.
!
	PROF_SPACING = dV_OBS_PROF/C_KMS	!For the main wind profile
	WING_SPACING = dV_OBS_WING/C_KMS	!For the e.s. wings
	BIG_SPACING  = dV_OBS_BIG/C_KMS    	!For the continuum
!                                                
	LN_INDX=1
	IF(.NOT. INCL_ALL_LINES)THEN
	  DO WHILE(LN_INDX .LE. N_LINES .AND. 
	1                       TRANS_TYPE(LN_INDX)(1:3) .NE. 'BLA')
	    LN_INDX=LN_INDX+1
	  END DO
	END IF
	INDX=1         
	FREQ(1)=NU_MAX
! 
! We can determine the frequency grid. We always chose the smallest needed
! frequency spacing.
!
	DO WHILE( FREQ(INDX) .GT. NU_MIN/MAX_RW_EXTENT )
!
! Continuum
!
	  dNU=FREQ(INDX)*BIG_SPACING
!
! Electron scattering wings. Since dV_OBS_KMS is the same for all lines,
! we can exit once e.s. wing region is found.
!
	  K=LN_INDX
	  DO WHILE( K .LE. N_LINES .AND.
	1               FREQ(INDX)-dNU .LT. 
	1               NU_STRT_LINE(MIN(K,N_LINES))*MAX_BW_EXTENT )
	    IF(INCL_ALL_LINES .OR. TRANS_TYPE(K) .EQ. 'BLANK')THEN
	      T1=FREQ(INDX)-dNU
	      IF( T1 .LE. NU_LINE(K)*MAX_BW_EXTENT .AND.
	1            FREQ(INDX) .GE. NU_LINE(K)*MAX_RW_EXTENT)THEN
	        T2=FREQ(INDX)*WING_SPACING
	        dNU=MIN(dNU,T2)
	        EXIT
	      END IF
	    END IF
	    K=K+1
	  END DO
!
! Profile. As NU_STRT_LINE monotonically decreases, we can stop looking at
!          line K when the next frequency is > NU_STRT_LINE(K). Since
!          dV_OBS_PROF is the same for all lines, we can exit immediately.
!
	  K=LN_INDX
	  DO WHILE( K .LE. N_LINES .AND. FREQ(INDX)-dNU .LT. 
	1               NU_STRT_LINE(MIN(K,N_LINES))*MAX_B_EXTENT )
	    IF(INCL_ALL_LINES .OR. TRANS_TYPE(K) .EQ. 'BLANK')THEN
	      T1=FREQ(INDX)-dNU
	      IF( T1 .LE. NU_LINE(K)*MAX_B_EXTENT .AND.
	1            FREQ(INDX) .GE. NU_LINE(K)*MAX_R_EXTENT)THEN
	        T2=FREQ(INDX)*PROF_SPACING
	        dNU=MIN(dNU,T2)
	        EXIT
	      END IF
	    END IF
	    K=K+1
	  END DO
!
! Doppler core. If FRAC_DOP_OBS is large, this section will be effectively
!               ignored. This sections is primarily needed when we see
!               photospheric lines (e.g. as for an O star).
!
	  K=LN_INDX
	  DO WHILE( K .LE. N_LINES .AND.
	1               FREQ(INDX)-dNU .LT. NU_STRT_LINE(MIN(K,N_LINES)))
	    IF(INCL_ALL_LINES .OR. TRANS_TYPE(K) .EQ. 'BLANK')THEN
	      T1=FREQ(INDX)-dNU
	      NU_END_LINE=NU_LINE(K)-(NU_STRT_LINE(K)-NU_LINE(K))
	      IF( (FREQ(INDX) .GT. NU_END_LINE .AND. 
	1          T1 .LE. NU_STRT_LINE(K)) )THEN
	        T2=0.5D0*FREQ(INDX)*FRAC_DOP_OBS*VEC_VMIN_VDOP(K)*
	1             SQRT( 1.0D0+ C_KMS*ABS(1.0D0-NU_LINE(K)/FREQ(INDX))/
	1                   VEC_VMIN_VDOP(K) )/C_KMS
	        dNU=MIN(dNU,T2)
	      END IF
	    END IF
	    K=K+1
	    IF(K .GT. N_LINES)EXIT
	  END DO
!
	  FREQ(INDX+1)=FREQ(INDX)-dNU
	  INDX=INDX+1
!
! Check whether finished with line.
!
	  K=LN_INDX
	  DO WHILE( FREQ(INDX) .LT. NU_LINE(K)*MAX_RW_EXTENT)
	    K=K+1
	    IF(K .GT. N_LINES)EXIT
	  END DO
!
	  IF(K .LE. N_LINES .AND. .NOT. INCL_ALL_LINES)THEN
	    DO WHILE(TRANS_TYPE(K) .NE. 'BLANK')
	      K=K+1
	      IF(K .GT. N_LINES)EXIT
	    END DO
	  END IF
	  LN_INDX=K
!
! INDX+1 since INDX+1 refers to index of the next frequency.
!
	  IF(INDX+1 .GT. NFREQ_MAX)THEN
	    WRITE(LU_ER,*)'Error in OBS_INS_LINE_V4'
	    WRITE(LU_ER,*)'Insufficent storage locations'
	    WRITE(LU_ER,*)'NFREQ_MAX=',NFREQ_MAX
	    WRITE(LU_ER,*)'Current frequency is',FREQ(INDX)
	    WRITE(LU_ER,*)'Current line is',LN_INDX
	    WRITE(LU_ER,*)'Maximum number of lines is',N_LINES
	    STOP
	  END IF
!
	END DO
!
	IF(FREQ(INDX) .GT. NU_MIN)THEN
	  INDX=INDX+1
	  FREQ(INDX)=NU_MIN
	END IF
!
! Set the number of frequencies
!
	NFREQ=INDX
!
! Test for monotocity of frequencies.
!
	DO J=1,NFREQ-1
	  IF(FREQ(J) .LE. FREQ(J+1))THEN
	    LU_ER=ERROR_LU()
	    WRITE(LU_ER,*)' Invalid frequency grid computed in INS_OBS_LINE_V4'
	    WRITE(LU_ER,*)' J=',J
	    DO I=MAX(1,J-20),MIN(NFREQ,J+20)
	      WRITE(LU_ER,*)I,FREQ(I)
	    END DO
	    STOP
	  END IF
	END DO
!
	OPEN(UNIT=77,FILE='OBS_FREQ')
	  DO I=1,NFREQ-1
	    WRITE(77,'(1X,ES12.6,3X,F12.3)')
	1              FREQ(I),C_KMS*(1.0D0-FREQ(I+1)/FREQ(I))
	  END DO
	CLOSE(UNIT=77)
!
	RETURN
	END
