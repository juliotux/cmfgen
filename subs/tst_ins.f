	PROGRAM TST_LINE
	IMPLICIT NONE
C
 	INTEGER NCF,NFREQ_MAX,N_LINES
	INTEGER NFREQ
	PARAMETER (NCF=2000)
	PARAMETER (N_LINES=1000)
	PARAMETER (NFREQ_MAX=NCF+N_LINES*21)
C
C Vecters returned by subroutine:
C
C Line+continuum frequencies
	REAL*8 FREQ(NFREQ_MAX)			!Continuum frequencies
	INTEGER LINES_THIS_FREQ(NFREQ_MAX) !Indicates that this 
						!  frequency has line 
						!  contriubutions,
C
	INTEGER LINE_ST_INDX(N_LINES)		!Start index for the line 
						!  in the NEW frequency array.
	INTEGER LINE_END_INDX(N_LINES)	!End index for the line 
						! in the NEW frequency array.
	CHARACTER*6 TRANS_TYPE(N_LINES)		!End index for the line 
C
C Passed vectors.
C
	REAL*8 NU_CONT(NCF)		!Continuum frequencies
	REAL*8 NU_LINE(N_LINES)		!Line frequencies
C
C Passed constants:
	REAL*8 VINF		!Terminal velocity of wind.
	REAL*8 V_DOP		!Doppler velocity (km/s).
	REAL*8 FRAC_DOP		!Indicates dNU across line in Doppler widths.
	REAL*8 MAX_DOP		!Half the extent of intrinsic profile
				!  in Doppler widths,
C
	REAL*8 dV_CMF_PROF
	REAL*8 dV_CMF_WING
	REAL*8 ES_WING_ExT
	REAL*8 R_CMF_WING_EXT
	INTEGER NCF1,N_LINES1
	INTEGER ML
	REAL*8 C_KMS
C
C External functions
C
	REAL*8 SPEED_OF_LIGHT
	EXTERNAL SPEED_OF_LIGHT
C
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
C
	OPEN(UNIT=10,FILE='CFDAT',STATUS='OLD',READONLY)
	DO ML=1,NCF
	  READ(10,*,END=100)NU_CONT(ML)
	END DO
100	NCF1=ML-1
C
	OPEN(UNIT=10,FILE='LINE',STATUS='OLD',READONLY)
	DO ML=1,N_LINES
	  READ(10,*,END=200)NU_LINE(ML)
	END DO
200	N_LINES1=ML-1
	TRANS_TYPE(1:N_LINES1)='BLANK'
C
	V_DOP=50.0		!km/s
	FRAC_DOP=1.0D0
	MAX_DOP=6.0D0
	VINF=840.0D0		!km/s
C
	dV_CMF_PROF=200.0  		!km/s
	dV_CMF_WING=300.0  		!km/s
	ES_WING_ExT=2500.0  		!km/s
	R_CMF_WING_EXT=3.0D0
C
	CALL INS_LINE_V4(FREQ,LINES_THIS_FREQ,NFREQ,NFREQ_MAX,
	1		NU_LINE,TRANS_TYPE,
	1               LINE_ST_INDX,LINE_END_INDX,N_LINES1,
	1		NU_CONT,NCF1,
	1		V_DOP,FRAC_DOP,MAX_DOP,VINF,dV_CMF_PROF,
	1               dV_CMF_WING,ES_WING_EXT,R_CMF_WING_EXT )
C
	DO ML=1,NFREQ-1
	  WRITE(17,'(I6,2X,1P2E15.6,3X,0P,F10.2,3X,I4)')
	1         ML,FREQ(ML), 0.01*C_KMS/FREQ(ML),
	1         C_KMS*(FREQ(ML)-FREQ(ML+1))/FREQ(ML),
	1         LINES_THIS_FREQ(ML)
	END DO
C
	DO ML=1,N_LINES
	  WRITE(18,'(1PE15.6,0P,3X,I6,3X,I6,F10.2)')
	1         NU_LINE(ML),LINE_ST_INDX(ML),LINE_END_INDX(ML),
	1         C_KMS*(FREQ(LINE_ST_INDX(ML))-FREQ(LINE_END_INDX(ML)))/
	1         NU_LINE(ML)
	END DO                           
C
	STOP
	END
