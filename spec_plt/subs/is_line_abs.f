C This program will calculate an absorption profile by computing an optical 
C depth at line center t0, with viogt profile.  It needs as input a 
C micro-turbulent parameter b(km/sec), a central wavelength lam0(ang),
C an ossillator (sp) strength f(dimensionless), and a line width (rather a 
C damping parameter) gam(sec^-1).  The temperature is also input but has little
C effect.  The voigt profile in IDL has a problem with
C negative independent variables so we do a onesided calculation and assume 
C symmetry.
C Started on 1-9-92 by SRM
C                      
C This program is a modification of absprof.pro.  This will calculate a
C set of optical depths as a function of wavelength for the Lyman series of 
C hydrogen.  When the optical depths get too large I truncate them at a 
C particular maximum finite value so that when I exponentiate and 
C convolve with the point spread function of the telescope I will (hopefully)
C get a well behaved transmission function which can be divided out of the
C EZ CMa spectrum to revel the ``true'' spectrum.  This modification started
C on 4-30-92 by SRM.
C
	SUBROUTINE IS_LINE_ABS(WAVE,FLUX,NLAM,V_TURB,LOG_NTOT,T_IN_K)
	IMPLICIT NONE
C
	INTEGER NLAM
	REAL*8 WAVE(NLAM)
	REAL*8 FLUX(NLAM)
	REAL*8 V_TURB		!km/s
	REAL*8 LOG_NTOT
	REAL*8 T_IN_K
C
C Local variables
C
	INTEGER NLINES_MAX
	PARAMETER (NLINES_MAX=200)
	REAL*8 LAM_ZERO(NLINES_MAX)
	REAL*8 NU_ZERO(NLINES_MAX)
	REAL*8 OSC(NLINES_MAX)
	REAL*8 GAM(NLINES_MAX)
	REAL*8 GL(NLINES_MAX)
	REAL*8 GUP(NLINES_MAX)
	REAL*8 VDOP(NLINES_MAX)
	REAL*8 VRAD(NLINES_MAX)
	REAL*8 NCOL(NLINES_MAX)
	REAL*8 NU_DOP(NLINES_MAX)
C
	REAL*8 CHIL(NLINES_MAX)
C
	REAL*8 C_KMS,PI
	REAL*8 OPLIN
	REAL*8 TAU
	REAL*8 NTOT
	REAL*8 PHI
	REAL*8 FREQ
	REAL*8 a
	REAL*8 v
	REAL*8 T1
C
	INTEGER I,J,IOS
	INTEGER NLINES
	CHARACTER(LEN=80) STRING
C
C Functions
C
	REAL*8 SPEED_OF_LIGHT
	REAL*8 FUN_PI
	REAL*8 VOIGT
	EXTERNAL SPEED_OF_LIGHT,FUN_PI,VOIGT
C
C First read in the atomic data, kindly provided by Chuck Bowers (CWB).
C
C This is for the first 49 lines of the Lyman series.
C
C 5th column: Wavelength (Angstroms)
C 6th column: Oscilator strength
C 7th column: G (ground state)
C 8th column: Damping parameter
C
	OPEN(UNIT=10,FILE='IS_LINE_LIST',ACTION='READ',STATUS='OLD',IOSTAT=IOS)
        IF(IOS .NE. 0)THEN
	   WRITE(6,*)'Unable to open IS_LINE_LIST FILE'
	   RETURN
	END IF
	I=0; NLINES=0
	DO WHILE(1 .EQ. 1)
	    STRING(1:1)='!'
	    DO WHILE(STRING(1:1) .EQ. '!' .OR. STRING .EQ. ' ')
	      READ(10,'(A)',END=10)STRING
	    END DO
	    STRING=ADJUSTL(STRING)
	    STRING=STRING(INDEX(STRING,' ')+1:)
	    I=I+1
	    READ(STRING,*,IOSTAT=IOS)LAM_ZERO(I),GL(I),GUP(I),OSC(I),NCOL(I),VDOP(I),VRAD(I)
	    IF(IOS .NE. 0)THEN
	      WRITE(6,*)'Error reading interstellar line data - erronous record follows'
	      WRITE(6,*)TRIM(STRING)
	      GOTO 10
	    END IF
	    NCOL(I)=10.0**NCOL(I)
	    NLINES=I
        END DO
10	CONTINUE
	CLOSE(UNIT=10)
	IF(NLINES .EQ. 0)THEN
	  WRITE(6,*)'Unable to read LINE DATA from IS_LINE_LIST in IS_LINE_ABS'
	  RETURN
	ELSE
	  WRITE(6,*)'Number of lines read in is',NLINES
	END IF
C
	C_KMS=1.0D-05*SPEED_OF_LIGHT()
	PI=FUN_PI()
C
	OPLIN=2.6540081E-02		!pi*e*e/m/c
	DO I=1,NLINES
	  NU_ZERO(I)=0.01*C_KMS*(1.0D0-VRAD(I)/C_KMS)/LAM_ZERO(I)   		!10^15 Hz
!	  NU_DOP(I)=12.85*NU_ZERO(I)*
!	1             SQRT( (T_IN_K/1.0D+04)+ (V_TURB/12.85)*2 )/C_KMS
	  NU_DOP(I)=NU_ZERO(I)*VDOP(I)/C_KMS
	  CHIL(I)=1.0D-15*OPLIN*OSC(I)*NCOL(I)/SQRT(PI)/NU_DOP(I)
	END DO
C
	DO J=1,NLAM
	  TAU=0.0D0
	  FREQ=0.01*C_KMS/WAVE(J)
	  DO I=1,NlINES
	      v=(FREQ-NU_ZERO(I))/NU_DOP(I)
	      IF(ABS(V) .LE. 100.0)THEN                  !Was 10
	        a=1.0D-15*GAM(I)/4/PI/NU_DOP(I)
	        PHI=VOIGT(a,v)    
	        TAU=TAU+CHIL(I)*PHI
	        FLUX(J)=FLUX(J)*EXP(-TAU)
	      END IF
	  END DO
	END DO
C
	RETURN
	END
