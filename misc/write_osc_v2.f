C
C General routine to write out Energy levels, Statistical Weights,
C Oscilator Strength, Transition wavelengths, and Einstein A values.
C The wavelength is in a vacuum for lambda < 2000 Ang, air otherwise.
C
C Output in a format that can be red by GENOSCIL.
C
	SUBROUTINE WRITE_OSC_V2(ACIV,TRANS,SECND,FEDGECIV,GCIV,
	1             CIVLEV,ARAD,GAM2,GAM4,KNOWN_LEVEL_ENERGY,
	1             IONCIV,ZSCR,LU,NCIV,
	1             BND_LEV,BRIEF,CIVOSC_DATE,DESCRIPTOR,FILENAME)
	IMPLICIT NONE
C
C Altered  17-Mar-1997 : Cleaned. FIlename inserted in CALL.
C                        Name changed form WRITE_OSCIL to WRITE_OSC.
C
C Modified 15-Aug-1990 : C now obtained by function call, and has exact value.
C                        Required change to GENOSCIL routine for consistency.
C                        Cleaned.
C Altered 29-Jun-1989 - NMAX was not beiing set if BND_LEV was true, and
C                       there were no levels above the ionization limit.
C                       Now subtract energy from ionization energy.
C                       Necessary if just handling triplets.
C Altered 19-Jul-1989 - FORMAT 120 altered. Check that no level names
C                       are identical.
C Altered 17-Mar-1989 - Extensively modifed.
C
	INTEGER LU,NCIV
C
C Constants for opacity etc.
C
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
	DOUBLE PRECISION CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
C
C Species arrays
C
	REAL*8 FEDGECIV(NCIV)			!Photoionization frequency
	REAL*8 ACIV(NCIV,NCIV)			!Oscilator strengths
	REAL*8 GCIV(NCIV)			!Statistical weights
	REAL*8 ARAD(NCIV)			!
	REAL*8 GAM2(NCIV)			!
	REAL*8 GAM4(NCIV)			!
	REAL*8 IONCIV,ZSCR
!
	LOGICAL KNOWN_LEVEL_ENERGY(NCIV)
C            
C On output, TRANS is an integer array containg the tranistion number
C numbered with all low i tranistions done first.
C
C SECND is not altered. IF secnd is TRUE, the oscilator strength is output
C as negative values.
C
	INTEGER TRANS(NCIV,NCIV)	
	LOGICAL  SECND(NCIV,NCIV)	
	CHARACTER*(*) CIVOSC_DATE
	CHARACTER*(*) CIVLEV(NCIV)
	CHARACTER*(*) DESCRIPTOR
	CHARACTER*(*) FILENAME
C
C Outputs only levels which are truely bound (i.e. those states with an
C ionization energy greater than 0).
C
	LOGICAL BND_LEV
C
C If BRIEF is false, the oscialtor strengths are not output in reverse
C order (i.e. listing of all transitions downward from level 2, downward 
C from level 3 etc.
C
	LOGICAL BRIEF
C
C External functions called.
C
	EXTERNAL SPEED_OF_LIGHT,LAMVACAIR
	REAL*8 SPEED_OF_LIGHT,LAMVACAIR
!
	REAL*8 NEW_ARAD(NCIV)			!
C
C Local Variables
C
	REAL*8 C_LIGHT,T1,T2
	REAL*8 NU_TO_EV
	INTEGER LNGTH,I,J,IOS
	INTEGER FMTGAP,NUM_TRAN,NMAX
	CHARACTER*1   FORMFEED
	CHARACTER*120 STARS
	CHARACTER*132 HEAD1
	CHARACTER*132 HEAD2
	CHARACTER*80 FMT_STR
!	DATA FORMFEED/Z'0C'/
C
	INTEGER, PARAMETER :: IZERO=0
	INTEGER, PARAMETER :: T_OUT=6
C
	FORMFEED=CHAR(12)                                
	C_LIGHT=SPEED_OF_LIGHT()		!cm s^-1
	NU_TO_EV=1.0D0/0.241798836		!Conversion factor from 10^15 Hz to eV
!
	DO I=1,LEN(STARS)
	  STARS(I:I)='*'
	END DO
C
C Compute the Einstein A coefficients, and the level lifetimes.
C
	T1=OPLIN/EMLIN*TWOHCSQ
	NEW_ARAD(:)=0.0D0
	DO I=1,NCIV-1
	  DO J=I+1,NCIV                
	    ACIV(J,I)=T1*ABS(ACIV(I,J))*GCIV(I)/GCIV(J)
	1     *( (FEDGECIV(I)-FEDGECIV(J))**2 )
	    NEW_ARAD(J)=NEW_ARAD(J)+ACIV(J,I)
	  END DO
	END DO
!
	CALL GEN_ASCI_OPEN(LU,'CHK_ARAD_','UNKNOWN',' ',' ',IZERO,IOS)
	  DO I=1,NCIV
	    IF(NEW_ARAD(I) .NE. 0)THEN
	      T1=100.0D0*(ARAD(I)/NEW_ARAD(I)-1.0D0)
	      WRITE(LU,'(2X,I5,3X,ES12.3,3X,ES12.3,F14.4)')I,NEW_ARAD(I),ARAD(I),T1
	    END IF
	  END DO
	CLOSE(UNIT=LU)
C
C If desired, only states which lie below the ground state
C ionization state will be considered. 
C
	NMAX=NCIV
	IF(BND_LEV)THEN
	  DO I=NCIV,1,-1
	    IF(FEDGECIV(I) .LT. 0)THEN
	       NMAX=I-1
	    END IF
	  END DO
	END IF
C
	DO I=1,NMAX-1
	  DO J=I+1,NMAX
	    IF( CIVLEV(I) .EQ. CIVLEV(J) )THEN
!	      WRITE(LU,1001)I,J
	      WRITE(T_OUT,1001)I,J,CIVLEV(I)
1001	      FORMAT(1X,'WARNING: Levels ',I3,' and ',I3,'have same name: ',A)
	    END IF
	  END DO
	END DO
C
C Determine number of transitions with non-zero f value.
C
	NUM_TRAN=0
	DO I=1,NMAX-1
	  DO J=I+1,NMAX
	    IF(ACIV(I,J) .NE. 0.0D0)THEN
	      NUM_TRAN=NUM_TRAN+1
	      TRANS(I,J)=NUM_TRAN
	    END IF
	  END DO
	END DO
!
! Find the maximum name length - for optimal formating.
!
	LNGTH=0
	DO I=1,NMAX
	  LNGTH=MAX( LEN_TRIM(CIVLEV(I)),LNGTH )
	END DO
	FMTGAP=2*(LNGTH-1)
	WRITE(6,*)NMAX,BND_LEV
	WRITE(6,*)NMAX,LNGTH
!
	HEAD1='      g        E(cm^-1)    10^15 Hz     eV       Lam(A)'
	HEAD2='       ID      ARAD       GAM2        GAM4'
	DO I=1,LNGTH+2
	  HEAD1=' '//HEAD1
	END DO
!
	CALL GEN_ASCI_OPEN(LU,FILENAME,'UNKNOWN',' ',' ',IZERO,IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(6,*)'Error in WRITE_OSC_V2'
	  WRITE(6,*)'Unable to open OSCILATOR file for output'
	  WRITE(6,*)'IOS=',IOS
	  WRITE(6,*)'FILE=',FILENAME
	  STOP
	END IF
C
	WRITE(LU,'(A)')STARS
	WRITE(LU,'(A)')' '
	WRITE(LU,'(17X,A,A)')'Energy Levels and statistical weights for ',
	1                      TRIM(DESCRIPTOR)
	WRITE(LU,'(A)')' '
        WRITE(LU,'(A,A)')TRIM(HEAD1),TRIM(HEAD2)
	WRITE(LU,'(A)')STARS
	WRITE(LU,'(A)')' '
!
	WRITE(LU,'(A,T50,''!Format date'')')'17-Oct-2000'
	WRITE(LU,'(A,T50,''!Date'')')TRIM(CIVOSC_DATE)
!
	WRITE(FMT_STR,'(I6)')NMAX
	DO WHILE(FMT_STR(1:1) .EQ. ' ')
	  FMT_STR(1:)=FMT_STR(2:)
	END DO
	WRITE(LU,'(A,T50,A)')TRIM(FMT_STR),'!Number of energy levels'
C
	WRITE(FMT_STR,'(F12.4)')IONCIV
	DO WHILE(FMT_STR(1:1) .EQ. ' ')
	  FMT_STR(1:)=FMT_STR(2:)
	END DO
	WRITE(LU,'(A,T50,A)')TRIM(FMT_STR),'!Ionization energy'

	WRITE(FMT_STR,'(F4.1)')ZSCR
	DO WHILE(FMT_STR(1:1) .EQ. ' ')
	  FMT_STR(1:)=FMT_STR(2:)
	END DO
	WRITE(LU,'(A,T50,A)')TRIM(FMT_STR),'!Screened nuclear charge'
C
	WRITE(FMT_STR,'(I10)')NUM_TRAN
	DO WHILE(FMT_STR(1:1) .EQ. ' ')
	  FMT_STR(1:)=FMT_STR(2:)
	END DO
	WRITE(LU,'(A,T50,A)')TRIM(FMT_STR),'!Number of transitions'
!
	WRITE(LU,'(A)')' '
	DO I=1,NMAX
	  J=I
	  IF(.NOT. KNOWN_LEVEL_ENERGY(I))J=-I
	  WRITE(LU,120)CIVLEV(I)(1:LNGTH),GCIV(I)
	1,  ( IONCIV-FEDGECIV(I)*1.0D+15/C_LIGHT )
	1,   FEDGECIV(I),FEDGECIV(I)*NU_TO_EV
	1,   ( 1.0D-07*C_LIGHT/FEDGECIV(I)),J 
	1,    NEW_ARAD(I),GAM2(I),GAM4(I)
	END DO
120	FORMAT(A,2X,F8.1,3X,F12.4,3X,F9.5,2X,F7.3,3X,
	1             ES10.3,2X,I6,3ES12.3)
C
C 
C
C Output Transitions to file. First Dtermine the Nuber of transitions.
C
	WRITE(FMT_STR,'(I2.2)')FMTGAP
	FMT_STR='(A12,'//FMT_STR(1:2)//'X,A1,11X,A1,10X,A6,6X,A3,4X,A8,/)'
	WRITE(LU,'(/,A1,/,A80,//,24X,A,A)')FORMFEED,STARS,
	1             'Oscillator strengths for ',DESCRIPTOR
	WRITE(LU,'(13X,A,//,A80,/)')
	1             'Wavelengths in air for lambda > 2000 Ang, else vacuum',STARS
	WRITE(LU,FMT_STR)'  Transition','f','A','Lam(A)','i-j','Trans. #'
C
	DO I=1,NMAX-1
	  DO J=I+1,NMAX
	    IF(ACIV(I,J) .NE. 0.0D0)THEN
	      T1=FEDGECIV(I)-FEDGECIV(J)		!Frequency (10^{15} Hz)
	      T1=LAMVACAIR(T1)				!Lam (Angstroms)
	      IF(.NOT. KNOWN_LEVEL_ENERGY(I) .OR.
	1          .NOT. KNOWN_LEVEL_ENERGY(J))T1=-T1
	      T2=ACIV(I,J)
	      IF(SECND(I,J))T2=-T2
              IF(ABS(T1) .LT. 100000.0)THEN
	        WRITE(LU,140)CIVLEV(I)(1:LNGTH),CIVLEV(J)(1:LNGTH),
	1         T2,ABS(ACIV(J,I)),T1,I,J,TRANS(I,J)
              ELSE
	        WRITE(LU,150)CIVLEV(I)(1:LNGTH),CIVLEV(J)(1:LNGTH),
	1         T2,ABS(ACIV(J,I)),T1,I,J,TRANS(I,J)
	      END IF
	    END IF
	  END DO
	END DO
C 
	IF(BRIEF)THEN
	  CLOSE(LU)
	  RETURN
	END IF
C
C Write transitions based on the upper level.
C
	WRITE(FMT_STR,'(I2.2)')FMTGAP
	FMT_STR='(A12,'//FMT_STR(1:2)//'X,A1,11X,A1,10X,A6,6X,A3,4X,A8,/)'
	WRITE(LU,'(/,A1,/,A80,//,24X,A,A)')FORMFEED,STARS,
	1             'Oscillator strengths for ',DESCRIPTOR
	WRITE(LU,'(13X,A,//,A80,/)')
	1             'Wavelengths in air for lambda > 2000 Ang, else vacuum',STARS
	WRITE(LU,FMT_STR)'  Transition','f','A','Lam(A)','i-j','Trans. #'
C
	DO I=2,NMAX
	  DO J=1,I-1
	    IF(ACIV(J,I) .NE. 0.0D0)THEN
	      T1=FEDGECIV(J)-FEDGECIV(I)		!Frequency (10^{15} Hz)
	      T1=LAMVACAIR(T1)				!Lam (Angstroms)
	      IF(.NOT. KNOWN_LEVEL_ENERGY(I) .OR.
	1          .NOT. KNOWN_LEVEL_ENERGY(J))T1=-T1
	      T2=ACIV(J,I)
	      IF( SECND(J,I) )T2=-T2
              IF(ABS(T1) .LT. 100000.0)THEN
	        WRITE(LU,140)CIVLEV(J)(1:LNGTH),CIVLEV(I)(1:LNGTH),
	1         T2,ABS(ACIV(I,J)),T1,J,I,TRANS(J,I)
              ELSE
	        WRITE(LU,150)CIVLEV(J)(1:LNGTH),CIVLEV(I)(1:LNGTH),
	1         T2,ABS(ACIV(I,J)),T1,J,I,TRANS(J,I)
	      END IF
	    END IF
	  END DO
	END DO
C
	CLOSE(LU)
C
	RETURN
C
140	FORMAT(A,'-',A,3X,1PE11.4,2X,E10.4,2X,0P,F10.3,4X,I4,
	1        '-',I4,3X,I6)
150	FORMAT(A,'-',A,3X,1PE11.4,2X,E10.4,2X,E10.3,4X,I4,
	1        '-',I4,3X,I6)

	END
