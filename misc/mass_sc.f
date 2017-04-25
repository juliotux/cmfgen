	PROGRAM MASS_SC
	USE GEN_IN_INTERFACE
	IMPLICIT NONE
!
! Altered: 29-May-2008: Revised input/output files
! Altered:  8-Feb-2008: More than one set of Solar abundaces can be
!                      read in. FILENAME read inserted on 27-Feb-2008.
!
! Routine to read in a set of abundances (generally solar) from a file.
!
! Abundances of individual species can be adjusted. Species not adjusted
! are held fixed at a constant mass-fraction, and a revised X/He abundace
! output.
!
	INTEGER, PARAMETER :: MAX_EL=92
	INTEGER, PARAMETER :: MAX_ABUND=5
!
	REAL*8 AT_NO(MAX_EL)			!Atomic number
	REAL*8 AT_MASS(MAX_EL)			!Atomic mass(amu)
	REAL*8 ABUND(MAX_EL,MAX_ABUND)		!Fractional abundace (relative)
	REAL*8 MASS_FRAC(MAX_EL,MAX_ABUND)	!Mass fractional abundace.
	CHARACTER*2 SYMB(MAX_EL)		!Element symbol
	CHARACTER*20 NAME(MAX_EL)		!Element name
!
	REAL*8 NEW_ABUND(MAX_EL,MAX_ABUND)
	REAL*8 NEW_MASS_FRAC(MAX_EL,MAX_ABUND)
	LOGICAL*1 ALTERED_ABUND(MAX_EL,MAX_ABUND)		!Indicates revised abundace
!
	CHARACTER*2 SPEC			!Used for IO
	LOGICAL*1 SYMB_OK
	CHARACTER*80 STRING
	CHARACTER*80 FILENAME
	REAL*4 VAL
!
	REAL*8 MASS
	REAL*8 OLD_MASS_SUM
	REAL*8 NEW_MASS_SUM
	REAL*8 NHE
	INTEGER I,J,N
	INTEGER IOS
	INTEGER N_ABUND
!
! Open output file.
!
	WRITE(6,*)'Data will be output to ABUNDANCE_SUM'
	OPEN(UNIT=10,FILE='ABUNDANCE_SUMMARY',STATUS='UNKNOWN',ACTION='WRITE')
!
! Open input file.
!
	IOS=10
	FILENAME='sol_abund.dat'
	DO WHILE(IOS .NE. 0)
	  CALL GEN_IN(FILENAME,'File with multiple sets of Solar abudances')
	  OPEN(UNIT=20,FILE=FILENAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	  IF(IOS .NE. 0)THEN
	    WRITE(6,*)'Error opening ',TRIM(FILENAME)
	    WRITE(6,*)'IOS=',IOS
	  END IF
	END DO
!
! Input abundance data from approproiately formated file. Strings must
! be enclosed in quotes to allow free format read.
!
! We skip comments by finding the first record with "Hydrogen" in it.
!

	 STRING=' '
	 DO WHILE(INDEX(STRING,'!Number of abundance columns') .EQ. 0)
	   READ(20,'(A)')STRING
	 END DO
	 READ(STRING,*)N_ABUND
!
	 DO WHILE(INDEX(STRING,'Hydrogen') .EQ. 0)
	   READ(20,'(A)')STRING
	 END DO
	 BACKSPACE(20)
!
	 I=0
	 DO WHILE(1 .EQ. 1)
	   READ(20,*,END=100)AT_NO(I+1),SYMB(I+1),
	1             NAME(I+1),AT_MASS(I+1),(ABUND(I+1,J),J=1,N_ABUND)
	   I=I+1
	 END DO
100	 CONTINUE
	 N=I
	CLOSE(UNIT=20)
!
! Decide on format for the abundances.
!
	DO J=1,N_ABUND
	  IF(ABS(ABUND(1,J)-12.0D0) .LT. 0.1)THEN
	    DO I=1,N
	      ABUND(I,J)=10.0D0**(ABUND(I,J)-12.0D0)
	    END DO
	  END IF
	END DO
!
! Compute the mass fractions of all species.
!
	DO J=1,N_ABUND
	  MASS=0.0D0
	  DO I=1,N
	    MASS_FRAC(I,J)=AT_MASS(I)*ABUND(I,J)
	    MASS=MASS+MASS_FRAC(I,J)
	  END DO
	  MASS_FRAC(1:N,J)=MASS_FRAC(1:N,J)/MASS
	END DO
!
! Output for check.
!
	DO I=1,N
	  WRITE(6,200)AT_NO(I),SYMB(I),NAME(I),AT_MASS(I),(12.0D0+LOG10(ABUND(I,J)/ABUND(1,J)),
	1              ABUND(I,J),MASS_FRAC(I,J),J=1,N_ABUND)
	  WRITE(10,200)AT_NO(I),SYMB(I),NAME(I),AT_MASS(I),(12.0D0+LOG10(ABUND(I,J)/ABUND(1,J)),
	1              ABUND(I,J),MASS_FRAC(I,J),J=1,N_ABUND)
	END DO
	WRITE(6,'(A)')' '
200	FORMAT(1X,F4.1,3X,A2,3X,A10,4X,F5.1,5(F8.2,ES10.2,ES10.2))
!
! We now begin the section to allow individual abundances to be varied.
! Each species can be separately changed by using it chenical abbreviation.
!
	NEW_ABUND(:,:)=0.0D0
	SPEC=' '
	DO WHILE(SPEC(1:1) .NE. 'E')
	  SPEC='E'
	  CALL GEN_IN(SPEC,'Chemical symbol')
	  SYMB_OK=.FALSE.
	  DO J=1,N_ABUND
	    DO I=1,N
	      IF(SPEC(1:2) .EQ. SYMB(I))THEN
	        VAL=ABUND(I,J)
	        CALL GEN_IN(VAL,'Abund for '//NAME(I))
	        NEW_ABUND(I,J)=VAL
	        ALTERED_ABUND(I,J)=.TRUE.
	        SYMB_OK=.TRUE.
	      END IF
	    END DO
	    IF(.NOT. SYMB_OK .AND. SPEC(1:1) .NE. 'E')
	1        WRITE(6,*)'Invalid chemical symbol'
	    END DO
	  END DO
!
! Bow adjust the new abundances so that their combined mass fraction is the 
! same as C in the solar data. As a consequence all other species will have 
! the same-mass fraction.
!
! Compute the scale factors.
!
	DO J=1,N_ABUND
	  OLD_MASS_SUM=0.0D0
	  NEW_MASS_SUM=0.0D0
	  DO I=1,N
	    IF(ALTERED_ABUND(I,J))THEN
	      OLD_MASS_SUM=OLD_MASS_SUM+MASS_FRAC(I,J)
	      NEW_MASS_FRAC(I,J)=AT_MASS(I)*NEW_ABUND(I,J)
	      NEW_MASS_SUM=NEW_MASS_SUM+NEW_MASS_FRAC(I,J)
	    END IF
	  END DO
	END DO
!
! Do the actual scaling.
!
	DO J=1,N_ABUND
	  MASS=0.0D0
	  DO I=1,N
	    IF(ALTERED_ABUND(I,J))THEN
	      NEW_MASS_FRAC(I,J)=NEW_MASS_FRAC(I,J)*OLD_MASS_SUM/NEW_MASS_SUM
	    ELSE
	      NEW_MASS_FRAC(I,J)=MASS_FRAC(I,J)              
	    END IF
	    MASS=MASS+NEW_MASS_FRAC(I,J)
	  END DO
	  WRITE(6,*)MASS
!
! Copute the relative fractional populations.
!
	  NHE=0.0D0
	  DO I=1,N
	    NEW_ABUND(I,J)=NEW_MASS_FRAC(I,J)*MASS/AT_MASS(I)
	    IF(SYMB(I) .EQ. 'He')NHE=NEW_ABUND(I,J)
	  END DO
	  IF(NHE .EQ. 0.0D0)NHE=MAXVAL(NEW_ABUND)
!
! Normalize the abundances so the the He abundace is 1.0
!
	NEW_ABUND(1:N,J)=NEW_ABUND(1:N,J)/NHE
!
! With the new fractional abundace we compute the revised mass-fractions. 
! Acts as a check.
!
	  MASS=0.0D0
	  DO I=1,N
	    NEW_MASS_FRAC(I,J)=AT_MASS(I)*NEW_ABUND(I,J)
	    MASS=MASS+NEW_MASS_FRAC(I,J)
	  END DO
	  NEW_MASS_FRAC(1:N,J)=NEW_MASS_FRAC(1:N,J)/MASS
!
	  WRITE(6,'(34X,A,A)')'  N(old)     M(old)  ',
	1                   '   N(new)     M(new)  '
	  DO I=1,N
	    WRITE(6,300)AT_NO(I),SYMB(I),NAME(I),AT_MASS(I),ABUND(I,J),
	1              MASS_FRAC(I,J),NEW_ABUND(I,J),NEW_MASS_FRAC(I,J)
	  END DO
300	  FORMAT(1X,F4.1,3X,A2,3X,A10,4X,F5.1,3X,1PE8.2,3X,E8.2,3X,E8.2,
	1       3X,E8.2)
!
	  WRITE(10,'(34X,A,A)')'  N(old)     M(old)  ',
	1                   '   N(new)     M(new)  '
	  DO I=1,N
	    WRITE(10,300)AT_NO(I),SYMB(I),NAME(I),AT_MASS(I),ABUND(I,J),
	1               MASS_FRAC(I,J),NEW_ABUND(I,J),NEW_MASS_FRAC(I,J)
	  END DO
	END DO
	CLOSE(UNIT=10)
!
	STOP
	END
