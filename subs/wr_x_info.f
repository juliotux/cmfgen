C
C The routines in this file have been created in order to produce
C a short formated summary of each model.
C
C Output atomic model information.
C
	SUBROUTINE WR_SL_INFO(STRING,NS,NF,ZION,DESC,LUOUT)
	IMPLICIT NONE
C
	INTEGER NF		!Number of level in full atom
	INTEGER NS		!Number of super levels
	INTEGER LUOUT		!Output unit for string if ``full''
	REAL*8 ZION		!Charge on ion (or core)
	CHARACTER*(*) STRING	!Output string
	CHARACTER*(*) DESC	!Species description (i.e. C2)
C
C Number of descriptors written per line, and the length each descriptor
C is assigned.
C
	INTEGER, PARAMETER :: NVAR_PER_LINE=4
	INTEGER, PARAMETER :: FIELD_LENGTH=20
C
	INTEGER NEXT_LOC
	INTEGER J
	INTEGER I_NF,I_NS
	REAL*8 T1
	CHARACTER*15 FMT
C
	T1=NVAR_PER_LINE
	NEXT_LOC=FIELD_LENGTH*MOD( (ZION-1), T1 )+1
	IF(STRING(NEXT_LOC:) .NE. ' ')THEN
	  WRITE(LUOUT,'(A)')TRIM(STRING)
	  STRING=' '
	END IF
	J=LEN_TRIM(DESC)
	STRING(NEXT_LOC:)=DESC(1:J)
	NEXT_LOC=NEXT_LOC+J
	I_NS=LOG10(FLOAT(NS))+1
	I_NF=LOG10(FLOAT(NF))+1
	WRITE(FMT,'(A,I1,A,I1,A)')'(A1,I',I_NF,',A,I',I_NS,',A)'
	WRITE(STRING(NEXT_LOC:),FMT)'[',NF,'/',NS,']'
C
	RETURN
	END
C
C Routine to output ND etc
C
	SUBROUTINE WR_INT_INFO(STRING,NEXT_LOC,DESC,IVAL)
	IMPLICIT NONE
	INTEGER NEXT_LOC	!Next location for output (updated)
	INTEGER IVAL		!Integer value to be output to string
	CHARACTER*(*) STRING	!Output string
	CHARACTER*(*) DESC	!Species description (i.e. C2)
C
C Number of descriptors written per line, and the length each descriptor
C is assigned.
C
!	INTEGER, PARAMETER :: NVAR_PER_LINE=4
	INTEGER, PARAMETER :: FIELD_LENGTH=20
C
	INTEGER I,J
	CHARACTER*15 FMT
C
C Output format is ND[60] etc
C
	J=LEN_TRIM(DESC)
	STRING(NEXT_LOC:)=DESC(1:J)
	NEXT_LOC=NEXT_LOC+J
	I=LOG10(FLOAT(ABS(IVAL)))+1
	IF(IVAL .LT. 0)I=I+1
	WRITE(FMT,'(A,I1,A)')'(A,I',I,',A)'
	WRITE(STRING(NEXT_LOC:),FMT)'[',IVAL,']'
	NEXT_LOC=NEXT_LOC+FIELD_LENGTH-J
C
	RETURN
	END
C
C Routine to output Rstar etc,
C
	SUBROUTINE WR_VAL_INFO(STRING,NEXT_LOC,DESC,VAL)
	IMPLICIT NONE
	INTEGER NEXT_LOC	!Next location for output (updated)
	REAL*8 VAL		!Real value to be output to string
	CHARACTER*(*) STRING	!Output string
	CHARACTER*(*) DESC	!Species description (i.e. C2)
C
C Number of descriptors written per line, and the length each descriptor
C is assigned.
C
!	INTEGER, PARAMETER :: NVAR_PER_LINE=4
	INTEGER, PARAMETER :: FIELD_LENGTH=20
C
	INTEGER I,J
	CHARACTER*15 FMT
C
C Ouput of the form R*=100.002 of Mdot=1.0E-05 
C
	J=LEN_TRIM(DESC)
	STRING(NEXT_LOC:)=DESC(1:J)//'='
	NEXT_LOC=NEXT_LOC+J+1
	IF(ABS(VAL) .LT. 1 .OR. ABS(VAL) .GT. 1.0E+04)THEN
	  FMT='(1P,E9.3)'
	  IF(VAL .LT. 0)FMT='(1P,E10.3)'
	ELSE
	  I=LOG10(ABS(VAL))+5
	  IF(VAL .LT. 0)I=I+1
	  WRITE(FMT,'(A,I1,A)')'(F',I,'.3)'
	END IF
	WRITE(STRING(NEXT_LOC:),FMT)VAL
	NEXT_LOC=NEXT_LOC+FIELD_LENGTH-J-1
C
	RETURN
	END	
C
C Routine to output abundance information: Both the relative number 
C abundance and the mass fractions are printed.
C
	SUBROUTINE WR_ABUND_INFO_V2(SPECIES,MASS,ABUND,ABUND_SUM,
	1           MEAN_ATOMIC_WEIGHT,SOL_MASS_FRAC,LUOUT)
C
C Altered 17-Dec-2007: Now output mean atomic mass.
C
	IMPLICIT NONE
	REAL*8 ABUND			!Relative abundance by number
	REAL*8 ABUND_SUM		!Sum of relative abundances
	REAL*8 MASS			!Mass in atomic mass units
	REAL*8 MEAN_ATOMIC_WEIGHT	!Mean ioic atomic weight (all)
	REAL*8 SOL_MASS_FRAC
	INTEGER LUOUT
	CHARACTER*(*) SPECIES		!e.g. HYD or CARB
C
	REAL*8 T1
	LOGICAL FIRST			!Indicate whether header should 
	DATA FIRST/.TRUE./		!    be output.
C
	IF(FIRST)THEN
	  WRITE(LUOUT,'(A)')' '
	  WRITE(LUOUT,'(A,F7.4)')'Mean atomic mass (amu) is: ',MEAN_ATOMIC_WEIGHT
	  WRITE(LUOUT,'(A)')' '
	  WRITE(LUOUT,'(3X,A,5X,A,5X,A,4X,A,6X,A)')
	1    'SPECIES','Rel. # Fraction','Mass Fraction','Z/Z(sun)','Z(sun)'
	  FIRST=.FALSE.    
	END IF
C
	T1=ABUND*MASS/MEAN_ATOMIC_WEIGHT/ABUND_SUM
	IF(ABUND .GT. 0.1)THEN
	  WRITE(LUOUT,'(5X,A,T16,F9.3,T38,1P,E9.3,6X,E8.2,5X,E8.2)')
	1    SPECIES,ABUND,T1,T1/SOL_MASS_FRAC,SOL_MASS_FRAC
	ELSE
	  WRITE(LUOUT,'(5X,A,T20,1P,E9.3,T38,E9.3,6X,E8.2,5X,E8.2)')
	1     SPECIES,ABUND,T1, T1/SOL_MASS_FRAC,SOL_MASS_FRAC
	END IF
C
	RETURN
	END
