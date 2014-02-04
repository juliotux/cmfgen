C
C Subroutine to read in Low temperature dielectronic data for
C an arbitrary species.
C
C The user indicates whether to include NORMAL (ie AUTO) and WI LTDR
C transitions. LTDR transition are returned only for those levels included
C in the model atom.
C
C NB : AUTO --- Levels permitted to autoionize in pure LS coupling.
C      WI   --- Levels not permitted to autoionize in LS coupling.
C
C LTDR transition are returned only for those levels included in the model atom.
C
	SUBROUTINE RDGENDIE_V4(C2NAME,INDX_C2,NC2,
	1        EDGEDIE,EINADIE,GUPDIE,
	1        LEVDIE,INDXDIE,SPECDIE,TRANSDIE,GION,
	1        DO_AUTO,DO_WI,
	1        DESC,LUIN,LUOUT,FILNAME,NMAX,NUM_DIE)
	IMPLICIT NONE
C
C Altered 25-May-1996 - GEN_SEQ_OPEN removed, and replaced by F90 version.
C Altered 18-Dec-1995 - DO_AUTO inserted. Can now specify whether to
C                       include:
C                         Dielectronic transitions to levels which
C                           autoionize in LS coupling (AUTO)
C                         Dielectronic transitions to levels which
C                           do NOT autoionize in LS coupling (WI)
C                        Either, or both may be included.
C
C Altered 27-Nov-1990 - Internal free format reads removed. Open changed.
C                       LUER installed : Cray compatibility. NDIEST was
C                       removed : bug on error message (30-Dec-91).
C Altered 10-Oct-1990 - Number of LTDR (non WI) transitions can now be zero.
C Altered 07-Sep-1989 - Extensive rewriting. WI option included in call.
C                       No error if lower level is not in model atom.
C                       Fraction of LTDR recombinations included is indicated.
C                       GION also included in call.
C Created 23-Mar-1989.
C
C
	INTEGER NC2,INDX_C2
	INTEGER NMAX,NUM_DIE,LUIN,LUOUT
C
	REAL*8 GION			!St. Weight of ion (eg CIII).
	REAL*8 EDGEDIE(NMAX)
	REAL*8 EINADIE(NMAX)
	REAL*8 GUPDIE(NMAX)
	INTEGER LEVDIE(NMAX)  	!Indicates NL of low state
	INTEGER INDXDIE(NMAX)		!Species identification.
C
	LOGICAL DO_WI
	LOGICAL DO_AUTO
C
	CHARACTER*(*) DESC,FILNAME
	CHARACTER*(*) SPECDIE(NMAX),TRANSDIE(NMAX),C2NAME(NC2)
C
	REAL*8 CHIBF,CHIFF,HDKT,TWOHCSQ,OPLIN,EMLIN
	COMMON/CONSTANTS/ CHIBF,CHIFF,HDKT,TWOHCSQ
	COMMON/LINE/ OPLIN,EMLIN
C
	INTEGER ERROR_LU
	EXTERNAL ERROR_LU
C
C Local variables
C
	INTEGER I,LOOP,IOS,L1,L2,LUER
	INTEGER MNL,ML,NUM_D_RD
	INTEGER INC,WIINC,MIS,WIMIS,INDX_HASH
	REAL*8 T1
	CHARACTER*132 STRING
C
	REAL*8 A10,A20,A30
	REAL*8 A1,A2,A3,M1,M2,M3		!Effective recombination rate.
	REAL*8 WIA1,WIA2,WIA3,WIM1,WIM2,WIM3    !(included and missing).
C
C Variables for free-format internal reads. NB. A string length greater
C than 80 is required for NIIIDIE, and OIVDIE. 132 is an absolute maximum.
C
	REAL*8 RD_FREE_VAL
	INTEGER NEXT,STR_LEN
	CHARACTER*132 ER_DESC
	DATA STR_LEN/132/
C
	LUER=ERROR_LU()
C
	IF(LEN(SPECDIE(1)) .LT. LEN(DESC))THEN
	  WRITE(LUER,*)'Error in RDGENDIE_V4'
	  WRITE(LUER,*)'Length of SPECDIE is too small for species name'
	  WRITE(LUER,*)'Species is ',DESC
	  STOP
	END IF
C
C LEVEL is the name of the lower bound state, FL is the frequency of the
C stabilizing transition, GuP is the statistical weight of the autoionizing
C level, and EINADIE is the einstein A of the stabilizing transition.
C
	OPEN(UNIT=LUIN,FILE=FILNAME,STATUS='OLD',ACTION='READ',IOSTAT=IOS)
	IF(IOS .NE. 0)THEN
	  WRITE(LUER,*)'Error in RDGENDIE - cant open '//FILNAME
	  STOP
	END IF
C
50	  READ(LUIN,'(A)',IOSTAT=IOS)STRING
	    L1=INDEX(STRING,'!Number of dielectronic transitions')
	    IF(IOS .NE. 0)THEN
	      WRITE(LUER,*)'Error reading # of levels from '//FILNAME
	      WRITE(LUER,*)'IOSTAT=',IOS
	      STOP
	    END IF
	    WRITE(LUOUT,'(A)')STRING
	  IF(L1 .EQ. 0)GOTO 50
C
	  ER_DESC='NUM_D_RD read in RDGENDIE-'//FILNAME
	  NUM_D_RD=RD_FREE_VAL(STRING,1,STR_LEN,NEXT,ER_DESC)
C
C Skip all records until we come across the header.
C
999	READ(LUIN,'(A)')STRING
	L1=INDEX(STRING,'Transition')
	IF(L1 .EQ. 0)GOTO 999
	L1=INDEX(STRING,'Lam(A)')
	IF(L1 .EQ. 0)THEN
	  WRITE(LUER,*)'Error reading oscillator header in '//FILNAME
	  STOP
	END IF
C
C Skip blank record before reading in dielectronic transitions.
C
	  READ(LUIN,'(A)')STRING
	  WRITE(LUOUT,'(A)')STRING
	  IF(STRING .NE. ' ')THEN
	    WRITE(LUER,*)'Error reading blank(1) from '//FILNAME
	    STOP
	  END IF	
C
C Zero summation scalers for computing recombination rates.
C
	  A1=0.0D0
	  A2=0.0D0
	  A3=0.0D0
	  M1=0.0D0
	  M2=0.0D0
	  M3=0.0D0
	  WIA1=0.0D0
	  WIA2=0.0D0
	  WIA3=0.0D0
	  WIM1=0.0D0
	  WIM2=0.0D0
	  WIM3=0.0D0
	  INC=0				!LTDR counters.
	  WIINC=0
	  MIS=0
	  WIMIS=0
C
	  ML=NUM_DIE
	  DO LOOP=1,NUM_D_RD
C
	    ML=ML+1
C
C Check whether this LTDR transition can be accommodated.
C
	    IF(ML .GT. NMAX)THEN
	      WRITE(LUER,*)
	1         'NMAX not large enough in RDGENDIE - '//FILNAME
	      WRITE(LUER,*)'NMAX=',NMAX,'NUM_DIE=',NUM_DIE+NUM_D_RD
	      STOP
	    END IF
C
C Read in file header.
C
	    READ(LUIN,'(A)')STRING
C
C Allows for gaps between first level name and '-'.
C
	    L1=INDEX(STRING,'-')
	    L1=INDEX(STRING(L1+1:),'  ')+L1-1
	    IF( L1 .LE. 0)THEN
	      WRITE(LUER,*)
	1       'Error reading in Transition Names from '//FILNAME
	      STOP
	    END IF
	    IF( LEN(TRANSDIE(ML)) .LT. L1+LEN(DESC)+2 )THEN
	      WRITE(LUER,*)'Error - transition name in RDGENIE too small'
	      STOP
	    END IF
	    TRANSDIE(ML)=DESC//'('//STRING(1:L1)//')'
C
C The first T1 is the oscillator strength, the second is the transition
C wavelength.
C
	    NEXT=L1+1
	    ER_DESC='Oscillator read in RDGENDIE-'//FILNAME
	    T1=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,ER_DESC)
	    ER_DESC='Einstein A read in RDGENDIE-'//FILNAME
	    EINADIE(ML)=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,ER_DESC)
	    ER_DESC='TRansition lambda read in RDGENDIE-'//FILNAME
	    T1=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,ER_DESC)
	    ER_DESC='GUPDIE read in RDGENDIE-'//FILNAME
	    GUPDIE(ML)=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,ER_DESC)
	    ER_DESC='EDGEDIE read in RDGENDIE-'//FILNAME
	    EDGEDIE(ML)=RD_FREE_VAL(STRING,NEXT,STR_LEN,NEXT,ER_DESC)
	    EDGEDIE(ML)=-EDGEDIE(ML)   	  !Since above ionization limit.
C
C Now compute the effective recombination coefficient (units of
C 10^{-12} for 10^4, 2 x 10^4, and 3 x 10^4 K.
C
	    T1=HDKT*EDGEDIE(ML)
	    A10=2.07D-10*GUPDIE(ML)*EINADIE(ML)/GION
	    A20=A10*EXP(0.5D0*T1)/(2.0D0**1.5)
	    A30=A10*EXP(T1/3.0D0)/(3.0D0**1.5)
	    A10=A10*EXP(T1)
C
C Is this lower level in model atom ?
C
	    L2=INDEX(STRING,'-')-1
	    IF( L2 .LE. 0 .OR. L2 .GE. L1)THEN
	      WRITE(LUER,*)'Error reading in Level Names from '//FILNAME
	      WRITE(LUER,'(A)')STRING
	      WRITE(LUER,*)'L2=',L2,'L1=',L1
	      STOP
	    END IF
	    MNL=0
	    DO I=1,NC2
	      IF(STRING(1:L2) .EQ. C2NAME(I))MNL=I
	    END DO
C
C Is it a WI transition ?
C
	    INDX_HASH=INDEX(STRING,'#')
	    IF(INDX_HASH .EQ. 0)THEN			!Not WI transition.
	      IF(DO_AUTO .AND. MNL .NE. 0)THEN
	        A1=A1+A10
	        A2=A2+A20
	        A3=A3+A30
	        INC=INC+1
	        LEVDIE(ML)=MNL
	        INDXDIE(ML)=INDX_C2
	        SPECDIE(ML)=DESC
	      ELSE
	        M1=M1+A10
	        M2=M2+A20
	        M3=M3+A30
	        ML=ML-1
	        MIS=MIS+1
	      END IF
	    ELSE
	      IF(DO_WI .AND. MNL .NE. 0)THEN
	        WIA1=WIA1+A10
	        WIA2=WIA2+A20
	        WIA3=WIA3+A30
	        WIINC=WIINC+1
	        LEVDIE(ML)=MNL
	        INDXDIE(ML)=INDX_C2
	        SPECDIE(ML)=DESC
	      ELSE
	        WIM1=WIM1+A10
	        WIM2=WIM2+A20
	        WIM3=WIM3+A30
	        ML=ML-1
	        WIMIS=WIMIS+1
	      END IF
	    END IF
C
	  END DO
C
C Return total number of dielectronic transitions.
C
	  NUM_DIE=ML
C
	  IF( (INC+WIINC+MIS+WIMIS) .NE. NUM_D_RD)THEN
	     WRITE(LUER,*)'Error in RDGENDIE -'//FILNAME
	     WRITE(LUER,*)'Invalid summation of included transitions'
	     STOP
	  END IF
C
	  WRITE(LUOUT,900)
900	  FORMAT(/,' Summary of dielectronic transitions ',
	1                  'included (LS : WI) ',/,
	1 '[Units 10^-12] ( ) denotes percentage of LTDR NOT included')
	  IF( (A1+M1) .NE. 0)THEN
	    M1=100.0*M1/(A1+M1)
	    M2=100.0*M2/(A2+M2)
	    M3=100.0*M3/(A3+M3)
	    WRITE(LUOUT,1000)INC,(INC+MIS),A1,M1,A2,M2,A3,M3
1000	    FORMAT( 1X,I4,'(',I4,')',3( 2X,1PE10.3,'(',0PF6.2,')' )  )
	  END IF
C
	  IF( (WIA1+WIM1) .NE. 0)THEN
	    WIM1=100.0*WIM1/(WIA1+WIM1)
	    WIM2=100.0*WIM2/(WIA2+WIM2)
	    WIM3=100.0*WIM3/(WIA3+WIM3)
	    WRITE(LUOUT,1000)WIINC,(WIINC+WIMIS),WIA1,WIM1,WIA2,WIM2,
	1                  WIA3,WIM3
	  END IF
C
	CLOSE(UNIT=LUIN)
C
	RETURN
	END
